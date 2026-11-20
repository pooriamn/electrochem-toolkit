"""Mixed-potential corrosion model and three ways of estimating corrosion rate.

A corroding electrode is the sum of independent partial reactions (Wagner and
Traud). Here: metal dissolution, oxygen reduction with a diffusion limit, and
hydrogen evolution, all as Tafel lines. The potentiostat sets the applied
potential; the interface sees that minus the uncompensated ohmic drop.

Units: potential in V, current density in A/cm2, resistance in ohm cm2,
Tafel slopes in V/decade. Corrosion rate follows ASTM G102.
"""
from dataclasses import dataclass, field
import math

import numpy as np
from scipy.optimize import least_squares

LN10 = math.log(10)
K1_MM_G_PER_UA_CM_Y = 3.27e-3          # ASTM G102, CR in mm/y with i in uA/cm2


@dataclass(frozen=True)
class Tafel:
    """i = i_ref * 10^(sign * (E - E_ref) / slope), optionally diffusion limited."""
    i_ref: float
    E_ref: float
    slope: float
    anodic: bool
    i_lim: float = math.inf

    def __post_init__(self):
        if not (self.i_ref > 0 and self.slope > 0 and self.i_lim > 0):
            raise ValueError("i_ref, slope and i_lim must be positive")

    def current(self, E):
        sign = 1 if self.anodic else -1
        ik = self.i_ref * 10 ** (sign * (np.asarray(E, float) - self.E_ref) / self.slope)
        if math.isfinite(self.i_lim):
            ik = ik * self.i_lim / (ik + self.i_lim)
        return ik if self.anodic else -ik


@dataclass(frozen=True)
class Electrode:
    reactions: tuple
    r_u: float = 0.0                       # uncompensated resistance, ohm cm2
    c_dl: float = 1e-4                     # double-layer capacitance, F/cm2
    name: str = ""

    def faradaic(self, E):
        return sum(r.current(E) for r in self.reactions)

    def anodic_metal(self, E):
        return sum(r.current(E) for r in self.reactions if r.anodic)

    def corrosion_point(self, lo=-1.5, hi=0.5):
        """True (Ecorr, icorr) where the net Faradaic current is zero."""
        E = _bisect(lambda e: self.faradaic(e), np.array([lo]), np.array([hi]))[0]
        return float(E), float(self.anodic_metal(E))


def _bisect(f, lo, hi, iters=80):
    """Vectorised bisection for monotonic increasing f."""
    lo, hi = lo.astype(float).copy(), hi.astype(float).copy()
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        up = f(mid) > 0
        hi = np.where(up, mid, hi)
        lo = np.where(up, lo, mid)
    return 0.5 * (lo + hi)


def measured_current(electrode, E_applied, scan_rate=0.0):
    """Current a potentiostat would record at the applied potentials.

    Solves E_int = E_app - i(E_int) R_u point by point, then adds the
    double-layer charging current for a linear sweep.
    """
    E_app = np.asarray(E_applied, float)
    if electrode.r_u == 0:
        i = electrode.faradaic(E_app)
    else:
        g = lambda e: e + electrode.faradaic(e) * electrode.r_u - E_app
        E_int = _bisect(g, E_app - 2.0, E_app + 2.0)
        i = electrode.faradaic(E_int)
    return i + electrode.c_dl * scan_rate


def potentiodynamic_scan(electrode, span=0.25, step=0.001, scan_rate=1.667e-4,
                         noise=0.02, floor=5e-9, seed=0):
    """Synthetic anodic-going scan from Ecorr - span to Ecorr + span.

    The default rate is 0.6 V/h. Noise is proportional to the current plus a
    small absolute floor.
    """
    ecorr, _ = electrode.corrosion_point()
    E = np.arange(ecorr - span, ecorr + span + step / 2, step)
    i = measured_current(electrode, E, scan_rate)
    rng = np.random.default_rng(seed)
    sigma = noise * np.abs(i) + floor
    return E, i + sigma * rng.normal(size=len(E)), sigma


def zero_crossing(E, i):
    k = np.flatnonzero(np.diff(np.sign(i)) > 0)
    if len(k) == 0:
        raise ValueError("measured current never changes sign")
    k = k[np.argmin(np.abs(k - len(E) // 2))]
    return float(E[k] - i[k] * (E[k + 1] - E[k]) / (i[k + 1] - i[k]))


@dataclass
class TafelExtrapolation:
    ecorr: float
    icorr: float
    beta_a: float
    beta_c: float
    anodic_line: tuple = field(repr=False)
    cathodic_line: tuple = field(repr=False)


def tafel_extrapolation(E, i, window=(0.05, 0.15)):
    """Fit straight lines to log|i| on both branches and intersect them.

    Window is the overpotential range from the measured zero-current potential,
    as commonly recommended (at least 50 mV away from Ecorr).
    """
    e0 = zero_crossing(E, i)
    lo, hi = window
    if not 0 < lo < hi:
        raise ValueError("window must satisfy 0 < low < high")
    a = (E - e0 >= lo) & (E - e0 <= hi) & (i > 0)
    c = (e0 - E >= lo) & (e0 - E <= hi) & (i < 0)
    if a.sum() < 5 or c.sum() < 5:
        raise ValueError("not enough points in the Tafel windows")
    ma, ba = np.polyfit(E[a], np.log10(i[a]), 1)
    mc, bc = np.polyfit(E[c], np.log10(-i[c]), 1)
    if ma <= 0 or mc >= 0:
        raise ValueError("Tafel lines have the wrong sign of slope")
    ec = (bc - ba) / (ma - mc)
    return TafelExtrapolation(ec, 10 ** (ma * ec + ba), 1 / ma, -1 / mc, (ma, ba), (mc, bc))


def stern_geary_B(beta_a, beta_c):
    """B = βa βc / (2.303 (βa + βc)); βc = inf for a diffusion-limited cathode."""
    if math.isinf(beta_c):
        return beta_a / LN10
    return beta_a * beta_c / (LN10 * (beta_a + beta_c))


def polarisation_resistance(E, i, half_width=0.010):
    """Rp = dE/di from a straight-line fit within ±half_width of Ecorr."""
    e0 = zero_crossing(E, i)
    m = np.abs(E - e0) <= half_width
    if m.sum() < 5:
        raise ValueError("not enough points inside the LPR window")
    slope, _ = np.polyfit(i[m], E[m], 1)
    return float(slope), e0


@dataclass
class MixedPotentialFit:
    model: str
    parameters: dict
    standard_errors: dict
    weighted_rss: float
    aicc: float
    n: int
    jacobian_condition: float
    predicted: np.ndarray = field(repr=False)


_NAMES = {"activation": ["Ecorr_V", "log10_icorr", "beta_a_V", "beta_c_V", "Ru_ohm_cm2"],
          "diffusion": ["Ecorr_V", "log10_icorr", "beta_a_V", "beta_c_V", "log10_ilim_over_icorr", "Ru_ohm_cm2"]}


def mixed_potential_current(theta, E_app, model):
    """Net current for the two-reaction model in Ecorr/icorr form, with IR drop."""
    if model == "activation":
        ec, lic, ba, bc, ru = theta
        ilim = math.inf
    elif model == "diffusion":
        ec, lic, ba, bc, lratio, ru = theta
        ilim = 10 ** (lic + lratio)
    else:
        raise ValueError("model must be 'activation' or 'diffusion'")
    icorr = 10 ** lic
    if ilim <= icorr:
        return np.full(np.shape(E_app), np.nan)
    ik0 = icorr if math.isinf(ilim) else icorr * ilim / (ilim - icorr)
    anodic = Tafel(icorr, ec, ba, True)
    cathodic = Tafel(ik0, ec, bc, False, ilim)
    return measured_current(Electrode((anodic, cathodic), r_u=max(ru, 0.0), c_dl=0.0), E_app)


def fit_mixed_potential(E, i, sigma, model="diffusion", start=None):
    """Weighted least squares on the full curve. Returns a MixedPotentialFit."""
    names = _NAMES[model]
    e0 = zero_crossing(E, i)
    if start is None:
        t = tafel_extrapolation(E, i)
        guess = [e0, math.log10(t.icorr), min(max(t.beta_a, 0.02), 0.3), min(max(t.beta_c, 0.02), 0.5)]
        if model == "diffusion":
            guess.append(math.log10(max(np.max(-i) / t.icorr, 1.5)))
        guess.append(0.0)
    else:
        guess = list(start)
    # The limiting current is fitted as log10(i_lim / i_corr) > 0, which keeps
    # every trial point physical (the cathodic plateau must exceed i_corr).
    lo = [e0 - 0.1, -9, 0.015, 0.015] + ([0.005] if model == "diffusion" else []) + [0.0]
    hi = [e0 + 0.1, -1, 0.4, 0.6] + ([6.0] if model == "diffusion" else []) + [5e3]
    guess = np.clip(guess, np.array(lo) + 1e-9, np.array(hi) - 1e-9)

    def resid(th):
        r = (mixed_potential_current(th, E, model) - i) / sigma
        return np.where(np.isfinite(r), r, 1e6)

    best = None
    for ru0 in (0.0, 30.0, 100.0, 300.0, 1000.0, 3000.0):
        g = guess.copy(); g[-1] = ru0
        sol = least_squares(resid, g, bounds=(lo, hi), x_scale="jac", method="trf", max_nfev=400)
        if best is None or sol.cost < best.cost:
            best = sol
    J = best.jac
    s = np.linalg.svd(J, compute_uv=False)
    cov = np.linalg.pinv(J.T @ J)
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    n, p = len(E), len(names)
    wrss = float(2 * best.cost)
    aicc = wrss + 2 * p + 2 * p * (p + 1) / (n - p - 1)   # known noise variance
    return MixedPotentialFit(model, dict(zip(names, map(float, best.x))), dict(zip(names, map(float, se))),
                             wrss, aicc, n, float(s[0] / s[-1]) if s[-1] > 0 else math.inf,
                             mixed_potential_current(best.x, E, model))


def corrosion_rate_mm_per_year(icorr_A_cm2, equivalent_weight, density_g_cm3):
    """ASTM G102: CR = K1 icorr EW / rho, with icorr in uA/cm2."""
    if icorr_A_cm2 < 0 or equivalent_weight <= 0 or density_g_cm3 <= 0:
        raise ValueError("icorr must be non-negative; EW and density positive")
    return K1_MM_G_PER_UA_CM_Y * icorr_A_cm2 * 1e6 * equivalent_weight / density_g_cm3


CARBON_STEEL = {"equivalent_weight": 27.92, "density_g_cm3": 7.87}   # Fe -> Fe2+
