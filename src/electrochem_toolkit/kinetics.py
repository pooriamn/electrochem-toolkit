"""Tafel analysis of an electrocatalytic sweep: what the slope depends on.

A cathodic linear sweep for the hydrogen evolution reaction is simulated from
Butler-Volmer kinetics, with uncompensated resistance, double-layer charging
and measurement noise. The functions here then estimate the Tafel slope and
exchange current density in the ways people commonly do, so the effect of each
choice can be measured against a known answer.

Sign convention: potentials vs RHE in V, so the HER equilibrium potential is 0
and cathodic overpotential is negative. Current density in A/cm2, cathodic
negative. Resistances in ohm cm2. Tafel slopes are reported as positive mV/dec.
"""
from dataclasses import dataclass, field
import math

import numpy as np
from scipy.optimize import least_squares

F = 96485.33212
R = 8.314462618
LN10 = math.log(10)


def butler_volmer(eta, j0, alpha_c, alpha_a=None, temperature_K=298.15):
    """Net current density; alpha_a defaults to 1 - alpha_c (single-step)."""
    if j0 <= 0 or not 0 < alpha_c < 2:
        raise ValueError("j0 must be positive and alpha_c in (0, 2)")
    aa = 1 - alpha_c if alpha_a is None else alpha_a
    f = F / (R * temperature_K)
    eta = np.asarray(eta, float)
    return j0 * (np.exp(aa * f * eta) - np.exp(-alpha_c * f * eta))


def tafel_slope_mV(alpha_c, temperature_K=298.15):
    return 1e3 * LN10 * R * temperature_K / (alpha_c * F)


@dataclass(frozen=True)
class HerElectrode:
    j0: float = 5e-6               # A/cm2
    alpha_c: float = 0.5
    r_u: float = 2.0               # ohm cm2
    c_dl: float = 1e-3             # F/cm2, rough electrode
    temperature_K: float = 298.15

    def kinetic(self, eta):
        return butler_volmer(eta, self.j0, self.alpha_c, temperature_K=self.temperature_K)


def _solve_interface(E_app, kinetic, r_u):
    """E_int with E_app = E_int + j(E_int) R_u, by vectorised bisection."""
    E_app = np.asarray(E_app, float)
    if r_u == 0:
        return E_app
    lo, hi = E_app - 3.0, E_app + 3.0
    for _ in range(90):
        mid = 0.5 * (lo + hi)
        up = mid + kinetic(mid) * r_u - E_app > 0
        hi = np.where(up, mid, hi)
        lo = np.where(up, lo, mid)
    return 0.5 * (lo + hi)


def sweep(el, E_start=0.05, E_end=-0.70, scan_rate=0.005, step=0.001, compensated=0.0,
          noise=0.01, floor=2e-7, seed=0):
    """Simulated cathodic LSV.

    compensated is the fraction of R_u removed by the potentiostat (positive
    feedback), so the data carry (1 - compensated) R_u. Scan rate in V/s; the
    sweep runs negative, so charging current is cathodic.
    Returns (E_applied, j_measured, sigma).
    """
    if not 0 <= compensated < 1:
        raise ValueError("compensated must be in [0, 1)")
    E = np.arange(E_start, E_end - step / 2, -step)
    r_eff = (1 - compensated) * el.r_u
    j = el.kinetic(_solve_interface(E, el.kinetic, r_eff)) - el.c_dl * scan_rate
    rng = np.random.default_rng(seed)
    sigma = noise * np.abs(j) + floor
    return E, j + sigma * rng.normal(size=len(E)), sigma


def ir_correct(E, j, r_u):
    return np.asarray(E, float) - np.asarray(j, float) * r_u


@dataclass
class LineFit:
    slope_mV: float
    j0: float
    r_squared: float
    n: int


def fit_tafel_line(eta, j, j_range=None, eta_range=None):
    """Least-squares line through eta against log10|j| on the cathodic branch.

    Select the region by current (j_range, A/cm2 magnitudes) or by overpotential
    magnitude (eta_range, V). Returns slope in mV/dec and j0 from the intercept.
    """
    eta, j = np.asarray(eta, float), np.asarray(j, float)
    m = j < 0
    if j_range is not None:
        m &= (-j >= j_range[0]) & (-j <= j_range[1])
    if eta_range is not None:
        m &= (-eta >= eta_range[0]) & (-eta <= eta_range[1])
    if m.sum() < 5:
        raise ValueError("fewer than five points in the selected region")
    x, y = np.log10(-j[m]), -eta[m]
    b, a = np.polyfit(x, y, 1)
    resid = y - (a + b * x)
    r2 = 1 - resid.var() / y.var() if y.var() > 0 else math.nan
    return LineFit(1e3 * b, 10 ** (-a / b), float(r2), int(m.sum()))


def local_tafel_slope(eta, j, width=0.020, min_points=8):
    """Tafel slope plot: slope over a sliding potential window of fixed width.

    Returns (mean log10|j|, slope in mV/dec) for each window centre. This is the
    construction recommended by van der Heijden et al., ACS Energy Lett. 2024.
    """
    eta, j = np.asarray(eta, float), np.asarray(j, float)
    keep = j < 0
    eta, j = eta[keep], j[keep]
    out_x, out_s = [], []
    for c in eta:
        m = np.abs(eta - c) <= width / 2
        if m.sum() < min_points:
            continue
        x = np.log10(-j[m])
        if np.ptp(x) == 0:
            continue
        b = np.polyfit(x, -eta[m], 1)[0]
        out_x.append(x.mean()); out_s.append(1e3 * b)
    return np.array(out_x), np.array(out_s)


def plateau_slope(log_j, slope, tolerance=0.05, min_decades=0.5):
    """Median slope over the longest run where the local slope stays within
    ±tolerance of its own running median. Returns (slope, log_j_low, log_j_high)."""
    order = np.argsort(log_j)
    x, s = log_j[order], slope[order]
    best = (math.nan, math.nan, math.nan, 0.0)
    for a in range(len(x)):
        for b in range(a + 5, len(x) + 1):
            med = np.median(s[a:b])
            if np.any(np.abs(s[a:b] - med) > tolerance * med):
                break
            span = x[b - 1] - x[a]
            if span > best[3]:
                best = (float(med), float(x[a]), float(x[b - 1]), float(span))
    if best[3] < min_decades:
        raise ValueError("no Tafel plateau spanning the required decades")
    return best[:3]


@dataclass
class BVFit:
    parameters: dict
    standard_errors: dict
    weighted_rss: float
    predicted: np.ndarray = field(repr=False)


def fit_butler_volmer(E_app, j, sigma, temperature_K=298.15, j_cap_known=None):
    """Fit log10 j0, alpha_c, R_u and a constant charging current to the sweep.

    Every point is used, including the region near equilibrium where the back
    reaction matters and the high-current region dominated by the ohmic drop.
    """
    E_app, j, sigma = (np.asarray(v, float) for v in (E_app, j, sigma))

    def model(th):
        lj0, ac, ru, jc = th
        kin = lambda e: butler_volmer(e, 10 ** lj0, ac, temperature_K=temperature_K)
        jcap = jc if j_cap_known is None else j_cap_known
        return kin(_solve_interface(E_app, kin, ru)) + jcap

    def resid(th):
        return (model(th) - j) / sigma

    start = [-5.0, 0.5, 1.0, 0.0]
    lo, hi = [-10, 0.1, 0.0, -1e-3], [-1, 1.5, 50.0, 1e-3]
    best = None
    for lj0 in (-7, -5, -3):
        for ru in (0.2, 2.0, 10.0):
            sol = least_squares(resid, [lj0, 0.5, ru, 0.0], bounds=(lo, hi), x_scale="jac", max_nfev=300)
            if best is None or sol.cost < best.cost:
                best = sol
    J = best.jac
    cov = np.linalg.pinv(J.T @ J)
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    names = ["log10_j0", "alpha_c", "Ru_ohm_cm2", "j_capacitive_A_cm2"]
    p = dict(zip(names, map(float, best.x)))
    s = dict(zip(names, map(float, se)))
    p["tafel_slope_mV"] = tafel_slope_mV(p["alpha_c"], temperature_K)
    s["tafel_slope_mV"] = p["tafel_slope_mV"] * s["alpha_c"] / p["alpha_c"]
    return BVFit(p, s, float(2 * best.cost), model(best.x))
