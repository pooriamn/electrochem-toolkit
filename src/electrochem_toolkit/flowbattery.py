"""Lumped redox flow battery model: one cell, two tanks, any two redox couples.

The state is four amounts of active species, in this order:

    0  negative electrolyte, charged form      (V(II) for all-vanadium)
    1  negative electrolyte, discharged form   (V(III))
    2  positive electrolyte, discharged form   (V(IV))
    3  positive electrolyte, charged form      (V(V))

What is specific to a chemistry (formal potential, rate constants, membrane
permeabilities and what happens to an ion that crosses the membrane) lives in a
`Chemistry`. Geometry, resistance, mass transfer and pumping live in `FlowCell`.
`VANADIUM` is provided as a worked example; other chemistries are described by
building a `Chemistry` of their own.

Charge transfer runs at constant current. Membrane diffusion and its follow-up
reactions are integrated with SciPy alongside it. Cell voltage combines a
Nernst term, symmetric Butler-Volmer activation at each electrode, a
mass-transfer term that depends on electrolyte velocity and single-pass
conversion, and an area-specific resistance. Optional hydrogen evolution takes
a fixed share of the negative-electrode charging current.

All quantities are SI unless the name says otherwise.
"""
from dataclasses import dataclass, replace
import math

import numpy as np
from scipy.integrate import solve_ivp

F = 96485.33212
R = 8.314462618


@dataclass(frozen=True)
class Chemistry:
    """The electrochemistry of one flow battery, independent of cell hardware.

    crossover[j] is the change in the four species when one mole of species j
    diffuses through the membrane. For a chemistry where both electrolytes
    share one active element (all-vanadium, for example) the crossing ion
    reacts with the other side and nothing is lost overall. For two different
    elements, use `lost_crossover()`: the ion leaves the active inventory.

    oxidation_states, when given, enables `remix()` and
    `average_oxidation_state()`; it only makes sense for a shared element.
    """
    name: str
    formal_potential_V: float             # cell potential at 50 % SOC, activity terms lumped
    k_neg_m_s: float
    k_pos_m_s: float
    diffusivity_m2_s: tuple               # membrane diffusivity of species 0-3
    crossover: tuple                      # 4 x 4, see above
    oxidation_states: tuple = None
    labels: tuple = ("neg charged", "neg discharged", "pos discharged", "pos charged")

    def __post_init__(self):
        for name in ("formal_potential_V", "k_neg_m_s", "k_pos_m_s"):
            value = getattr(self, name)
            if not (math.isfinite(value) and value > 0):
                raise ValueError(f"{name} must be positive and finite")
        if len(self.diffusivity_m2_s) != 4 or any(d < 0 for d in self.diffusivity_m2_s):
            raise ValueError("diffusivity_m2_s needs four non-negative values")
        m = np.asarray(self.crossover, dtype=float)
        if m.shape != (4, 4):
            raise ValueError("crossover must be a 4 x 4 table")
        for j in range(4):
            if m[j][j] >= 0:
                raise ValueError("a crossing species must leave its own side")
        if self.oxidation_states is not None and len(self.oxidation_states) != 4:
            raise ValueError("oxidation_states needs four values")


def lost_crossover():
    """Crossover table for two different active elements: crossing ions are lost."""
    return tuple(tuple(-1.0 if i == j else 0.0 for i in range(4)) for j in range(4))


# All-vanadium: a crossing ion reacts at once with the other electrolyte
# (Tang, Bao and Skyllas-Kazacos, J. Power Sources 196 (2011) 10737):
#   positive side  V2 + 2 V5 -> 3 V4        V3 + V5 -> 2 V4
#   negative side  V4 + V2   -> 2 V3        V5 + 2 V2 -> 3 V3
# Diffusivities are illustrative; V(II) and V(IV) cross fastest, which gives the
# net drift to the positive side usually reported for cation-exchange membranes.
VANADIUM = Chemistry(
    name="all-vanadium",
    formal_potential_V=1.40,
    k_neg_m_s=1.7e-7,
    k_pos_m_s=6.8e-7,
    diffusivity_m2_s=(1.2e-11, 0.6e-11, 1.0e-11, 0.5e-11),
    crossover=((-1, 0, 3, -2), (0, -1, 2, -1), (-1, 2, -1, 0), (-2, 3, 0, -1)),
    oxidation_states=(2, 3, 4, 5),
    labels=("V(II)", "V(III)", "V(IV)", "V(V)"),
)


@dataclass(frozen=True)
class FlowCell:
    chemistry: Chemistry = VANADIUM
    area_m2: float = 0.01                 # 100 cm2 active area
    electrode_height_m: float = 0.10      # flow length through the felt
    electrode_thickness_m: float = 4e-3
    specific_area_m2_m3: float = 3.5e4    # felt surface per electrode volume
    mass_transfer_b: float = 2.5e-5       # k_m = b * u**0.4, u in m/s
    asr_ohm_m2: float = 1.0e-4            # 1.0 ohm cm2
    membrane_thickness_m: float = 125e-6
    viscosity_Pa_s: float = 5e-3
    permeability_m2: float = 1e-10
    pump_efficiency: float = 0.7
    temperature_K: float = 298.15

    def __post_init__(self):
        positive = ["area_m2", "electrode_height_m", "electrode_thickness_m",
                    "specific_area_m2_m3", "mass_transfer_b", "membrane_thickness_m", "viscosity_Pa_s",
                    "permeability_m2", "pump_efficiency", "temperature_K"]
        for name in positive:
            value = getattr(self, name)
            if not (math.isfinite(value) and value > 0):
                raise ValueError(f"{name} must be positive and finite")
        if self.asr_ohm_m2 < 0 or not math.isfinite(self.asr_ohm_m2):
            raise ValueError("asr_ohm_m2 must be non-negative")
        if self.pump_efficiency > 1:
            raise ValueError("pump_efficiency cannot exceed 1")

    @property
    def flow_section_m2(self):
        width = self.area_m2 / self.electrode_height_m
        return width * self.electrode_thickness_m

    @property
    def area_ratio(self):
        """Electrochemically active area per unit geometric area (a L)."""
        return self.specific_area_m2_m3 * self.electrode_thickness_m

    def flow_rate_m3_s(self, velocity_m_s):
        return velocity_m_s * self.flow_section_m2

    def pressure_drop_Pa(self, velocity_m_s):
        """Darcy flow through the felt only; manifolds and pipes are excluded."""
        return self.viscosity_Pa_s * velocity_m_s * self.electrode_height_m / self.permeability_m2

    def pump_power_W(self, velocity_m_s):
        """Hydraulic power for both half-cells divided by pump efficiency."""
        return 2 * self.pressure_drop_Pa(velocity_m_s) * self.flow_rate_m3_s(velocity_m_s) / self.pump_efficiency


@dataclass(frozen=True)
class Electrolyte:
    volume_m3: float = 1e-4               # 100 mL per side
    active_mol_m3: float = 1600.0         # total active species per side, 1.6 M
    initial_soc: float = 0.02

    def __post_init__(self):
        if not (self.volume_m3 > 0 and self.active_mol_m3 > 0):
            raise ValueError("volume and concentration must be positive")
        if not 0 < self.initial_soc < 1:
            raise ValueError("initial_soc must lie strictly between 0 and 1")

    def initial_state(self):
        n = self.volume_m3 * self.active_mol_m3
        s = self.initial_soc
        return np.array([s * n, (1 - s) * n, (1 - s) * n, s * n])


def crossover_rates(cell, state, volume_m3):
    """Rate of change of the four species from membrane crossover (mol/s)."""
    chem = cell.chemistry
    c = np.clip(state, 0, None) / volume_m3
    flux = np.array([d * cell.area_m2 * ci / cell.membrane_thickness_m
                     for d, ci in zip(chem.diffusivity_m2_s, c)])
    return flux @ np.asarray(chem.crossover, dtype=float)


def concentrations(state, volume_m3):
    return np.asarray(state, dtype=float) / volume_m3


def soc(state):
    """State of charge of the negative and positive electrolytes."""
    a, b, c, d = state
    return a / (a + b), d / (c + d)


def average_oxidation_state(state, chemistry=None):
    """Mean oxidation state of the active element; shared-element chemistries only."""
    ox = _oxidation_states(chemistry or VANADIUM)
    n = np.asarray(state, dtype=float)
    return float(np.dot(ox, n) / n.sum())


def _oxidation_states(chemistry):
    if chemistry.oxidation_states is None:
        raise ValueError(f"{chemistry.name}: oxidation states are only defined for a shared active element")
    return np.asarray(chemistry.oxidation_states, dtype=float)


def cell_voltage(cell, state, volume_m3, current_A, velocity_m_s):
    """Return (voltage, breakdown). Positive current charges the cell.

    Returns nan when the requested current exceeds what the electrolyte can
    supply in a single pass, which the integrator treats as a cut-off.
    """
    c2, c3, c4, c5 = concentrations(state, volume_m3)   # neg charged, neg discharged, pos discharged, pos charged
    if min(c2, c3, c4, c5) <= 0:
        return math.nan, {}
    rt_f = R * cell.temperature_K / F
    chem = cell.chemistry
    ocv = chem.formal_potential_V + rt_f * math.log((c2 * c5) / (c3 * c4))
    if current_A == 0:
        return ocv, {"ocv": ocv, "activation": 0.0, "mass_transfer": 0.0, "ohmic": 0.0}
    i = abs(current_A) / cell.area_m2
    al = cell.area_ratio
    i0_neg = F * chem.k_neg_m_s * al * math.sqrt(c2 * c3)
    i0_pos = F * chem.k_pos_m_s * al * math.sqrt(c4 * c5)
    eta_act = 2 * rt_f * (math.asinh(i / (2 * i0_neg)) + math.asinh(i / (2 * i0_pos)))
    # Charging consumes the discharged forms; discharging consumes the charged forms.
    reactants = (c3, c4) if current_A > 0 else (c2, c5)
    q = cell.flow_rate_m3_s(velocity_m_s)
    km = cell.mass_transfer_b * velocity_m_s ** 0.4
    eta_mt = 0.0
    for c_in in reactants:
        c_out = c_in - abs(current_A) / (F * q)
        if c_out <= 0:
            return math.nan, {}
        c_mean = 0.5 * (c_in + c_out)
        i_lim = F * km * al * c_mean
        if i >= i_lim:
            return math.nan, {}
        eta_mt += -rt_f * math.log(1 - i / i_lim)
    eta_ohm = i * cell.asr_ohm_m2
    sign = 1 if current_A > 0 else -1
    v = ocv + sign * (eta_act + eta_mt + eta_ohm)
    return v, {"ocv": ocv, "activation": eta_act, "mass_transfer": eta_mt, "ohmic": eta_ohm}


def _rhs(cell, volume_m3, current_A, h2_fraction):
    def f(_, y):
        dy = crossover_rates(cell, y, volume_m3)
        faradaic = current_A / F
        to_v2 = faradaic * (1 - h2_fraction) if current_A > 0 else faradaic
        dy += np.array([to_v2, -to_v2, -faradaic, faradaic])
        return dy
    return f


def half_cycle(cell, electrolyte, state, current_A, velocity_m_s, v_limit,
               h2_fraction=0.0, sample_s=20.0, max_time_s=5e4):
    """Run at constant current until the voltage reaches v_limit.

    Returns a dict of time series and the final state. The half-cycle also stops
    if any concentration is exhausted or the single-pass supply is exceeded.
    """
    vol = electrolyte.volume_m3
    charging = current_A > 0

    def limit(_, y):
        v, _ = cell_voltage(cell, y, vol, current_A, velocity_m_s)
        if not math.isfinite(v):
            return -1.0
        return (v_limit - v) if charging else (v - v_limit)
    limit.terminal = True
    limit.direction = -1

    v0, _ = cell_voltage(cell, state, vol, current_A, velocity_m_s)
    if not math.isfinite(v0) or limit(0, state) <= 0:
        return {"t": np.array([0.0]), "V": np.array([v0]), "state": np.asarray(state, float),
                "charge_Ah": 0.0, "energy_Wh": 0.0, "duration_s": 0.0}
    sol = solve_ivp(_rhs(cell, vol, current_A, h2_fraction), (0, max_time_s), state,
                    events=limit, rtol=1e-8, atol=1e-10, max_step=sample_s, dense_output=True)
    t_end = sol.t[-1]
    t = np.append(np.arange(0, t_end, sample_s), t_end)
    y = sol.sol(t)
    v = np.array([cell_voltage(cell, y[:, k], vol, current_A, velocity_m_s)[0] for k in range(len(t))])
    v[-1] = v_limit if not math.isfinite(v[-1]) else v[-1]
    energy = _trapz(np.abs(current_A) * v, t) / 3600 if len(t) > 1 else 0.0
    return {"t": t, "V": v, "state": sol.y[:, -1], "charge_Ah": abs(current_A) * t_end / 3600,
            "energy_Wh": float(energy), "duration_s": float(t_end)}


def cycle(cell, electrolyte, state, current_A, velocity_m_s, v_max=1.65, v_min=1.00,
          h2_fraction=0.0, sample_s=20.0):
    """One charge followed by one discharge. Returns (charge, discharge, metrics)."""
    if current_A <= 0:
        raise ValueError("current_A is the charging current and must be positive")
    if not v_min < v_max:
        raise ValueError("v_min must be below v_max")
    ch = half_cycle(cell, electrolyte, state, current_A, velocity_m_s, v_max, h2_fraction, sample_s)
    dis = half_cycle(cell, electrolyte, ch["state"], -current_A, velocity_m_s, v_min, h2_fraction, sample_s)
    pump_W = cell.pump_power_W(velocity_m_s)
    pump_ch, pump_dis = pump_W * ch["duration_s"] / 3600, pump_W * dis["duration_s"] / 3600
    ce = dis["charge_Ah"] / ch["charge_Ah"] if ch["charge_Ah"] else math.nan
    ee = dis["energy_Wh"] / ch["energy_Wh"] if ch["energy_Wh"] else math.nan
    metrics = {
        "charge_capacity_Ah": ch["charge_Ah"],
        "discharge_capacity_Ah": dis["charge_Ah"],
        "coulombic_efficiency": ce,
        "voltage_efficiency": ee / ce if ce else math.nan,
        "energy_efficiency": ee,
        "pump_energy_Wh": pump_ch + pump_dis,
        "system_energy_efficiency": (dis["energy_Wh"] - pump_dis) / (ch["energy_Wh"] + pump_ch)
                                    if ch["energy_Wh"] else math.nan,
    }
    return ch, dis, metrics


def remix(state, chemistry=None):
    """Mix both tanks, split the inventory equally and restore a discharged state.

    Only meaningful when both electrolytes share one active element. Crossover
    then only moves material between tanks, so remixing recovers it. A net
    change in average oxidation state (from hydrogen evolution, for example)
    stays behind as charged positive species, or as charged negative species if
    the shift is the other way, and keeps reducing usable capacity.
    """
    ox = _oxidation_states(chemistry or VANADIUM)
    n = np.asarray(state, dtype=float)
    total = n.sum()
    mid = 0.5 * (ox[1] + ox[2])                       # fully discharged and balanced
    excess = float(np.dot(ox, n)) - mid * total       # electrons relative to that state
    half = total / 2
    if excess >= 0:
        if excess >= half:
            raise ValueError("electrolyte is too oxidised to rebalance by remixing")
        return np.array([0.0, half, half - excess, excess]) + _floor(total)
    if -excess >= half:
        raise ValueError("electrolyte is too reduced to rebalance by remixing")
    return np.array([-excess, half + excess, half, 0.0]) + _floor(total)


def _floor(total):
    # Keep every species slightly positive so the Nernst term stays defined.
    eps = 1e-4 * total
    return np.array([eps, -eps, -eps, eps])


def _trapz(y, x):
    y, x = np.asarray(y, float), np.asarray(x, float)
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))
