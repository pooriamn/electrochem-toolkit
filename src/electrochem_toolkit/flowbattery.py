"""Lumped all-vanadium redox flow battery (VRFB) model.

One cell, two well-mixed tanks. Four vanadium species are tracked as moles:

    negative tank   n2 = V(II),  n3 = V(III)
    positive tank   n4 = V(IV),  n5 = V(V)

Charge transfer at constant current, diffusion of each species through the
membrane and the self-discharge reactions that follow (Tang, Bao and
Skyllas-Kazacos, J. Power Sources 196 (2011) 10737) are integrated with SciPy.
Cell voltage combines a Nernst term, symmetric Butler-Volmer activation at each
electrode, a mass-transfer term that depends on electrolyte velocity and
single-pass conversion, and an area-specific resistance.

Optional hydrogen evolution takes a fixed fraction of the negative-electrode
charging current. Unlike crossover, it changes the average oxidation state of
the electrolyte and cannot be undone by remixing.

All quantities are SI unless the name says otherwise.
"""
from dataclasses import dataclass, replace
import math

import numpy as np
from scipy.integrate import solve_ivp

F = 96485.33212
R = 8.314462618


@dataclass(frozen=True)
class FlowCell:
    area_m2: float = 0.01                 # 100 cm2 active area
    electrode_height_m: float = 0.10      # flow length through the felt
    electrode_thickness_m: float = 4e-3
    specific_area_m2_m3: float = 3.5e4    # felt surface per electrode volume
    formal_potential_V: float = 1.40      # cell potential at 50 % SOC, acid terms lumped
    k_neg_m_s: float = 1.7e-7             # V(II)/V(III) rate constant
    k_pos_m_s: float = 6.8e-7             # V(IV)/V(V) rate constant
    mass_transfer_b: float = 2.5e-5       # k_m = b * u**0.4, u in m/s
    asr_ohm_m2: float = 1.0e-4            # 1.0 ohm cm2
    membrane_thickness_m: float = 125e-6
    diffusivity_m2_s: tuple = (1.2e-11, 0.6e-11, 1.0e-11, 0.5e-11)  # V2, V3, V4, V5
    viscosity_Pa_s: float = 5e-3
    permeability_m2: float = 1e-10
    pump_efficiency: float = 0.7
    temperature_K: float = 298.15

    def __post_init__(self):
        positive = ["area_m2", "electrode_height_m", "electrode_thickness_m",
                    "specific_area_m2_m3", "formal_potential_V", "k_neg_m_s", "k_pos_m_s",
                    "mass_transfer_b", "membrane_thickness_m", "viscosity_Pa_s",
                    "permeability_m2", "pump_efficiency", "temperature_K"]
        for name in positive:
            value = getattr(self, name)
            if not (math.isfinite(value) and value > 0):
                raise ValueError(f"{name} must be positive and finite")
        if self.asr_ohm_m2 < 0 or not math.isfinite(self.asr_ohm_m2):
            raise ValueError("asr_ohm_m2 must be non-negative")
        if self.pump_efficiency > 1:
            raise ValueError("pump_efficiency cannot exceed 1")
        if len(self.diffusivity_m2_s) != 4 or any(d < 0 for d in self.diffusivity_m2_s):
            raise ValueError("diffusivity_m2_s needs four non-negative values")

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
    vanadium_mol_m3: float = 1600.0       # 1.6 M
    initial_soc: float = 0.02

    def __post_init__(self):
        if not (self.volume_m3 > 0 and self.vanadium_mol_m3 > 0):
            raise ValueError("volume and concentration must be positive")
        if not 0 < self.initial_soc < 1:
            raise ValueError("initial_soc must lie strictly between 0 and 1")

    def initial_state(self):
        n = self.volume_m3 * self.vanadium_mol_m3
        s = self.initial_soc
        return np.array([s * n, (1 - s) * n, (1 - s) * n, s * n])


def crossover_rates(cell, state, volume_m3):
    """Rate of change of the four species from membrane diffusion (mol/s).

    Each crossing ion reacts immediately with the species on the other side:
        V2 + 2 V5 -> 3 V4      V3 + V5 -> 2 V4      (in the positive tank)
        V4 + V2   -> 2 V3      V5 + 2 V2 -> 3 V3    (in the negative tank)
    """
    c = np.clip(state, 0, None) / volume_m3
    j2, j3, j4, j5 = (d * cell.area_m2 * ci / cell.membrane_thickness_m
                      for d, ci in zip(cell.diffusivity_m2_s, c))
    return np.array([
        -j2 - j4 - 2 * j5,          # V2
        -j3 + 2 * j4 + 3 * j5,      # V3
        3 * j2 + 2 * j3 - j4,       # V4
        -2 * j2 - j3 - j5,          # V5
    ])


def concentrations(state, volume_m3):
    return np.asarray(state, dtype=float) / volume_m3


def soc(state):
    n2, n3, n4, n5 = state
    return n2 / (n2 + n3), n5 / (n4 + n5)


def average_oxidation_state(state):
    n = np.asarray(state, dtype=float)
    return float(np.dot([2, 3, 4, 5], n) / n.sum())


def cell_voltage(cell, state, volume_m3, current_A, velocity_m_s):
    """Return (voltage, breakdown). Positive current charges the cell.

    Returns nan when the requested current exceeds what the electrolyte can
    supply in a single pass, which the integrator treats as a cut-off.
    """
    c2, c3, c4, c5 = concentrations(state, volume_m3)
    if min(c2, c3, c4, c5) <= 0:
        return math.nan, {}
    rt_f = R * cell.temperature_K / F
    ocv = cell.formal_potential_V + rt_f * math.log((c2 * c5) / (c3 * c4))
    if current_A == 0:
        return ocv, {"ocv": ocv, "activation": 0.0, "mass_transfer": 0.0, "ohmic": 0.0}
    i = abs(current_A) / cell.area_m2
    al = cell.area_ratio
    i0_neg = F * cell.k_neg_m_s * al * math.sqrt(c2 * c3)
    i0_pos = F * cell.k_pos_m_s * al * math.sqrt(c4 * c5)
    eta_act = 2 * rt_f * (math.asinh(i / (2 * i0_neg)) + math.asinh(i / (2 * i0_pos)))
    # Reactants: charging consumes V3 and V4, discharging consumes V2 and V5.
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


def remix(state):
    """Mix both tanks, split the vanadium equally and restore a discharged state.

    Crossover only moves vanadium between tanks, so remixing recovers it. A net
    change in average oxidation state (from hydrogen evolution) remains as V(V)
    in the positive tank, or as V(II) in the negative tank if the shift is the
    other way, and keeps reducing usable capacity.
    """
    n = np.asarray(state, dtype=float)
    total = n.sum()
    excess = float(np.dot([2, 3, 4, 5], n)) - 3.5 * total   # electrons relative to V3/V4
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
