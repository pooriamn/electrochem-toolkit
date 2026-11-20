# Vanadium flow battery: flow rate, pumping and capacity fade

Two questions come up whenever a flow cell goes on test. How much electrolyte should we pump? And why has the capacity dropped by cycle 50, and will remixing the tanks bring it back? This study answers both for a lumped single-cell model of an all-vanadium battery, with the physics kept simple enough that every term can be checked by hand.

## Model

The cell has 100 cm² of carbon felt on each side and two 100 mL tanks of 1.6 M vanadium. Four species are tracked in moles: V(II) and V(III) in the negative tank, V(IV) and V(V) in the positive tank. Each tank and its half-cell are treated as one well-mixed volume.

**Cell voltage** at current density *i* (positive on charge, *s* = ±1):

```text
E_ocv = E0' + (RT/F) ln( c_V2 c_V5 / (c_V3 c_V4) )
η_act = (2RT/F) [ asinh(i / 2i0,−) + asinh(i / 2i0,+) ],   i0 = F k aL √(c_ox c_red)
η_mt  = −(RT/F) Σ ln(1 − i / i_lim),   i_lim = F k_m aL c̄,   k_m = b u^0.4
V     = E_ocv + s (η_act + η_mt + i ASR)
```

E0' = 1.40 V is a formal cell potential that absorbs the proton and activity terms, so it holds only for a fixed supporting-acid concentration. *aL* is the felt surface area per unit geometric area. For the mass-transfer term, c̄ is the mean of the inlet and outlet reactant concentrations, with the outlet set by single-pass conversion, c_out = c_in − I/(FQ). If c_out reaches zero the electrolyte simply can't carry enough reactant through the cell, and the half-cycle ends. At low flow rates that is what actually limits capacity, more than the kinetics.

**Pumping.** Pressure drop is Darcy flow through the felt, ΔP = μuL/K. Pump power for both half-cells is 2ΔPQ/η_pump. Pipes, manifolds and shunt currents are left out.

**Crossover.** Each ion diffuses through a 125 µm membrane at the rate D_i A c_i / δ, and reacts as soon as it arrives on the other side:

```text
positive side:  V²⁺ + 2VO₂⁺ → 3VO²⁺        V³⁺ + VO₂⁺ → 2VO²⁺
negative side:  VO²⁺ + V²⁺ → 2V³⁺           VO₂⁺ + 2V²⁺ → 3V³⁺
```

This is the mass-balance scheme of Tang, Bao and Skyllas-Kazacos [1]. It conserves both vanadium and electrons, so crossover alone moves vanadium between the tanks without changing the average oxidation state. I chose diffusion coefficients that put V(II) and V(IV) through the membrane faster than V(III) and V(V). That gives a net drift of vanadium to the positive side, which is the direction usually reported for cation-exchange membranes. The values are illustrative, not fitted to a particular membrane.

**Hydrogen evolution** is optional and takes a fixed share of the negative-electrode charging current. It is the one process here that changes the average oxidation state.

## Numerical approach

- Each half-cycle runs at constant current, integrated with `solve_ivp` (rtol 1e-8). It terminates on a voltage event at 1.65 V (charge) or 1.00 V (discharge), or when the single-pass supply runs out.
- Efficiencies come from the third cycle, once the start-up transient has gone. Coulombic efficiency is discharge over charge Ah, energy efficiency is discharge over charge Wh, and voltage efficiency is their ratio. System efficiency charges pump energy against each half-cycle.
- Remixing combines both tanks, splits the vanadium equally and returns it to a discharged V(III)/V(IV) state. Any net shift in oxidation state stays behind as V(V) in the positive tank.

## Results

![Vanadium flow cell study figure](results.svg)

At 100 mA cm⁻² and 120 mL/min per half-cell, the settled cycle gives 96.0 % coulombic, 83.5 % voltage and 80.2 % energy efficiency. Once pumping is counted, energy efficiency is 78.5 %.

**Flow rate (b).** At the lowest flow, around 30 mL/min, the charge stops before the voltage limit because a single pass can no longer deliver the reactant the current needs. Between there and about 100 mL/min, mass-transfer overpotential at high SOC pushes the cell to its voltage limit early, so utilisation climbs steadily with flow. Above about 250 mL/min the pump losses grow quickly, since they scale with u². Taking net energy per cycle (discharge energy minus all pump energy) as the measure, the best flow is close to 115 mL/min for this geometry. Stack energy efficiency, the dashed line, hardly moves across the whole range. That's the trap: judged on stack efficiency alone, there's no reason not to pump harder.

**Capacity fade (c, d).** With crossover alone, capacity falls by 0.11 % per cycle over the first 80 cycles, and the negative tank ends up holding 46.9 % of the vanadium. Remixing restores the first-cycle capacity exactly, because nothing was lost; it was only in the wrong tank.

Add 0.2 % hydrogen evolution and the story changes. The fade before remixing is actually slower, at 0.06 % per cycle, because the oxidation-state drift (average 3.57 by cycle 80) partly offsets the vanadium imbalance. Remixing then removes the vanadium imbalance but not the drift, and capacity drops further, from 3.45 to 3.26 Ah. A cell like this needs an oxidation-state correction, chemical or electrochemical, as well as a remix. That's worth knowing before deciding remixing "didn't work".

## Reproduce

```sh
python -m pip install ".[advanced]"
python studies/flow-battery/run.py
```

The script writes `summary.json`, `flow_sweep.csv` and the figure to `study-output/flow-battery/`. The 2 × 120-cycle runs take a few seconds.

## Where this model stops

The tanks and half-cells are treated as one well-mixed volume, so there are no concentration gradients along the felt and no spatial current distribution. Electro-osmotic and osmotic water transfer are excluded, and in real cells they are often as important as diffusion for long-term imbalance. Migration of vanadium under the cell field, temperature effects on diffusion, air oxidation of V(II) and V(V) precipitation are also excluded. The mass-transfer coefficient *b* is set so transport matters at high SOC and low flow; it is not a correlation for any particular felt. Use the model to reason about trends and to design an experiment, not to predict the cycle life of a real stack.

## References

1. A. Tang, J. Bao and M. Skyllas-Kazacos, "Dynamic modelling of the effects of ion diffusion and side reactions on the capacity loss for vanadium redox flow battery", *J. Power Sources* 196 (2011) 10737–10747. [doi:10.1016/j.jpowsour.2011.09.003](https://doi.org/10.1016/j.jpowsour.2011.09.003)

[Implementation](../../src/electrochem_toolkit/flowbattery.py) · [Verification tests](../../tests/test_flowbattery.py)
