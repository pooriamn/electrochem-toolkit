"""Flow battery: flow rate, pumping power and crossover-driven capacity fade.

The cell, flow and pumping parts apply to any two-couple flow battery. The
worked example uses the all-vanadium chemistry (fb.VANADIUM); swap in another
fb.Chemistry to run the same study on a different system.
"""
import json
from pathlib import Path
import numpy as np
from electrochem_toolkit import plotstyle
plotstyle.apply()
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from electrochem_toolkit import flowbattery as fb

CURRENT_A = 10.0          # 100 mA cm-2 on 100 cm2
NOMINAL_U = 5e-3          # m/s superficial velocity in the felt
H2_FRACTION = 0.002       # share of negative charging current lost to H2
REMIX_AT = 80
CYCLES = 120


def ml_min(cell, u):
    return cell.flow_rate_m3_s(u) * 6e7


def settled_cycle(cell, el, u, n=3):
    state = el.initial_state()
    for _ in range(n):
        ch, dis, m = fb.cycle(cell, el, state, CURRENT_A, u)
        state = dis["state"]
    return ch, dis, m


def long_run(cell, el, h2):
    state = el.initial_state()
    rows = []
    for k in range(1, CYCLES + 1):
        if k == REMIX_AT + 1:
            state = fb.remix(state, cell.chemistry)
        _, dis, m = fb.cycle(cell, el, state, CURRENT_A, NOMINAL_U, h2_fraction=h2)
        state = dis["state"]
        rows.append({"cycle": k, "capacity_Ah": m["discharge_capacity_Ah"],
                     "coulombic_efficiency": m["coulombic_efficiency"],
                     "negative_share": float(state[:2].sum() / state.sum()),
                     "aos": fb.average_oxidation_state(state, cell.chemistry)})
    return rows


def main():
    out = Path("study-output/flow-battery"); out.mkdir(parents=True, exist_ok=True)
    cell, el = fb.FlowCell(chemistry=fb.VANADIUM), fb.Electrolyte()
    theoretical_Ah = el.volume_m3 * el.active_mol_m3 * fb.F / 3600

    curves = {u: settled_cycle(cell, el, u) for u in (2e-3, NOMINAL_U, 2e-2)}
    sweep_u = np.geomspace(1.2e-3, 3e-2, 15)
    sweep = []
    for u in sweep_u:
        _, dis, m = settled_cycle(cell, el, u)
        sweep.append({"velocity_m_s": float(u), "flow_mL_min": ml_min(cell, u),
                      "utilisation": m["discharge_capacity_Ah"] / theoretical_Ah,
                      "energy_efficiency": m["energy_efficiency"],
                      "system_energy_efficiency": m["system_energy_efficiency"],
                      "net_energy_Wh": dis["energy_Wh"] - m["pump_energy_Wh"]})
    best = max(sweep, key=lambda r: r["net_energy_Wh"])
    runs = {"crossover only": long_run(cell, el, 0.0), "crossover + 0.2 % H2": long_run(cell, el, H2_FRACTION)}

    nominal = curves[NOMINAL_U][2]
    summary = {
        "nominal": {"flow_mL_min": ml_min(cell, NOMINAL_U), **{k: float(v) for k, v in nominal.items()}},
        "best_net_energy_flow_mL_min": best["flow_mL_min"],
        "theoretical_capacity_Ah": theoretical_Ah,
        "chemistry": cell.chemistry.name,
    }
    for name, rows in runs.items():
        before = rows[REMIX_AT - 1]["capacity_Ah"]
        summary[name] = {"cycle1_Ah": rows[0]["capacity_Ah"], f"cycle{REMIX_AT}_Ah": before,
                         "fade_percent_per_cycle": 100 * (1 - before / rows[0]["capacity_Ah"]) / (REMIX_AT - 1),
                         "after_remix_Ah": rows[REMIX_AT]["capacity_Ah"],
                         "aos_before_remix": rows[REMIX_AT - 1]["aos"],
                         "negative_share_before_remix": rows[REMIX_AT - 1]["negative_share"]}
    summary = json.loads(json.dumps(summary, default=float))
    (out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    with (out / "flow_sweep.csv").open("w") as f:
        f.write(",".join(sweep[0]) + "\n")
        for r in sweep: f.write(",".join(f"{v:.6g}" for v in r.values()) + "\n")
    print(json.dumps(summary, indent=2))

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.4), layout="constrained")
    ax = axes[0, 0]
    colours = {2e-3: plotstyle.AMBER, NOMINAL_U: plotstyle.TEAL, 2e-2: plotstyle.SLATE}
    for u, (ch, dis, _) in curves.items():
        c = colours[u]; lab = f"{ml_min(cell, u):.0f} mL/min"
        ax.plot(ch["t"] * CURRENT_A / 3600, ch["V"], color=c, label=lab)
        ax.plot(dis["t"] * CURRENT_A / 3600, dis["V"], color=c, ls="--")
    ax.set(xlabel="charge passed / Ah", ylabel="cell voltage / V", title="Charge (solid) and discharge (dashed)")
    ax.legend(loc="upper left", fontsize=8)

    ax = axes[0, 1]
    q = [r["flow_mL_min"] for r in sweep]
    ax.semilogx(q, [100 * r["utilisation"] for r in sweep], color=plotstyle.SLATE, label="electrolyte utilisation")
    ax.semilogx(q, [100 * r["energy_efficiency"] for r in sweep], color=plotstyle.TEAL, ls="--", label="stack energy efficiency")
    ax.semilogx(q, [100 * r["system_energy_efficiency"] for r in sweep], color=plotstyle.TEAL, label="after pumping")
    ax.axvline(best["flow_mL_min"], color=plotstyle.AMBER, lw=1.2)
    ax.text(best["flow_mL_min"] * 1.08, 34, f"most net energy\n{best['flow_mL_min']:.0f} mL/min",
            fontsize=8.5, color=plotstyle.INK, va="bottom")
    ax.set_xticks([30, 50, 100, 200, 500])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(FuncFormatter(lambda v, _: ""))
    ax.set(xlabel="electrolyte flow per half-cell / mL min⁻¹", ylabel="%", ylim=(30, 100),
           title="Flow rate: capacity against pumping")
    ax.legend(loc="lower left", fontsize=8)

    styles = {"crossover only": dict(color=plotstyle.TEAL), "crossover + 0.2 % H2": dict(color=plotstyle.AMBER)}
    ax = axes[1, 0]
    for name, rows in runs.items():
        ax.plot([r["cycle"] for r in rows], [r["capacity_Ah"] for r in rows], label=name.replace("H2", "H₂"), **styles[name])
    ax.axvline(REMIX_AT + 0.5, color=plotstyle.NICKEL, ls=":", lw=1.2)
    ax.text(REMIX_AT + 2, ax.get_ylim()[0] + 0.02 * np.ptp(ax.get_ylim()), "tanks remixed", fontsize=8.5, color=plotstyle.NICKEL)
    ax.set(xlabel="cycle", ylabel="discharge capacity / Ah", title="Capacity fade and what remixing recovers")
    ax.legend(loc="lower left", fontsize=8)

    ax = axes[1, 1]
    for name, rows in runs.items():
        ax.plot([r["cycle"] for r in rows], [100 * r["negative_share"] for r in rows], **styles[name])
    ax.axhline(50, color=plotstyle.INK, lw=0.6)
    ax.set_ylim(top=50.4)
    ax.axvline(REMIX_AT + 0.5, color=plotstyle.NICKEL, ls=":", lw=1.2)
    h2 = runs["crossover + 0.2 % H2"][REMIX_AT - 1]["aos"]
    ax.text(0.03, 0.06, f"average oxidation state before remixing\ncrossover only 3.50, with H₂ {h2:.2f}",
            transform=ax.transAxes, fontsize=8.5, color=plotstyle.INK)
    ax.set(xlabel="cycle", ylabel="active species in negative tank / %", title="Where the electrolyte went")
    for ax, l in zip(axes.flat, "abcd"): plotstyle.panel_label(ax, l)
    fig.suptitle("Flow battery at 100 mA cm⁻²: flow, pumping and crossover (all-vanadium example)", x=.01, ha="left", fontsize=13, fontweight="semibold")
    fig.savefig(out / "results.svg", metadata={"Date": None})
    return summary


if __name__ == "__main__":
    main()
