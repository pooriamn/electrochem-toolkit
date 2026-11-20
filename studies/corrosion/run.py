"""Corrosion rate from a polarisation scan: Tafel extrapolation, LPR and a full mixed-potential fit."""
import json
import math
from pathlib import Path
import numpy as np
from electrochem_toolkit import plotstyle
plotstyle.apply()
import matplotlib.pyplot as plt
from electrochem_toolkit import corrosion as co

T = co.Tafel
FE_NEUTRAL = T(1e-6, -0.75, 0.070, True)
O2 = T(1e-6, -0.30, 0.120, False, i_lim=3e-5)
HER_NEUTRAL = T(1e-6, -0.95, 0.120, False)
SCENARIOS = {
    "A": co.Electrode((T(5e-5, -0.50, 0.060, True), T(5e-5, -0.50, 0.120, False)), name="deaerated acid"),
    "B": co.Electrode((FE_NEUTRAL, O2, HER_NEUTRAL), name="aerated chloride"),
    "C": co.Electrode((FE_NEUTRAL, O2, HER_NEUTRAL), r_u=800.0, name="aerated, low conductivity"),
}
METHODS = ["Tafel extrapolation", "LPR, B = 26 mV", "LPR, B from Tafel fit", "mixed-potential fit"]


def analyse(el, seed):
    ecorr, icorr = el.corrosion_point()
    E, i, s = co.potentiodynamic_scan(el, seed=seed)
    taf = co.tafel_extrapolation(E, i)
    rp, _ = co.polarisation_resistance(E, i)
    fits = {m: co.fit_mixed_potential(E, i, s, m) for m in ("activation", "diffusion")}
    best = min(fits.values(), key=lambda f: f.aicc)
    est = {
        METHODS[0]: taf.icorr,
        METHODS[1]: 0.026 / rp,
        METHODS[2]: co.stern_geary_B(taf.beta_a, taf.beta_c) / rp,
        METHODS[3]: 10 ** best.parameters["log10_icorr"],
    }
    return {"E": E, "i": i, "sigma": s, "ecorr": ecorr, "icorr": icorr, "tafel": taf, "rp": rp,
            "fits": fits, "best": best, "estimates": est}


def main():
    out = Path("study-output/corrosion"); out.mkdir(parents=True, exist_ok=True)
    res = {k: analyse(el, seed) for (k, el), seed in zip(SCENARIOS.items(), (11, 12, 13))}
    summary = {}
    for k, r in res.items():
        b = r["best"]
        summary[k] = {
            "description": SCENARIOS[k].name,
            "true": {"Ecorr_V": r["ecorr"], "icorr_uA_cm2": r["icorr"] * 1e6,
                     "corrosion_rate_mm_y": co.corrosion_rate_mm_per_year(r["icorr"], **co.CARBON_STEEL)},
            "tafel_slopes_from_extrapolation_V": {"beta_a": r["tafel"].beta_a, "beta_c": r["tafel"].beta_c},
            "Rp_ohm_cm2": r["rp"],
            "estimates_uA_cm2": {m: v * 1e6 for m, v in r["estimates"].items()},
            "ratio_to_true": {m: v / r["icorr"] for m, v in r["estimates"].items()},
            "fit": {"model_selected": b.model, "parameters": b.parameters, "standard_errors": b.standard_errors,
                    "aicc": {m: f.aicc for m, f in r["fits"].items()}, "n": b.n},
        }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False, default=float) + "\n")
    print(json.dumps({k: {m: round(v, 3) for m, v in s["ratio_to_true"].items()} for k, s in summary.items()}, indent=2))

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.6), layout="constrained")
    for ax, key, title in [(axes[0, 0], "A", "A  Deaerated acid: both branches activation controlled"),
                           (axes[0, 1], "B", "B  Aerated chloride: oxygen reduction diffusion limited")]:
        r, el = res[key], SCENARIOS[key]
        E, i = r["E"], r["i"]
        ax.scatter(np.abs(i) * 1e6, E, s=5, color=plotstyle.INK, alpha=.55, lw=0, label="measured |i|")
        grid = np.linspace(E[0], E[-1], 400)
        ax.plot(el.anodic_metal(grid) * 1e6, grid, color=plotstyle.TEAL, ls="--", lw=1.3, label="Fe dissolution")
        cath = -(el.faradaic(grid) - el.anodic_metal(grid))
        ax.plot(cath * 1e6, grid, color=plotstyle.AMBER, ls="--", lw=1.3, label="cathodic reaction(s)")
        t = r["tafel"]
        for (m, b), lim in [(t.anodic_line, (t.ecorr, E[-1])), (t.cathodic_line, (E[0], t.ecorr))]:
            ee = np.linspace(*lim, 50)
            ax.plot(10 ** (m * ee + b) * 1e6, ee, color=plotstyle.SLATE, lw=1.2)
        for lo, hi in [(0.05, 0.15), (-0.15, -0.05)]:
            ax.axhspan(r["ecorr"] + lo, r["ecorr"] + hi, color=plotstyle.SLATE, alpha=.06, lw=0)
        ax.plot([r["icorr"] * 1e6], [r["ecorr"]], "o", ms=7, mfc="white", mec=plotstyle.INK, mew=1.5, label="true corrosion point")
        ax.plot([t.icorr * 1e6], [t.ecorr], "x", ms=8, color=plotstyle.SLATE, mew=2, label="Tafel extrapolation")
        ax.set_xscale("log"); ax.set_xlim(0.5, 3e3 if key == "A" else 500)
        ax.set(xlabel="|i| / µA cm⁻²", ylabel="E / V vs SCE", title=title)
        ax.title.set_fontsize(10)
    axes[0, 1].legend(loc="upper right", fontsize=7.5)

    ax = axes[1, 0]
    for key, c in (("B", plotstyle.TEAL), ("C", plotstyle.AMBER)):
        r = res[key]
        m = np.abs(r["E"] - r["ecorr"]) <= 0.02
        ax.scatter(r["i"][m] * 1e6, (r["E"][m] - r["ecorr"]) * 1e3, s=9, color=c, lw=0, alpha=.8)
        x = np.linspace(-12, 12, 10)
        ax.plot(x, x * 1e-6 * r["rp"] * 1e3, color=c,
                label=f"{key}: Rp = {r['rp']:.0f} Ω cm²" + ("  (includes Ru = 800)" if key == "C" else ""))
    ax.axhspan(-10, 10, color=plotstyle.NICKEL, alpha=.07, lw=0)
    ax.axhline(0, color=plotstyle.INK, lw=.6); ax.axvline(0, color=plotstyle.INK, lw=.6)
    ax.set(xlabel="i / µA cm⁻²", ylabel="E − Ecorr / mV", title="Linear polarisation: the ohmic drop hides in Rp", xlim=(-12, 12))
    ax.legend(loc="upper left", fontsize=8)

    ax = axes[1, 1]
    ax.axhspan(0.8, 1.25, color=plotstyle.TEAL, alpha=.08, lw=0)
    ax.axhline(1, color=plotstyle.INK, lw=.7)
    colours = {"A": plotstyle.SLATE, "B": plotstyle.TEAL, "C": plotstyle.AMBER}
    for j, (key, r) in enumerate(res.items()):
        xs = np.arange(len(METHODS)) + (j - 1) * 0.2
        ys = [r["estimates"][m] / r["icorr"] for m in METHODS]
        ax.scatter(xs, ys, s=42, color=colours[key], zorder=3, label=f"{key}  {SCENARIOS[key].name}")
        se = r["best"].standard_errors["log10_icorr"]
        ax.errorbar(xs[-1], ys[-1], yerr=[[ys[-1] * (1 - 10 ** (-2 * se))], [ys[-1] * (10 ** (2 * se) - 1)]],
                    color=colours[key], capsize=3, lw=1.2)
    ax.set_yscale("log"); ax.set_ylim(0.4, 2.5)
    ax.set_yticks([0.5, 0.8, 1, 1.25, 2]); ax.set_yticklabels(["0.5", "0.8", "1", "1.25", "2"])
    ax.minorticks_off()
    ax.set_xticks(range(len(METHODS)), [m.replace(", ", ",\n") for m in METHODS], fontsize=8.5)
    ax.set(ylabel="estimated ÷ true corrosion current", title="Which estimate lands within 20 %?")
    ax.legend(loc="upper left", fontsize=7.5)
    for ax, l in zip(axes.flat, "abcd"): plotstyle.panel_label(ax, l)
    fig.suptitle("Carbon steel corrosion rate from one polarisation scan, three ways", x=.01, ha="left", fontsize=13, fontweight="semibold")
    fig.savefig(out / "results.svg", metadata={"Date": None})
    return summary


if __name__ == "__main__":
    main()
