"""Tafel analysis of a synthetic HER sweep: how iR, window choice and method change the slope."""
import json
import math
from pathlib import Path
import numpy as np
from electrochem_toolkit import plotstyle
plotstyle.apply()
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from electrochem_toolkit import kinetics as k

EL = k.HerElectrode()
TRUE_SLOPE = k.tafel_slope_mV(EL.alpha_c)
LOW, HIGH = (1e-3, 1e-2), (1e-2, 5e-2)


def main():
    out = Path("study-output/tafel-analysis"); out.mkdir(parents=True, exist_ok=True)
    E, j, s = k.sweep(EL, seed=1)
    E85, j85, _ = k.sweep(EL, compensated=0.85, seed=1)
    corr = k.ir_correct(E, j, EL.r_u)
    hfr = k.ir_correct(E, j, 1.05 * EL.r_u)
    lx, ls = k.local_tafel_slope(corr, j)
    plateau, p_lo, p_hi = k.plateau_slope(lx, ls)
    bv = k.fit_butler_volmer(E, j, s)

    methods = {
        "raw, 10–50 mA cm⁻²": k.fit_tafel_line(E, j, j_range=HIGH),
        "raw, 1–10 mA cm⁻²": k.fit_tafel_line(E, j, j_range=LOW),
        "85 % iR compensated": k.fit_tafel_line(E85, j85, j_range=LOW),
        "iR corrected, Ru 5 % high": k.fit_tafel_line(hfr, j, j_range=LOW),
    }
    estimates = {m: {"slope_mV": f.slope_mV, "j0_A_cm2": f.j0, "r_squared": f.r_squared} for m, f in methods.items()}
    estimates["Tafel slope plot plateau"] = {"slope_mV": plateau, "j0_A_cm2": None, "r_squared": None}
    estimates["Butler–Volmer fit, Ru free"] = {"slope_mV": bv.parameters["tafel_slope_mV"],
                                               "j0_A_cm2": 10 ** bv.parameters["log10_j0"], "r_squared": None}

    starts, widths = np.arange(0.02, 0.401, 0.01), np.arange(0.03, 0.201, 0.01)
    err = np.full((len(widths), len(starts)), np.nan)
    for a, st in enumerate(starts):
        for b, w in enumerate(widths):
            try:
                err[b, a] = 100 * (k.fit_tafel_line(E, j, eta_range=(st, st + w)).slope_mV / TRUE_SLOPE - 1)
            except ValueError:
                pass
    good = np.abs(err) <= 5

    summary = {"truth": {"slope_mV": TRUE_SLOPE, "j0_A_cm2": EL.j0, "Ru_ohm_cm2": EL.r_u,
                         "j_capacitive_A_cm2": -EL.c_dl * 0.005},
               "estimates": estimates,
               "tafel_plot_plateau_log10_j": [p_lo, p_hi],
               "bv_fit": {"parameters": bv.parameters, "standard_errors": bv.standard_errors,
                          "weighted_rss": bv.weighted_rss, "n": len(E)},
               "raw_windows_within_5_percent": int(good.sum()), "raw_windows_total": int(np.isfinite(err).sum())}
    (out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False, default=float) + "\n")
    print(json.dumps({m: round(v["slope_mV"], 1) for m, v in estimates.items()}, indent=2))

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.6), layout="constrained")
    ax = axes[0, 0]
    m = j < 0
    ax.plot(np.log10(-j[m]), -E[m] * 1e3, color=plotstyle.AMBER, lw=1.6, label="as measured")
    ax.plot(np.log10(-j85[j85 < 0]), -E85[j85 < 0] * 1e3, color=plotstyle.SLATE, lw=1.6, label="85 % compensated")
    ax.plot(np.log10(-j[m]), -corr[m] * 1e3, color=plotstyle.TEAL, lw=1.6, label="iR corrected")
    xl = np.linspace(-6, -1, 50)
    ax.plot(xl, TRUE_SLOPE * (xl - math.log10(EL.j0)), color=plotstyle.INK, ls=":", lw=1.2, label="true Tafel line")
    for (lo, hi), c, lab in ((LOW, plotstyle.NICKEL, "1–10\nmA cm⁻²"), (HIGH, plotstyle.AMBER, "10–50\nmA cm⁻²")):
        ax.axvspan(math.log10(lo), math.log10(hi), color=c, alpha=.08, lw=0)
        ax.text((math.log10(lo) + math.log10(hi)) / 2, 18, lab, ha="center", fontsize=7.5, color=plotstyle.NICKEL)
    ax.set(xlim=(-6, -1), ylim=(0, 720), xlabel="log₁₀ |j| / A cm⁻²", ylabel="−η / mV",
           title="Tafel plot: the same sweep three ways")
    ax.legend(loc="upper left", fontsize=8)

    ax = axes[0, 1]
    for eta, jj, c, lab in ((E, j, plotstyle.AMBER, "as measured"), (E85, j85, plotstyle.SLATE, "85 % compensated"),
                            (corr, j, plotstyle.TEAL, "iR corrected")):
        x, sl = k.local_tafel_slope(eta, jj)
        ax.plot(x, sl, color=c, lw=1.4, label=lab)
    ax.axhline(TRUE_SLOPE, color=plotstyle.INK, ls=":", lw=1.2)
    ax.axvspan(p_lo, p_hi, color=plotstyle.TEAL, alpha=.08, lw=0)
    ax.text((p_lo + p_hi) / 2, 330, f"plateau {plateau:.1f} mV/dec", ha="center", fontsize=8.5)
    ax.annotate("back reaction and\ncharging current", xy=(-5.55, 45), xytext=(-5.9, 215), fontsize=8,
                color=plotstyle.NICKEL, arrowprops=dict(arrowstyle="-", color=plotstyle.NICKEL, lw=.8))
    ax.set(xlim=(-6, -1), ylim=(0, 400), xlabel="log₁₀ |j| / A cm⁻²", ylabel="local slope / mV dec⁻¹",
           title="Tafel slope plot (20 mV windows)")
    ax.legend(loc="upper left", fontsize=8)

    ax = axes[1, 0]
    norm = TwoSlopeNorm(vmin=-20, vcenter=0, vmax=60)
    pc = ax.pcolormesh(starts * 1e3, widths * 1e3, np.clip(err, -20, 60), cmap="BrBG_r", norm=norm, shading="nearest")
    ax.contour(starts * 1e3, widths * 1e3, np.abs(err), levels=[5], colors=[plotstyle.INK], linewidths=1.2)
    cb = fig.colorbar(pc, ax=ax, shrink=.9); cb.set_label("slope error / %")
    cb.outline.set_visible(False)
    ax.grid(False)
    ax.text(150, 34, "between the lines:\nwithin 5 %", fontsize=8.5, color=plotstyle.INK)
    ax.set(xlabel="window starts at −η / mV (uncorrected)", ylabel="window width / mV",
           title="Fitted slope against window choice, raw data")

    ax = axes[1, 1]
    names = list(estimates)
    vals = [estimates[n]["slope_mV"] for n in names]
    ypos = np.arange(len(names))[::-1]
    cols = [plotstyle.AMBER, plotstyle.AMBER, plotstyle.SLATE, plotstyle.TEAL, plotstyle.TEAL, plotstyle.TEAL]
    ax.barh(ypos, vals, color=cols, height=.6)
    ax.axvline(TRUE_SLOPE, color=plotstyle.INK, ls=":", lw=1.2)
    for y, n, v in zip(ypos, names, vals):
        j0 = estimates[n]["j0_A_cm2"]
        txt = f"{v:.0f}" + (f"   j₀ ×{j0 / EL.j0:.2g}" if j0 else "")
        ax.text(v + 4, y, txt, va="center", fontsize=8.5)
    ax.set_yticks(ypos, names, fontsize=8.5)
    ax.set_xlim(0, 300); ax.grid(axis="y", visible=False)
    ax.set(xlabel="Tafel slope / mV dec⁻¹", title=f"Reported slope (true {TRUE_SLOPE:.0f} mV/dec)")
    for ax, l in zip(axes.flat, "abcd"): plotstyle.panel_label(ax, l)
    fig.suptitle("One hydrogen evolution sweep, six Tafel slopes", x=.01, ha="left", fontsize=13, fontweight="semibold")
    fig.savefig(out / "results.svg", metadata={"Date": None})
    return summary


if __name__ == "__main__":
    main()
