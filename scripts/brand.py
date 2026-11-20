"""Generate the README hero, profile banner and social preview images.

    python scripts/brand.py            # writes into docs/assets/

The motif is a family of impedance arcs for a resistor in parallel with a
constant-phase element, Z = R / (1 + (jωτ)^α), drawn for α from 1 down to 0.55.
As α falls the arc flattens: the familiar "depressed semicircle" of real
electrodes. Text is converted to paths so the SVGs render identically on GitHub.
"""
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from electrochem_toolkit import plotstyle  # noqa: E402

plotstyle.apply()
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "docs" / "assets"

THEMES = {
    "light": dict(bg="#F3F6F5", ink="#12242C", muted="#5E6F77", rule="#D5DEDD",
                  arc="#0B7A75", hot="#C98A16"),
    "dark": dict(bg="#0C171C", ink="#E6EEEC", muted="#8FA3A8", rule="#20343B",
                 arc="#4FC1B8", hot="#E9B44C"),
}


def arc(alpha, rs=0.0, r=1.0, n=400):
    w = np.logspace(-4, 4, n)
    z = rs + r / (1 + (1j * w) ** alpha)
    return z.real, -z.imag


def mix(c, bg, t):
    a, b = np.array(to_rgb(c)), np.array(to_rgb(bg))
    return tuple(b + (a - b) * t)


def draw_arcs(ax, th, x0, y0, scale, alphas=(1.0, 0.9, 0.8, 0.7, 0.62, 0.55)):
    ax.plot([x0 - 0.04 * scale, x0 + 1.08 * scale], [y0, y0], color=th["rule"], lw=1.2)
    ax.plot([x0, x0], [y0 - 0.02 * scale, y0 + 0.6 * scale], color=th["rule"], lw=1.2)
    for k, a in enumerate(alphas):
        x, y = arc(a)
        t = 1.0 - k / (len(alphas) + 0.5)
        ax.plot(x0 + x * scale, y0 + y * scale, color=mix(th["arc"], th["bg"], t),
                lw=2.6 if k == 2 else 1.5, solid_capstyle="round")
    # sampled "measurement" points on the α = 0.8 arc
    w = np.logspace(-2.2, 2.2, 15)
    z = 1 / (1 + (1j * w) ** 0.8)
    ax.scatter(x0 + z.real * scale, y0 - z.imag * scale, s=34, facecolor=th["bg"],
               edgecolor=th["hot"], linewidth=1.6, zorder=5)


def canvas(w, h, th):
    fig = plt.figure(figsize=(w / 100, h / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, w)
    ax.set_ylim(0, h)
    ax.axis("off")
    ax.grid(False)
    fig.patch.set_facecolor(th["bg"])
    return fig, ax


def save(fig, name):
    fig.savefig(OUT / name, metadata={"Date": None}, facecolor=fig.get_facecolor())
    plt.close(fig)


def toolkit_hero(theme):
    th = THEMES[theme]
    fig, ax = canvas(1600, 520, th)
    ax.text(80, 330, "Electrochemical", fontsize=58, fontweight="semibold", color=th["ink"])
    ax.text(80, 255, "Engineering Toolkit", fontsize=58, fontweight="semibold", color=th["ink"])
    ax.text(82, 185, "Impedance fitting, electrolyser stacks and measurement",
            fontsize=21, color=th["muted"])
    ax.text(82, 152, "reports, with the physics written down and tested.",
            fontsize=21, color=th["muted"])
    ax.text(82, 82, "Python 3.10+   /   MIT licence   /   tested on every commit", fontsize=16,
            color=th["arc"])
    draw_arcs(ax, th, x0=1000, y0=110, scale=520)
    ax.text(1540, 70, "R with a constant-phase element, α from 1 to 0.55", fontsize=13, color=th["muted"], ha="right")
    save(fig, f"hero-{theme}.svg")


def profile_banner(theme):
    th = THEMES[theme]
    fig, ax = canvas(1600, 400, th)
    ax.text(80, 236, "Pooria Moozarm Nia", fontsize=56, fontweight="semibold", color=th["ink"])
    ax.text(83, 170, "Electrochemist and Principal Scientist, London", fontsize=23, color=th["muted"])
    ax.text(83, 112, "Water electrolysis  /  impedance  /  scientific software",
            fontsize=18, color=th["arc"])
    draw_arcs(ax, th, x0=1060, y0=80, scale=440)
    save(fig, f"profile-banner-{theme}.svg")


def social_preview():
    th = THEMES["dark"]
    fig, ax = canvas(1280, 640, th)
    ax.text(70, 420, "Electrochemical", fontsize=52, fontweight="semibold", color=th["ink"])
    ax.text(70, 350, "Engineering Toolkit", fontsize=52, fontweight="semibold", color=th["ink"])
    ax.text(72, 280, "EIS inverse modelling, electrolyser stacks,", fontsize=21, color=th["muted"])
    ax.text(72, 246, "electrothermal simulation, all tested in CI", fontsize=21, color=th["muted"])
    ax.text(72, 90, "github.com/pooriamn/electrochem-toolkit", fontsize=17, color=th["arc"])
    draw_arcs(ax, th, x0=760, y0=150, scale=460)
    fig.savefig(OUT / "social-preview.png", dpi=100, facecolor=th["bg"])
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    for t in THEMES:
        toolkit_hero(t)
        profile_banner(t)
    social_preview()
    print("wrote", *sorted(p.name for p in OUT.glob("*-*.svg")), "social-preview.png")
