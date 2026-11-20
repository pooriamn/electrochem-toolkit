"""Shared Matplotlib style for the study figures.

Kept in one place so every published figure uses the same palette, type and
axis conventions. Import is lazy: the core toolkit never needs Matplotlib.
"""

TEAL = "#0B7A75"
AMBER = "#C98A16"
SLATE = "#34508A"
NICKEL = "#7A8A91"
INK = "#12242C"
RULE = "#D5DEDD"
PAPER = "#FFFFFF"


def apply():
    import matplotlib
    from matplotlib import font_manager

    matplotlib.use("Agg")
    families = {f.name for f in font_manager.fontManager.ttflist}
    sans = ["Fira Sans"] if "Fira Sans" in families else []
    matplotlib.rcParams.update({
        "svg.hashsalt": "electrochem-toolkit-study",
        "svg.fonttype": "path",
        "font.family": "sans-serif",
        "font.sans-serif": sans + ["DejaVu Sans"],
        "font.size": 10,
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.edgecolor": NICKEL,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.titlesize": 11,
        "axes.titleweight": "medium",
        "axes.titlelocation": "left",
        "axes.titlepad": 10,
        "axes.grid": True,
        "axes.axisbelow": True,
        "axes.prop_cycle": matplotlib.cycler(color=[TEAL, AMBER, SLATE, NICKEL]),
        "grid.color": RULE,
        "grid.linewidth": 0.6,
        "xtick.color": NICKEL,
        "ytick.color": NICKEL,
        "xtick.labelcolor": INK,
        "ytick.labelcolor": INK,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.minor.size": 1.5,
        "ytick.minor.size": 1.5,
        "lines.linewidth": 1.8,
        "legend.frameon": False,
        "legend.fontsize": 9,
        "figure.facecolor": PAPER,
        "savefig.facecolor": PAPER,
    })


def panel_label(ax, letter):
    """Put a bold (a), (b)... marker at the top-left of a panel."""
    ax.text(-0.02, 1.02, letter, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=11, fontweight="bold", color=INK)
