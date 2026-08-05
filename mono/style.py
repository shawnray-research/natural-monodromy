"""
Shared publication style: real Computer Modern via LaTeX, so the figures match
the typography of the papers they sit beside.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, AutoMinorLocator

# ---------------------------------------------------------------- palette ---
SHAPE = "#1b2733"      # the manifold M
BB = "#1f6fb2"         # birth-birth: minima, inner bitangent circle, births
DD = "#c0392b"         # death-death: maxima, outer bitangent circle, deaths
LOOP = "#6d28d9"       # the observation loop gamma
MUTED = "#8b98a5"
INK = "#1b2733"


def use_paper_style(usetex=True):
    rc = {
        "font.family": "serif",
        "font.size": 9.5,
        "axes.titlesize": 10,
        "axes.labelsize": 10,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.linewidth": 0.7,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 3.2,
        "ytick.major.size": 3.2,
        "xtick.minor.size": 1.8,
        "ytick.minor.size": 1.8,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "legend.frameon": True,
        "legend.framealpha": 1.0,
        "legend.edgecolor": "#c8d0d8",
        "legend.borderpad": 0.45,
        "legend.handlelength": 1.6,
        "legend.labelspacing": 0.42,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "lines.solid_capstyle": "round",
    }
    if usetex:
        rc.update({
            "text.usetex": True,
            "text.latex.preamble": r"\usepackage{amsmath}\usepackage{amssymb}",
        })
    else:
        rc.update({"text.usetex": False, "mathtext.fontset": "stix",
                   "font.serif": ["STIX Two Text", "STIXGeneral", "DejaVu Serif"]})
    plt.rcParams.update(rc)


def tidy(ax, nx=4, ny=4, minor=True):
    """Sparse, non-overlapping ticks."""
    ax.xaxis.set_major_locator(MaxNLocator(nx, prune="both"))
    ax.yaxis.set_major_locator(MaxNLocator(ny, prune="both"))
    if minor:
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which="both", top=True, right=True)


def panel(ax, letter, dx=-0.02, dy=1.04, size=11):
    ax.text(dx, dy, rf"\textbf{{({letter})}}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=size)


def panel3d(ax, letter, x=0.02, y=0.94, size=11):
    ax.text2D(x, y, rf"\textbf{{({letter})}}", transform=ax.transAxes,
              ha="left", va="top", fontsize=size)
