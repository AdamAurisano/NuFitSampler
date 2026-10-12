"""Shared plot styling. Colours follow a fixed role assignment, never cycled:
exact/direct sampler = blue, PT = orange, legacy = aqua, NuFIT published = ink."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8985"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"
DIRECT = "#2a78d6"
PT = "#eb6834"
LEGACY = "#1baf7a"
EXTRA = "#6250d6"
CMAP = "Blues"
DIVERGING = "RdBu_r"

# 2-dof Delta chi^2 levels for 1, 2, 3 sigma and the matching credible masses
LEVELS_2D = [2.30, 6.18, 11.83]
MASSES = [0.6827, 0.9545, 0.9973]

plt.rcParams.update({
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "axes.edgecolor": MUTED,
    "axes.labelcolor": INK,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "legend.frameon": False,
    "lines.linewidth": 1.6,
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def page(nrows, ncols, title, w=11, h=8.5):
    fig, axs = plt.subplots(nrows, ncols, figsize=(w, h), squeeze=False)
    fig.suptitle(title, fontsize=12, color=INK, x=0.01, ha="left")
    return fig, axs


def note(fig, text):
    fig.text(0.01, 0.005, text, fontsize=7.5, color=INK2, ha="left", va="bottom", wrap=True)
