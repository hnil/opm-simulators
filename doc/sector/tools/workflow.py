#!/usr/bin/env python3
"""Draw the sector modelling workflow diagram.

Usage: workflow.py <output.png>
"""

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

fig, ax = plt.subplots(figsize=(10, 4.2))
ax.set_xlim(0, 10)
ax.set_ylim(0, 4.2)
ax.axis("off")


def box(x, y, w, h, text, color):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05",
                                fc=color, ec="0.3", lw=1.2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=9)


def arrow(x0, y0, x1, y1, text=""):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                                 mutation_scale=14, lw=1.2, color="0.25"))
    if text:
        ax.text((x0 + x1) / 2, (y0 + y1) / 2 + 0.12, text, ha="center",
                va="bottom", fontsize=8, style="italic")


parent = "#d9e8f5"
tool = "#fbe5c8"
data = "#e3f1d8"
reduced = "#f3d9e6"

box(0.1, 2.75, 2.2, 1.1, "Parent run\nDUMPFLUX, FLUXNUM\nFLUXALL", parent)
box(0.1, 0.35, 2.2, 1.1, "Parent run\nFLUXALL\n(restarts)", parent)
box(3.0, 0.35, 1.6, 1.1, "make_flux", tool)
box(5.1, 1.55, 1.6, 1.1, "FLUX file\nboundary,\nsummary, wells", data)
box(7.4, 2.75, 2.5, 1.1, "Reduced run on parent deck\nUSEFLUX, FLUXREG\n(inactive outside)", reduced)
box(7.4, 0.35, 2.5, 1.1, "Reduced run on\ncut-out sector deck\nUSEFLUX, FLUXREG", reduced)

arrow(2.35, 3.3, 5.05, 2.45, "during the run")
arrow(2.35, 0.9, 2.95, 0.9, "output")
arrow(4.65, 0.9, 5.05, 1.75)
arrow(6.75, 2.3, 7.35, 3.2)
arrow(6.75, 1.9, 7.35, 0.95)

fig.tight_layout()
fig.savefig(sys.argv[1], dpi=130)
