#!/usr/bin/env python3
"""Compare a full-field run with reduced sector runs of it.

Writes PNG plots of field and well vectors and a markdown table of the
largest deviations from the full run.

Usage:
  compare.py --pydir=<opm-common python dir> --out=<dir> --prefix=<name>
             --full=<case> [--wells=A3,A4,A6] [--extra=KEY,...]
             <label>=<case> [<label>=<case> ...]
"""

import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--pydir", required=True)
parser.add_argument("--out", required=True)
parser.add_argument("--prefix", required=True)
parser.add_argument("--full", required=True)
parser.add_argument("--wells", default="A3,A4,A6")
parser.add_argument("--extra", default="")
parser.add_argument("cases", nargs="+")
args = parser.parse_args()

sys.path.insert(0, args.pydir)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from opm.io.ecl import ESmry

full = ESmry(args.full + ".SMSPEC")
cases = [(lab, ESmry(case + ".SMSPEC"))
         for lab, case in (c.split("=", 1) for c in args.cases)]

WELLS = [w for w in args.wells.split(",") if w]
WELL_VECTORS = ["WOPR", "WOPT", "WWIR", "WWIT", "WGIR", "WGIT", "WBHP", "WTHP"]
UNITS = {"OPR": "Sm3/day", "OPT": "Sm3", "WIR": "Sm3/day", "WIT": "Sm3",
         "GIR": "Sm3/day", "GIT": "Sm3", "BHP": "bar", "THP": "bar",
         "PR": "bar", "PRFLUX": "bar", "UOUTOPR": "Sm3/day", "USECOPR": "Sm3/day"}
STYLES = ["C0--", "C1-.", "C3:", "C2--"]


def days(s):
    return np.asarray(s["TIME"])


def series(s, key):
    return np.asarray(s[key]) if key in s else None


def on_full_times(s, key):
    v = series(s, key)
    return None if v is None else np.interp(days(full), days(s), v)


def unit(key):
    return UNITS.get(key.split(":")[0][1:], "")


well_keys = [f"{v}:{w}" for w in WELLS for v in WELL_VECTORS
             if series(full, f"{v}:{w}") is not None
             and np.max(np.abs(series(full, f"{v}:{w}"))) > 0]
field_keys = ["FOPR", "FOPT", "FPR"]
extra_keys = [k for k in args.extra.split(",") if k]

rows = []
header = ["Vector", "Full (end)"]
for lab, _ in cases:
    header += [lab]
rows.append("| " + " | ".join(header) + " |")
rows.append("|" + "|".join([":--"] + ["--:"] * (len(header) - 1)) + "|")


def table_row(key, ref_key=None):
    ref = series(full, ref_key or key)
    label = f"`{key}`" + (f" (full: `{ref_key}`)" if ref_key else "")
    cells = [label, f"{ref[-1]:.5g}"]
    for _, s in cases:
        val = on_full_times(s, key)
        if val is None:
            cells += ["–"]
            continue
        scale = np.max(np.abs(ref))
        rel = np.max(np.abs(val - ref)) / scale if scale > 0 else 0.0
        cells += [f"{100 * rel:.2f} %"]
    rows.append("| " + " | ".join(cells) + " |")


for k in field_keys:
    table_row(k)
table_row("FPRFLUX", "FPR")
for k in extra_keys:
    table_row(k)
for k in well_keys:
    table_row(k)

with open(f"{args.out}/{args.prefix}_table.md", "w") as f:
    f.write("\n".join(rows) + "\n")


def plot(keys, fname, ref_keys=None, ncols=2):
    nrows = (len(keys) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 2.9 * nrows),
                             squeeze=False)
    for ax, key in zip(axes.flat, keys):
        ref_key = (ref_keys or {}).get(key, key)
        ax.plot(days(full), series(full, ref_key), "k-", lw=2.2,
                label="Full" + (f" ({ref_key})" if ref_key != key else ""))
        for (lab, s), style in zip(cases, STYLES):
            if key in s:
                ax.plot(days(s), series(s, key), style, lw=1.5, label=lab)
        ax.set_title(key, fontsize=10)
        ax.set_xlabel("Days", fontsize=8)
        ax.set_ylabel(unit(key), fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    for ax in list(axes.flat)[len(keys):]:
        ax.set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{args.out}/{fname}", dpi=100)
    plt.close(fig)


plot(["FOPR", "FOPT"], f"{args.prefix}_field_oil.png")
plot(["FPR", "FPRFLUX"], f"{args.prefix}_field_pressure.png",
     ref_keys={"FPRFLUX": "FPR"})
if extra_keys:
    plot(extra_keys, f"{args.prefix}_extra.png")
for w in WELLS:
    keys = [k for k in well_keys if k.endswith(":" + w)]
    if keys:
        plot(keys, f"{args.prefix}_well_{w}.png")

print("\n".join(rows))
