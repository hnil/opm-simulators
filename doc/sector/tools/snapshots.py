#!/usr/bin/env python3
"""Grid snapshots for the sector modelling documentation, made with ResInsight.

Usage: snapshots.py <run dir> <output dir>

The run directory is the one run-examples.sh fills in. Set
RESINSIGHT_EXECUTABLE if ResInsight is not on the PATH.
"""

import glob
import os
import sys

import rips
from PIL import Image, ImageChops

run = os.path.abspath(sys.argv[1])
out = os.path.abspath(sys.argv[2])
os.makedirs(out, exist_ok=True)

resinsight = rips.Instance.launch(
    resinsight_executable=os.environ.get("RESINSIGHT_EXECUTABLE", "resinsight"))
if resinsight is None:
    sys.exit("could not launch ResInsight")


def trim(path):
    """Crop the white margin ResInsight leaves around the model."""
    img = Image.open(path).convert("RGB")
    bg = Image.new("RGB", img.size, (255, 255, 255))
    box = ImageChops.difference(img, bg).getbbox()
    if box:
        pad = 10
        box = (max(box[0] - pad, 0), max(box[1] - pad, 0),
               min(box[2] + pad, img.size[0]), min(box[3] + pad, img.size[1]))
        img.crop(box).save(path)


def snapshot(case, prop, ptype, step, name):
    view = case.create_view()
    view.apply_cell_result(result_type=ptype, result_variable=prop)
    view.set_time_step(step)
    tmp = os.path.join(out, "tmp")
    os.makedirs(tmp, exist_ok=True)
    view.export_snapshot(prefix=name + "_", export_folder=tmp)
    produced = glob.glob(os.path.join(tmp, name + "_*.png"))
    target = os.path.join(out, name + ".png")
    os.replace(produced[0], target)
    trim(target)


try:
    project = resinsight.project
    full = project.load_case(f"{run}/DUMP-SECTOR.EGRID")
    same_grid = project.load_case(f"{run}/USE-SECTOR.EGRID")
    cutout = project.load_case(f"{run}/sector/job/SECTOR-UF-R2.EGRID")

    def last(case):
        return len(case.time_steps()) - 1

    snapshot(full, "FLUXNUM", "STATIC_NATIVE", 0, "full_fluxnum")
    snapshot(full, "PRESSURE", "DYNAMIC_NATIVE", last(full), "full_pressure")
    snapshot(same_grid, "PRESSURE", "DYNAMIC_NATIVE", last(same_grid), "usesector_pressure")
    snapshot(cutout, "PRESSURE", "DYNAMIC_NATIVE", last(cutout), "cutout_pressure")
finally:
    resinsight.exit()

for f in glob.glob(os.path.join(out, "tmp", "*")):
    os.remove(f)
print(sorted(p for p in os.listdir(out) if p.endswith(".png")))
