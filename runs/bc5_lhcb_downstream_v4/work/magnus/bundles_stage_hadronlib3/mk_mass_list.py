"""Build the mass list used by the Pythia hadronisation-library job.

Run from anywhere on the host: the scan grid is imported from bc5.analysis.scan, so the list is
always the same as the production grid.  The output is written next to this script.
"""
import math
import sys
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(RUN))

from bc5.analysis import scan as R

grid = np.exp(np.linspace(math.log(R.MASS_MIN), math.log(R.MASS_MAX), 52))
masses = sorted(set(round(float(m), 6) for m in R.merge_mass_grid(grid)))
print("n =", len(masses))
out = Path(__file__).resolve().parent / "mass_list.txt"
with open(out, "w") as fh:
    for m in masses:
        fh.write("%.6f\n" % m)
print(open(out).read().split("\n")[:5], "...", open(out).read().strip().split("\n")[-3:])
