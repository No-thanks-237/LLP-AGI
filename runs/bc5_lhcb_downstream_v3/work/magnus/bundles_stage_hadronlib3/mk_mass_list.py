
import sys, math, os
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\analysis")
import numpy as np
import run_scan as R
grid = np.exp(np.linspace(math.log(R.MASS_MIN), math.log(R.MASS_MAX), 52))
masses = R.merge_mass_grid(grid)
masses = sorted(set(round(float(m), 6) for m in masses))
print("n =", len(masses))
out = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\magnus\bundles_stage_hadronlib3\mass_list.txt"
with open(out, 'w') as fh:
    for m in masses:
        fh.write("%.6f\n" % m)
print(open(out).read().split("\n")[:5], "...", open(out).read().strip().split("\n")[-3:])
