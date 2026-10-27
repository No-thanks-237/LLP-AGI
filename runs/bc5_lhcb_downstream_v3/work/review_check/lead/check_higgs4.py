
import sys, os, math
import numpy as np
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\analysis")
import run_scan as R
smp = R.load_sample(R.SAMPLE_H)
sel = np.abs(smp['status']) == 62
px, py, pz = smp['px'][sel], smp['py'][sel], smp['pz'][sel]
pt = np.hypot(px, py)
eta = np.abs(np.arcsinh(pz/np.maximum(pt, 1e-12)))
print("rows %d  status-62 events %d  |eta| max %.4f" % (len(smp['ev']), len(np.unique(smp['ev'][sel])), eta.max()))
for lo in (5.5, 6.0, 6.5, 7.0):
    print("  fraction |eta| > %.1f : %.5f" % (lo, (eta > lo).mean()))
