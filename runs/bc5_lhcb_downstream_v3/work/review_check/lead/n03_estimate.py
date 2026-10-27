
import sys, os, math
import numpy as np
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\analysis")
import run_scan as R
smp = R.load_sample(R.SAMPLE_H)
sel = np.abs(smp['status']) == 62
px, py, pz = smp['px'][sel], smp['py'][sel], smp['pz'][sel]
pt = np.hypot(px, py)
eta = np.arcsinh(pz/np.maximum(pt, 1e-12))
a = np.abs(eta)
print("n = %d, |eta| max = %.5f" % (len(a), a.max()))
for lo in (4.0, 4.5, 5.0, 5.5, 5.7, 5.9):
    print("  fraction |eta| > %.1f : %.4f" % (lo, (a > lo).mean()))
# exponential tail fit above 5.0 in |eta|
x = a[a > 5.0]
print("tail above 5.0: n=%d  mean |eta| = %.3f" % (len(x), x.mean()))
# fit lambda from the survival: N(>eta) ~ A exp(-eta/lambda) for eta in [5,6]
bins = np.linspace(5.0, 6.0, 11)
h, _ = np.histogram(a, bins=bins)
cent = 0.5*(bins[1:]+bins[:-1])
m = h > 0
coef = np.polyfit(cent[m], np.log(h[m]), 1)
lam = -1.0/coef[0]
print("fitted exponential slope lambda = %.3f (dN/deta ~ exp(-eta/lambda))" % lam)
frac_above6 = (h[-1]/ (h.sum())) * (lam/(bins[1]-bins[0])) * 0  # placeholder
# integrate the fitted tail from 6 to infinity, normalised to the total
d_eta = bins[1]-bins[0]
n_tail_meas = h.sum()
A = h[0]/math.exp(-cent[0]/lam)
n_above6 = A*lam*math.exp(-6.0/lam)
print("estimated fraction |eta|>6: %.4f  (from the fitted tail, %d events observed in [5,6])" % (
    n_above6/len(a), n_tail_meas))
