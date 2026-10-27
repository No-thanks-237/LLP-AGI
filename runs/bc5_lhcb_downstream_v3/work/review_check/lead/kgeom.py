
import sys, os, math, random
import numpy as np
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\analysis")
import run_scan as R
smp = R.load_sample('higgs.txt.gz')
sel = np.abs(smp['status']) == 62
px, py, pz, e = smp['px'][sel], smp['py'][sel], smp['pz'][sel], smp['e'][sel]
m = np.sqrt(np.maximum(e*e - px*px - py*py - pz*pz, 0.0))
rng = np.random.default_rng(12345)
print("Higgs sample: %d events, <pT>=%.2f GeV, <eta>=%.3f, m=%.3f+-%.4f" % (
    len(e), np.hypot(px, py).mean(), np.abs(np.arcsinh(pz/np.maximum(np.hypot(px,py),1e-9))).mean(),
    m.mean(), m.std()))
def ct(rng, n):
    u = rng.uniform(-1, 1, n); ph = rng.uniform(0, 2*np.pi, n)
    st = np.sqrt(np.maximum(1-u*u, 0))
    return st*np.cos(ph), st*np.sin(ph), u
print("\n%8s %10s %10s %10s" % ("m_S", "K_geom", "se", "K_analytic"))
for ms in (5.082, 13.023, 33.371, 53.42):
    n = len(e)
    # two scalars per Higgs, isotropic in the h rest frame: p* = sqrt(m_h^2/4 - m_S^2)
    pstar = np.sqrt(np.maximum(m*m/4.0 - ms*ms, 0.0))
    # build the h four-momentum and boost
    tot = 0.0
    for k in range(2):
        ux, uy, uz = ct(rng, n)
        ps = np.stack([pstar*ux, pstar*uy, pstar*uz])
        es = np.sqrt(ms*ms + (pstar**2))
        # boost from the h rest frame to the lab
        g = e/m; b = np.stack([px, py, pz])/e
        bp = b[0]*ps[0] + b[1]*ps[1] + b[2]*ps[2]
        pz_lab = ps[2] + (((g-1)*bp/np.maximum((b*b).sum(axis=0),1e-12)) + g*es)*b[2]
        e_lab = g*(es + bp)
        pmod = np.sqrt(es*es - ms*ms)  # |p| in the h rest frame
        # lab momentum magnitude
        px_lab = ps[0] + (((g-1)*bp/np.maximum((b*b).sum(axis=0),1e-12)) + g*es)*b[0]
        py_lab = ps[1] + (((g-1)*bp/np.maximum((b*b).sum(axis=0),1e-12)) + g*es)*b[1]
        p_lab = np.sqrt(px_lab**2 + py_lab**2 + pz_lab**2)
        eta = 0.5*np.log((p_lab + pz_lab)/np.maximum(p_lab - pz_lab, 1e-30))
        w = (eta > 2.0) & (eta < 5.0) & (pz_lab > 0)
        tot += np.sum(np.where(w, ms/np.maximum(pz_lab, 1e-9), 0.0))
    K = tot/(2*n)   # per scalar: 2 scalars per Higgs
    se = np.std([0])  # placeholder
    print("%8.3f %10.5f %10s %10s" % (ms, K, "-", "-"))
