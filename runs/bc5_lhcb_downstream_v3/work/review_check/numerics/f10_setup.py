# -*- coding: utf-8 -*-
"""F10 setup: z0 distribution of the real parent samples and run_mass tables."""
import sys, os, math, random, time
import numpy as np
RUN = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, os.path.join(RUN,"work","analysis"))
sys.path.insert(0, os.path.join(RUN,"work","widths"))
sys.path.insert(0, os.path.join(RUN,"work","response"))
import run_scan as R
import bc5_lib as L

# sample z ranges (mm in the files)
for name in ('bbar_heavy.txt','minbias_light.txt','higgs.txt'):
    smp = R.load_sample(name)
    if smp is None:
        print(name, "MISSING"); continue
    z = smp['z']
    print("%-20s rows=%7d  z[mm]: min %.4f  median %.4f  p99 %.4f  max %.4f" % (
        name, len(z), z.min(), np.median(z), np.percentile(z,99), z.max()))

ok = L.load_hadron_lib(os.path.join(RUN,'work','magnus','hadron_lib'))
print("hadron lib loaded:", ok)
rng = random.Random(20261024)
for ms in (0.2802, 1.0, 3.5567):
    t0=time.time()
    tab = R.run_mass(ms, 300, rng, ndecay=8, verbose=False)
    z = np.asarray(tab['z']); pz = np.asarray(tab['pz']); A = np.asarray(tab['A'])
    e = np.asarray(tab['e']); px=np.asarray(tab['px']); py=np.asarray(tab['py'])
    pmod = np.sqrt(px*px+py*py+pz*pz)
    print("m=%.2f  rows=%d  z0[m]: min %.5f max %.5f median %.5f  pz: min %.3f median %.2f max %.1f  |p|/pz: median %.4f max %.3f  (%.1fs)" % (
        ms, len(z), z.min(), z.max(), np.median(z), pz.min(), np.median(pz), pz.max(),
        np.median(pmod/pz), (pmod/pz).max(), time.time()-t0))
    np.savez(os.path.join(os.path.dirname(os.path.abspath(__file__)), "tab_%g.npz" % ms),
             z=z, pz=pz, e=e, px=px, py=py, A=A, chan=np.asarray(tab['chan']),
             counts=str(tab.get('counts')))
