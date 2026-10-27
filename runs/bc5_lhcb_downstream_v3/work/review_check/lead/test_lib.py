
import sys, os, math, random
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\analysis")
import numpy as np
import bc5_lib as L
rng = random.Random(1)
# F03 tests
arr = np.array([L.nbody(1.0, [0.0]*4, rng) for _ in range(20000)])
print("massless 4-body <E/M>:", np.round(arr[:,:,1].mean(axis=0), 4))
print("  se:", np.round(arr[:,:,1].std(axis=0)/math.sqrt(20000), 4))
M, m3 = 1.0, [0.13957039]*3
arr3 = np.array([L.nbody(M, m3, rng) for _ in range(20000)])
s23 = (M-arr3[:,0,1])**2 - np.sum(arr3[:,0,2:5]**2, axis=1)
print("massive 3-body <s23>:", round(s23.mean(),5), "+-", round(s23.std()/math.sqrt(20000),5), "(exact 0.352813)")
print("  <E>:", np.round(arr3[:,:,1].mean(axis=0),4), "(exact 0.3333)")
# 4-momentum conservation
for n,ms in ((3,[0.14]*3),(4,[0.14]*4)):
    d = L.nbody(2.0, ms, rng)
    tot = np.sum([a[1] for a in d]), np.sum([a[2] for a in d]), np.sum([a[3] for a in d]), np.sum([a[4] for a in d])
    print("conservation n=%d: E=%.12f p=(%.2e,%.2e,%.2e)" % (n, tot[0]-2.0, tot[1], tot[2], tot[3]))
# F10 test
for z0 in (0.0, 0.5, 1.01, 1.249, 2.4):
    for lam in (0.1, 0.3, 1.0, 10.0):
        ms, pz = 1.0, 10.0
        ctau = lam*ms/pz
        got = L.lifetime_probability(ms, (0.0, 0.0, 0.0, pz), z0, np.ones(len(L.ZGRID)), ctau)
        lo = max(1.0, z0)
        exact = math.exp(-(lo-z0)/lam) - math.exp(-(2.5-z0)/lam)
        print("z0=%.3f L=%.2f  code=%.10f exact=%.10f ratio=%.10f" % (z0, lam, got, exact, got/exact))
