
import sys, math
import numpy as np
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\analysis")
import bc5_lib as L
ok = True
for z0 in (0.0, 0.5, 1.01, 1.249, 2.4):
    for lam in (0.1, 0.3, 1.0, 10.0):
        ms, pz = 1.0, 10.0
        got = L.lifetime_probability(ms, (0.0, 0.0, 0.0, pz), z0, np.ones(len(L.ZGRID)), lam*ms/pz)
        ex = math.exp(-(max(1.0, z0)-z0)/lam) - math.exp(-(2.5-z0)/lam)
        if abs(got/ex - 1) > 1e-9:
            ok = False
            print("FAIL z0=%.3f L=%.2f code=%.10f exact=%.10f" % (z0, lam, got, ex))
print("analytic lifetime test:", "PASS" if ok else "FAIL")
# four-body + rotation checks after the edit
import random
rng = random.Random(3)
arr = np.array([L.nbody(1.0, [0.0]*4, rng) for _ in range(4000)])
print("massless 4-body <E/M>:", np.round(arr[:,:,1].mean(axis=0), 3))
L.load_hadron_lib(r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\magnus\hadron_lib")
vals = []
for _ in range(1500):
    for a in L.s_decay_products(9.516888, 'gg', rng):
        p2 = a[2]**2+a[3]**2+a[4]**2
        if p2 > 0: vals.append(a[4]**2/p2)
print("<cos2th_z> gg@9.52: %.4f (isotropic 0.3333)" % np.mean(vals))
