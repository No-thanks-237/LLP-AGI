
import sys, os, math, random
import numpy as np
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\analysis")
import bc5_production as P, bc5_lib as L
rng = random.Random(3)
ms, mX, mXp, ch = 1.0, 5.27934, 0.493677, 'K'
qs = np.array([P.sample_dipion_mass(ms, mX, mXp, ch, rng) for _ in range(4000)])
lo, hi = 4*ms*ms, (mX-mXp)**2
grid = np.linspace(lo, hi, 400)
w = np.array([P._q2_weight(ms, mX, mXp, ch, q) for q in grid])
w /= np.trapezoid(w, grid)
# compare quantiles of the sampled sqrt(q2)^2 distribution with the exact density
print("sampled mean q2 = %.4f   exact %.4f" % ((qs**2).mean(), np.trapezoid(grid*w, grid)))
for frac in (0.1, 0.5, 0.9):
    cw = np.cumsum(w*(grid[1]-grid[0])); cw /= cw[-1]
    print("  %.1f quantile: sampled %.4f  exact %.4f" % (frac, np.quantile(qs**2, frac),
          float(np.interp(frac, cw, grid))))
# multi composition
from collections import Counter
cnt = Counter()
for _ in range(3000):
    d = L.s_decay_products(2.0, 'multi', rng)
    cnt[len([a for a in d if a[5] != 0])] += 1
print("multi charged multiplicity at 2 GeV:", dict(sorted(cnt.items())), "(expect 4: %.3f, 2: %.3f, 0: %.3f)" % (4/9, 4/9, 1/9))
c2 = Counter()
for _ in range(3000):
    d = L.s_decay_products(0.2, 'pipi', rng)
    c2[len(d)] += 1
print("pipi at 0.2 GeV (below threshold -> empty):", dict(c2))
