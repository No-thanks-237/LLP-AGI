
import sys, os, math, random
from collections import defaultdict
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\analysis")
import bc5_lib as L
rng = random.Random(5)
s_p4 = (100.0, 0.0, 0.0, 99.99)   # a boosted S so the kaon flight is large
shared = 0; pairs = 0
for _ in range(4000):
    d = L.s_decay_products(1.5, 'KK', rng, s_p4=s_p4)
    if not d: continue
    offs = [a[6] for a in d]
    # daughters of one kaon were appended consecutively; group by identical offset object/value
    if len(d) >= 2:
        for i in range(len(d)):
            for j in range(i+1, len(d)):
                if offs[i] is not None and offs[i] == offs[j]:
                    shared += 1; pairs += 1
print("pairs sharing the same birth offset: %d / %d" % (shared, pairs))
# also: how large is the displacement?
big = [a[6][2] for _ in range(2000) for a in (L.s_decay_products(1.5, 'KK', rng, s_p4=s_p4) or []) if a[6]]
import numpy as np
v = np.array(big)
print("delta-z of kaon daughters: median |dz| = %.3f m, 90%% = %.2f m, max = %.1f m" % (
    np.median(np.abs(v)), np.quantile(np.abs(v), 0.9), np.abs(v).max()))
