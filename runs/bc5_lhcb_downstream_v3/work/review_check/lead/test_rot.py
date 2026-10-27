
import sys, os, math, random
import numpy as np
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\analysis")
import bc5_lib as L
L.load_hadron_lib(r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\magnus\hadron_lib")
rng = random.Random(7)
print("%6s %6s %10s %10s" % ("channel","m_S","<cos2th_z>","se"))
for ch, ms in (("gg", 9.516888), ("bb", 28.527411), ("tau", 28.527411), ("tau", 62.49)):
    vals = []
    for _ in range(4000):
        d = L.s_decay_products(ms, ch, rng)
        for a in d:
            p2 = a[2]**2 + a[3]**2 + a[4]**2
            if p2 > 0:
                vals.append(a[4]**2/p2)
    v = np.array(vals)
    print("%6s %6.2f %10.4f %10.4f  (isotropic: 0.3333)" % (ch, ms, v.mean(), v.std()/math.sqrt(len(v))))
