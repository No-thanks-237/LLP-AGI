
import sys, math, random
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\analysis")
import bc5_lib as L
rng = random.Random(9)
thr = 4*L.M_PI
print("4 m_pi = %.8f" % thr)
for eps in (1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2):
    ms = thr*(1+eps)
    ok = fail = nones = 0
    for _ in range(300):
        try:
            d = L.nbody(ms, [L.M_PI]*4, rng)
            if d: ok += 1
            else: nones += 1
        except RuntimeError:
            fail += 1
    print("m/4mpi-1 = %.0e  success %3d  None %3d  raise %3d  wmax=%.3e" % (
        eps, ok, nones, fail, L._g4_envelope(ms, [L.M_PI]*4)[4]))
