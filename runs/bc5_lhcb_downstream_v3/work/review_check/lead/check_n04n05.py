
import sys, os, math, random
import numpy as np
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\analysis")
import bc5_lib as L
rng = random.Random(11)
# N04: near-threshold four-body success rate
print("N04: four-body failure rate near threshold (1000 draws each)")
for ms in (0.558392, 0.56, 0.565616, 0.58, 0.6, 0.661654):
    ok = 0
    for _ in range(1000):
        try:
            d = L.nbody(ms, [L.M_PI]*4, rng)
            ok += 1 if d else 0
        except RuntimeError:
            ok += 0
    print("   m=%.6f  success %.1f%%" % (ms, ok/10.0))
# conservation and symmetry still hold
arr = np.array([L.nbody(2.0, [0.14]*4, rng) for _ in range(20000)])
print("   four-body <E/M> at 2 GeV:", np.round(arr[:,:,1].mean(axis=0)/2.0, 4))
# N05: lifetime slope term vs high-precision
import decimal
def exact_pwl(alpha, z0, A, z):
    # 70-digit reference with mpmath-free Decimal series on each cell
    decimal.getcontext().prec = 50
    tot = decimal.Decimal(0)
    K = len(z)
    for k in range(K-1):
        z1, z2 = decimal.Decimal(z[k]), decimal.Decimal(z[k+1])
        lo = max(z1, decimal.Decimal(z0)); d = z2 - lo
        if d <= 0: continue
        s = (decimal.Decimal(A[k+1])-decimal.Decimal(A[k]))/(z2-z1)
        a = decimal.Decimal(alpha)
        f1 = (-a*(lo-decimal.Decimal(z0))).exp()
        A_lo = decimal.Decimal(A[k]) + s*(lo-z1)
        w0 = (1-(-a*d).exp())/a
        w1 = (1-(1+a*d)*(-a*d).exp())/(a*a)
        tot += f1*(A_lo*w0 + s*w1)
    return float(tot)
z = np.linspace(1.0, 2.5, 13)
for alpha in (8.6438826e-8, 1e-6, 1e-3):
    for name, A in (("A=(z-1)/1.5", [(x-1)/1.5 for x in z]),
                    ("ramp last cell", [0.0]*11 + [0.0, 1.0])):
        got = L.lifetime_probability(1.0, (0.0, 0.0, 0.0, 1.0/8.6438826e-8*1.0), 0.0, A, 1.0/(1.0*8.6438826e-8)) if False else None
        p = L._p_window(np.array([alpha]), np.array([0.0]), np.array([A]), z)[0]
        ex = exact_pwl(alpha, 0.0, A, list(z))
        print("   alpha=%.3e %-16s code=%.12e exact=%.12e ratio=%.9f" % (alpha, name, p, ex, p/ex))
