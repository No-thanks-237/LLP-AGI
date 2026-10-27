
import sys, os, math, time, random
import numpy as np
RUN = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, os.path.join(RUN,"work","analysis"))
sys.path.insert(0, os.path.join(RUN,"work","widths"))
sys.path.insert(0, os.path.join(RUN,"work","response"))
import bc5_lib as L

M_PI = 0.13957039
M = 1.0
rng = random.Random(12345)

def fourmomenta(parts):
    E = sum(p[1] for p in parts); P = [sum(p[i] for p in parts) for i in (2,3,4)]
    return E, P

# ---------- (a) massless 4-body: mean E_i / M ----------
N = 30000
Es = np.zeros((N,4))
bad = 0
for k in range(N):
    d = L.nbody(M, [0.0,0.0,0.0,0.0], rng)
    if d is None:
        bad += 1; continue
    for i in range(4): Es[k,i] = d[i][1]
print("massless 4-body current nbody: N=%d None=%d" % (N, bad))
for i in range(4):
    print("  <E_%d/M> = %.6f +- %.6f   (1/4 = 0.25)" % (i+1, Es[:,i].mean()/M, Es[:,i].std(ddof=1)/math.sqrt(N)/M))
# momentum conservation and massless check
d = L.nbody(M, [0.0]*4, rng)
E,P = fourmomenta(d)
print("  one event: E=%.12f |P|=%.3e" % (E, math.sqrt(sum(x*x for x in P))))

# ---------- (c) permutation test massless 4-body ----------
# relabel identical particles: feed masses in a permuted order; sampling must be symmetric
# test: for identical masses, the marginal of E_1 must equal that of E_4 (already above);
# stronger: mean E_2 == mean E_3 etc. Also check a *different* ordering of identical masses
# gives the same joint law -> compare <E_1 E_2> vs <E_1 E_4>
print("  <E1 E2>/M^2 = %.6f  <E1 E4>/M^2 = %.6f  <E2 E3>/M^2 = %.6f" % (
    (Es[:,0]*Es[:,1]).mean()/M**2, (Es[:,0]*Es[:,3]).mean()/M**2, (Es[:,1]*Es[:,2]).mean()/M**2))

# ---------- (b) massive 3-body: <s23> vs exact ----------
# exact: dPhi3/ds ~ sqrt(lam(M^2,m1^2,s))*sqrt(lam(s,m2^2,m3^2))/s,  s = s23
def lam(a,b,c):
    return (a-(math.sqrt(b)+math.sqrt(c))**2)*(a-(math.sqrt(b)-math.sqrt(c))**2)
ms = M; m = M_PI
smin = (2*m)**2; smax = (ms-m)**2
sg = np.linspace(smin, smax, 20001)
f = np.sqrt(np.maximum(lam(ms*ms, m*m, sg),0.0))*np.sqrt(np.maximum(lam(sg, m*m, m*m),0.0))/sg
I0 = np.trapezoid(f, sg) if hasattr(np,'trapezoid') else np.trapz(f,sg)
I1 = np.trapezoid(f*sg, sg) if hasattr(np,'trapezoid') else np.trapz(f*sg,sg)
print("exact <s23> (3 pions, M=1) = %.6f GeV^2" % (I1/I0))

Ns3 = 200000
s23 = []
none3 = 0
for k in range(Ns3):
    d = L.nbody(ms, [m,m,m], rng)
    if d is None:
        none3 += 1; continue
    p2 = d[1]; p3 = d[2]
    s = (p2[1]+p3[1])**2 - (p2[2]+p3[2])**2 - (p2[3]+p3[3])**2 - (p2[4]+p3[4])**2
    s23.append(s)
s23 = np.array(s23)
print("code  <s23> = %.6f +- %.6f  (None=%d/%d)" % (s23.mean(), s23.std(ddof=1)/math.sqrt(len(s23)), none3, Ns3))
print("ratio code/exact = %.5f" % (s23.mean()/(I1/I0)))
