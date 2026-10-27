# -*- coding: utf-8 -*-
"""Validate correct_ps.nbody for n=2,3,4 and compare cost with the production nbody."""
import sys, os, math, time, random
import numpy as np
RUN = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, os.path.join(RUN,"work","analysis"))
sys.path.insert(0, os.path.join(RUN,"work","widths"))
sys.path.insert(0, os.path.join(RUN,"work","response"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bc5_lib as L
import correct_ps as C
import correct_ps
M_PI = 0.13957039
tr = np.trapezoid if hasattr(np,'trapezoid') else np.trapz

def lam_arr(s, a2, b2):
    return (s-(np.sqrt(a2)+np.sqrt(b2))**2)*(s-(np.sqrt(a2)-np.sqrt(b2))**2)

def mom_check(d):
    E = sum(p[1] for p in d); P = [sum(p[i] for p in d) for i in (2,3,4)]
    m2 = [p[1]**2-p[2]**2-p[3]**2-p[4]**2 for p in d]
    return E, math.sqrt(sum(x*x for x in P)), m2

rng = random.Random(24680)

print("="*80)
print("(1) n=4 massless: <E_i/M>, permutation symmetry, 4-momentum conservation")
print("="*80)
N = 200000
Es = np.zeros((N,4)); Emax=0.0; Econs=0.0; Pcons=0.0
for k in range(N):
    d = C.nbody(1.0, [0.0,0.0,0.0,0.0], rng)
    for i in range(4): Es[k,i] = d[i][1]
    E,P,m2 = mom_check(d)
    Econs = max(Econs, abs(E-1.0)); Pcons = max(Pcons, P)
print("  new sampler  N=%d" % N)
for i in range(4):
    print("    <E_%d/M> = %.6f +- %.6f  (exact 0.25)" % (i+1, Es[:,i].mean(), Es[:,i].std(ddof=1)/math.sqrt(N)))
print("    max |sum E - M| = %.3e ; max |sum p| = %.3e" % (Econs, Pcons))
print("    pair products <E_i E_j>/M^2 (all must be equal by permutation symmetry):")
for i in range(4):
    print("      " + "  ".join("%.4f" % (Es[:,i]*Es[:,j]).mean() for j in range(4)))

print()
print("  production nbody same test, N=30000:")
Np=30000; Ep=np.zeros((Np,4))
for k in range(Np):
    d = L.nbody(1.0, [0.0,0.0,0.0,0.0], rng)
    for i in range(4): Ep[k,i]=d[i][1]
for i in range(4):
    print("    <E_%d/M> = %.6f +- %.6f" % (i+1, Ep[:,i].mean(), Ep[:,i].std(ddof=1)/math.sqrt(Np)))

print()
print("="*80)
print("(2) n=3 massive: s23 density vs exact, 3 pions M=1")
print("="*80)
ms=1.0; m=M_PI
smin=(2*m)**2; smax=(ms-m)**2
sg=np.linspace(smin,smax,40001)
f=np.sqrt(np.maximum(lam_arr(ms*ms,m*m,sg),0.0))*np.sqrt(np.maximum(lam_arr(sg,m*m,m*m),0.0))/sg
exact_mean=tr(f*sg,sg)/tr(f,sg)
print("  exact <s23> = %.6f GeV^2" % exact_mean)
for tag, sampler in (("new", C.nbody), ("production", L.nbody)):
    Ns=150000; vals=[]; none=0
    for k in range(Ns):
        d = sampler(ms,[m,m,m],rng)
        if d is None: none+=1; continue
        p2,p3=d[1],d[2]
        vals.append((p2[1]+p3[1])**2-sum((p2[j]+p3[j])**2 for j in (2,3,4)))
    v=np.array(vals)
    print("  %-10s <s23> = %.6f +- %.6f   ratio to exact %.5f  (None %d)" % (
        tag, v.mean(), v.std(ddof=1)/math.sqrt(len(v)), v.mean()/exact_mean, none))
    # chi2 against exact density in 20 equal-probability bins of the exact density
    cdf = np.cumsum(f); cdf = cdf/cdf[-1]
    edges = np.interp(np.linspace(0,1,21), cdf, sg)
    cnt,_ = np.histogram(v, bins=edges)
    exp = len(v)/20.0
    chi2 = ((cnt-exp)**2/exp).sum()
    print("               chi2/ndf vs exact density (20 bins) = %.2f / 19" % chi2)
print()
print("  massless 3-body: exact <x>=1/3, <x^2>=1/6, x=s23/M^2")
for tag, sampler in (("new", C.nbody), ("production", L.nbody)):
    Ns=150000; xs=[]
    for k in range(Ns):
        d=sampler(1.0,[0.0,0.0,0.0],rng)
        p2,p3=d[1],d[2]
        s=(p2[1]+p3[1])**2-sum((p2[j]+p3[j])**2 for j in (2,3,4))
        xs.append(s)
    x=np.array(xs)
    print("  %-10s <x>=%.6f+-%.6f  <x^2>=%.6f+-%.6f" % (tag, x.mean(), x.std(ddof=1)/math.sqrt(len(x)),
          (x*x).mean(), (x*x).std(ddof=1)/math.sqrt(len(x))))

print()
print("="*80)
print("(3) n=4 massive: conservation, on-shell, s12 marginal vs exact")
print("="*80)
ms=2.0; masses=[0.5,0.4,0.3,0.2]
N=200000
e1s=[]; s12s=[]; Emax=0.0; Pmax=0.0; m2err=0.0
for k in range(N):
    d=C.nbody(ms,masses,rng)
    E,P,m2=mom_check(d)
    Emax=max(Emax,abs(E-ms)); Pmax=max(Pmax,P)
    for i in range(4): m2err=max(m2err,abs(m2[i]-masses[i]**2))
    s12=(d[0][1]+d[1][1])**2-sum((d[0][j]+d[1][j])**2 for j in (2,3,4))
    s12s.append(s12); e1s.append(d[0][1])
print("  max|sumE-ms|=%.2e max|sum p|=%.2e max|m^2-m_nominal^2|=%.2e" % (Emax,Pmax,m2err))
# exact marginal of s12
def marg4(s12, ms=ms, masses=masses):
    m1,m2,m3,m4=masses
    lo34=(m3+m4)**2; hi34=max(lo34,(ms-math.sqrt(s12))**2)
    if hi34<=lo34: return 0.0
    s34=np.linspace(lo34,hi34,6001)
    fo=np.sqrt(np.maximum(lam_arr(ms*ms,s12,s34),0.0))/ms**2
    f12=np.sqrt(np.maximum(lam_arr(s12,m1*m1,m2*m2),0.0))/s12
    f34=np.sqrt(np.maximum(lam_arr(s34,m3*m3,m4*m4),0.0))/s34
    return tr(fo*f12*f34,s34)
s12s=np.array(s12s)
s12min=(0.5+0.4)**2; s12max=(ms-0.3-0.2)**2
grid=np.linspace(s12min+1e-4, s12max-1e-4, 201)
fv=np.array([marg4(s) for s in grid])
cdf=np.cumsum(fv); cdf=cdf/cdf[-1]
edges=np.interp(np.linspace(0,1,21),cdf,grid)
cnt,_=np.histogram(s12s,bins=edges); exp=N/20.0
print("  s12 chi2/ndf vs exact marginal (20 bins) = %.2f / 19" % (((cnt-exp)**2/exp).sum()))
# also s34 marginal (clusters 3,4) using the same formula by symmetry of the code (not of masses)
def marg34(s34, ms=ms, masses=masses):
    m1,m2,m3,m4=masses
    lo12=(m1+m2)**2; hi12=max(lo12,(ms-math.sqrt(s34))**2)
    if hi12<=lo12: return 0.0
    s12=np.linspace(lo12,hi12,6001)
    fo=np.sqrt(np.maximum(lam_arr(ms*ms,s12,s34),0.0))/ms**2
    f12=np.sqrt(np.maximum(lam_arr(s12,m1*m1,m2*m2),0.0))/s12
    f34=np.sqrt(np.maximum(lam_arr(s34,m3*m3,m4*m4),0.0))/s34
    return tr(fo*f12*f34,s12)
s34s=[]
for k in range(20000):
    d=C.nbody(ms,masses,rng)
    s34s.append((d[2][1]+d[3][1])**2-sum((d[2][j]+d[3][j])**2 for j in (2,3,4)))
s34s=np.array(s34s)
s34min=(0.3+0.2)**2; s34max=(ms-0.5-0.4)**2
grid2=np.linspace(s34min+1e-4,s34max-1e-4,201)
fv2=np.array([marg34(s) for s in grid2]); cdf2=np.cumsum(fv2); cdf2=cdf2/cdf2[-1]
edges2=np.interp(np.linspace(0,1,21),cdf2,grid2)
cnt2,_=np.histogram(s34s,bins=edges2); exp2=len(s34s)/20.0
print("  s34 chi2/ndf vs exact marginal (20 bins) = %.2f / 19" % (((cnt2-exp2)**2/exp2).sum()))

print()
print("  n=4 equal massive 0.3 each, M=2: permutation symmetry <E_i>")
Es=np.zeros((50000,4))
for k in range(50000):
    d=C.nbody(2.0,[0.3]*4,rng)
    for i in range(4): Es[k,i]=d[i][1]
for i in range(4):
    print("    <E_%d> = %.6f +- %.6f" % (i+1, Es[:,i].mean(), Es[:,i].std(ddof=1)/math.sqrt(50000)))

print()
print("="*80)
print("(4) acceptance and cost per call")
print("="*80)
# acceptance of the 4-body rejection (same proposal, count trials)
def acc4(ms,masses,rng,n=200000):
    m1,m2,m3,m4=masses
    lo12=(m1+m2)**2; hi12=(ms-m3-m4)**2
    lo34=(m3+m4)**2; hi34=(ms-m1-m2)**2
    ms2=ms*ms; tot=0
    for _ in range(n):
        while True:
            tot+=1
            s12=lo12+(hi12-lo12)*rng.random(); s34=lo34+(hi34-lo34)*rng.random()
            l12=C._lam(s12,m1*m1,m2*m2)
            if l12<=0: continue
            l34=C._lam(s34,m3*m3,m4*m4)
            if l34<=0: continue
            lout=C._lam(ms2,s12,s34)
            if lout<=0: continue
            if (math.sqrt(l12)/s12)*(math.sqrt(l34)/s34)*(math.sqrt(lout)/ms2) >= rng.random(): break
    return n/tot
for (ms,masses) in ((1.0,[0.0]*4),(2.0,[0.5,0.4,0.3,0.2]),(2.0,[0.3]*4),(5.0,[1.0,1.0,1.0,1.0])):
    a=acc4(ms,masses,rng,30000)
    print("   4-body M=%.1f masses=%s acceptance=%.3f" % (ms,masses,a))

def bench(fn, ms, masses, n=60000):
    for _ in range(2000): fn(ms,masses,rng)
    t0=time.perf_counter()
    for _ in range(n): fn(ms,masses,rng)
    return 1e6*(time.perf_counter()-t0)/n
print("  cost per call (us):")
print("    n=3 pions M=1     production %.2f   new %.2f" % (bench(L.nbody,1.0,[M_PI]*3), bench(C.nbody,1.0,[M_PI]*3)))
print("    n=4 massless M=1  production %.2f   new %.2f" % (bench(L.nbody,1.0,[0.0]*4), bench(C.nbody,1.0,[0.0]*4)))
print("    n=4 massive  M=2  production %.2f   new %.2f" % (bench(L.nbody,2.0,[0.5,0.4,0.3,0.2]), bench(C.nbody,2.0,[0.5,0.4,0.3,0.2])))
