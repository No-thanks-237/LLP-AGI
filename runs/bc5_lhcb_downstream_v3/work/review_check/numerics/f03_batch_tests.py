# -*- coding: utf-8 -*-
"""Validate and benchmark the vectorised exact sampler (batch_ps)."""
import sys, os, math, time, random
import numpy as np
sys.path.insert(0, r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3/work/review_check/numerics")
import correct_ps as C, batch_ps as B
M_PI=0.13957039
tr = np.trapezoid if hasattr(np,'trapezoid') else np.trapz
def lam_arr(s,a2,b2):
    return (s-(np.sqrt(a2)+np.sqrt(b2))**2)*(s-(np.sqrt(a2)-np.sqrt(b2))**2)
rng = np.random.default_rng(20261024)

def mom(d):
    E,PX,PY,PZ = d
    return E, PX, PY, PZ

print("="*78)
print("(1) batch n=3, 3 pions M=1: <s23> vs exact 0.352813")
print("="*78)
d = B.nbody_batch(1.0, [M_PI]*3, 400000, rng)
E,X,Y,Z = d
s23 = (E[:,1]+E[:,2])**2 - (X[:,1]+X[:,2])**2 - (Y[:,1]+Y[:,2])**2 - (Z[:,1]+Z[:,2])**2
print("  <s23> = %.6f +- %.6f  ratio %.5f  (exact 0.352813)" % (
    s23.mean(), s23.std(ddof=1)/math.sqrt(len(s23)), s23.mean()/0.352813))
# 4-momentum conservation
print("  max|sum E - M| = %.2e  max|sum p| = %.2e" % (
    np.abs(E.sum(1)-1.0).max(), np.sqrt((X.sum(1))**2+(Y.sum(1))**2+(Z.sum(1))**2).max()))
ms=1.0; m=M_PI; sg=np.linspace((2*m)**2,(ms-m)**2,40001)
f=np.sqrt(np.maximum(lam_arr(ms*ms,m*m,sg),0))*np.sqrt(np.maximum(lam_arr(sg,m*m,m*m),0))/sg
cdf=np.cumsum(f); cdf/=cdf[-1]
edges=np.interp(np.linspace(0,1,21),cdf,sg)
cnt,_=np.histogram(s23,bins=edges); exp=len(s23)/20.
print("  chi2/ndf vs exact density = %.2f/19" % (((cnt-exp)**2/exp).sum()))

print()
print("="*78)
print("(2) batch n=4 massless: <E_i/M> = 1/4, permutation symmetry")
print("="*78)
d = B.nbody_batch(1.0, [0.0]*4, 400000, rng)
E,X,Y,Z = d
for i in range(4):
    print("  <E_%d/M> = %.6f +- %.6f" % (i+1, E[:,i].mean(), E[:,i].std(ddof=1)/math.sqrt(len(E))))
print("  max|sum E - M| = %.2e  max|sum p| = %.2e" % (
    np.abs(E.sum(1)-1.0).max(), np.sqrt((X.sum(1))**2+(Y.sum(1))**2+(Z.sum(1))**2).max()))
print("  <E1 E4> = %.5f  <E2 E3> = %.5f  <E1 E2> = %.5f" % (
    (E[:,0]*E[:,3]).mean(), (E[:,1]*E[:,2]).mean(), (E[:,0]*E[:,1]).mean()))

print()
print("="*78)
print("(3) batch n=4 massive M=2 (0.5,0.4,0.3,0.2): conservation + s12 marginal")
print("="*78)
d = B.nbody_batch(2.0, [0.5,0.4,0.3,0.2], 400000, rng)
E,X,Y,Z = d
print("  max|sum E - 2| = %.2e  max|sum p| = %.2e" % (
    np.abs(E.sum(1)-2.0).max(), np.sqrt((X.sum(1))**2+(Y.sum(1))**2+(Z.sum(1))**2).max()))
s12 = (E[:,0]+E[:,1])**2-(X[:,0]+X[:,1])**2-(Y[:,0]+Y[:,1])**2-(Z[:,0]+Z[:,1])**2
m2 = E**2-X**2-Y**2-Z**2
print("  max|m_i^2 - nominal| = %.2e" % np.abs(m2-np.array([0.25,0.16,0.09,0.04])).max())
def marg4(s12, ngrid=8001):
    m1,m2_,m3,m4=0.5,0.4,0.3,0.2
    lo34=(m3+m4)**2; hi34=(2.0-math.sqrt(s12))**2
    if hi34<=lo34: return 0.0
    t=np.linspace(0,1,ngrid); s34=lo34+(hi34-lo34)*t*t; jac=2*(hi34-lo34)*t
    fo=np.sqrt(np.maximum(lam_arr(4.0,s12,s34),0.0))/4.0
    f12=np.sqrt(np.maximum(lam_arr(s12,m1*m1,m2_*m2_),0.0))/s12
    f34=np.sqrt(np.maximum(lam_arr(s34,m3*m3,m4*m4),0.0))/s34
    return tr(fo*f12*f34*jac,t)
grid=np.linspace((0.5+0.4)**2+1e-6,(2.0-0.3-0.2)**2-1e-6,1501)
fv=np.array([marg4(s) for s in grid]); cdf=np.cumsum(fv); cdf/=cdf[-1]
edges=np.interp(np.linspace(0,1,21),cdf,grid)
cnt,_=np.histogram(s12,bins=edges); exp=len(s12)/20.
print("  s12 chi2/ndf vs exact marginal = %.2f/19" % (((cnt-exp)**2/exp).sum()))
print("  exact <s12> = %.6f  MC = %.6f +- %.6f" % (
    tr(fv*grid,grid)/tr(fv,grid), s12.mean(), s12.std(ddof=1)/math.sqrt(len(s12))))

print()
print("="*78)
print("(4) speed: batch cost per event vs the production per-call sampler")
print("="*78)
def bench_batch(fn, ms, masses, n, rng):
    fn(ms,masses,2000,rng)
    t0=time.perf_counter(); fn(ms,masses,n,rng); dt=time.perf_counter()-t0
    return 1e6*dt/n
for (tag, ms, masses) in (("n=3 pions M=1",1.0,[M_PI]*3),
                          ("n=4 massless M=1",1.0,[0.0]*4),
                          ("n=4 massive M=2",2.0,[0.5,0.4,0.3,0.2])):
    tb = bench_batch(B.nbody_batch, ms, masses, 1000000, rng)
    print("  batch %-18s %.3f us/event  (1e7 decays: %.1f s)" % (tag, tb, tb*1e7/1e6))
print()
# production cost (pure python per call)
import importlib.util
sys.path.insert(0, r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3/work/analysis")
import bc5_lib as L
rng2 = random.Random(7)
def bench_py(fn, ms, masses, n=50000):
    for _ in range(2000): fn(ms,masses,rng2)
    t0=time.perf_counter()
    for _ in range(n): fn(ms,masses,rng2)
    return 1e6*(time.perf_counter()-t0)/n
for (tag, ms, masses) in (("n=3 pions M=1",1.0,[M_PI]*3),
                          ("n=4 massless M=1",1.0,[0.0]*4),
                          ("n=4 massive M=2",2.0,[0.5,0.4,0.3,0.2])):
    prod = bench_py(L.nbody, ms, masses); new = bench_py(C.nbody, ms, masses)
    print("  per-call %-18s production %6.2f us   new scalar %6.2f us  (ratio %.2f)" % (tag, prod, new, prod/new))
