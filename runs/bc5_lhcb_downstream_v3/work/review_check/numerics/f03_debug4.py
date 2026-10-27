
import sys, os, math, random
import numpy as np
sys.path.insert(0, r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3/work/review_check/numerics")
import correct_ps as C
tr = np.trapezoid if hasattr(np,'trapezoid') else np.trapz
def lam_arr(s,a2,b2):
    return (s-(np.sqrt(a2)+np.sqrt(b2))**2)*(s-(np.sqrt(a2)-np.sqrt(b2))**2)
rng=random.Random(111)
ms=2.0; masses=[0.5,0.4,0.3,0.2]; m1,m2,m3,m4=masses
def marg4(s12, ngrid=20001):
    lo34=(m3+m4)**2; hi34=(ms-math.sqrt(s12))**2
    if hi34<=lo34: return 0.0
    # substitution s34 = lo34 + (hi34-lo34)*t^2 handles sqrt endpoints
    t=np.linspace(0,1,ngrid)
    s34=lo34+(hi34-lo34)*t*t
    jac=2*(hi34-lo34)*t
    fo=np.sqrt(np.maximum(lam_arr(ms*ms,s12,s34),0.0))/ms**2
    f12=np.sqrt(np.maximum(lam_arr(s12,m1*m1,m2*m2),0.0))/s12
    f34=np.sqrt(np.maximum(lam_arr(s34,m3*m3,m4*m4),0.0))/s34
    return tr(fo*f12*f34*jac,t)
s12min=(m1+m2)**2; s12max=(ms-m3-m4)**2
N=200000
s12s=[]
for k in range(N):
    d=C.nbody(ms,masses,rng)
    s12s.append((d[0][1]+d[1][1])**2-sum((d[0][j]+d[1][j])**2 for j in (2,3,4)))
s12s=np.array(s12s)
grid=np.linspace(s12min+1e-6,s12max-1e-6,5001)
fv=np.array([marg4(s) for s in grid])
cdf=np.cumsum(fv); cdf=cdf/cdf[-1]
edges=np.interp(np.linspace(0,1,21),cdf,grid)
cnt,_=np.histogram(s12s,bins=edges); exp=N/20.
print("s12 bin counts / expected (finer marg grid):", np.round(cnt/exp,3))
print("chi2=%.2f/19" % (((cnt-exp)**2/exp).sum()))
print("exact <s12> from marg = %.6f ; MC = %.6f +- %.6f" % (
   tr(fv*grid,grid)/tr(fv,grid), s12s.mean(), s12s.std(ddof=1)/math.sqrt(N)))
# cross-check: exact marginal computed the other way -- integrate over s12 for s34
def marg34(s34, ngrid=20001):
    lo12=(m1+m2)**2; hi12=(ms-math.sqrt(s34))**2
    if hi12<=lo12: return 0.0
    t=np.linspace(0,1,ngrid); s12=lo12+(hi12-lo12)*t*t; jac=2*(hi12-lo12)*t
    fo=np.sqrt(np.maximum(lam_arr(ms*ms,s12,s34),0.0))/ms**2
    f12=np.sqrt(np.maximum(lam_arr(s12,m1*m1,m2*m2),0.0))/s12
    f34=np.sqrt(np.maximum(lam_arr(s34,m3*m3,m4*m4),0.0))/s34
    return tr(fo*f12*f34*jac,t)
s34min=(m3+m4)**2; s34max=(ms-m1-m2)**2
N2=200000
s34s=np.array([ (lambda d: (d[2][1]+d[3][1])**2-sum((d[2][j]+d[3][j])**2 for j in (2,3,4)))(C.nbody(ms,masses,rng)) for _ in range(N2)])
grid2=np.linspace(s34min+1e-6,s34max-1e-6,5001)
fv2=np.array([marg34(s) for s in grid2]); cdf2=np.cumsum(fv2); cdf2=cdf2/cdf2[-1]
edges2=np.interp(np.linspace(0,1,21),cdf2,grid2)
cnt2,_=np.histogram(s34s,bins=edges2); exp2=N2/20.
print("s34 bin counts / expected:", np.round(cnt2/exp2,3))
print("s34 chi2=%.2f/19" % (((cnt2-exp2)**2/exp2).sum()))
print("exact <s34> = %.6f ; MC = %.6f +- %.6f" % (tr(fv2*grid2,grid2)/tr(fv2,grid2), s34s.mean(), s34s.std(ddof=1)/math.sqrt(N2)))
