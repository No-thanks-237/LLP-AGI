# -*- coding: utf-8 -*-
import sys, os, math, random
import numpy as np
RUN = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, os.path.join(RUN,"work","analysis")); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import correct_ps as C
M_PI=0.13957039
tr = np.trapezoid if hasattr(np,'trapezoid') else np.trapz
def lam_arr(s,a2,b2):
    return (s-(np.sqrt(a2)+np.sqrt(b2))**2)*(s-(np.sqrt(a2)-np.sqrt(b2))**2)
rng=random.Random(999)

ms=1.0; m=M_PI; m1=m2=m3=m
smin=(m2+m3)**2; smax=(ms-m1)**2
# fine gmax
sg=np.linspace(smin,smax,200001)
g=np.sqrt(np.maximum(lam_arr(ms*ms,m1*m1,sg),0.0))*np.sqrt(np.maximum(lam_arr(sg,m2*m2,m3*m3),0.0))/sg
print("gmax fine=%.10e  at s=%.5f" % (g.max(), sg[np.argmax(g)]))
sminc,smaxc,gmaxc=C._g3max(ms,m1,m2,m3)
print("gmax cached=%.10e  ratio cached/fine=%.6f" % (gmaxc, gmaxc/g.max()))
# exact mean and chi2 with production-density-free test: weighted integration check
I0=tr(g,sg); I1=tr(g*sg,sg); exact=I1/I0
print("exact <s>=%.6f" % exact)
N=400000
vals=[]
for k in range(N):
    d=C.nbody(ms,[m,m,m],rng); p2,p3=d[1],d[2]
    vals.append((p2[1]+p3[1])**2-sum((p2[j]+p3[j])**2 for j in (2,3,4)))
v=np.array(vals)
print("new <s>=%.6f +- %.6f  ratio %.5f" % (v.mean(), v.std(ddof=1)/math.sqrt(N), v.mean()/exact))
# KS test against exact CDF
cdf=np.cumsum(g); cdf=cdf/cdf[-1]
vs=np.sort(v)
F=np.interp(vs, sg, cdf)
n=len(vs)
D=max(np.max(np.arange(1,n+1)/n - F), np.max(F - np.arange(n)/n))
print("KS D=%.5f  (95%% crit ~ %.5f)" % (D, 1.36/math.sqrt(n)))
# equal-probability bins chi2 with more events
edges=np.interp(np.linspace(0,1,21),cdf,sg)
cnt,_=np.histogram(v,bins=edges); exp=N/20.
print("chi2/ndf=%.2f/19" % (((cnt-exp)**2/exp).sum()))
print("bin counts/expected:", np.round(cnt/exp,3))

# massless 3-body
for N2 in (400000,):
    xs=[]
    for k in range(N2):
        d=C.nbody(1.0,[0.,0.,0.],rng); p2,p3=d[1],d[2]
        xs.append((p2[1]+p3[1])**2-sum((p2[j]+p3[j])**2 for j in (2,3,4)))
    x=np.array(xs)
    print("massless: <x>=%.6f+-%.6f (1/3)  <x^2>=%.6f+-%.6f (1/6)" % (
      x.mean(), x.std(ddof=1)/math.sqrt(N2), (x*x).mean(), (x*x).std(ddof=1)/math.sqrt(N2)))
    # chi2 against 2(1-x)
    ed=np.linspace(0,1,21); cnt,_=np.histogram(x,bins=ed)
    exp=N2*(ed[1:]**2-ed[:-1]**2) - N2*2*(ed[1:]**3-ed[:-1]**3)/3  # integral of 2(1-x)
    exp=np.diff(2*ed-ed**2)*N2/1.0
    print("  massless chi2/ndf=%.2f/19" % (((cnt-exp)**2/exp).sum()))
