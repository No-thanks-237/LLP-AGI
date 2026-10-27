
import sys, os, math, random
import numpy as np
RUN = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, os.path.join(RUN,"work","analysis"))
sys.path.insert(0, os.path.join(RUN,"work","widths"))
sys.path.insert(0, os.path.join(RUN,"work","response"))
import bc5_lib as L
M_PI = 0.13957039
rng = random.Random(54321)

def lam_arr(s, a2, b2):
    # lambda(s, a^2, b^2) with a=sqrt(a2), b=sqrt(b2)
    return (s-(np.sqrt(a2)+np.sqrt(b2))**2)*(s-(np.sqrt(a2)-np.sqrt(b2))**2)

def exact_s23(ms, m1, m2, m3):
    smin=(m2+m3)**2; smax=(ms-m1)**2
    sg=np.linspace(smin, smax, 40001)
    f=np.sqrt(np.maximum(lam_arr(ms*ms, m1*m1, sg),0.0))*np.sqrt(np.maximum(lam_arr(sg,m2*m2,m3*m3),0.0))/sg
    tr = np.trapezoid if hasattr(np,'trapezoid') else np.trapz
    return tr(f*sg,sg)/tr(f,sg)

for (ms,m1,m2,m3,tag) in [(1.0,M_PI,M_PI,M_PI,"3 pions M=1"),
                          (1.0,0.4,0.3,0.2,"mixed massive M=1"),
                          (2.0,0.5,0.5,0.5,"3 x 0.5 GeV M=2")]:
    ex = exact_s23(ms,m1,m2,m3)
    N=100000; vals=[]; none=0
    for k in range(N):
        d=L.nbody(ms,[m1,m2,m3],rng)
        if d is None: none+=1; continue
        p2=d[1]; p3=d[2]
        s=(p2[1]+p3[1])**2-sum((p2[j]+p3[j])**2 for j in (2,3,4))
        vals.append(s)
    v=np.array(vals)
    print("%-20s exact <s23>=%.6f  code=%.6f +- %.6f  ratio=%.5f  None=%d" % (
        tag, ex, v.mean(), v.std(ddof=1)/math.sqrt(len(v)), v.mean()/ex, none))
