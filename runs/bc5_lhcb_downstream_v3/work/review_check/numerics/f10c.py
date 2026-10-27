# -*- coding: utf-8 -*-
"""F10 detail: how many real scalars are produced inside [1,2.5] m, and the FY12 visible-rate effect."""
import sys, os, math, json
import numpy as np
RUN = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, os.path.join(RUN,"work","analysis")); sys.path.insert(0, os.path.join(RUN,"work","widths")); sys.path.insert(0, os.path.join(RUN,"work","response"))
import run_scan as R, bc5_lib as L, scalar_widths as sw
HERE = os.path.dirname(os.path.abspath(__file__))
Z=np.asarray(L.ZGRID); H=Z[1]-Z[0]
scan=json.load(open(os.path.join(RUN,"outputs","data","signal_scan.json")))
th=np.array(scan['theta2']); sig=scan['signal_per_fb_per_channel']
def crossings_scan(m, level):
    row=sig[str(m)]; s=np.array([r['B']+r['K']+r['H'] for r in row])
    lx=np.log10(th); ly=np.log10(np.maximum(s,1e-300)); lv=math.log10(level); out=[]
    for i in range(len(th)-1):
        f0=ly[i]-lv; f1=ly[i+1]-lv
        if f0*f1<0: t=f0/(f0-f1); out.append(10**(lx[i]+t*(lx[i+1]-lx[i])))
    return out
def trap_p(alpha, z0, A):
    dz=Z-z0; w=np.where(dz>0, alpha*np.exp(-np.clip(dz,0,None)*alpha),0.0)
    return float(np.trapezoid(w*A,Z) if hasattr(np,'trapezoid') else np.trapz(w*A,Z))
print("  rows with 1 < z0 < 2.5 m (birth point inside the window) and their weight")
for ms_ in (0.2802, 1.0, 3.5567):
    d=np.load(os.path.join(HERE,"tab_%g.npz"%ms_)); z0=d['z']; pz=d['pz']; A=d['A']; chan=d['chan']
    inside=(z0>1.0)&(z0<2.5)
    print("   m=%.4f: inside-window rows %d/%d (%.3f%%)  z0 range %s" % (
        ms_, inside.sum(), len(z0), 100*inside.mean(),
        ("%.3f..%.3f" % (z0[inside].min(), z0[inside].max())) if inside.any() else "-"))
    for L in (25.0, 300.0):
        cr=crossings_scan(ms_, 2.3/L)
        t2=min(cr); alpha=ms_*sw.gamma_total(ms_)*t2/(sw.HBAR_C*pz)
        p=np.array([trap_p(alpha[i],z0[i],A[i]) for i in range(len(z0))])
        print("     L=%3.0f long-lived: sum p = %.6e ; sum p(inside window) = %.3e ; ratio inside share = %.2e ; mean code/exact for inside rows = %.4f" % (
            L, p.sum(), p[inside].sum(), p[inside].sum()/max(p.sum(),1e-300), float(np.mean(p[inside]/np.maximum([1e-300]*inside.sum(),1e-300))) if inside.any() else float('nan')))
    # per channel fraction of rows outside window
    for c in sorted(set(chan.tolist())):
        sel=(chan==c)
        print("     chan %s: n=%4d  z0<1m %4d  inside %3d  z0>2.5m %4d" % (c, sel.sum(), int((sel&(z0<1)).sum()), int((sel&inside).sum()), int((sel&(z0>2.5)).sum())))
