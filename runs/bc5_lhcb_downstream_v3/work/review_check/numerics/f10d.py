# -*- coding: utf-8 -*-
import sys, os, math, json
import numpy as np
RUN=r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, os.path.join(RUN,"work","analysis")); sys.path.insert(0, os.path.join(RUN,"work","widths")); sys.path.insert(0, os.path.join(RUN,"work","response"))
import run_scan as R, bc5_lib as L, scalar_widths as sw
HERE=os.path.dirname(os.path.abspath(__file__)); Z=np.asarray(L.ZGRID)
def _seg(zk,zk1,z1,z2,Ak,Ak1,alpha,z0):
    d=z2-z1
    if d<=0: return 0.0
    x=alpha*d; f1=math.exp(-alpha*(z1-z0))
    if x<1e-8: w0=d*(1-x/2+x*x/6); w1=d*d*(0.5-x/3+x*x/8)
    else: w0=-math.expm1(-x)/alpha; w1=(1-(1+x)*math.exp(-x))/(alpha*alpha)
    sk=(Ak1-Ak)/(zk1-zk)
    return alpha*f1*(Ak*w0+sk*((z1-zk)*w0+w1))
def exact_pwl(alpha,z0,A):
    if alpha<=0 or z0>=Z[-1]: return 0.0
    zstart=max(z0,Z[0]); k0=int(np.searchsorted(Z,zstart,side='right'))-1; k0=max(0,min(k0,len(Z)-2))
    out=_seg(Z[k0],Z[k0+1],zstart,Z[k0+1],A[k0],A[k0+1],alpha,z0)
    for k in range(k0+1,len(Z)-1): out+=_seg(Z[k],Z[k+1],Z[k],Z[k+1],A[k],A[k+1],alpha,z0)
    return out
scan=json.load(open(os.path.join(RUN,"outputs","data","signal_scan.json")))
th=np.array(scan['theta2']); sig=scan['signal_per_fb_per_channel']
def crossings_scan(m,level):
    s=np.array([r['B']+r['K']+r['H'] for r in sig[str(m)]])
    lx=np.log10(th); ly=np.log10(np.maximum(s,1e-300)); lv=math.log10(level); out=[]
    for i in range(len(th)-1):
        f0=ly[i]-lv; f1=ly[i+1]-lv
        if f0*f1<0: t=f0/(f0-f1); out.append(10**(lx[i]+t*(lx[i+1]-lx[i])))
    return out
for ms_ in (0.2802,1.0,3.5567):
    d=np.load(os.path.join(HERE,"tab_%g.npz"%ms_)); z0=d['z']; pz=d['pz']; A=d['A']
    print("m=%.4f" % ms_)
    for L in (25.0,300.0):
        cr=crossings_scan(ms_,2.3/L)
        for tag,t2 in (("long",min(cr)),("short",max(cr))):
            if tag=="short" and len(cr)<2: continue
            alpha=ms_*sw.gamma_total(ms_)*t2/(sw.HBAR_C*pz)
            pc=R.per_scalar_p(ms_,None,t2,pz,z0,A)
            pe=np.array([exact_pwl(alpha[i],z0[i],A[i]) for i in range(len(z0))])
            print("  L=%3.0f %-5s theta2=%.3e  sum_pc=%.6e sum_pe=%.6e  ratio=%.8f" % (L,tag,t2,pc.sum(),pe.sum(),pc.sum()/pe.sum()))
