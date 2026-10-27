# -*- coding: utf-8 -*-
"""F10 part 2: realistic signal-weighted bias, exact-formula verification and convergence."""
import sys, os, math, json
import numpy as np
RUN = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, os.path.join(RUN,"work","analysis")); sys.path.insert(0, os.path.join(RUN,"work","widths"))
sys.path.insert(0, os.path.join(RUN,"work","response"))
import run_scan as R, bc5_lib as L, scalar_widths as sw
HERE = os.path.dirname(os.path.abspath(__file__))
Z = np.asarray(L.ZGRID, dtype=float); H = Z[1]-Z[0]

def _seg(zk, zk1, z1, z2, Ak, Ak1, alpha, z0):
    d = z2-z1
    if d <= 0: return 0.0
    x = alpha*d
    f1 = math.exp(-alpha*(z1-z0))
    if x < 1e-8:
        w0 = d*(1 - x/2 + x*x/6); w1 = d*d*(0.5 - x/3 + x*x/8)
    else:
        w0 = -math.expm1(-x)/alpha; w1 = (1-(1+x)*math.exp(-x))/(alpha*alpha)
    sk = (Ak1-Ak)/(zk1-zk)
    return alpha*f1*(Ak*w0 + sk*((z1-zk)*w0 + w1))

def lifetime_exact_pwl(alpha, z0, A, zgrid=Z):
    z=zgrid; n=len(z)
    if alpha <= 0 or z0 >= z[-1]: return 0.0
    zstart=max(z0,z[0]); k0=int(np.searchsorted(z,zstart,side='right'))-1
    k0=max(0,min(k0,n-2)); out=_seg(z[k0],z[k0+1],zstart,z[k0+1],A[k0],A[k0+1],alpha,z0)
    for k in range(k0+1,n-1):
        out += _seg(z[k],z[k+1],z[k],z[k+1],A[k],A[k+1],alpha,z0)
    return out

def lifetime_A1(alpha,z0,lo=1.0,hi=2.5):
    zs=max(z0,lo)
    if zs>=hi: return 0.0
    return math.exp(-(zs-z0)*alpha)-math.exp(-(hi-z0)*alpha)

# -------- (a) exact PWL formula == analytic for A=1 (machine precision) ----------
err=0.0
for z0 in (-0.1, 0.0, 0.5, 0.999, 1.0, 1.01, 1.0625, 1.249, 1.9, 2.49, 2.5, 3.0):
    for alpha in (0.0, 1e-6, 1e-3, 0.05, 0.5, 1.0, 2.0, 8.0, 50.0, 500.0):
        e=lifetime_exact_pwl(alpha,z0,np.ones(len(Z))); a=lifetime_A1(alpha,z0)
        err=max(err,abs(e-a))
print("(a) exact PWL formula vs analytic A=1: max |diff| over (alpha,z0) grid = %.3e" % err)

# -------- (b) trapezoid convergence O(h^2) vs the exact integral -----------------
print()
print("(b) 13-point trapezoid (the code) vs exact, A=1, z0=0:")
for alpha in (0.1, 0.5, 1.0, 2.0, 5.0):
    exact=lifetime_A1(alpha,0.0)
    def trap(n):
        z=np.linspace(1.0,2.5,n); h=z[1]-z[0]
        w=alpha*np.exp(-(z-0.0)*alpha)
        return float(np.trapezoid(w,z) if hasattr(np,'trapezoid') else np.trapz(w,z))
    r13=trap(13)/exact
    # grid refinement
    errs=[abs(trap(13*2**k)/exact-1) for k in range(5)]
    print("   alpha=%5.1f  code(13 nodes)/exact=%.5f   trapezooid error vs 1/h: %s" % (
        alpha, r13, " ".join("%.2e" % e for e in errs)))
print("   (error halves by ~4 per doubling of the node count -> O(h^2); the exact formula is exact)")

# -------- (c) realistic signal-weighted ratio at the actual crossings -----------
print()
print("(c) realistic bias with the run_mass (z0, pz, A) tables")
scan=json.load(open(os.path.join(RUN,"outputs","data","signal_scan.json")))
th=np.array(scan['theta2']); sig=scan['signal_per_fb_per_channel']
def crossings_scan(m, level):
    row=sig[str(m)]; s=np.array([r['B']+r['K']+r['H'] for r in row])
    lx=np.log10(th); ly=np.log10(np.maximum(s,1e-300)); lv=math.log10(level); out=[]
    for i in range(len(th)-1):
        f0=ly[i]-lv; f1=ly[i+1]-lv
        if f0*f1<0: t=f0/(f0-f1); out.append(10**(lx[i]+t*(lx[i+1]-lx[i])))
    return out
for ms_ in (0.2802, 1.0, 3.5567):
    d=np.load(os.path.join(HERE,"tab_%g.npz"%ms_)); z0=d['z']; pz=d['pz']; A=d['A']; w=d.get('w') if 'w' in d else None
    ctau1=sw.HBAR_C/sw.gamma_total(ms_)
    print("\n m_S=%.4f GeV  rows=%d  z0<1m: %d/%d  max(z0<1m)=%.4f m  frac z0<0: %.3f" % (
        ms_, len(z0), int((z0<1).sum()), len(z0), z0[z0<1].max() if (z0<1).any() else 0, float((z0<0).mean())))
    for L in (25.0, 300.0):
        cr=crossings_scan(ms_, 2.3/L)
        if not cr: continue
        for tag,t2 in (("long-lived", min(cr)), ("short-lived", max(cr))):
            if len(cr)<2 and tag=="short-lived": continue
            alpha=ms_*sw.gamma_total(ms_)*t2/(sw.HBAR_C*pz)
            pc=R.per_scalar_p(ms_, None, t2, pz, z0, A)
            pe=np.array([lifetime_exact_pwl(alpha[i],z0[i],A[i]) for i in range(len(z0))])
            ratio=float(pc.sum()/pe.sum())
            print("   L=%3.0f %-12s th2=%.3e  <alpha>=%.2e  <p_code>=%.4f <p_exact>=%.4f  sum ratio=%.5f" % (
                L, tag, t2, alpha.mean(), pc.mean(), pe.mean(), ratio))
