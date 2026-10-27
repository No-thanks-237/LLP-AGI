# -*- coding: utf-8 -*-
"""F10: verify the lifetime-integral claim and quantify the realistic bias."""
import sys, os, math, json, random
import numpy as np
RUN = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, os.path.join(RUN,"work","analysis"))
sys.path.insert(0, os.path.join(RUN,"work","widths"))
sys.path.insert(0, os.path.join(RUN,"work","response"))
import run_scan as R
import bc5_lib as L
import scalar_widths as sw
HERE = os.path.dirname(os.path.abspath(__file__))
Z = np.asarray(L.ZGRID, dtype=float)
H = Z[1]-Z[0]

# -------- exact piecewise-linear-A lifetime integral (the recommended fix) --------
def lifetime_exact_pwl(alpha, z0, A, zgrid=Z):
    """int_{max(z0,zgrid[0])}^{zgrid[-1]} dz alpha exp(-(z-z0)alpha) * A_lin(z), A linear per cell."""
    z = zgrid; n = len(z)
    if alpha <= 0.0 or z0 >= z[-1] or z0 < -1e9:
        return 0.0
    zstart = max(z0, z[0])
    k0 = int(np.searchsorted(z, zstart, side='right')) - 1
    k0 = max(0, min(k0, n-2))
    out = 0.0
    # partial first cell [zstart, z[k0+1]]
    out += _seg(z[k0], z[k0+1], zstart, z[k0+1], A[k0], A[k0+1], alpha, z0)
    for k in range(k0+1, n-1):
        out += _seg(z[k], z[k+1], z[k], z[k+1], A[k], A[k+1], alpha, z0)
    return out

def _seg(zk, zk1, z1, z2, Ak, Ak1, alpha, z0):
    d = z2 - z1
    if d <= 0: return 0.0
    x = alpha * d
    f1 = math.exp(-alpha * (z1 - z0))
    if x < 1e-8:
        w0 = d * (1.0 - x/2.0 + x*x/6.0)
        w1 = d*d * (0.5 - x/3.0 + x*x/8.0)
    else:
        w0 = -math.expm1(-x) / alpha
        w1 = (1.0 - (1.0 + x) * math.exp(-x)) / (alpha*alpha)
    sk = (Ak1 - Ak) / (zk1 - zk)
    return alpha * f1 * (Ak * w0 + sk * ((z1 - zk) * w0 + w1))

def lifetime_A1_analytic(alpha, z0, lo=1.0, hi=2.5):
    zs = max(z0, lo)
    if zs >= hi: return 0.0
    return math.exp(-(zs - z0)*alpha) - math.exp(-(hi - z0)*alpha)

# -------- 1. reproduce the review's A=1 benchmark table --------------
print("="*86)
print("(1) A(z)=1 benchmark: review's table, re-derived with the code's own per_scalar_p")
print("="*86)
ms = 1.0; pz0 = 1.0
def code_p_for_alpha(alpha, z0):
    theta2 = alpha * sw.HBAR_C / (ms * sw.gamma_total(ms) * pz0)
    p = R.per_scalar_p(ms, None, theta2, np.array([pz0]), np.array([z0]), np.ones((1,len(Z))))
    return float(p[0])
print("  %-8s %-10s %-14s %-14s %-8s" % ("z0[m]","Lz=1/a[m]","code","exact A=1","code/exact"))
for z0, Lz in ((0.0,0.1),(0.0,0.3),(1.01,0.1),(1.01,1.0),(1.249,1.0)):
    a = 1.0/Lz
    c = code_p_for_alpha(a, z0); e = lifetime_A1_analytic(a, z0)
    print("  %-8.3f %-10.3f %-14.8f %-14.8f %-8.5f" % (z0, Lz, c, e, c/e))
print("  (review table quotes 1.12694, 1.01443, 0.55473, 0.93280, 1.08737)")

# -------- 2. realistic (z0, pz) distributions from run_mass --------------
print()
print("="*86)
print("(2) realistic bias from run_scan.run_mass tables (nsamp=300, ndecay=8)")
print("="*86)
scan = json.load(open(os.path.join(RUN,"outputs","data","signal_scan.json")))
th = np.array(scan['theta2'], dtype=float)
sig = scan['signal_per_fb_per_channel']
masses = np.array(scan['masses'], dtype=float)
def _row(m):
    for k in (str(m), "%.6f" % m, "%.4f" % m):
        if k in sig:
            return sig[k]
    raise KeyError(m)   # not a scan mass point
def crossings_scan(m, level):
    row = _row(m)
    s = np.array([(r['B']+r['K']+r['H']) for r in row], dtype=float)
    lx=np.log10(th); ly=np.log10(np.maximum(s,1e-300)); lv=math.log10(level)
    out=[]
    for i in range(len(th)-1):
        f0=ly[i]-lv; f1=ly[i+1]-lv
        if f0*f1<0:
            t=f0/(f0-f1); out.append(10**(lx[i]+t*(lx[i+1]-lx[i])))
    return out
for ms_ in (0.2802, 1.0, 3.5567):
    d = np.load(os.path.join(HERE, "tab_%g.npz" % ms_))
    z0 = d['z']; pz = d['pz']; A = d['A']
    ctau1 = sw.HBAR_C / sw.gamma_total(ms_)     # ctau at theta^2=1
    cr25 = crossings_scan(ms_, 2.3/25.0); cr300 = crossings_scan(ms_, 2.3/300.0)
    print("\n m_S=%.2f GeV  rows=%d  z0 max=%.4f m  (all z0 << 1 m)" % (ms_, len(z0), z0.max()))
    print("   crossings 2.3/L:  25 fb^-1: %s" % ["%.3e"%x for x in cr25[:1]+cr25[-1:]])
    print("                     300 fb^-1: %s" % ["%.3e"%x for x in cr300[:1]+cr300[-1:]])
    cases = [("long-lived edge, 25 fb", (min(cr25) if cr25 else float('nan'))),
             ("long-lived edge, 300 fb", (min(cr300) if cr300 else float('nan'))),
             ("short-lived edge, 25 fb", (max(cr25) if cr25 else float('nan'))),
             ("short-lived edge, 300 fb", (max(cr300) if cr300 else float('nan')))]
    for tag, t2 in cases:
        if not np.isfinite(t2):
            print("   %-26s  no crossing" % tag); continue
        alpha = ms_*sw.gamma_total(ms_)*t2/(sw.HBAR_C*pz)
        pc = R.per_scalar_p(ms_, None, t2, pz, z0, A)
        pe = np.array([lifetime_exact_pwl(alpha[i], z0[i], A[i]) for i in range(len(z0))])
        pa1 = np.array([lifetime_A1_analytic(alpha[i], z0[i]) for i in range(len(z0))])
        ratio = float(np.mean(pc)/np.mean(pe))
        rat1 = float(np.mean(pc)/np.mean(pa1))
        print("   %-26s theta2=%.3e  <alpha>=%.3e m^-1  <alpha>*dz=%.3f  code/exact(PWL)=%.5f  code/A=1 exact=%.5f"
              % (tag, t2, alpha.mean(), alpha.mean()*H, ratio, rat1))
