# -*- coding: utf-8 -*-
"""Final review, checks 3-4: lifetime/decay-length consistency and the low-mass end."""
import json, csv, math, sys, gzip
import numpy as np
R = r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream"
sys.path.insert(0, R + "/work/review_check")
import audit_2hdm as AU
S90 = -math.log(0.1)
d = json.load(open(R + "/outputs/data/signal_scan.json"))
mg, tg = d['mass_grid'], d['tanb_grid']
idx = lambda v, g: min(range(len(g)), key=lambda i: abs(g[i]-v))
P = {}
for p in d['points']:
    P[(idx(p['mH'], mg), idx(p['tanb'], tg))] = p
key = lambda L: str(int(L))

# ---- H kinematics from the B sample (B -> K H two-body, isotropic, correct Lorentz boost)
mB = 5.27934; mK = 0.493677
data = []
with gzip.open(R + "/work/magnus/samples/bbar_heavy.txt.gz", 'rt') as f:
    for i, ln in enumerate(f):
        if i >= 400000: break
        a = ln.split()
        if len(a) >= 10 and a[1] in ('521','-521','511','-511'):
            data.append((float(a[2]), float(a[3]), float(a[4]), float(a[5])))
D = np.array(data); px,py,pz,E = D.T
pm = np.linalg.norm(D[:,:3],axis=1); u = D[:,:3]/pm[:,None]
b = pm/E; g = E/mB
rng = np.random.default_rng(11)
def bgz_percentiles(mH, qs=(10,50,90)):
    ps = math.sqrt(max((mB**2-(mK+mH)**2)*(mB**2-(mK-mH)**2),0.0))/(2*mB)
    Es = math.sqrt(ps*ps+mH*mH)
    ct = rng.uniform(-1,1,len(D)); ph = rng.uniform(0,2*np.pi,len(D)); st = np.sqrt(1-ct*ct)
    h = np.stack([ps*st*np.cos(ph), ps*st*np.sin(ph), ps*ct], axis=1)
    bp = np.einsum('ij,ij->i', h, u)
    k = (g-1.0)*bp + (g*b)*Es
    pl = h + k[:,None]*u
    eta = 0.5*np.log((np.linalg.norm(pl,axis=1)+pl[:,2])/(np.linalg.norm(pl,axis=1)-pl[:,2]+1e-300))
    ok = (eta>2)&(eta<5)&(pl[:,2]>0)
    v = pl[ok,2]/mH
    return {q: float(np.percentile(v,q)) for q in qs}, float(ok.mean())

def ctau_interp(im, tcross):
    """linear interpolation of log ctau vs log tanb on the scan grid"""
    lts = [math.log(tg[i]) for i in range(len(tg))]
    lcs = [math.log(P[(im,i)]['ctau_m']) for i in range(len(tg))]
    lt = math.log(tcross)
    for i in range(len(tg)-1):
        if lts[i] <= lt <= lts[i+1]:
            f = (lt-lts[i])/(lts[i+1]-lts[i])
            return math.exp(lcs[i]+f*(lcs[i+1]-lcs[i]))
    return float('nan')

print("="*104)
print("CHECK 3: decay length at the long-lived edges vs the detector decay window z = 1.0 - 2.5 m")
print("="*104)
cases = [(1.0, 25.0, 357.029), (1.0, 300.0, 487.491), (3.50112, 25.0, 273.366), (3.50112, 300.0, 410.768),
         (0.289535, 25.0, 179.391), (0.021, 25.0, 21.8753)]
print(f"{'mH':>8} {'L':>5} {'edge tanb':>10} {'ctau [m]':>10} {'b*g_z p10/50/90':>22} {'ctau*b*g_z [m] (med)':>22} {'x window(1.5m)':>14}")
W = {"B": None}
for (m, L, tc) in cases:
    im = idx(m, mg)
    ctm = ctau_interp(im, tc)
    pct,frac = bgz_percentiles(m)
    med = pct[50]
    print(f"{m:8.5f} {L:5.0f} {tc:10.4f} {ctm:10.4g}   {pct[10]:7.2f}/{med:7.2f}/{pct[90]:7.2f}"
          f"    {ctm*med:10.3g} ({ctm*pct[10]:.3g}-{ctm*pct[90]:.3g})   {ctm*med/1.5:8.2f}")
print()
print("Where the acceptance peaks (decay length ~ window): tanb such that ctau*b*g_z ~ 1.5 m")
for (m, L, _) in cases[:4]:
    im = idx(m, mg)
    pct,_ = bgz_percentiles(m); med = pct[50]
    for it in range(len(tg)-1):
        c0, c1 = P[(im,it)]['ctau_m'], P[(im,it+1)]['ctau_m']
        l0, l1 = tg[it], tg[it+1]
        if (c0*med-1.5)*(c1*med-1.5) <= 0:
            f = (1.5-c0*med)/(c1*med-c0*med)
            print(f"  mH={m:8.5f}: ctau*b*g_z(median) crosses 1.5 m at tanb = {math.exp(math.log(l0)+f*(math.log(l1)-math.log(l0))):.2f}")
            break
print()
print("="*104)
print("CHECK 4: low-mass end mH = 0.021 GeV -- where is the short-lived edge?")
print("="*104)
im = idx(0.021, mg)
print(f"{'tanb':>10} {'ctau [m]':>11} {'b*g_z(med)':>11} {'ctau*b*g_z [m]':>15} {'s_25':>12} {'s_300':>12}")
pct,_ = bgz_percentiles(0.021); med = pct[50]
for it in (0,1,2,4,8,12,16):
    p = P[(im,it)]
    print(f"{p['tanb']:10.4f} {p['ctau_m']:11.4g} {med:11.2f} {p['ctau_m']*med:15.4g} {p['s']['25']:12.5g} {p['s']['300']:12.5g}")
# where ctau*bgz = 1 m and 2.5 m (interpolating log ctau vs log tanb)
lts = [math.log(tg[i]) for i in range(len(tg))]; lcs = [math.log(P[(im,i)]['ctau_m']) for i in range(len(tg))]
for target in (1.0, 2.5):
    for i in range(len(tg)-1):
        c0, c1 = P[(im,i)]['ctau_m']*med, P[(im,i+1)]['ctau_m']*med
        if (c0-target)*(c1-target) <= 0:
            f = (target-c0)/(c1-c0)
            print(f"  ctau*b*g_z = {target} m at tanb = {math.exp(lts[i]+f*(lts[i+1]-lts[i])):.3f}")
            break
print()
print("mH=0.021 s(tanb=5) at both lumis vs s90 = 2.3026:",
      P[(im,0)]['s']['25'], P[(im,0)]['s']['300'], "-> already far above threshold;")
print("the only crossing in [5,1e5] is the long-lived one at", 21.8753, "(25 fb^-1) /", 29.8414, "(300 fb^-1)")
print()
print("validity check for the short-lived crossings:")
for (m, L, tc) in [(1.0, 25.0, 6.0601), (1.0, 300.0, 5.8437), (3.50112, 25.0, 24.8203), (3.50112, 300.0, 21.5034)]:
    print(f"  mH={m:8.5f} L={L:3.0f} short-lived crossing tanb={tc:8.4f}: Br(B->XsH) = {AU.br_B_to_XsH(m, tc)[0]:.4g}")
