
import json, numpy as np, csv, math
RUN = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
scan = json.load(open(RUN+"/outputs/data/signal_scan.json"))
masses = np.array(scan['masses']); thetas = np.array(scan['theta2'])
sig = scan['signal_per_fb_per_channel']
print("theta2 grid: %.3e .. %.3e  n=%d"%(thetas[0],thetas[-1],len(thetas)))
print("masses:", np.round(masses,4))
# boundary table
rows = list(csv.DictReader(open(RUN+"/outputs/data/exclusion_boundary.csv")))
def bnd(m,L):
    for r in rows:
        if abs(float(r['mass_GeV'])-m)<1e-9 and abs(float(r['lumi_fb'])-L)<1e-9:
            return r
    return None
def chan_at(m, t2):
    row = sig[str(m)]
    # nearest theta2
    i = int(np.argmin(np.abs(np.log(thetas/t2))))
    r = row[i]
    return r, thetas[i]
def interp_frac(m, t2):
    # log-log interpolate each channel
    row = sig[str(m)]
    out={}
    for c in ('B','K','H'):
        y = np.array([r.get(c,0.0) for r in row])
        lt = np.log(thetas); 
        y = np.maximum(y,1e-300)
        out[c]=float(np.exp(np.interp(np.log(t2), lt, np.log(y))))
    s = sum(out.values())
    return out, s
print("\n--- K fraction vs mass at long-lived edge (25 fb-1) ---")
for m in masses:
    b = bnd(m,25.0)
    if b is None: continue
    lo = float(b['theta2_long_lived_edge'])
    hi = float(b['theta2_short_lived_edge']) if b['theta2_short_lived_edge'] not in ('','nan') else float('nan')
    out, s = interp_frac(m, lo) if np.isfinite(lo) else ({'B':0,'K':0,'H':0},0)
    if 0.05 < m < 2.2:
        print("m=%8.4f lo=%.3e K/tot=%8.3e  B=%9.3e K=%9.3e H=%9.3e  s25=%.3e"%(
            m, lo, out['K']/s if s>0 else 0, out['B']*25, out['K']*25, out['H']*25, s*25))
print("\n--- K fraction vs mass at short-lived edge (25 fb-1) ---")
for m in masses:
    b = bnd(m,25.0)
    if b is None: continue
    hi = b['theta2_short_lived_edge']
    if hi in ('','nan'): continue
    hi=float(hi)
    out, s = interp_frac(m, hi)
    if 0.05 < m < 2.2:
        print("m=%8.4f hi=%.3e K/tot=%8.3e  B=%9.3e K=%9.3e H=%9.3e  s25=%.3e"%(
            m, hi, out['K']/s if s>0 else 0, out['B']*25, out['K']*25, out['H']*25, s*25))
print("\n--- K fraction through low-mass region at fixed theta2 = long-lived edge values ---")
for m in masses:
    if m>0.4: break
    b = bnd(m,25.0)
    lo = float(b['theta2_long_lived_edge'])
    out, s = interp_frac(m, lo)
    print("m=%8.4f lo=%.3e  frac B=%.3f K=%.4f H=%.3f  (s25=%.3e)"%(m,lo,out['B']/s,out['K']/s,out['H']/s,s*25))
