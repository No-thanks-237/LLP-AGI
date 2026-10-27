# -*- coding: utf-8 -*-
"""F11 addendum: nsamp needed to reach 5% rel err at the short-lived edge, and cap check."""
import sys, os, math, json
import numpy as np
RUN=r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
scan=json.load(open(os.path.join(RUN,"outputs","data","signal_scan.json")))
th=np.array(scan['theta2']); sig=scan['signal_per_fb_per_channel']
masses=np.array(scan['masses']); nsamp={float(k):float(v) for k,v in scan['nsamp_by_mass'].items()}
def s_rel(row):
    s=np.array([r['B']+r['K']+r['H'] for r in row])
    var=np.array([(r['B']*r.get('relerr_B',0))**2+(r['K']*r.get('relerr_K',0))**2+(r['H']*r.get('relerr_H',0))**2 for r in row])
    return s, np.where(s>0,np.sqrt(var)/np.maximum(s,1e-300),0.0)
def crossings(x,y,level):
    lx=np.log10(x); ly=np.log10(np.maximum(y,1e-300)); lv=math.log10(level); out=[]
    for i in range(len(x)-1):
        f0=ly[i]-lv; f1=ly[i+1]-lv
        if f0*f1<0: t=f0/(f0-f1); out.append(10**(lx[i]+t*(lx[i+1]-lx[i])))
    return out
def rel_at(x,rel): return float(np.interp(math.log10(x),np.log10(th),rel))
for L in (25.0,300.0):
    print("L=%.0f fb^-1" % L)
    for tag, which in (("long-lived",'lo'),("short-lived",'hi')):
        needs=[]; bad=[]
        for m in masses:
            s,rel=s_rel(sig[str(m)]); cr=crossings(th,s,2.3/L)
            if not cr: continue
            if which=='hi' and len(cr)<2: continue
            x=min(cr) if which=='lo' else max(cr)
            r=rel_at(x,rel); n=nsamp.get(m,float('nan'))
            need=n*(r/0.05)**2
            needs.append(need)
            if r>0.05: bad.append((m,r,n,need))
        needs=np.array(needs)
        print("  %-13s >5%%: %2d/%2d  median rel %.2f%%  needed nsamp: median %.0f  max %.0f  >25000 cap: %d" % (
            tag, len(bad), len(needs), 100*np.median([b[1] for b in bad]) if bad else 0,
            np.median(needs), needs.max(), int((needs>25000).sum())))
        if which=='hi':
            print("     worst 5 short-lived masses (m, rel, nsamp_now, nsamp_needed):")
            for m,r,n,need in sorted(bad,key=lambda t:-t[1])[:5]:
                print("       m=%9.5f rel=%5.1f%%  now=%.0f need=%.0f" % (m,100*r,n,need))
