
import json, os, csv, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__))
def bnd(fn):
    out = {}
    with open(fn) as f:
        for row in csv.DictReader(f):
            if float(row['lumi_fb']) != 25.0: continue
            out[round(float(row['mass_GeV']),6)] = float(row['theta2_long_lived_edge'])
    return out
old = bnd(os.path.join(D, 'bnd_v34.csv'))
new = json.load(open(os.path.join(D, 'scan_v35_base.json')))
th = np.array(new['theta2'])
def cross(s, level):
    f = np.log10(np.maximum(s,1e-300))-math.log10(level); out=[]
    for i in range(len(th)-1):
        if f[i]*f[i+1]<0:
            out.append(10**(np.log10(th[i])+f[i]/(f[i]-f[i+1])*(np.log10(th[i+1])-np.log10(th[i]))))
    return out
print("%9s %12s %12s %8s" % ("m","v3.4","v3.5(base)","ratio"))
rs=[]
for m, row in new['signal_per_fb_per_channel'].items():
    s = np.array([r['B']+r['K']+r['H'] for r in row])
    cr = cross(s, 2.3/25.0)
    if not cr: continue
    mm = round(float(m),6); lo = min(cr); o = old.get(mm)
    if o: rs.append(lo/o)
    print("%9.5f %12.3e %12.3e %8.3f" % (mm, o if o else float('nan'), lo, lo/o if o else float('nan')))
rs=np.array(rs)
print("\nmedian ratio %.4f  min %.3f max %.3f" % (np.median(rs), rs.min(), rs.max()))
