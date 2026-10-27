
import json, os, math
import numpy as np
RUN = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3"
scan = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "scan_v35_base.json")))
print("model_tag", scan.get("model_tag"), " masses", len(scan["masses"]))
th = np.array(scan["theta2"])
def cross(s, level):
    f = np.log10(np.maximum(s,1e-300))-math.log10(level); out=[]
    for i in range(len(th)-1):
        if f[i]*f[i+1]<0:
            out.append(10**(np.log10(th[i])+f[i]/(f[i]-f[i+1])*(np.log10(th[i+1])-np.log10(th[i]))))
    return out
for lum in (25.0, 300.0):
    rels=[]
    for m, row in scan["signal_per_fb_per_channel"].items():
        s = np.array([r['B']+r['K']+r['H'] for r in row])
        var = np.array([(r['B']*r.get('relerr_B',0))**2+(r['K']*r.get('relerr_K',0))**2+(r['H']*r.get('relerr_H',0))**2 for r in row])
        rel = np.where(s>0, np.sqrt(var)/np.maximum(s,1e-300), 0.0)
        cr = cross(s, 2.3/lum)
        if not cr: continue
        lo = min(cr)
        i = int(np.argmin(np.abs(np.log10(th)-np.log10(lo))))
        rels.append((rel[i], float(m)))
    rels.sort(reverse=True)
    print("L=%.0f: median %.2f%%  max %.2f%%  above 5%%: %d/%d" % (
        lum, 100*np.median([r for r,_ in rels]), 100*rels[0][0],
        sum(1 for r,_ in rels if r>0.05), len(rels)))
    print("   worst:", ", ".join("%.5g: %.2f%%" % (m, 100*r) for r,m in rels[:8]))