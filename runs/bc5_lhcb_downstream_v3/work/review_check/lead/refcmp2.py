
import json, os, csv, math
import numpy as np
RUN = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3"
curves = {}
with open(os.path.join(RUN, "work", "widths", "reference_figures", "bc5_reference_curves.csv")) as f:
    for row in csv.DictReader(f):
        curves.setdefault(row["curve_id"], []).append((float(row["mass_GeV"]), float(row["theta2"])))
for k in curves: curves[k].sort()
def ref_at(cid, m, which="lower"):
    pts = curves.get(cid, [])
    if not pts: return None
    xs = np.array([p[0] for p in pts]); ys = np.array([p[1] for p in pts])
    return float(np.interp(m, xs, ys))
new = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "scan_peek.json")))
th = np.array(new["theta2"])
def cross(s, level):
    f = np.log10(np.maximum(s,1e-300))-math.log10(level); out=[]
    for i in range(len(th)-1):
        if f[i]*f[i+1]<0:
            out.append(10**(np.log10(th[i])+f[i]/(f[i]-f[i+1])*(np.log10(th[i+1])-np.log10(th[i]))))
    return out
import csv as _csv
oldb = {}
with open(os.path.join(RUN, "outputs", "data", "exclusion_boundary.csv")) as f:
    for row in _csv.DictReader(f):
        if float(row["lumi_fb"]) == 25.0:
            oldb[round(float(row["mass_GeV"]),6)] = float(row["theta2_long_lived_edge"])
res = []
for m, row in new["signal_per_fb_per_channel"].items():
    s = np.array([r['B']+r['K']+r['H'] for r in row])
    cr = cross(s, 2.3/25.0)
    if not cr: continue
    mm = float(m); lo = min(cr)
    r0 = ref_at("downstream_25fb", mm)
    if r0 is None or r0 <= 0: continue
    res.append((mm, oldb.get(round(mm,6)), lo, r0))
print("%9s %11s %11s %11s %8s %8s" % ("m","v3.3","v3.4","reference","old/ref","new/ref"))
rat_old, rat_new = [], []
for mm, o, lo, r0 in sorted(res):
    ro = o/r0 if o else float('nan'); rn = lo/r0
    rat_old.append(ro); rat_new.append(rn)
    print("%9.4f %11.3e %11.3e %11.3e %8.3f %8.3f" % (mm, o if o else float('nan'), lo, r0, ro, rn))
print()
print("median |ln(old/ref)| = %.3f,  median |ln(new/ref)| = %.3f" % (
    np.median(np.abs(np.log(np.array(rat_old)))), np.median(np.abs(np.log(np.array(rat_new))))))
print("median old/ref %.3f  new/ref %.3f" % (np.median(rat_old), np.median(rat_new)))
