
import json, os, csv, math
import numpy as np
RUN = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3"
REF = os.path.join(RUN, "work", "widths", "reference_figures", "bc5_reference_curves.csv")
REPORT = os.path.join(RUN, "work", "widths", "reference_figures", "bc5_digitisation_report.json")
rep = json.load(open(REPORT))
rows = []
with open(REF) as f:
    for row in csv.DictReader(f):
        rows.append((row['curve_id'], float(row['mass_GeV']), float(row['theta2'])))
branches = {}
for cid, info in rep['curves'].items():
    branches[cid] = []
    for seg in info['segments']:
        a = rows[seg['row_start']-2: seg['row_end']-1]
        branches[cid].append((np.array([r[1] for r in a]), np.array([r[2] for r in a]),
                              seg.get('branch','?')))
for cid, brs in branches.items():
    print(cid, [(b[2], len(b[0]), round(float(b[0][0]),3), round(float(b[0][-1]),3)) for b in brs])
new = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "scan_peek.json")))
th = np.array(new["theta2"])
def cross(s, level):
    f = np.log10(np.maximum(s,1e-300))-math.log10(level); out=[]
    for i in range(len(th)-1):
        if f[i]*f[i+1]<0:
            out.append(10**(np.log10(th[i])+f[i]/(f[i]-f[i+1])*(np.log10(th[i+1])-np.log10(th[i]))))
    return out
oldb = {}
with open(os.path.join(RUN, "outputs", "data", "exclusion_boundary.csv")) as f:
    for row in csv.DictReader(f):
        if float(row["lumi_fb"]) == 25.0:
            oldb[round(float(row["mass_GeV"]),6)] = float(row["theta2_long_lived_edge"])
lower = branches["downstream_25fb"][0]
xr, yr = lower[0], lower[1]
print("\n%9s %11s %11s %11s %8s %8s" % ("m","v3.3","v3.4","ref(low)","old/ref","new/ref"))
ro_, rn_ = [], []
for m, row in new["signal_per_fb_per_channel"].items():
    s = np.array([r['B']+r['K']+r['H'] for r in row])
    cr = cross(s, 2.3/25.0)
    if not cr: continue
    mm = float(m); lo = min(cr)
    if mm < xr[0] or mm > xr[-1]: continue
    r0 = float(np.interp(mm, xr, yr))
    o = oldb.get(round(mm,6))
    if not o: continue
    ro, rnn = o/r0, lo/r0
    ro_.append(ro); rn_.append(rnn)
    print("%9.4f %11.3e %11.3e %11.3e %8.3f %8.3f" % (mm, o, lo, r0, ro, rnn))
print("\nmedian old/ref %.3f   new/ref %.3f" % (np.median(ro_), np.median(rn_)))
print("median |ln|  old %.3f   new %.3f" % (np.median(np.abs(np.log(ro_))), np.median(np.abs(np.log(rn_)))))
