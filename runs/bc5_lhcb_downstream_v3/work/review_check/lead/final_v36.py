
import json, os, csv, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__))
RUN = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3"
REF = os.path.join(RUN, "work", "widths", "reference_figures", "bc5_reference_curves.csv")
REPORT = os.path.join(RUN, "work", "widths", "reference_figures", "bc5_digitisation_report.json")
rep = json.load(open(REPORT))
rows = []
with open(REF) as f:
    for row in csv.DictReader(f):
        rows.append((row['curve_id'], float(row['mass_GeV']), float(row['theta2'])))
br = {}
for cid, info in rep['curves'].items():
    br[cid] = []
    for seg in info['segments']:
        a = rows[seg['row_start']-2: seg['row_end']-1]
        br[cid].append((np.array([r[1] for r in a]), np.array([r[2] for r in a]), seg.get('branch','?')))
def bnd(fn, lumi=25.0):
    out = {}
    with open(fn) as f:
        for row in csv.DictReader(f):
            if float(row['lumi_fb']) != lumi: continue
            out[round(float(row['mass_GeV']),6)] = float(row['theta2_long_lived_edge'])
    return out
new = bnd(os.path.join(RUN, "outputs", "data", "exclusion_boundary.csv"))
old = bnd(os.path.join(D, "bnd_v35.csv"))
print("%9s %12s %12s %8s   %12s %8s" % ("m","v3.5","v3.6","36/35","ref(low)","36/ref"))
xr, yr, _ = br["downstream_25fb"][0]
rs_old, rs_new = [], []
for m in sorted(new):
    if not (np.isfinite(new[m]) and new[m] > 0): continue
    if m < xr[0] or m > xr[-1]: continue
    r0 = float(np.interp(m, xr, yr)); o = old.get(m)
    if o: rs_old.append(o/r0)
    rs_new.append(new[m]/r0)
    print("%9.4f %12.3e %12.3e %8.3f   %12.3e %8.3f" % (m, o if o else float('nan'), new[m], new[m]/o if o else float('nan'), r0, new[m]/r0))
print()
print("v3.5 vs ref: median %.3f  |ln| median %.3f" % (np.median(rs_old), np.median(np.abs(np.log(rs_old)))))
print("v3.6 vs ref: median %.3f  |ln| median %.3f" % (np.median(rs_new), np.median(np.abs(np.log(rs_new)))))
# validation scan
newr = bnd(os.path.join(RUN, "outputs", "data", "exclusion_boundary_refbr.csv"))
x3, y3, _ = br["downstream_25fb"][0]
rr = [newr[m]/float(np.interp(m, x3, y3)) for m in sorted(newr) if newr.get(m) and np.isfinite(newr[m]) and newr[m]>0 and x3[0] <= m <= x3[-1]]
print("v3.5 REFBR vs ref: median %.3f  |ln| median %.3f" % (np.median(rr), np.median(np.abs(np.log(rr)))))