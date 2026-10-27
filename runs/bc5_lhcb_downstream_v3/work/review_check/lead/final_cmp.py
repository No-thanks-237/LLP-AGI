
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
br = {}
for cid, info in rep['curves'].items():
    br[cid] = []
    for seg in info['segments']:
        a = rows[seg['row_start']-2: seg['row_end']-1]
        br[cid].append((np.array([r[1] for r in a]), np.array([r[2] for r in a]), seg.get('branch','?')))
def boundary(fn, lumi=25.0):
    out = {}
    with open(fn) as f:
        for row in csv.DictReader(f):
            if float(row['lumi_fb']) != lumi: continue
            out[round(float(row['mass_GeV']),6)] = (float(row['theta2_long_lived_edge']),
                                                    float(row['theta2_short_lived_edge']))
    return out
new = boundary(os.path.join(RUN, "outputs", "data", "exclusion_boundary.csv"))
xr, yr, _ = br["downstream_25fb"][0]
rat = []
print("%9s %11s %11s %8s   %11s %11s" % ("m","ours(long)","ref(lower)","ratio","ours(short)","ref(upper)"))
for m in sorted(new):
    lo, hi = new[m]
    if not np.isfinite(lo) or lo <= 0: continue
    if m < xr[0] or m > xr[-1]: continue
    r0 = float(np.interp(m, xr, yr))
    rat.append(lo/r0)
    print("%9.4f %11.3e %11.3e %8.3f" % (m, lo, r0, lo/r0))
rat = np.array(rat)
print("\n25 fb-1 long-lived edge vs reference lower branch: n=%d  median %.3f  |ln| median %.3f  range %.2f-%.2f" % (
    len(rat), np.median(rat), np.median(np.abs(np.log(rat))), rat.min(), rat.max()))
# and 300 fb-1
b300 = boundary(os.path.join(RUN, "outputs", "data", "exclusion_boundary.csv"), 300.0)
x3, y3, _ = br["downstream_300fb"][0]
r3 = []
for m in sorted(b300):
    lo, hi = b300[m]
    if not np.isfinite(lo) or lo <= 0 or m < x3[0] or m > x3[-1]: continue
    r3.append(lo/float(np.interp(m, x3, y3)))
r3 = np.array(r3)
print("300 fb-1 long-lived edge vs reference: n=%d  median %.3f  |ln| median %.3f  range %.2f-%.2f" % (
    len(r3), np.median(r3), np.median(np.abs(np.log(r3))), r3.min(), r3.max()))
