
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
xr, yr, _ = br["downstream_25fb"][0]
def bnd(fn, lumi=25.0):
    out = {}
    with open(fn) as f:
        for row in csv.DictReader(f):
            if float(row['lumi_fb']) != lumi: continue
            out[round(float(row['mass_GeV']),6)] = float(row['theta2_long_lived_edge'])
    return out
prod = bnd(os.path.join(RUN, "outputs", "data", "exclusion_boundary.csv"))
val  = bnd(os.path.join(RUN, "outputs", "data", "exclusion_boundary_refbr.csv"))
print("high-mass comparison: ours/ref at the 25 fb-1 long-lived edge")
print("%9s %12s %8s %12s %8s" % ("m_S","prod","p/ref","refbr","v/ref"))
for m in sorted(prod):
    if m < 3.0 or m > 61.0 or m < xr[0] or m > xr[-1]: continue
    r0 = float(np.interp(m, xr, yr))
    print("%9.4f %12.3e %8.3f %12.3e %8.3f" % (m, prod[m], prod[m]/r0, val.get(m, float('nan')), val.get(m, float('nan'))/r0))
