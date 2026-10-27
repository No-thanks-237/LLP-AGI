
import sys, os, math
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\analysis")
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\widths")
import numpy as np, scalar_widths as sw, reference_br_model as rb
sw.load_hadronic_table(os.path.join(os.path.dirname(sw.__file__), 'hadronic_table.csv'))
# review's numbers (high_mass_reference_diagnosis_20261025.md)
K_ratio = {5.082255: 0.830, 13.023095: 0.810, 33.371213: 0.796, 53.419632: 0.761, 62.49: 0.796}
resp_ratio = {5.082255: 1.041, 13.023095: 1.013, 33.371213: 0.990, 53.419632: 0.955, 62.49: 0.951}
import csv, json
RUN = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3"
REF = os.path.join(RUN, "work", "widths", "reference_figures", "bc5_reference_curves.csv")
REPORT = os.path.join(RUN, "work", "widths", "reference_figures", "bc5_digitisation_report.json")
rep = json.load(open(REPORT)); rows = []
with open(REF) as f:
    for row in csv.DictReader(f):
        rows.append((row['curve_id'], float(row['mass_GeV']), float(row['theta2'])))
br = {}
for cid, info in rep['curves'].items():
    br[cid] = []
    for seg in info['segments']:
        a = rows[seg['row_start']-2: seg['row_end']-1]
        br[cid].append((np.array([r[1] for r in a]), np.array([r[2] for r in a])))
xr, yr = br["downstream_25fb"][0]
def bnd(fn):
    out = {}
    with open(fn) as f:
        for row in csv.DictReader(f):
            if float(row['lumi_fb']) != 25.0: continue
            out[round(float(row['mass_GeV']),6)] = float(row['theta2_long_lived_edge'])
    return out
prod = bnd(os.path.join(RUN, "outputs", "data", "exclusion_boundary.csv"))
val  = bnd(os.path.join(RUN, "outputs", "data", "exclusion_boundary_refbr.csv"))
print("%8s %8s %8s | %7s | %8s %8s %8s %8s" % ("m_S","G_ref/G","K_ours/K_ref","resp","p/ref","v/ref","pred(prod)","pred(val)"))
for m in sorted(K_ratio):
    g_ours = sum(sw.widths(m).values()); g_ref = rb.gamma_total(m)
    w = g_ref/g_ours                 # boundary ratio from the width alone
    k = K_ratio[m]; r = resp_ratio[m]
    p = prod.get(round(m,6)); v = val.get(round(m,6))
    r0 = float(np.interp(m, xr, yr))
    print("%8.3f %8.3f %8.3f | %7.3f | %8.3f %8.3f %8.3f %8.3f" % (
        m, w, k, r, p/r0 if p else float('nan'), v/r0 if v else float('nan'), w*k*r, k*r))
