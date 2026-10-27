
import json, os, math
import numpy as np
RUN = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3"
new = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "scan_peek.json")))
old = json.load(open(os.path.join(RUN, "outputs", "data", "exclusion_boundary.csv_all_backup.json")))["x"] if False else None
import csv
oldb = {}
with open(os.path.join(RUN, "outputs", "data", "exclusion_boundary.csv")) as f:
    for row in csv.DictReader(f):
        if float(row["lumi_fb"]) == 25.0:
            oldb[round(float(row["mass_GeV"]),6)] = float(row["theta2_long_lived_edge"])
th = np.array(new["theta2"])
def cross(s, level):
    f = np.log10(np.maximum(s,1e-300))-math.log10(level); out=[]
    for i in range(len(th)-1):
        if f[i]*f[i+1]<0:
            out.append(10**(np.log10(th[i])+f[i]/(f[i]-f[i+1])*(np.log10(th[i+1])-np.log10(th[i]))))
    return out
print("25 fb-1 long-lived edge, v3.3 -> v3.4 (base pass, nsamp 1200)")
print("%9s %12s %12s %7s" % ("m", "v3.3", "v3.4", "ratio"))
rows = []
for m, row in new["signal_per_fb_per_channel"].items():
    s = np.array([r['B']+r['K']+r['H'] for r in row])
    cr = cross(s, 2.3/25.0)
    if not cr: continue
    lo = min(cr); mm = round(float(m),6)
    o = oldb.get(mm)
    rows.append((mm, o, lo))
for mm, o, lo in sorted(rows):
    print("%9.5f %12s %12.3e %7s" % (mm, ("%.3e" % o) if o else "-", lo, ("%.3f" % (lo/o)) if o else "-"))
