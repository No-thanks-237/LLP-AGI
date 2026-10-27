
import json, os, sys, math
import numpy as np
RUN = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3"
sys.path.insert(0, os.path.join(RUN, "work", "analysis"))
scan = json.load(open(os.path.join(RUN, "outputs", "data", "signal_scan.json")))
for m in ["5.082255", "28.527411"]:
    row = scan["signal_per_fb_per_channel"][m]
    print("m", m, "n_theta", len(row), "keys", sorted(row[0].keys()))
    for i in [0, 60, 120, 180, 239]:
        t2 = scan["theta2"][i]
        print("  t2=%.3e B=%.4g K=%.4g H=%.4g sum=%.4g" % (t2, row[i]['B'], row[i]['K'], row[i]['H'],
              row[i]['B']+row[i]['K']+row[i]['H']))
    print("  nsamp", scan["nsamp_by_mass"].get(m))
