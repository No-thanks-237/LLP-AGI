
import json, os
import numpy as np
RUN = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3"
a = json.load(open(os.path.join(RUN, "outputs", "data", "signal_scan.json")))
b = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "scan_check.json")))
print("theta grids identical:", np.allclose(a["theta2"], b["theta2"]))
for m in ["5.082255", "28.527411"]:
    ra = a["signal_per_fb_per_channel"][m]; rb = b["signal_per_fb_per_channel"][m]
    for i in [0, 60, 120, 180, 200]:
        pa = ra[i]; pb = rb[i]
        print("m=%s i=%3d  stored B=%.4g K=%.4g H=%.4g | rerun B=%.4g K=%.4g H=%.4g  ratioH=%.4f" % (
            m, i, pa['B'], pa['K'], pa['H'], pb['B'], pb['K'], pb['H'],
            pb['H']/pa['H'] if pa['H'] else float('nan')))
