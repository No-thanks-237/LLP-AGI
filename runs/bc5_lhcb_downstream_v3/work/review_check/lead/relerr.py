
import json, os, math
import numpy as np
RUN = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3"
scan = json.load(open(os.path.join(RUN, "outputs", "data", "signal_scan.json")))
th = np.array(scan["theta2"])
L = 25.0
print("%10s %6s %10s %10s %10s %10s" % ("m", "nsamp", "lo", "relB", "relK", "relH"))
for m in ["1.0", "2.0", "3.174881", "5.082255", "7.5", "10.5834", "15.234347", "28.527411", "45.66582", "62.49"]:
    row = scan["signal_per_fb_per_channel"][m]
    s = np.array([r['B']+r['K']+r['H'] for r in row])
    f = np.log10(np.maximum(s,1e-300)) - math.log10(2.3/L)
    lo = None
    for i in range(len(th)-1):
        if f[i]*f[i+1] < 0:
            lo = 10**(np.log10(th[i]) + f[i]/(f[i]-f[i+1])*(np.log10(th[i+1])-np.log10(th[i])))
            break
    if lo is None:
        print("%10s %6s  no crossing" % (m, scan['nsamp_by_mass'].get(m))); continue
    i = int(np.argmin(np.abs(np.log10(th)-np.log10(lo))))
    r = row[i]
    print("%10s %6.0f %10.3e %10.2f %10.2f %10.2f   (B=%.3g K=%.3g H=%.3g)" % (
        m, scan['nsamp_by_mass'].get(m, 0), lo, r['relerr_B'], r['relerr_K'], r['relerr_H'],
        r['B'], r['K'], r['H']))
