
import json, os, math
import numpy as np
RUN = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3"
scan = json.load(open(os.path.join(RUN, "outputs", "data", "signal_scan.json")))
th = np.array(scan["theta2"]); L = 25.0
print("%10s %8s | %10s %10s %10s   (share of s at the long-lived crossing)" % ("m","share","B","K","H"))
for m in ["0.021","0.1","0.2118","0.5","1.0","2.0","3.174881","5.082255","10.5834","28.527411","62.49"]:
    if m not in scan["signal_per_fb_per_channel"]: 
        print(m, "missing"); continue
    row = scan["signal_per_fb_per_channel"][m]
    s = np.array([r['B']+r['K']+r['H'] for r in row])
    f = np.log10(np.maximum(s,1e-300))-math.log10(2.3/L)
    lo = None
    for i in range(len(th)-1):
        if f[i]*f[i+1] < 0:
            lo = 10**(np.log10(th[i])+f[i]/(f[i]-f[i+1])*(np.log10(th[i+1])-np.log10(th[i]))); break
    if lo is None: 
        print("%10s no crossing" % m); continue
    i = int(np.argmin(np.abs(np.log10(th)-np.log10(lo))))
    r = row[i]; tot = r['B']+r['K']+r['H']
    print("%10s %8.3f | %10.3g %10.3g %10.3g  B%%=%.2f K%%=%.2f H%%=%.2f" % (
        m, lo, r['B'], r['K'], r['H'], 100*r['B']/tot, 100*r['K']/tot, 100*r['H']/tot))
