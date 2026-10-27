
import json, numpy as np
RUN = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
scan = json.load(open(RUN+"/outputs/data/signal_scan.json"))
sig = scan['signal_per_fb_per_channel']
ks = sorted(float(k) for k in sig)
print(ks[:5], '...', len(ks))
thetas = np.array(scan['theta2'])
def get(m):
    k = min(sig.keys(), key=lambda x: abs(float(x)-m))
    return sig[k], float(k)
for m in [0.1008, 0.3021, 0.3541, 0.4835, 1.0, 1.5, 2.0]:
    row, mk = get(m)
    print("\n### m=%s"%mk)
    for i in [0, 40, 80, 120, 160, 200, 239]:
        r = row[i]
        print("  t2=%.2e  B=%.4e K=%.4e H=%.4e"%(thetas[i], r.get('B',0), r.get('K',0), r.get('H',0)))
    Ks = np.array([r.get('K',0) for r in row]); Bs=np.array([r.get('B',0) for r in row])
    imax = int(np.argmax(Ks))
    print("  K max = %.4e at t2=%.3e ; B there=%.4e  K/(B+K+H)=%.3e"%(Ks[imax], thetas[imax], Bs[imax], Ks[imax]/max(Ks[imax]+Bs[imax]+[r.get('H',0) for r in row][imax],1e-300)))
