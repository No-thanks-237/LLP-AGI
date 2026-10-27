
import json, numpy as np, math
R = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
scan = json.load(open(R+'/outputs/data/signal_scan.json'))
thetas = np.array(scan['theta2']); sc = scan['signal_per_fb_per_channel']
def get(m):
    k = min(sc.keys(), key=lambda x: abs(float(x)-m)); return float(k), sc[k]
def crossings(x,y,level):
    out=[]; lx=np.log10(x); ly=np.log10(np.maximum(y,1e-300))
    for i in range(len(x)-1):
        f0=ly[i]-math.log10(level); f1=ly[i+1]-math.log10(level)
        if f0==0: out.append(x[i])
        elif f0*f1<0:
            t=f0/(f0-f1); out.append(10**(lx[i]+t*(lx[i+1]-lx[i])))
    return out
for m in [1.0,1.2,1.4494,1.6955,1.9833,2.0]:
    mk,row = get(m)
    s1 = np.array([r.get('B',0)+r.get('K',0)+r.get('H',0) for r in row])
    cr = crossings(thetas, s1*25.0, 2.3)
    print('m=%.4f crossings(25fb, s=2.3): %s'%(mk, ['%.3e'%c for c in cr]))
    print('   max s1*25 = %.3e at t2=%.3e   s1 at t2=1e-5: %.3e'%(np.max(s1*25), thetas[np.argmax(s1)], s1[np.argmin(abs(thetas-1e-5))]*25))
    print('   first/last 6 theta and s:', ['%.2e:%.2e'%(thetas[i], s1[i]*25) for i in range(0,6)], '...', ['%.2e:%.2e'%(thetas[i], s1[i]*25) for i in range(len(thetas)-4, len(thetas))])
