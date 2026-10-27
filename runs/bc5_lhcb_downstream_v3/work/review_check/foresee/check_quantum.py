
import numpy as np, math
FS = 'temp/foresee/files/hadrons/14TeV/Pythia8/Pythia8_14TeV_25.txt'
raw=[]
for line in open(FS):
    if line.startswith('#'): continue
    raw.append([float(x) for x in line.split()])
raw=np.array(raw); lth,lp,w = raw[:,0],raw[:,1],raw[:,2]
nz = w>1e-6
ww = w[nz]
print('n nonzero', len(ww), 'sum', ww.sum()) 
q = ww.min()
print('min weight', q)
print('w/q for the 20 smallest:', np.round(np.sort(ww)[:20]/q, 4))
print('w/q for 10 random bins  :', np.round(np.sort(ww)[-10:]/q, 4))
r = ww/q
dev = np.abs(r-np.round(r))
print('max |w/q - round| / 1 :', dev.max(), ' frac with dev<0.01:', np.mean(dev<0.01))
# try candidate quanta: min, and w.min*? 
for k in [1,2,0.5]:
    qq=q*k; rr=ww/qq; d=np.abs(rr-np.round(rr))
    print('quantum %.4e: max dev %.4f  frac<0.01 %.3f' % (qq, d.max(), np.mean(d<0.01)))
# implied number of events if quantum = sigma_forward/N
sig = w.sum()
print('implied N (sigma_forward/quantum) = %.4e' % (sig/q))
print('implied N if 2 events per quantum = %.4e' % (sig/(2*q)))
# Check: are weights multiples of min over a wide dynamic range?
print('weights sorted desc, w/min for a few:', np.round(np.sort(ww)[::400]/q,3))
