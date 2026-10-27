
import numpy as np, math
f = r'temp/foresee/files/hadrons/14TeV/Pythia8/Pythia8_14TeV_25.txt'
hdr = {}
data=[]
for line in open(f):
    if line.startswith('#'):
        continue
    data.append([float(x) for x in line.split()])
d = np.array(data)
print('rows', d.shape)
lth, lp, w = d[:,0], d[:,1], d[:,2]
print('unique log10theta:', len(np.unique(lth)), 'range', lth.min(), lth.max())
print('unique log10p    :', len(np.unique(lp)), 'range', lp.min(), lp.max())
lt = np.unique(lth); lpv = np.unique(lp)
print('dlogth', np.diff(lt)[:3], '... min/max', np.diff(lt).min(), np.diff(lt).max())
print('dlogp ', np.diff(lpv)[:3], '... min/max', np.diff(lpv).min(), np.diff(lpv).max())
print('sum weights = %.6f pb' % w.sum())
print('nonzero bins:', int((w>0).sum()), '/', len(w))
print('max weight %.4e at logth=%.4f logp=%.4f' % (w.max(), lth[w.argmax()], lp[w.argmax()]))
# momentum range where nonzero
nz = w>0
print('nonzero logp range %.3f..%.3f  -> p %.3f..%.1f' % (lp[nz].min(), lp[nz].max(), 10**lp[nz].min(), 10**lp[nz].max()))
print('nonzero theta range %.4f..%.4f -> th %.3e..%.3f' % (lth[nz].min(), lth[nz].max(), 10**lth[nz].min(), 10**lth[nz].max()))
# nonzero pT
th = 10**lth[nz]; p=10**lp[nz]
pt = p*np.sin(th)
print('nonzero pT range %.3f .. %.1f GeV' % (pt.min(), pt.max()))
# check if zero for p<60
for cut in (10,30,50,60,62):
    m = lp < math.log10(cut)
    print('  sum w for p<%.0f GeV: %.4e   (nonzero %d)' % (cut, w[m].sum(), (w[m]>0).sum()))
# eta coverage
eta = -np.log(np.tan(th/2))
print('nonzero eta range: %.3f .. %.3f' % (eta.min(), eta.max()))
print('sum w with eta in 2..5: %.4f' % w[nz & (eta>2) & (eta<5)].sum())
