
import numpy as np, math
f = r'temp/foresee/files/hadrons/14TeV/Pythia8/Pythia8_14TeV_25.txt'
data=[]
for line in open(f):
    if line.startswith('#'): continue
    data.append([float(x) for x in line.split()])
d = np.array(data); lth,lp,w = d[:,0],d[:,1],d[:,2]
nz = w>0
print('nonzero bins with p < 30 GeV:')
idx = np.where(nz & (lp<math.log10(30)))[0]
order = idx[np.argsort(-w[idx])]
print('%-10s %-10s %-12s %-12s %-10s' % ('log10th','log10p','p','theta','w'))
for i in order[:25]:
    print('%-10.4f %-10.4f %-12.4g %-12.4g %-10.4g' % (lth[i], lp[i], 10**lp[i], 10**lth[i], w[i]))
print()
# total weight in these
print('count p<30 nonzero:', len(idx), ' sum w = %.4e' % w[idx].sum())
# what about require eta in 2..5 for all bins
th = 10**lth
eta = -np.log(np.tan(th/2))
print('sum w with eta in (2,5): %.5f pb' % w[nz & (eta>2) & (eta<5)].sum())
print('sum w with eta > 2      : %.5f pb' % w[nz & (eta>2)].sum())
print('fraction of total forward cross section with 2<eta<5: %.4f' % (w[nz & (eta>2) & (eta<5)].sum()/w.sum()))
# check the pT spectrum: for each bin pT=p*sin(th)
p = 10**lp; pt = p*np.sin(th)
print()
print('cross section vs pT (sum w):')
for lo,hi in [(0,5),(5,10),(10,20),(20,30),(30,50),(50,80),(80,120),(120,200),(200,400),(400,1e9)]:
    m = nz & (pt>=lo) & (pt<hi)
    print('  pT %6.1f-%7.1f : %10.5f pb' % (lo,hi,w[m].sum()))
print()
print('cross section vs eta:')
for lo,hi in [(-13,0),(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,13)]:
    m = nz & (eta>=lo) & (eta<hi)
    print('  eta %6.1f..%5.1f : %10.5f pb' % (lo,hi,w[m].sum()))
print()
print('cross section vs p:')
for lo,hi in [(0,62),(62,100),(100,150),(150,300),(300,600),(600,1e9)]:
    m = nz & (p>=lo) & (p<hi)
    print('  p %6.1f..%7.1f : %10.5f pb' % (lo,hi,w[m].sum()))
