
import sys, os, json, math
import numpy as np
sys.path.insert(0, 'runs/bc5_lhcb_downstream_v3/work/analysis')
import run_scan as R
smp = R.load_sample('higgs.txt.gz')
print('keys:', list(smp.keys()))
print('rows:', len(smp['ev']))
for k in ('ev','pdg','px','py','pz','e','x','y','z','t','status'):
    if k in smp:
        a=smp[k]
        print(k, 'len', len(a), 'min', a.min() if len(a) else None, 'max', a.max() if len(a) else None, 'dtype', a.dtype)
print()
print('unique events:', len(np.unique(smp['ev'])))
print('rows with |status|==62:', int(np.sum(np.abs(smp['status'])==62)))
print('status value counts (top 20):')
u,c = np.unique(smp['status'], return_counts=True)
for uu,cc in sorted(zip(u.tolist(),c.tolist()), key=lambda z:-z[1])[:20]:
    print('  status %6d : %d' % (uu,cc))
print()
print('unique pdg:', np.unique(smp['pdg'])[:20], '... n=', len(np.unique(smp['pdg'])))
cand = np.where(np.abs(smp['status'])==62)[0]
evc = smp['ev'][cand]
uev, cnt = np.unique(evc, return_counts=True)
print('events with a 62 row:', len(uev), ' per-event 62-row count: min %d max %d mean %.4f' % (cnt.min(), cnt.max(), cnt.mean()))
print('total events in file:', len(np.unique(smp['ev'])))
print('events WITHOUT a 62 row:', len(np.unique(smp['ev']))-len(uev))
print('pdg of 62 rows:')
uu,cc = np.unique(smp['pdg'][cand], return_counts=True)
for a,b in zip(uu,cc): print('   pdg %6d: %d' % (a,b))
