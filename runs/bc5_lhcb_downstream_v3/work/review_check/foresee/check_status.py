
import sys, math
import numpy as np
sys.path.insert(0, 'runs/bc5_lhcb_downstream_v3/work/analysis')
import run_scan as R
smp = R.load_sample('higgs.txt.gz')
st = smp['status']; ev = smp['ev']
px,py,pz,E = smp['px'],smp['py'],smp['pz'],smp['e']
pt = np.hypot(px,py); p3=np.sqrt(px*px+py*py+pz*pz)
m = np.sqrt(np.maximum(E*E-p3*p3,0))
eta = np.arcsinh(pz/np.maximum(pt,1e-300))
print('status classes:')
for s in np.unique(st):
    sel = st==s
    print('  %4d : N=%8d  <pT>=%8.2f  <eta>=%7.3f  <m>=%8.4f  pt=0 frac %.4f' % (
        s, sel.sum(), pt[sel].mean(), eta[sel].mean(), m[sel].mean(), np.mean(pt[sel]<1e-9)))
# per-event: compare the 62 copy with the last -44 copy
order = np.lexsort((st, ev))
ev_s = ev[order]; st_s = st[order]; pt_s = pt[order]; pz_s=pz[order]; eta_s=eta[order]
idx62 = np.where(st_s==62)[0]
idx44 = np.where(st_s==-44)[0]
ev62 = ev_s[idx62]
# last -44 per event (event ids are contiguous blocks in the sorted array)
last44 = np.ones(len(idx44), bool); last44[:-1] = ev_s[idx44][1:] != ev_s[idx44][:-1]
idx44last = idx44[last44]
print()
print('62 copy : N=%d  <pT>=%.2f  <eta>=%.3f' % (len(idx62), pt_s[idx62].mean(), eta_s[idx62].mean()))
print('last -44: N=%d  <pT>=%.2f  <eta>=%.3f' % (len(idx44last), pt_s[idx44last].mean(), eta_s[idx44last].mean()))
# match by event id
d = dict(zip(ev_s[idx44last].tolist(), np.arange(len(idx44last))))
common = [i for i in idx62 if ev_s[i] in d]
a = np.array([pt_s[i] for i in common]); b = np.array([pt_s[idx44last[d[ev_s[i]]]] for i in common])
print('events with both: %d   max|dpT(62-last44)| = %.3e   <pT> %.3f vs %.3f' % (len(common), np.abs(a-b).max(), a.mean(), b.mean()))
# first row per event (what the OLD code used)
first = np.ones(len(ev), bool); first[1:] = ev[1:] != ev[:-1]
fi = np.where(first)[0]
print()
print('first record row of each event: statuses', dict(zip(*[x.tolist() for x in np.unique(st[fi], return_counts=True)])))
print('  <pT>=%.2f  <eta>=%.3f' % (pt[fi].mean(), eta[fi].mean()))
# 62 copy: momentum conservation check - the sum of all pT in the event should be ~0?
print()
print('62-copy pT quantiles:', np.round(np.percentile(pt[idx62],[5,25,50,75,95]),2))
