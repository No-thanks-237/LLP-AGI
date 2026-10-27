
import sys
import numpy as np
sys.path.insert(0, 'runs/bc5_lhcb_downstream_v3/work/analysis')
import run_scan as R
smp = R.load_sample('higgs.txt.gz')
cand = np.where(np.abs(smp['status'])==62)[0]
E,px,py,pz = smp['e'][cand],smp['px'][cand],smp['py'][cand],smp['pz'][cand]
m = np.sqrt(np.maximum(E*E-px*px-py*py-pz*pz,0))
print('mH: min %.4f max %.4f mean %.5f  median %.5f' % (m.min(), m.max(), m.mean(), np.median(m)))
print('percentiles [0.1,1,5,50,95,99,99.9]:', np.round(np.percentile(m,[0.1,1,5,50,95,99,99.9]),4))
for ms in (5.08,13.02,33.37,53.42):
    n = int(np.sum(m < 2*ms))
    print('m_S=%.2f: events with mH < 2mS = %d  (fraction %.2e)' % (ms, n, n/len(m)))
