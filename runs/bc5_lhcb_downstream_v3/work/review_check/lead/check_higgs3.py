
import sys, os
import numpy as np
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\analysis")
import run_scan as R
smp = R.load_sample(R.SAMPLE_H)
print("rows", len(smp['ev']), "events", len(np.unique(smp['ev'])), "has status:", 'status' in smp)
if 'status' in smp:
    ids, cnt = np.unique(smp['status'], return_counts=True)
    print("status codes:", dict(zip(ids.tolist(), cnt.tolist())))
    sel = np.abs(smp['status']) == 62
    print("rows with |status|=62:", int(sel.sum()), "events:", len(np.unique(smp['ev'][sel])))
    p = np.hypot(smp['px'], smp['py'])
    first = np.unique(smp['ev'], return_index=True)[1]
    print("<pT> first copy %.2f   <pT> |62| copy %.2f" % (p[first].mean(), p[sel].mean()))
