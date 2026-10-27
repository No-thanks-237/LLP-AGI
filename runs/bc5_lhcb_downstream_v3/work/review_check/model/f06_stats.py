import numpy as np, collections, os
RUN=r'E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3'
d=np.loadtxt(os.path.join(RUN,'work','magnus','samples','higgs.txt'),comments='#')
ev=d[:,0].astype(int); code=d[:,1].astype(int)
px=d[:,2];py=d[:,3];pz=d[:,4];e=d[:,5]
u,first,counts=np.unique(ev,return_index=True,return_counts=True)
last=first+counts-1
pt=np.hypot(px,py)
for nm,idx in (('first',first),('last',last)):
    q=np.percentile(pt[idx],[1,5,10,25,50,75,90,95,99])
    print(nm,'pT percentiles',np.round(q,2),'mean %.2f'%pt[idx].mean())
# copies distribution
cnt=collections.Counter(counts.tolist())
print('copy count histogram',sorted(cnt.items()))
# by subprocess code: count of events, mean pT first
for c in sorted(set(code.tolist())):
    evsel=u[np.array([code[f]==c for f in first])]
    print('code',c,'events',len(evsel))
# per code mean first pT and copies
m=collections.defaultdict(list)
for i,f in enumerate(first):
    m[code[f]].append((pt[f], counts[i]))
for c in sorted(m):
    a=np.array(m[c]); print('code %d n=%d mean pT_first %.2f mean copies %.2f'%(c,len(a),a[:,0].mean(),a[:,1].mean()))
# correlation copies vs pT
print('corr copies, pT_first', np.corrcoef(counts, pt[first])[0,1])
# check last copy relative momentum diff
print('mean |dpt| first->last %.2f  mean |dp| %.2f'%(np.abs(pt[last]-pt[first]).mean(),
    np.linalg.norm(d[last,2:5]-d[first,2:5],axis=1).mean()))
# distribution of ncopies for events with 1 copy
print('events with 1 copy',int((counts==1).sum()))
