
import numpy as np, sys
R = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, R+'/work/analysis'); sys.path.insert(0, R+'/work/response'); sys.path.insert(0,R+'/work/widths')
import bc5_lib as L
d = np.load(R+'/work/review_check/kaon/f04_m1.5000.npz')
for k in ('BAc','BAn','BAa','Bmeta'):
    a = d[k]; print(k, a.shape, 'min',a.min(),'max',a.max(), 'mean',a.mean())
print('meta sample rows:')
print(d['Bmeta'][:3])
print('Z0 range', d['Bmeta'][:,2].min(), d['Bmeta'][:,2].max())
print('pz range', d['Bmeta'][:,1].min(), d['Bmeta'][:,1].max())
