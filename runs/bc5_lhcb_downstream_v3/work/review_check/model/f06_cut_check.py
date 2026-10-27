import numpy as np, os, json
RUN=r'E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3'
d=np.loadtxt(os.path.join(RUN,'work','magnus','samples','higgs.txt'),comments='#')
e=d[:,5]; px=d[:,2]; py=d[:,3]; pz=d[:,4]
pt=np.hypot(px,py)
print('rows',len(e),'min pT %.6g'%pt.min(),'rows with pT<1e-4:',int((pt<1e-4).sum()))
print('rows with |E-125.25|<0.01:',int((abs(e-125.25)<0.01).sum()),' min E %.4f'%e.min())
ev=d[:,0].astype(int); u,first,counts=np.unique(ev,return_index=True,return_counts=True)
last=first+counts-1
print('first-copy pT<1e-4 frac',float((pt[first]<1e-4).mean()),'last-copy pT<1e-4 frac',float((pt[last]<1e-4).mean()))
print('first-copy E<126 frac',float((e[first]<126).mean()))
# the hard-process copy would have E=125.25, pT=0; both would be cut by |eta|<6 -> show none
