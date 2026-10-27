import numpy as np, os
RUN=r'E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3'
d=np.loadtxt(os.path.join(RUN,'work','magnus','samples','higgs.txt'),comments='#')
ev=d[:,0].astype(int); cd=d[:,1].astype(int)
px=d[:,2];py=d[:,3];pz=d[:,4];e=d[:,5]
u,first,counts=np.unique(ev,return_index=True,return_counts=True)
pt=np.hypot(px,py); p=np.sqrt(px**2+py**2+pz**2)
eta=0.5*np.log((p+pz)/np.maximum(p-pz,1e-300))
idx=first[cd[first]==902]
print('code902 first-copy pT percentiles',np.round(np.percentile(pt[idx],[1,5,10,25,50,75,90,99]),1))
print('code902 first-copy |p| percentiles',np.round(np.percentile(p[idx],[1,5,10,25,50,75,90,99]),1))
print('code902 first-copy eta percentiles',np.round(np.percentile(eta[idx],[1,5,10,25,50,75,90,99]),2))
print('fraction pT<20 of code902 first',np.mean(pt[idx]<20))
# distribution of copy multiplicity vs pT for code 902
cc=counts[cd[first]==902]
print('mean copies code902 by pT quartile:')
qs=np.percentile(pt[idx],[25,50,75])
for lo,hi in zip([-1]+list(qs),list(qs)+[1e9]):
    s=(pt[idx]>lo)&(pt[idx]<=hi); print('  %.1f-%.1f n=%d mean copies %.2f'%(lo,hi,s.sum(),cc[s].mean()))
# check the first copy pT vs last for code 902
last=first+counts-1
print('code902 pT first %.2f last %.2f'%(pt[first[cd[first]==902]].mean(),pt[last[cd[first]==902]].mean()))
# consecutive differences
sel=cd[first]==902
f=first[sel]; l=last[sel]
dpt=pt[l]-pt[f]
print('code902 dpt last-first mean %.2f median %.2f frac increase %.3f'%(dpt.mean(),np.median(dpt),(dpt>0).mean()))
