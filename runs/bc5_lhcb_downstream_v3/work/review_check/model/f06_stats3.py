import numpy as np, os, collections
RUN=r'E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3'
d=np.loadtxt(os.path.join(RUN,'work','magnus','samples','higgs.txt'),comments='#')
ev=d[:,0].astype(int); cd=d[:,1].astype(int)
px=d[:,2];py=d[:,3];pz=d[:,4];e=d[:,5]
m=np.sqrt(np.maximum(e**2-px**2-py**2-pz**2,0))
u,first,counts=np.unique(ev,return_index=True,return_counts=True)
pt=np.hypot(px,py)
print('mass quantiles',np.round(np.percentile(m,[0.1,1,5,50,95,99,99.9]),5),'std %.5f'%m.std())
for n in [1,2,3,5,6,8]:
    sel=counts==n
    f=first[sel]
    print('ncopies',n,'nev',int(sel.sum()),'pT pos0 mean %.1f median %.1f'%(pt[f].mean(),np.median(pt[f])),
          'pos-last mean %.1f'%pt[first[sel]+n-1].mean())
# Is the sequence monotone in pT?
mono=0; tot=0
for i in range(0,len(u),7):
    s=first[i]; n=counts[i]
    if n>2:
        tot+=1
        p=pt[s:s+n]
        if np.all(np.diff(p)>0): mono+=1
print('strictly increasing pT seq fraction',mono/max(tot,1),'of',tot)
# print 8 events of n=6 all copies pT and E
shown=0
for i in range(len(u)):
    if counts[i]==6 and shown<6:
        s=first[i]
        print('event',u[i],'pT',np.round(pt[s:s+6],2),'E',np.round(e[s:s+6],2),'code',cd[s])
        shown+=1
