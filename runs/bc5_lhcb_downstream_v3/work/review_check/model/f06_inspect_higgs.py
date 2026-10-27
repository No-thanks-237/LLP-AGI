import numpy as np, collections, math, os
RUN=r'E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3'
p=os.path.join(RUN,'work','magnus','samples','higgs.txt')
ev=[];code=[];px=[];py=[];pz=[];e=[]
with open(p) as f:
    for line in f:
        if line.startswith('#'):
            print('HEADER:',line.strip()); continue
        v=line.split()
        ev.append(int(v[0])); code.append(int(v[1]))
        px.append(float(v[2]));py.append(float(v[3]));pz.append(float(v[4]));e.append(float(v[5]))
ev=np.array(ev);code=np.array(code);px=np.array(px);py=np.array(py);pz=np.array(pz);e=np.array(e)
print('rows',len(ev),'events',len(np.unique(ev)))
u,start,counts=np.unique(ev,return_index=True,return_counts=True)
print('unique events',len(u),'min/max copies',counts.min(),counts.max(),'mean',counts.mean())
m=np.sqrt(np.maximum(e**2-px**2-py**2-pz**2,0))
pt=np.hypot(px,py)
first=start; last=start+counts-1
for name,idx in (('first',first),('last',last)):
    print(name,'E mean %.4f  mass mean %.6f std %.6f  pT mean %.4f  pz mean %.4f'%(
        e[idx].mean(),m[idx].mean(),m[idx].std(),pt[idx].mean(),pz[idx].mean()))
d=(e[last]-e[first]); print('E last-first mean %.4f  maxabs %.4f'%(d.mean(),np.abs(d).max()))
same=(np.abs(e[first]-e[last])<1e-6)&(np.abs(px[first]-px[last])<1e-6)
print('first==last fourmom frac',same.mean())
nmax=counts.max()
for pos in range(0,min(int(nmax),9)):
    sel=start[counts>pos]+pos
    print('pos',pos,'n',len(sel),'meanE %.3f mass %.5f pT %.3f'%(e[sel].mean(),m[sel].mean(),pt[sel].mean()),'codes',sorted(set(code[sel].tolist()))[:8])
# distribution of copy count vs subprocess code
print('subprocess codes',collections.Counter(code.tolist()).most_common())
# for a few events print all copies
for k in range(3):
    i=start[k]; n=counts[k]
    print('event',ev[i],'ncopies',n)
    for j in range(i,i+n):
        print('  code',code[j],'E %.4f px %.4f py %.4f pz %.4f m %.5f'%(e[j],px[j],py[j],pz[j],m[j]))
