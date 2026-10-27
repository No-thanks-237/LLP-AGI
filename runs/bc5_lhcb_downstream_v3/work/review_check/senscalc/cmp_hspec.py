
import numpy as np, gzip
px=[];py=[];pz=[];E=[]
with gzip.open(r'E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3/work/magnus/samples/higgs.txt.gz','rt') as f:
    for line in f:
        if line.startswith('#'): continue
        a=line.split()
        if int(a[10])!=62: continue
        px.append(float(a[2])); py.append(float(a[3])); pz.append(float(a[4])); E.append(float(a[5]))
px=np.array(px);py=np.array(py);pz=np.array(pz);E=np.array(E)
print('n status62 =',len(E))
pt=np.hypot(px,py)
eta=0.5*np.log((E+pz)/(E-pz))
print('E range %.1f-%.1f, med %.1f'%(E.min(),E.max(),np.median(E)))
print('eta range %.2f-%.2f, median %.2f'%(eta.min(),eta.max(),np.median(eta)))
print('pT median %.1f mean %.1f'%(np.median(pt),pt.mean()))
bins=np.linspace(-6,8,57); h,_=np.histogram(eta,bins=bins,density=True); c=0.5*(bins[1:]+bins[:-1])
p=r'E:/AAAwork/LLP/code/SensCalc-main/spectra/SM particles/DoubleDistr_LHC_h.txt'
d=np.loadtxt(p); th=np.unique(d[:,0]); EE=np.unique(d[:,1]); V=d[:,2].reshape(len(th),len(EE))
etaT=-np.log(np.tan(th/2)); order=np.argsort(etaT); etaT=etaT[order]; V=V[order,:]
g=np.array([np.trapezoid(V[i,:],EE) for i in range(len(th))])
dthdeta=2*np.exp(-etaT)/(1+np.exp(-2*etaT)); dNdeta=g*np.abs(dthdeta)
dNdeta=dNdeta/np.trapezoid(dNdeta,etaT)
print('\n eta   ours_dN/deta   SensCalc_dN/deta   ratio')
for e0 in [1.,2.,3.,4.,5.,6.]:
    io=np.argmin(np.abs(c-e0)); it=np.argmin(np.abs(etaT-e0))
    print('%5.2f  %12.4e  %12.4e   %.3f'%(e0,h[io],dNdeta[it],h[io]/dNdeta[it]))
for lo,hi in [(2,5),(2,3),(3,5)]:
    mo=(c>=lo)&(c<=hi); mt=(etaT>=lo)&(etaT<=hi)
    fo=np.trapz(h[mo],c[mo]); ft=np.trapezoid(dNdeta[mt],etaT[mt])
    print('fraction eta in [%g,%g]: ours %.4f  SensCalc %.4f'%(lo,hi,fo,ft))
