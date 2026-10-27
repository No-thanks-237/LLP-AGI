
import numpy as np
base=r'E:/AAAwork/LLP/code/SensCalc-main/spectra/SM particles/'
def etadist(fn):
    d=np.loadtxt(base+fn); th=np.unique(d[:,0]); E=np.unique(d[:,1]); V=d[:,2].reshape(len(th),len(E))
    eta=-np.log(np.tan(th/2)); o=np.argsort(eta); eta=eta[o]; V=V[o,:]
    g=np.array([np.trapezoid(V[i,:],E) for i in range(len(th))])
    dthdeta=np.abs(2*np.exp(-eta)/(1+np.exp(-2*eta)))
    dn=g*dthdeta; dn=dn/np.trapezoid(dn,eta)
    return eta,dn,E,V
for fn in ['DoubleDistr_LHC_h.txt','DoubleDistr_LHC_Bcharged.txt','DoubleDistr_LHC_Upsilon.txt']:
    eta,dn,E,V=etadist(fn)
    s=''
    for e0 in [0,1,2,3,4,5,6]:
        i=np.argmin(np.abs(eta-e0)); s+=' %.3f'%dn[i]
    fr=[]
    for lo,hi in [(2,5),(2,3),(3,5)]:
        m=(eta>=lo)&(eta<=hi); fr.append(np.trapezoid(dn[m],eta[m]))
    print('%-42s dN/deta(0,1,2,3,4,5,6) =%s ; frac[2,5]=%.3f frac[2,3]=%.3f frac[3,5]=%.3f ; Egrid %.0f-%.0f'%(fn,s,fr[0],fr[1],fr[2],E[0],E[-1]))
