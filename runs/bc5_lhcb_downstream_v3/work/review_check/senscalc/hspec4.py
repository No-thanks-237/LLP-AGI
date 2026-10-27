
import numpy as np
p=r'E:/AAAwork/LLP/code/SensCalc-main/spectra/SM particles/DoubleDistr_LHC_h.txt'
d=np.loadtxt(p); th=np.unique(d[:,0]); E=np.unique(d[:,1]); V=d[:,2].reshape(len(th),len(E))
i,j=np.unravel_index(np.argmax(V),V.shape)
print('max %.4g at theta=%.4g rad (eta=%.2f) E=%.1f'%(V[i,j],th[i],-np.log(np.tan(th[i]/2)),E[j]))
print('\nE-slice at several theta (weight vs E):')
for it in [0,10,20,30,40,55,70,90]:
    row=V[it,:]
    k=np.argmax(row)
    print(' theta=%.4g (eta=%6.2f): peak at E=%.1f w=%.4g ; w at E=200,1000,7000: %.3g %.3g %.3g'%(th[it],-np.log(np.tan(th[it]/2)),E[k],row[k],row[np.argmin(abs(E-200))],row[np.argmin(abs(E-1000))],row[np.argmin(abs(E-7000))]))
print('\ntheta-slice at several E:')
for je in [0,10,20,30,40,50,59]:
    col=V[:,je]
    k=np.argmax(col)
    print(' E=%.1f: peak at theta=%.4g (eta=%.2f) w=%.4g'%(E[je],th[k],-np.log(np.tan(th[k]/2)),col[k]))
