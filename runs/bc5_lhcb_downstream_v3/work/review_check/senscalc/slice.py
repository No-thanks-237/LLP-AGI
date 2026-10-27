
import numpy as np, os
base=r'E:/AAAwork/LLP/code/SensCalc-main/Tabulated Nevents/Scalar/LHCb-downstream-D'
fn=[f for f in os.listdir(base) if f.startswith('Nevents_Scalar_h-quartic')][0]
d=np.loadtxt(os.path.join(base,fn))
ms=np.unique(d[:,0]); th=np.unique(d[:,1])
V=d[:,2].reshape(len(ms),len(th))
print('shape',V.shape,'max',V.max(),'min',V.min())
im=np.argmin(np.abs(ms-13.0))
print('m=',ms[im])
for j in range(0,len(th),100):
    print('  th2=%.3e N=%.4g'%(th[j],V[im,j]))
print(' last couplings:')
for j in range(len(th)-6,len(th)):
    print('  th2=%.3e N=%.4g'%(th[j],V[im,j]))
