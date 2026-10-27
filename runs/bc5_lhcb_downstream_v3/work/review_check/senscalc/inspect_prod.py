
import numpy as np, os
base=r'E:/AAAwork/LLP/code/SensCalc-main/phenomenology/Scalar/Production probabilities'
for fn in ['BtoXsStotal.dat','BtoXsSinclusive.dat','BtoKplusS.dat','BtoKstar892S.dat']:
    d=np.loadtxt(os.path.join(base,fn))
    print('==',fn,'n=',len(d),'m_min=%.3f m_max=%.3f'%(d[0,0],d[-1,0]),'step=%.4f'%(d[1,0]-d[0,0]))
    print('   start:', d[0], ' end:', d[-1])
d1=np.loadtxt(os.path.join(base,'BtoXsStotal.dat'))
d2=np.loadtxt(os.path.join(base,'BtoXsSinclusive.dat'))
# compare on common grid
m1,v1=d1[:,0],d1[:,1]; m2,v2=d2[:,0],d2[:,1]
print('\nratio inclusive/exclusive (interp on exclusive grid):')
for m in [0.0,0.5,1.0,2.0,3.0,4.0,4.3,4.5,4.7,4.8,4.9,5.0,5.1]:
    a=np.interp(m,m1,v1); b=np.interp(m,m2,v2)
    print('  m=%.2f  excl=%.5f  incl=%.5f  ratio=%.4f'%(m,a,b,b/a if a else float('nan')))
print('\nexclusive values near 4.3-4.4:')
for m in np.arange(4.25,4.46,0.02):
    print('  m=%.3f  excl=%.6f  incl=%.6f'%(m,np.interp(m,m1,v1),np.interp(m,m2,v2)))
print('\nexclusive/inclusive near 4.6-4.9:')
for m in np.arange(4.55,4.96,0.05):
    print('  m=%.3f  excl=%.6f  incl=%.6f'%(m,np.interp(m,m1,v1),np.interp(m,m2,v2)))
