
import numpy as np, os
base=r'E:/AAAwork/LLP/code/SensCalc-main/phenomenology/Scalar/Production probabilities'
def load(f):
    d=np.loadtxt(os.path.join(base,f)); return d[:,0],d[:,1]
mS,eS=load('BtoXsStotal.dat'); mK,eK=load('BtoKplusS.dat'); mKs,eKs=load('BtoKstar892S.dat')
print(' m      sum    B->K    B->K*   sum-K-K*')
for m in [0.,1.,2.,3.,3.5,4.0,4.3,4.38,4.4,4.5,4.7,4.75,4.78,4.79,4.9,5.0,5.1]:
    print(' %5.2f  %7.4f %7.4f %7.4f  %8.4f'%(m,np.interp(m,mS,eS),np.interp(m,mK,eK),np.interp(m,mKs,eKs),np.interp(m,mS,eS)-np.interp(m,mK,eK)-np.interp(m,mKs,eKs)))
