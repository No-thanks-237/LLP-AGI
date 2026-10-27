
import numpy as np, os, sys
base=r'E:/AAAwork/LLP/code/SensCalc-main/phenomenology/Scalar/Production probabilities'
def load(f): 
    d=np.loadtxt(os.path.join(base,f)); return d[:,0],d[:,1]
m1,e1=load('BtoXsStotal.dat'); m2,e2=load('BtoXsSinclusive.dat')
def last_nonzero(m,v):
    nz=np.nonzero(v>0)[0]
    return m[nz[-1]] if len(nz) else None
print('exclusive last nonzero m =',last_nonzero(m1,e1), ' first zero-run from', m1[np.argmax(np.cumsum(e1[::-1]>0)[::-1]>0)] if False else '')
print('inclusive last nonzero m =',last_nonzero(m2,e2))
print()
print('exclusive near 4.76-4.82:')
for m in np.arange(4.76,4.83,0.01):
    print('  m=%.2f excl=%.6f incl=%.6f'%(m,np.interp(m,m1,e1),np.interp(m,m2,e2)))
print()
print('inclusive near 4.10-4.25:')
for m in np.arange(4.10,4.26,0.02):
    print('  m=%.2f incl=%.6e excl=%.6f'%(m,np.interp(m,m2,e2),np.interp(m,m1,e1)))
# compare shape of inclusive with our formula
sys.path.insert(0, r'E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3/work/analysis')
import bc5_production as bp
ms=[0.,0.5,1.,2.,3.,4.]
print()
print(' m   incl_table  incl/ours  excl_table  excl/ours')
for m in ms:
    o=bp.br_B_to_S(m,'B+')['B_incl']
    print(' %.2f  %.4f  %.4f   %.4f  %.3f'%(m,np.interp(m,m2,e2),np.interp(m,m2,e2)/o,np.interp(m,m1,e1),np.interp(m,m1,e1)/o))
