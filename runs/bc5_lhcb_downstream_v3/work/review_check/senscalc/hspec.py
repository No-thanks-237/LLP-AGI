
import numpy as np
p=r'E:/AAAwork/LLP/code/SensCalc-main/spectra/SM particles/DoubleDistr_LHC_h.txt'
d=np.loadtxt(p)
print('shape',d.shape)
c1=np.unique(d[:,0]); c2=np.unique(d[:,1])
print('col1 n=%d min=%.6g max=%.6g'%(len(c1),c1[0],c1[-1]))
print('col1 sample:', c1[:6], c1[-4:])
print('col2 n=%d min=%.6g max=%.6g'%(len(c2),c2[0],c2[-1]))
print('col2 sample:', c2[:6], c2[-4:])
print('col3 sum=%.4f min=%.3g max=%.4g'%(d[:,2].sum(), d[:,2].min(), d[:,2].max()))
# is col1 log10(theta)? check spacing
print('col1 diffs:', np.diff(c1)[:8])
print('col2 diffs:', np.diff(c2)[:8])
# rows per col1?
import collections
cnt=collections.Counter(d[:,0])
print('rows per col1 value (first 5):', [cnt[x] for x in c1[:5]])
