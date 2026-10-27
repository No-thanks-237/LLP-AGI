import numpy as np, json
d=r'E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3/work/physics_review_20261025/'
a=np.loadtxt(d+'reference_2023_Higgs.txt'); b=np.loadtxt(r'E:/AAAwork/LLP/code/SensCalc-main/spectra/SM particles/DoubleDistr_LHC_h.txt')
print('2023 h spectrum',a.shape,'cur',b.shape)
print('shape equal:',a.shape==b.shape)
if a.shape==b.shape:
    print('theta col identical:',np.array_equal(a[:,0],b[:,0]))
    print('E col identical:',np.array_equal(a[:,1],b[:,1]), ' first E 2023 %.4f cur %.4f'%(a[0,1],b[0,1]))
    diff=np.abs(a[:,2]-b[:,2])
    print('max |diff| weights %.3e'%diff.max())
c=np.loadtxt(d+'reference_2023_ctau.dat'); ms=c[:,0]
print('ctau table: %d rows, m %.4f-%.4f'%(len(ms),ms[0],ms[-1]))
print(' unique steps:', np.unique(np.round(np.diff(ms),6))[:12])
for m in [5.082255,13.023095,33.371213,53.419632,62.49]:
    i=np.argmin(np.abs(ms-m)); ct=c[i,1]; G=1.973269804593025e-16/ct
    print('  m=%.4f ctau=%.6e m  Gamma=%.6e GeV'%(c[i,0],ct,G))