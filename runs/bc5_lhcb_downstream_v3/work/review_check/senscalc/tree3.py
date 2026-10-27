import json,os
d=r'E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3/work/physics_review_20261025/'
L=r'E:/AAAwork/LLP/code/SensCalc-main/'
t=json.load(open(d+'reference_2023_tree.json'))
for e in t['tree']:
    p=e['path']
    if ('Scalar' in p or 'FIP' in p or 'DoubleDistr_LHC_h' in p or 'ctauS' in p or 'BrRatioScalar' in p or 'LLP distribution' in p):
        print('%-95s %10d %s'%(p,e['size'],e['sha'][:10]))
print()
print('--- local sizes ---')
for p in ['1. Acceptances.nb','2. LLP distribution.nb','3. Scalar sensitivity.nb','codes/LLP distribution/prod-pheno-Scalar.nb','phenomenology/Scalar/decay widths/BrRatioScalar.mx','phenomenology/Scalar/decay widths/ctauS.mx','spectra/SM particles/DoubleDistr_LHC_h.txt']:
    f=os.path.join(L,p)
    print('%-95s %10d'%(p,os.path.getsize(f)))
