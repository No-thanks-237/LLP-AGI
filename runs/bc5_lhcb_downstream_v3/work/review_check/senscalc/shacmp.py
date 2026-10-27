import json,os,hashlib
d=r'E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3/work/physics_review_20261025/'
L=r'E:/AAAwork/LLP/code/SensCalc-main/'
t=json.load(open(d+'reference_2023_tree.json'))
blobs={e['path']:e['sha'] for e in t['tree'] if e['type']=='blob'}
def gitsha(path):
    b=open(path,'rb').read()
    h=hashlib.sha1(); h.update(b'blob %d\0'%len(b)); h.update(b); return h.hexdigest()
pairs=[('phenomenology/Scalar/branching ratios/BtoXsStotal.dat','phenomenology/Scalar/Production probabilities/BtoXsStotal.dat'),
('phenomenology/Scalar/branching ratios/BtoXsSinclusive.dat','phenomenology/Scalar/Production probabilities/BtoXsSinclusive.dat'),
('phenomenology/Scalar/branching ratios/BplustoXsSStotal.dat','phenomenology/Scalar/Production probabilities/BplustoXsSStotal.dat'),
('phenomenology/Scalar/branching ratios/BstoSS.dat','phenomenology/Scalar/Production probabilities/BstoSS.dat'),
('phenomenology/Scalar/branching ratios/KplustoPiS.dat','phenomenology/Scalar/Production probabilities/KplustoPiS.dat'),
('phenomenology/Scalar/branching ratios/K0LtoPiS.dat','phenomenology/Scalar/Production probabilities/K0LtoPiS.dat'),
('phenomenology/Scalar/decay widths/ctauS.dat','phenomenology/Scalar/decay widths/ctauS.mx'),
('phenomenology/Scalar/decay widths/BrRatios.m','phenomenology/Scalar/decay widths/BrRatioScalar.mx')]
for old,new in pairs:
    s2023=blobs.get(old,'MISSING')
    ln=os.path.join(L,new)
    sl=gitsha(ln) if os.path.exists(ln) else 'MISSING'
    print('%-70s 2023=%s local=%s  %s'%(os.path.basename(new), s2023[:10], sl[:10], 'SAME' if s2023==sl else 'DIFFERENT'))
print()
for p in sorted(blobs):
    if 'spectra/SM particles' in p and ('_h.' in p or 'DoubleDistr_LHC_h' in p):
        print('2023 spectra entry:',p, blobs[p][:10])
