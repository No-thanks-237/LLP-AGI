
# -*- coding: utf-8 -*-
"""Validate the F04 reimplementation of acceptance_profile (current treatment) and inspect r(z)."""
import numpy as np, sys, math, random
R = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, R+'/work/analysis'); sys.path.insert(0, R+'/work/response'); sys.path.insert(0,R+'/work/widths')
import bc5_lib as L
import importlib.util
spec = importlib.util.spec_from_file_location('f04', R+'/work/review_check/kaon/f04_mc.py')
f04 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f04)

def ref_KK(ms, s, rng, ndecay):
    tag,e,px,py,pz,x,y,z = s
    pmod=math.sqrt(px*px+py*py+pz*pz)
    x0,y0,z0=x*1e-3,y*1e-3,z*1e-3
    dx,dy = px/pz, py/pz
    A=np.zeros(13); n=0
    for _ in range(ndecay):
        rest = L.s_decay_products(ms,'KK',rng)
        if not rest: continue
        n+=1
        chg=[a[5] for a in rest]
        lab5=L.boost([(a[0],a[1],a[2],a[3],a[4]) for a in rest], e,(px,py,pz))
        lab=[(b[0],b[1],b[2],b[3],b[4],c) for b,c in zip(lab5,chg)]
        for k,zd in enumerate(f04.ZG):
            A[k]+=L._passes_pair(lab, x0+(zd-z0)*dx, y0+(zd-z0)*dy, zd)
    return A/n if n else A

rng = random.Random(4321)
f04._load()
S,nd = f04.sample_S_B(1.5, 6, rng)
print('ndraws', nd)
for s in S[:4]:
    e=s[1]; p=math.sqrt(s[2]**2+s[3]**2+s[4]**2)
    print('\nS eta=%.2f p=%.1f'%(0.5*math.log((p+s[4])/(p-s[4])), p))
    a1 = f04.s_arrays(1.5, s, random.Random(11), ndecay=300, ndecay_all=0)
    a2 = ref_KK(1.5, s, random.Random(22), 300)
    r1 = a1[0]; r2 = a2
    print('  my Ac   ', np.round(r1,3))
    print('  ref KK  ', np.round(r2,3))
    print('  ratio(new/cur) per z', np.round(a1[1]/np.maximum(r1,1e-12),2))
