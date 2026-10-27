
import numpy as np, sys, math, random
R = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, R+'/work/analysis'); sys.path.insert(0, R+'/work/response'); sys.path.insert(0,R+'/work/widths')
import bc5_lib as L, bc5_production as P, run_scan as RS, downstream_ctypes as dk
import importlib.util
spec = importlib.util.spec_from_file_location('f04', R+'/work/review_check/kaon/f04_mc.py')
f04 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f04)
rng = random.Random(999)
f04._load()
S, nd = f04.sample_S_B(1.5, 5, rng)
for s in S:
    tag,e,px,py,pz,x,y,z = s
    m = math.sqrt(max(e*e-px*px-py*py-pz*pz,0))
    pmod = math.sqrt(px*px+py*py+pz*pz)
    eta = 0.5*math.log((pmod+pz)/(pmod-pz))
    print('S m=%.4f e=%.2f p=%.2f eta=%.2f z0=%.4f'%(m,e,pmod,eta,z*1e-3))
    res = f04.s_arrays(1.5, s, rng, ndecay=4, ndecay_all=2)
    if res is None: print('   rejected by eta'); continue
    Ac, An, Aa, eta2, pz2, z02, nkk, nall, fl = res
    print('   Ac', np.round(Ac,2), ' An', np.round(An,2), ' Aa', np.round(Aa,2))
    # manual look at one KK decay
    kk = f04._split_kk(1.5, rng)
    if kk is not None:
        kind, drest, kinfo = kk
        lab5 = L.boost([(a[0],a[1],a[2],a[3],a[4]) for a in drest], e, (px,py,pz))
        for b in lab5:
            print('     kind',kind,'dau e=%.2f pT=%.2f pz=%.1f'%(b[1], math.hypot(b[2],b[3]), b[4]))
