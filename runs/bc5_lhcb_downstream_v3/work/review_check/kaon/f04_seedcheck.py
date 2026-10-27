
import numpy as np, sys, math, random
R = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, R+'/work/analysis'); sys.path.insert(0, R+'/work/response'); sys.path.insert(0,R+'/work/widths')
import bc5_lib as L
import importlib.util
spec = importlib.util.spec_from_file_location('f04', R+'/work/review_check/kaon/f04_mc.py')
f04 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f04)
ms=1.5
# same initial rng seed -> the KK final state streams must be identical
r1=random.Random(77); r2=random.Random(77)
ok=True
for t in range(200):
    kk = f04._split_kk(ms, r1)
    ref = L.s_decay_products(ms,'KK', r2)
    if kk is None:
        same = (ref == [])
    else:
        kind, drest, kinfo = kk
        same = len(drest)==len(ref) and all(abs(a[i]-b[i])<1e-12 for a,b in zip(drest,ref) for i in range(6))
    if not same:
        ok=False; print('mismatch at draw', t); print(kk); print(ref); break
print('same random stream -> identical final states:', ok)
