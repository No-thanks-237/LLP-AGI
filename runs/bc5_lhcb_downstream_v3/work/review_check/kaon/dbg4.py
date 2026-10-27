
import numpy as np, sys, math, random
R = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, R+'/work/analysis'); sys.path.insert(0, R+'/work/response'); sys.path.insert(0,R+'/work/widths')
import bc5_lib as L, downstream_ctypes as dk
import importlib.util
spec = importlib.util.spec_from_file_location('f04', R+'/work/review_check/kaon/f04_mc.py')
f04 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f04)
rng = random.Random(999)
f04._load()
S, nd = f04.sample_S_B(1.5, 5, rng)
s = S[1]
tag,e,px,py,pz,x,y,z = s
print('S', tag, e, px, py, pz, x, y, z)
kk = f04._split_kk(1.5, rng)
kind, drest, kinfo = kk
print('kind', kind, 'ndau', len(drest))
lab5 = L.boost([(a[0],a[1],a[2],a[3],a[4]) for a in drest], e, (px,py,pz))
dau = [(b[1],b[2],b[3],b[4],c[5],c[5]) for b,c in zip(lab5,drest)]
x0,y0,z0 = x*1e-3,y*1e-3,z*1e-3
for zd in (1.0,1.5,2.0,2.5):
    bx = x0+(zd-z0)*px/pz; by=y0+(zd-z0)*py/pz
    print('zd',zd,'birth',bx,by)
    for (dpx,dpy,dpz,de,ch,q) in dau:
        bc = dk.passes_birth_kinematics(dpx,dpy,dpz,de,5.0,0.5)
        st = dk.transport_track(bx,by,zd,dpx,dpy,dpz,de,ch)
        print('   dau pT=%.2f e=%.2f ch=%d birth=%s -> %s'%(math.hypot(dpx,dpy),de,ch,bc,st[1]))
