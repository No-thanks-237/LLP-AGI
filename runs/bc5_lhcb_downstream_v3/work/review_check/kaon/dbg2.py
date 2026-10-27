
import numpy as np, sys, math, random
R = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, R+'/work/analysis'); sys.path.insert(0, R+'/work/response'); sys.path.insert(0,R+'/work/widths')
import bc5_lib as L, downstream_ctypes as dk
print('DLL', dk.DLL_PATH)
# 1) direct kernel test
for (x,y,z,px,py,pz,e,ch) in [
    (0.1,0.05,1.5, 0.5,0.1,3.0, math.hypot(3.0,math.hypot(0.5,0.1)), 1),
    (0.1,0.05,1.5, 0.5,0.1,50.0, math.hypot(50.0,math.hypot(0.5,0.1)), 1),
    (0.0,0.0,1.5, 0.0,0.0,50.0, 50.0, 1),
]:
    print(dk.transport_track(x,y,z,px,py,pz,e,ch))
# 2) single S -> KK check through real code path
rng = random.Random(7)
e,px,py,pz = 100.0, 6.0, 0.0, 95.0
ms=1.5
pmod=math.sqrt(px*px+py*py+pz*pz); print('S eta', 0.5*math.log((pmod+pz)/(pmod-pz)))
rest = L.s_decay_products(ms,'KK',rng)
print('rest', len(rest), rest)
lab5 = L.boost([(a[0],a[1],a[2],a[3],a[4]) for a in rest], e, (px,py,pz))
chg=[a[5] for a in rest]
lab=[(b[0],b[1],b[2],b[3],b[4],c) for b,c in zip(lab5,chg)]
for b in lab: print('  dau m=%.4f e=%.2f p=(%.2f,%.2f,%.2f) ch=%d pT=%.2f'%(b[0],b[1],b[2],b[3],b[4],b[5],math.hypot(b[2],b[3])))
print('birth', [dk.passes_birth_kinematics(b[2],b[3],b[4],b[1],dk.BIRTH_CUTS['minE'],dk.BIRTH_CUTS['minPT']) for b in lab])
for zd in (1.0,1.5,2.0):
    print(' zd',zd, 'passes_pair=', L._passes_pair(lab, 0.0+(zd-0.0)*px/pz, 0.0+(zd-0.0)*py/pz, zd))
