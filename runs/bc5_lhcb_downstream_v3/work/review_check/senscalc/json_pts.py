
import json, numpy as np, os
p=r'E:/AAAwork/LLP/code/SensCalc-main/Sensitivity domains/Scalar/LHCb-downstream-D/Sensitivity_Scalar_at_LHCb-downstream-D_BrhToSS=0.01_Nev=2.3_Npot=1.75e15.json'
d=json.load(open(p))
pts=np.array(d[0]['Sensitivity domains'][0], dtype=float)
print('npts',len(pts))
# sort by mass
o=np.argsort(pts[:,0])
pts=pts[o]
print('m range %.4f - %.4f'%(pts[0,0],pts[-1,0]))
for mq in [2.0,3.0,4.5,5.0,5.08,10.,13.02,20.,30.,33.37,40.,50.,53.42,60.,62.]:
    i=np.argmin(np.abs(pts[:,0]-mq))
    print('  m=%.3f -> th2min=%.6e'%(pts[i,0],pts[i,1]))
# print a few points where mass > 4.9
sel=pts[pts[:,0]>4.9]
print('n high-mass pts', len(sel))
print(sel[:10])
