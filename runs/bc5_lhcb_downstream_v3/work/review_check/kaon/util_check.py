
import sys, os, time, math, random
RUN = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, RUN+'/work/analysis'); sys.path.insert(0, RUN+'/work/response'); sys.path.insert(0,RUN+'/work/widths')
import bc5_lib as L, bc5_production as P, run_scan as RS
import downstream_ctypes as dk
import numpy as np
print('HBAR_C', L.sw.HBAR_C)
print('M_K0', L.M_K0L, 'M_PI', L.M_PI, 'M_PI0', L.M_PI0)
print('ctau_KL = HBAR_C/GAM_KL =', L.sw.HBAR_C/P.GAM_KL, 'm')
for m in [1.0,1.2,1.5,1.8,2.0,0.1,0.3]:
    br,_ = L.channel_branching(m)
    top = sorted(br.items(), key=lambda kv:-kv[1])[:5]
    print('m=%.2f Gam(1)=%.4e  ctau(1)=%.3e m  BR: %s'%(m, L.sw.gamma_total(m), L.sw.HBAR_C/L.sw.gamma_total(m), ['%s:%.3f'%kv for kv in top]))
# benchmark transport
rng=random.Random(1)
t0=time.time(); N=200000
for i in range(N):
    dk.transport_track(0.1,0.05,1.5, 0.5,0.1,3.0, 3.1, 1)
print('transport calls/s: %.0f  (%.1f us/call)'%(N/(time.time()-t0), 1e6*(time.time()-t0)/N))
t0=time.time()
for i in range(N):
    dk.passes_birth_kinematics(0.5,0.1,3.0,3.1,5.0,0.5)
print('birth calls/s: %.0f'%(N/(time.time()-t0)))
