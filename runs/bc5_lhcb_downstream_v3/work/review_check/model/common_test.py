import sys, os, time, math, random
import numpy as np
RUN=r'E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3'
sys.path.insert(0, os.path.join(RUN,'work','analysis'))
sys.path.insert(0, os.path.join(RUN,'work','widths'))
import bc5_lib as L, bc5_production as P, run_scan as R
print('import ok, kernel', L.dk.__file__ if hasattr(L.dk,'__file__') else L.dk)
print('hadronic table loaded', L.sw.HADRONIC_TABLE is not None)
t0=time.time(); ok=L.load_hadron_lib(os.path.join(RUN,'work','magnus','hadron_lib')); print('hadron lib',ok,'%.1fs'%(time.time()-t0))
# quick acceptance test
rng=random.Random(1)
p4=(60.1,2.0,5.7,59.7)
a=L.acceptance_profile(1.99,p4,(0,-2,1500),rng,ndecay=100)
print('acceptance profile',np.round(a,3))
print('ZGRID',L.ZGRID)
