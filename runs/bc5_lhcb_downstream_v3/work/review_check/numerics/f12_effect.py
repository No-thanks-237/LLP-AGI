# -*- coding: utf-8 -*-
"""F12: quantify the visible-rate effect of dropping the gamma channel (assigned to 'gg')."""
import sys, os, math, random, json
import numpy as np
RUN = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, os.path.join(RUN,"work","analysis")); sys.path.insert(0, os.path.join(RUN,"work","widths")); sys.path.insert(0, os.path.join(RUN,"work","response"))
import run_scan as R, bc5_lib as L, scalar_widths as sw
HERE=os.path.dirname(os.path.abspath(__file__))
ok=L.load_hadron_lib(os.path.join(RUN,'work','magnus','hadron_lib'))
print("hadron lib:", ok)
rng=random.Random(4242)
# check the inside-window rows of m=0.2802
d=np.load(os.path.join(HERE,"tab_0.2802.npz")); z0=d['z']; A=d['A']
ins=(z0>1.0)&(z0<2.5)
print("\nF10 check: m=0.2802 inside-window rows: A.sum(z) = %s" % np.round(A[ins].sum(1),3))

# scalar kinematics: eta_S = 3.0, pz = 50 GeV, pt = 5 GeV
pz=50.0; pt=5.0; px=pt; py=0.0
print("\nF12: per-decay visible acceptance at eta_S=3, pz=50 GeV, z_dec=1.5 m, x=y=0")
print("  %-10s %-10s %-10s %-10s %-10s %-10s" % ("m_S","BR_gamma","P_gg","A_mix","dA/A","n_ev"))
scan=json.load(open(os.path.join(RUN,"outputs","data","signal_scan.json")))
masses=np.array(scan['masses'])
def br_table(m):
    w=sw.widths(m); tot=sum(w.values())
    br={c: w.get(c,0.0)/tot for c in L.CHANNELS}
    return br, w.get('gamma',0.0)/tot, tot
for m in (0.188692, 0.2118, 0.25821, 0.2802, 0.3541, 0.5, 0.7, 0.9884, 1.5, 5.0):
    br, brg, tot = br_table(m)
    pmod=math.sqrt(px*px+py*py+pz*pz); E=math.sqrt(m*m+pmod*pmod)
    N=600
    Pc={}; nch={}
    for c in L.CHANNELS:
        if br[c] <= 0: 
            Pc[c]=0.0; continue
        good=0; n=0
        for _ in range(N):
            rest=L.s_decay_products(m, c, rng)
            if not rest: continue
            chg=[a[5] for a in rest]
            lab5=L.boost([(a[0],a[1],a[2],a[3],a[4]) for a in rest], E, (px,py,pz))
            lab=[(b[0],b[1],b[2],b[3],b[4],cc) for b,cc in zip(lab5,chg)]
            n+=1; good+=int(L._passes_pair(lab, 0.0, 0.0, 1.5))
        Pc[c]=good/N if N else 0.0
    Amix=sum(br[c]*Pc[c] for c in L.CHANNELS)       # what the code effectively does (remainder -> gg)
    Atrue=Amix + 0.0                                # correct: gamma contributes 0
    # dA/A = brg*P_gg / Amix  (code = Amix + brg*P_gg because remainder folded into gg)
    Pgg=Pc.get('gg',0.0)
    dA = brg*Pgg/max(Amix,1e-300)
    print("  %-10.5f %-10.5f %-10.4f %-10.5f %-10.5f" % (m, brg, Pgg, Amix, dA))
    print("           per-channel P:", " ".join("%s=%.3f"%(c,Pc[c]) for c in L.CHANNELS if br[c]>0))
