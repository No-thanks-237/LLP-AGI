
import sys, os, math
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\analysis")
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\widths")
import scalar_widths as sw
# old gamma-gamma (with Q_b^2 = 4/9) reconstructed by patching the loop charge
import types
orig = sw.gamma_gammagamma
def old_gg(ms):
    A = complex(0.0, 0.0)
    for m, q2, nc in ((sw.M_E,1.0,1),(sw.M_MU,1.0,1),(sw.M_TAU,1.0,1),
                      (sw.M_U_MSBAR,4.0/9.0,3),(sw.M_D_MSBAR,1.0/9.0,3),
                      (sw.MSB_S_2GEV,1.0/9.0,3),(sw.MSB_C_2GEV,4.0/9.0,3),
                      (sw.MSB_B_MB,4.0/9.0,3),(sw.M_T_POLE,4.0/9.0,3)):
        A += nc*q2*sw.F_half(4.0*m*m/ms**2)
    A += sw.F_one(4.0*sw.M_W**2/ms**2)
    return sw.G_F*sw.ALPHA_EM_0**2*ms**3*abs(A)**2/(128.0*math.sqrt(2.0)*math.pi**3)
print("%9s %12s %12s %10s %10s" % ("m_S","gam_gg_new","gam_gg_old","new/old","BR_gam"))
for m in (0.05,0.1,0.188692,0.2118,0.28,0.5,1.0,1.5,2.0,3.0,10.0,62.0):
    wn = sw.widths(m); gn = wn['gamma']; go = old_gg(m)
    tot = sum(wn.values())
    print("%9.4f %12.4e %12.4e %10.3f %10.4f" % (m, gn, go, gn/go, gn/tot))
print()
print("%9s %12s %12s %8s" % ("m_S","Gamma_tot_new","Gamma_tot_old","ratio"))
for m in (0.05,0.1,0.188692,0.2118,0.28,0.5,1.0,1.5,2.0):
    wn = sw.widths(m); totn = sum(wn.values())
    toto = totn - wn['gamma'] + old_gg(m)
    print("%9.4f %12.5e %12.5e %8.4f" % (m, totn, toto, totn/toto))
