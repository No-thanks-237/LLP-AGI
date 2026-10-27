
import sys, os, math
RUN = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3"
sys.path.insert(0, os.path.join(RUN, "work", "analysis"))
sys.path.insert(0, os.path.join(RUN, "work", "widths"))
import bc5_production as P, scalar_widths as sw

GF = 1.1663787e-5
MT = 173.0
VTS, VTB = 0.0404, 0.9991
MB, MS_ = 4.18, 0.0924
MBP, MK = P.M_BP, P.M_KP
V = sw.V_HIGGS

def xi_from_CKM():
    return 3*math.sqrt(2)*GF*MT*MT*abs(VTS*VTB)/(16*math.pi**2)

def lam_sqrt(mx, my, mz):
    a = (mx*mx-(my-mz)**2)*(mx*mx-(my+mz)**2)
    return math.sqrt(max(a,0.0))/(mx*mx)

def br_winkler_KS(ms):
    """BR(B->K S)/theta^2 from Winkler appendix (eq:BKphi + |<K|sL bR|B>|^2)."""
    if MBP < MK + ms: return 0.0
    g = (MB/V)*xi_from_CKM()                      # g_phi sb / theta
    fK = 0.33/(1.0-ms*ms/37.5)
    ME2 = 0.25*(MBP**2-MK**2)**2/(MB-MS_)**2 * fK*fK
    Gam = g*g*ME2*lam_sqrt(MBP, MK, ms)/(16*math.pi*MBP)
    return Gam/P.GAM_BP

print("xi from Winkler's CKM/G_F formula: %.5g   (code XI_SB = %.5g)" % (xi_from_CKM(), P.XI_SB))
print()
print("%8s %14s %14s %8s" % ("m_S", "our 'K' ch.", "Winkler app.", "ratio"))
for ms in [0.0, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 4.7, 4.78]:
    ours = P.br_B_to_S(ms, 'B+')['K']
    w = br_winkler_KS(ms)
    print("%8.3f %14.6g %14.6g %8.4f" % (ms, ours, w, ours/w if w else float('nan')))
