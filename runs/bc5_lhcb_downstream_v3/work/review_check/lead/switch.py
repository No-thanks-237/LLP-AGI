
import sys, os, math
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\analysis")
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\widths")
import scalar_widths as sw
GF, MT, VTS, VTB = 1.1663787e-5, 173.0, 0.0404, 0.9991
XI = 3*math.sqrt(2)*GF*MT*MT*abs(VTS*VTB)/(16*math.pi**2)
MBQ, MSQ, V = 4.18, 0.0924, sw.V_HIGGS
MBP, MB0, MBS = 5.27934, 5.27965, 5.36688
MK, MPHI = 0.493677, 1.019461
GAM_BP, GAM_B0, GAM_BS = 6.582119569e-25/1.638e-12, 6.582119569e-25/1.519e-12, 6.582119569e-25/1.520e-12
print("XI_sb(Winkler) = %.6g  (code 3.6e-4)" % XI)

def incl(ms, mX, GX):
    return (MBQ/V)**2*XI**2*(mX*mX-ms*ms)**2/(32*math.pi*mX**3*GX)

def AK(ms):  return 1.36/(1-ms*ms/27.9) - 0.99/(1-ms*ms/36.8)
def fK(ms):  return 0.33/(1-ms*ms/37.5)

def br_ps(ms, mX, GX, mY):
    if mX < mY+ms: return 0.0
    ME2 = 0.25*(mX*mX-mY*mY)**2/(MBQ-MSQ)**2*fK(ms)**2
    lam = math.sqrt(max((mX*mX-(mY-ms)**2)*(mX*mX-(mY+ms)**2),1e-300))/(mX*mX)
    return (MBQ/V)**2*XI**2*ME2*lam/(16*math.pi*mX*GX)

def br_vec(ms, mX, GX, mY):
    if mX < mY+ms: return 0.0
    lam_dimless = ((mX*mX-(mY-ms)**2)*(mX*mX-(mY+ms)**2))/mX**4
    ME2 = 0.25*mX**4*lam_dimless/(MBQ+MSQ)**2*AK(ms)**2
    lam = math.sqrt(max((mX*mX-(mY-ms)**2)*(mX*mX-(mY+ms)**2),1e-300))/(mX*mX)
    return (MBQ/V)**2*XI**2*ME2*lam/(16*math.pi*mX*GX)

print("\n B+ :  m    inclusive    B->K S     ratio")
for m in [0,1,2,3,4,4.4,4.5,4.6,4.7,4.72,4.75,4.78]:
    a, b = incl(m, MBP, GAM_BP), br_ps(m, MBP, GAM_BP, MK)
    print("  %5.2f  %10.4g %10.4g  %8.3f" % (m, a, b, a/b if b else float('nan')))
lo, hi = 3.0, MBP-MK
for _ in range(60):
    mid = 0.5*(lo+hi)
    if incl(mid, MBP, GAM_BP) > br_ps(mid, MBP, GAM_BP, MK): lo = mid
    else: hi = mid
print(" B+ crossing point: %.4f GeV  (Winkler quotes 4.7)" % mid)
print("\n Bs :  m    inclusive    Bs->phi S  ratio")
for m in [0,1,2,3,3.5,3.8,4.0,4.2,4.3,4.34,4.347,4.4]:
    a, b = incl(m, MBS, GAM_BS), br_vec(m, MBS, GAM_BS, MPHI)
    print("  %6.3f %10.4g %10.4g  %8.3f" % (m, a, b, a/b if b else float('nan')))
lo, hi = 2.0, MBS-MPHI
for _ in range(60):
    mid = 0.5*(lo+hi)
    if incl(mid, MBS, GAM_BS) > br_vec(mid, MBS, GAM_BS, MPHI): lo = mid
    else: hi = mid
print(" Bs crossing point: %.4f GeV" % mid)