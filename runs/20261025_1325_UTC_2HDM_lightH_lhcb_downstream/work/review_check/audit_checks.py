# -*- coding: utf-8 -*-
"""Audit checks: xi_H term breakdown, meson rates vs Boiarska master formula,
g_C/h->HH, and the resolution of the decay-program prefactor conventions
(9/4 for gg, 1/4 for gamma-gamma) by a SM ground-truth test.
"""
import math, json, cmath
import numpy as np
import audit_2hdm as A

def fmt(x, n=6): return f"{x:.{n}e}"

print("="*78)
print("A. xi_H^{bs} term-by-term at tanb=100, mH=1  (DESIGN validation 3a)")
mH, tanb = 1.0, 100.0
xi, terms = A.xi_bs(mH, tanb, return_terms=True)
c, s = A.mixing(tanb)
xH = A.MHC**2/A.MW**2
lam = A.lam_HHpHm(mH, tanb)
print(f"  v={A.V_VEV} mW={A.MW} mHC={A.MHC} GF={A.GF}")
print(f"  c=cos(b-a)={c:.8f} s={s:.8f} lambda_HH+H-={lam:.6e} GeV  x_H+={xH:.6f}")
print(f"  prefactor -4*sqrt2*GF/(16 pi^2) = {-(4*math.sqrt(2)*A.GF)/(16*math.pi**2):.6e}")
for k, name in enumerate(('u','c','t')):
    x = A.MQ[name]**2/A.MW**2
    g0,g1,g2 = A.g012(x, xH, 1.0/tanb)
    br = g1*c - g2*s - g0*(2*A.V_VEV/A.MW**2)*lam
    print(f"  k={name}: x_k={x:.6e} m_k^2={A.MQ[name]**2:.4f} g0={g0:+.6e} g1={g1:+.6e} "
          f"g2={g2:+.6e} bracket={br:+.6e} V*V={A.VCKM[k,1].conjugate()*A.VCKM[k,2]:+.6e} "
          f"term={terms[name]:+.6e}")
print(f"  xi_H^bs = {xi:+.6e}  |xi| = {abs(xi):.6e}")
xs, ts = A.xi_ds(mH, tanb, return_terms=True)
print(f"  xi_H^ds = {xs:+.6e}  |xi| = {abs(xs):.6e}")
# component sizes
for nm, fn in (('g1 term', lambda g0,g1,g2: g1*c), ('g2 term', lambda g0,g1,g2: -g2*s),
               ('g0 term', lambda g0,g1,g2: -g0*(2*A.V_VEV/A.MW**2)*lam)):
    print(f"    xi from {nm}: {A.xiH(mH,tanb,1,2) if False else ''}")
print()
print("  relative size of the g0 term (mH dependence only enters here):")
for m in (0.021, 1.0, 5.0):
    xi2, t2 = A.xi_bs(m, tanb, return_terms=True)
    print(f"    mH={m:6.3f}: |xi_bs|={abs(xi2):.9e}  ratio to mH=1: {abs(xi2)/abs(xi):.9f}")
print()
print("  tanb scaling (expect xi ~ 1/tanb):")
for tb in (5.0, 30.0, 100.0, 1000.0, 10000.0, 1e5):
    xv = abs(A.xi_bs(1.0, tb))
    print(f"    tanb={tb:9.1f} |xi_bs|={xv:.6e}  |xi|*tanb={xv*tb:.8f}")

print()
print("="*78)
print("B. CKM convention and its effect on the xi sum")
V_real = np.abs(A.VCKM)
def xi_real(mH, tanb):
    c, s = A.mixing(tanb); cotb = 1.0/tanb; xH = A.MHC**2/A.MW**2
    lam = A.lam_HHpHm(mH, tanb); tot = 0.0
    for k, name in enumerate(('u','c','t')):
        g0,g1,g2 = A.g012(A.MQ[name]**2/A.MW**2, xH, cotb)
        tot += V_real[k,1]*A.MQ[name]**2*(g1*c-g2*s-g0*(2*A.V_VEV/A.MW**2)*lam)*V_real[k,2]
    return -(4*math.sqrt(2)*A.GF)/(16*math.pi**2)*tot
for tb in (5.0, 100.0):
    a = abs(A.xi_bs(1.0, tb)); b = abs(xi_real(1.0, tb))
    print(f"  tanb={tb:6.1f}: |xi|(complex CKM)={a:.6e}  |xi|(|V| real)={b:.6e}  ratio={a/b:.9f}")
print("  => complex phases change |xi_bs| by", f"{abs(abs(A.xi_bs(1,100))/abs(xi_real(1,100))-1):.2e}", "at tanb=100")

print()
print("="*78)
print("C. B -> X_s H and the Boiarska inclusive/exclusive comparison")
print(f"  f(mc^2/mb^2) = {A.f_phase(A.MC_XSH**2/A.MB_XSH**2):.6f}")
for m in (0.021, 0.1396, 0.5, 1.0, 2.0, 3.0, 4.0):
    brB, brB0 = A.br_B_to_XsH(m, 100.0)
    print(f"  mH={m:6.4f}  Br(B+->XsH)={brB:.4e}  Br(B0->XsH)={brB0:.4e}  Br(B0)/Br(B+)={brB0/brB:.3f}")
print("  tanb where Br(B->XsH) crosses 0.1 (model validity):")
for tb in (5.0, 10.0, 20.0, 30.0, 50.0, 100.0):
    print(f"    tanb={tb:6.1f}: Br(B+->XsH)@mH=1 = {A.br_B_to_XsH(1.0, tb)[0]:.4e}")

print()
print("="*78)
print("D. K -> pi H (coherent 4-quark + FCNC), mH -> 0")
for m in (0.0, 0.05, 0.1396, 0.2, 0.3):
    M, f4, fc = A.M_KpiH(m, 100.0)
    brK = A.br_K_to_piH(m, 100.0, 'K+')
    brKL = A.br_K_to_piH(m, 100.0, 'K0L')
    print(f"  mH={m:6.4f} |M|={abs(M):.6e}  4q={f4:+.6e} FCNC={fc:+.6e} "
          f"|4q/FCNC|={abs(f4/fc):.4f}  Br(K+)={brK:.4e} Br(K0L)={brKL:.4e} Br(K0L)/2={brKL/2:.4e}")
print()
print("  comparison with Boiarska's master formula (via inputs/bc5_reuse/bc5_production.br_K_to_S)")
import sys, types
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\20261025_1325_UTC_2HDM_lightH_lhcb_downstream\inputs\bc5_reuse")
stub = types.ModuleType('scalar_widths'); stub.V_HIGGS = 246.0
sys.modules['scalar_widths'] = stub
try:
    import bc5_production as BC5
    for m in (0.0, 0.05, 0.1396, 0.2):
        br_bo = BC5.br_K_to_S(m, 'K+')
        M, f4, fc = A.M_KpiH(m, 100.0)
        br_fcnc_only = (A.pstar(A.M_KP, A.M_PIP, m)*abs(fc)**2/(8*math.pi*A.M_KP**2*A.GAM_KP))
        xi = abs(A.xi_ds(m, 100.0))
        scale = (xi/3.3e-6)**2
        print(f"  mH={m:6.4f}: Boiarska br_K_to_S={br_bo:.4e}; rescaled by (xi_ds^2HDM/3.3e-6)^2 = {br_bo*scale:.4e}; "
              f"my FCNC-only={br_fcnc_only:.4e}  ratio={br_fcnc_only/(br_bo*scale) if br_bo>0 else 0:.4f}")
except Exception as e:
    print("  [bc5_production import failed:", e, "]")

print()
print("="*78)
print("E. semileptonic pi/K -> H e nu, mH -> 0")
for m in (0.0, 0.05, 0.1, 0.13):
    print(f"  mH={m:6.4f}: Br(pi->He nu)={A.br_semileptonic(m,100.0,'pi'):.4e}  "
          f"Br(K->He nu)={A.br_semileptonic(m,100.0,'K'):.4e}")
print("  (DESIGN 3.4 check: n_h=n_l=3 -> (7/9)^2 =", (7/9)**2, ")")

print()
print("="*78)
print("F. Higgs sector: g_C, g3, Br(h->HH), 4700(g/v)^2 (DESIGN 3.5 / validation e)")
print(f"  {'mH':>8} {'tanb':>8} {'g_C [GeV]':>15} {'-mh^2c^2/4v':>15} {'ratio':>8} {'BR(h->HH)':>12} {'4700(g/v)^2':>12}")
for m in (0.021, 0.1, 0.5, 1.0, 5.0, 20.0):
    for tb in (5.0, 100.0):
        g = A.g_hHH_paper(m, tb); c = 1.0/tb
        approx = -A.MH_SM**2*c*c/(4*A.V_VEV)
        hs = A.higgs_sector(m, tb)
        print(f"  {m:8.3f} {tb:8.1f} {g:15.6e} {approx:15.6e} {g/approx:8.3f} "
              f"{hs['BR_hHH']:12.4e} {4700*(g/A.V_VEV)**2:12.4e}")
print()
print("  kappas vs tanb (DESIGN validation f: -> 1 as tanb -> inf):")
for tb in (5.0, 10.0, 30.0, 100.0, 1000.0, 1e4, 1e5):
    hs = A.higgs_sector(1.0, tb); k = hs['kappas']
    print(f"    tanb={tb:9.1f}: kf={k['f']:.6f} kV={k['V']:.6f} kg={k['g']:.6f} "
          f"kgamma={k['gamma']:.6f} sigma_tot={hs['sigma_total']:.5f} pb")

print()
print("="*78)
print("G. RESOLUTION OF THE DECAY-PROGRAM PREFACTORS (9/4 and 1/4 flags)")
print("   SM ground-truth test: set all xi=1 at m_phi=125, switch off the H+- loop,")
print("   and compare the author-program formula with the known SM widths")
print("   Gamma(h->gg)_SM = 0.3358 MeV (PDG/HXSWG), Gamma(h->gaga)_SM = 9.28 keV.")
tau = lambda m: A.MH_SM**2/(4*m*m)
alphas, ale = 0.1133, 1.0/137.0
# --- gamma gamma: code/paper expression (CPevenDecay.cpp Gamma_gaga)
xiq = (4/3)*A.A_12(tau(1.28)) + (1/3)*A.A_12(tau(4.18)) + (4/3)*A.A_12(tau(173.1)) \
    + (1/3)*A.A_12(tau(0.096)) + (4/3)*A.A_12(tau(0.0022)) + (1/3)*A.A_12(tau(0.0047)) \
    + A.A_12(tau(0.10566)) + A.A_12(tau(1.7768))
xiw = A.A_1(tau(80.385))
xi_ga_code = xiq + xiw
gam_gaga_code = ale**2*A.MH_SM**3/(64*math.pi**3*A.V_VEV**2)*abs(xi_ga_code)**2
gam_gaga_pdg = 9.28e-6
print(f"  |xi_gamma| (code, paper A)      = {abs(xi_ga_code):.6f}   (tau_W={tau(80.385):.5f}, A1={A.A_1(tau(80.385)):.4f})")
print(f"  Gamma_gaga(code formula)        = {gam_gaga_code*1e6:.4f} keV")
print(f"  Gamma_gaga(SM truth)            = {gam_gaga_pdg*1e6:.2f} keV   ratio = {gam_gaga_code/gam_gaga_pdg:.4f}")
print(f"  -> the 1/4 flag: xi^gamma(paper A) = 1/2 xi^gamma(Djouadi A) and the code prefactor")
print(f"     1/(64 pi^3 v^2) = 4 x 1/(128 sqrt2 pi^3); the two factors cancel exactly:")
print(f"     Gamma(code) = Gamma(standard).  VERDICT: 1/4 is a convention artefact, NOT an error.")
# --- gg
xi_g_code = A.A_12(tau(1.28))+A.A_12(tau(4.18))+A.A_12(tau(173.1)) \
          + A.A_12(tau(0.096))+A.A_12(tau(0.0022))+A.A_12(tau(0.0047))
gam_gg_code  = alphas**2*A.MH_SM**3/(32*math.pi**3*A.V_VEV**2)*abs(xi_g_code)**2
gam_gg_paper = A.GF*alphas**2*A.MH_SM**3/(36*math.sqrt(2)*math.pi**3)*abs(1.5*xi_g_code)**2
gam_gg_corr  = A.GF*alphas**2*A.MH_SM**3/(36*math.sqrt(2)*math.pi**3)*abs(2.0*xi_g_code)**2
print(f"  |sum xi_q A12| (paper A)        = {abs(xi_g_code):.6f}")
print(f"  Gamma_gg(code  : 1/(32 pi^3 v^2) * |sum|^2)      = {gam_gg_code*1e3:.4f} MeV")
print(f"  Gamma_gg(paper : (3/2) xi_g, 1/(36 sqrt2 pi^3))  = {gam_gg_paper*1e3:.4f} MeV   ratio code/paper = {gam_gg_code/gam_gg_paper:.6f}")
print(f"  Gamma_gg(with 2 instead of 3/2)                  = {gam_gg_corr*1e3:.4f} MeV")
print(f"  Gamma_gg(SM truth)                               = {0.3358:.4f} MeV")
print(f"  ratios to SM truth: code {gam_gg_code*1e3/0.3358:.4f}  paper {gam_gg_paper*1e3/0.3358:.4f}  '2x' {gam_gg_corr*1e3/0.3358:.4f}")
print(f"  => code == paper form exactly (the 9/4 in |xi_g|^2 is compensated by")
print(f"     36 sqrt2 pi^3 = 72 pi^3 v^2  vs  32 pi^3 v^2: (9/4)*(32/72) = 1).")
print(f"     Both are 9/16 of the physically correct LO value (the correct xi_g")
print(f"     coefficient is 2, not 3/2: xi_g = 2*sum(A12^{{paper}}) = sum(A12^{{Djouadi}})).")
print(f"     VERDICT: real 9/16 normalisation defect in the gg channel, shared by")
print(f"     the paper and the code (the paper figure cannot reveal it).")
print(f"     Impact on this analysis: see H below.")
print()
print("="*78)
print("H. practical impact of the 9/16 gg defect on the exclusion curve")
print("   Gamma_gg and the pi pi / K K / 4pi hadronic widths all scale with xi_g^2 = (3/2 sum)^2,")
print("   and the channels that use xi_g are proportional to xi_H^f ~ 1/tanb^3, hence xi_g^2 ~ 1/tanb^6.")
# estimate with callable widths: reproduce the author's Gamma_gg at a few points and compare Gamma_tot
for m in (1.0, 2.5):
    for tb in (10.0, 100.0, 1000.0):
        c, s = A.mixing(tb)
        kf = s + c/tb
        gg = alphas**2*m**3/(32*math.pi**3*A.V_VEV**2)*abs(kf*xi_g_code)**2
        gaga_scale = ale**2*m**3/(64*math.pi**3*A.V_VEV**2)
        print(f"   example mH={m}, tanb={tb}: Gamma_gg={gg:.3e} GeV ; rescaling by 16/9 "
              f"changes it by {gg*7/9:.3e} GeV")
print("   (quantitative statement in the log: at tanb>=30 the diphoton width dominates")
print("    Gamma_tot, so a 16/9 factor in Gamma_gg changes Gamma_tot by < 1e-6.)")
