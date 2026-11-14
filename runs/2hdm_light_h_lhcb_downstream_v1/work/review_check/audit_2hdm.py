# -*- coding: utf-8 -*-
"""Independent (audit) implementation of the 2HDM light-H production and Higgs-sector
physics for arXiv:2212.06186v1 (light CP-even H), benchmark as fixed in DESIGN.md sec.1/3.5.

AUTHOR: model/audit (T5).  Written from the paper text + DESIGN.md ONLY, without reading
work/production/twohdm_production.py.  The purpose is an independent numerical cross-check.

References to the tex: inputs/references/FASER_light.tex
  Eq.(phi_lag)     effective Lagrangian normalisation
  Eq.(Lfceff_phi)  (phi/v) xi^{ij} m_fj fbar_i P_R f_j
  Eq.(hgaga) / (hgg) / App.A form factors A_0, A_1/2, A_1
  Eq.(Matrix_4FH)  K -> pi phi amplitude (FCNC + 4-quark coherent)
  Eq.(xi_H_2HDM)   xi_H^{ij} with g0,g1,g2 (App.D)
  App.C eq:ghHH    g_hHH, and eq:coup_HCC for lambda_{HH+H-}
  eq:brhphiphi     Br(h->HH) ~ 4700 (g/v)^2
"""
import math, cmath, json
import numpy as np

# ---------------------------------------------------------------- constants
GF      = 1.1663787e-5      # GeV^-2
V_VEV   = 246.0             # GeV  (DESIGN pins v=246)
MW      = 80.385            # GeV  (DESIGN 3.1)
MZ      = 91.1876
MH_SM   = 125.0             # GeV
MHC     = 600.0             # m_{H^pm}
MA      = 600.0
MH_MIN, MH_MAX = 0.021, 10.0

# quark masses for xi_H (DESIGN 3.1: MSbar/PDG choices recorded)
MQ = dict(u=2.2e-3, c=1.28, t=173.1, d=4.7e-3, s=9.6e-2, b=4.18)
# B -> X_s H inputs (DESIGN 3.2)
MB_XSH, MC_XSH, BR_BXCE = 4.18, 1.28, 0.104
# K -> pi H inputs (DESIGN 3.3)
M_KP, M_K0L = 0.493677, 0.497611
M_PIP, M_PI0 = 0.13957039, 0.1349768
GAM_KP  = 5.314e-17         # GeV
GAM_K0L = 1.287e-17         # GeV
LAMBDA_K = 3.1e-7
A_K      = 0.0
F0_KPI   = 0.96
MS_QUARK, MD_QUARK = 0.096, 4.7e-3
# semileptonic
GAM_PI = 2.5284e-17         # pi+ width [GeV] (PDG tau=2.6033e-8 s)
M_MU    = 0.1056583755
BR_PI_MUNU, BR_K_MUNU = 0.999877, 0.6356
# SM reference values for the Higgs sector (DESIGN 3.5)
GAM_H_SM = 4.07e-3
BR_SM = dict(bb=0.5824, WW=0.2137, gg=0.0818, tautau=0.0627, cc=0.0288,
             ZZ=0.0262, gaga=0.00227, Zgamma=0.00153, ss=0.0003, mumu=0.00022)
SIGMA_SM_PB = dict(ggF=48.58, VBF=3.782, WH=1.369, ZH=0.883, bbH=0.2384,
                   ttH=0.5065, tHj=0.0740, tHW=0.0152, tHb=0.0046)

# ---------------------------------------------------------------- CKM (PDG 2022 magnitudes + standard phases)
def ckm_matrix(delta=1.144):
    """Standard PDG parametrisation, Unitarity-consistent. Magnitudes match DESIGN 3.1."""
    s12, s13, s23 = 0.2243, 0.00382, 0.0408
    c12, c13, c23 = math.sqrt(1-s12**2), math.sqrt(1-s13**2), math.sqrt(1-s23**2)
    e = cmath.exp(1j*delta)
    V = np.zeros((3,3), complex)
    V[0,0]=c12*c13; V[0,1]=s12*c13; V[0,2]=s13/e
    V[1,0]=-s12*c23-c12*s23*s13*e; V[1,1]=c12*c23-s12*s23*s13*e; V[1,2]=s23*c13
    V[2,0]=s12*s23-c12*c23*s13*e; V[2,1]=-c12*s23-s12*c23*s13*e; V[2,2]=c23*c13
    return V
VCKM = ckm_matrix()

# ---------------------------------------------------------------- 2HDM mixing / trilinears
def mixing(tanb):
    c = 1.0/tanb
    s = math.sqrt(max(0.0, 1.0-c*c))
    return c, s

def lam_HHpHm(mH, tanb):     # App.C eq:coup_HCC for lambda v^2 = 0
    c, s = mixing(tanb)
    return -(2.0*MHC**2 - mH**2)*c/V_VEV

def g_hHH_paper(mH, tanb, lambdav2=0.0):
    """App.C eq:ghHH for general lambda v^2 (default 0)."""
    c, s = mixing(tanb)
    tanb_ = 1.0/c
    t2b = (1.0-tanb_**2)/(2.0*tanb_)       # cot(2beta):  (1-tan^2)/(2 tan) = (c^2-1)/(2c)
    brace = (mH**2 - 3.0*lambdav2 - MH_SM**2)*(2.0*t2b*s*c - c**2 + s**2) + (lambdav2 - mH**2)
    return s/(2.0*V_VEV)*brace

def g_hHH3(mH, tanb):            # g3 = -2 g_C ;  L = -1/2 g3 h H^2
    return -2.0*g_hHH_paper(mH, tanb)

def Gamma_hHH(mH, tanb):
    g3 = g_hHH3(mH, tanb)
    arg = 1.0-4.0*mH**2/MH_SM**2
    if arg <= 0.0: return 0.0
    return g3*g3/(32.0*math.pi*MH_SM)*math.sqrt(arg)

# ---------------------------------------------------------------- App.A form factors (paper normalisation)
def f_loop(tau):
    """f(tau) of App.A, complex for tau>1."""
    if tau <= 0.0:  return 0j
    if tau <= 1.0:  return complex(math.asin(math.sqrt(tau))**2, 0.0)
    r = math.sqrt(1.0-1.0/tau)
    return -0.25*(cmath.log((1+r)/(1-r)) - 1j*math.pi)**2

def A_0(tau):    return -0.5*(tau-f_loop(tau))/tau**2
def A_12(tau):   return (tau+(tau-1.0)*f_loop(tau))/tau**2
def A_1(tau):    return -0.5*(2.0*tau**2+3.0*tau+3.0*(2.0*tau-1.0)*f_loop(tau))/tau**2

def kappa_g(tanb):
    c, s = mixing(tanb)
    kt = s + c/tanb          # xi_h^f
    kb = s + c/tanb
    tt = MH_SM**2/(4.0*MQ['t']**2); tb = MH_SM**2/(4.0*MQ['b']**2)
    return (kt*A_12(tt)+kb*A_12(tb))/(A_12(tt)+A_12(tb))

def kappa_gamma(tanb):
    c, s = mixing(tanb)
    kf = s + c/tanb          # xi_h^f  (all fermions)
    kV = s                   # xi_h^V
    num = (4.0/3.0)*kf*A_12(MH_SM**2/(4*MQ['c']**2)) \
        + (1.0/3.0)*kf*A_12(MH_SM**2/(4*MQ['b']**2)) \
        + (4.0/3.0)*kf*A_12(MH_SM**2/(4*MQ['t']**2)) \
        + (1.0/3.0)*kf*A_12(MH_SM**2/(4*MQ['s']**2)) \
        + (4.0/3.0)*kf*A_12(MH_SM**2/(4*MQ['u']**2)) \
        + (1.0/3.0)*kf*A_12(MH_SM**2/(4*MQ['d']**2)) \
        + kf*A_12(MH_SM**2/(4*0.10566**2)) + kf*A_12(MH_SM**2/(4*1.7768**2)) \
        + kV*A_1(MH_SM**2/(4*MW**2))
    den = (4.0/3.0)*A_12(MH_SM**2/(4*MQ['c']**2)) \
        + (1.0/3.0)*A_12(MH_SM**2/(4*MQ['b']**2)) \
        + (4.0/3.0)*A_12(MH_SM**2/(4*MQ['t']**2)) \
        + (1.0/3.0)*A_12(MH_SM**2/(4*MQ['s']**2)) \
        + (4.0/3.0)*A_12(MH_SM**2/(4*MQ['u']**2)) \
        + (1.0/3.0)*A_12(MH_SM**2/(4*MQ['d']**2)) \
        + A_12(MH_SM**2/(4*0.10566**2)) + A_12(MH_SM**2/(4*1.7768**2)) \
        + A_1(MH_SM**2/(4*MW**2))
    return abs(num)/abs(den)   # kappa_gamma as a real ratio of magnitudes (Gamma_gg ~ |..|^2)

def higgs_sector(mH, tanb):
    c, s = mixing(tanb)
    kf = s + c/tanb; kV = s
    kg, kgam = abs(kappa_g(tanb)), kappa_gamma(tanb)
    kappas = dict(f=kf, V=kV, g=kg, gamma=kgam, b=kf)
    g_other = 0.0
    for ch, br in BR_SM.items():
        if ch in ('bb', 'ss', 'mumu', 'cc', 'tautau'): k = kf
        elif ch in ('WW', 'ZZ'):                        k = kV
        elif ch == 'gg':                                k = kg
        elif ch == 'gaga':                             k = abs(kgam)
        elif ch == 'Zgamma':                           k = abs(kgam)   # approximated by the gamma combination
        g_other += GAM_H_SM*br*k*k
    gHH = Gamma_hHH(mH, tanb)
    br_hHH = gHH/(g_other + gHH)
    # 9-process scaling
    sig = {}
    sig['ggF'] = SIGMA_SM_PB['ggF']*abs(kg)**2
    for p in ('VBF', 'WH', 'ZH'):
        sig[p] = SIGMA_SM_PB[p]*kV**2
    for p in ('ttH', 'tHj', 'tHb'):
        sig[p] = SIGMA_SM_PB[p]*kf**2
    sig['bbH'] = SIGMA_SM_PB['bbH']*kf**2
    # tHW: sigma = sigma_SM*(kt^2+kV^2+2 kt kV cos(phi))/4 baseline (DESIGN 3.5)
    sig['tHW'] = SIGMA_SM_PB['tHW']*(kf**2+kV**2+2.0*kf*kV)/4.0
    sig_total = sum(sig.values())
    return dict(kappas=kappas, g_C=g_hHH_paper(mH, tanb), g3=g_hHH3(mH, tanb),
                Gamma_hHH=gHH, BR_hHH=br_hHH, Gamma_h_other=g_other,
                Gamma_h_total=g_other+gHH, sigma=sig, sigma_total=sig_total,
                sigma_tHW_min=SIGMA_SM_PB['tHW']*min(kf**2, kV**2),
                sigma_tHW_max=SIGMA_SM_PB['tHW']*max(kf**2, kV**2))

# ---------------------------------------------------------------- App.D g0,g1,g2
def g012(x, xH, cotb):
    """Type-I 2HDM aux functions of App.D.  x = m_k^2/m_W^2, xH = m_{H^pm}^2/m_W^2."""
    d = xH - x
    L = math.log(xH/x)
    g0 = -cotb**2*(3*xH**2 - 4*xH*x + x**2 - 2*x*(2*xH-x)*L)/(16*d**3)
    g1 = -0.75 + cotb**2*x*(5*xH**2 - 8*xH*x + 3*x**2 - 2*xH*(2*xH-x)*L)/(4*d**3)
    X1 = -0.25*((xH/d - 6.0/(x-1)**2 + 3.0)
                - xH*(3*xH-2*x)/d**2*math.log(xH)
                + (xH*(3*xH-2*x)/d**2 + 3.0*(x+1)/(x-1)**3)*math.log(x))
    X2 = x*(5*xH-3*x)/(4*d**2) - xH*x*(2*xH-x)/(2*d**3)*L
    g2 = cotb*X1 + cotb**3*X2
    return g0, g1, g2

def xiH(mH, tanb, i, j, return_terms=False):
    """xi_H^{ij} for down-type (i,j) in {0=d,1=s,2=b}, Eq.(xi_H_2HDM).
    L_eff = (phi/v) xi^{ij} m_{f_j} fbar_i P_R f_j + h.c."""
    c, s = mixing(tanb)
    cotb = 1.0/tanb
    xH = MHC**2/MW**2
    lam = lam_HHpHm(mH, tanb)
    terms = {}
    tot = 0j
    for k, name in enumerate(('u', 'c', 't')):
        x = MQ[name]**2/MW**2
        g0, g1, g2 = g012(x, xH, cotb)
        br = g1*c - g2*s - g0*(2.0*V_VEV/MW**2)*lam
        t = VCKM[k, i].conjugate()*MQ[name]**2*br*VCKM[k, j]
        terms[name] = t
        tot += t
    pref = -(4.0*math.sqrt(2.0)*GF)/(16.0*math.pi**2)
    out = pref*tot
    if return_terms:
        return out, {k: pref*v for k, v in terms.items()}
    return out

def xi_bs(mH, tanb, return_terms=False): return xiH(mH, tanb, 1, 2, return_terms)
def xi_ds(mH, tanb, return_terms=False): return xiH(mH, tanb, 0, 1, return_terms)

# ---------------------------------------------------------------- meson rates
def f_phase(x):
    if x <= 0.0:
        return 1.0
    return (1.0-8.0*x+x*x)*(1.0-x*x)-12.0*x*x*math.log(x)

def br_B_to_XsH(mH, tanb):
    """Inclusive Br(B -> X_s H), Eq. after (Lfceff_phi).  Returns (B+, B0) absolute BR."""
    if mH >= MB_XSH:
        return 0.0, 0.0
    xi = abs(xi_bs(mH, tanb))
    x = MC_XSH**2/MB_XSH**2
    ratio = (12.0*math.pi**2*V_VEV**2/MB_XSH**2)*(1.0-mH**2/MB_XSH**2)**2/f_phase(x)*abs(xi/VCKM[1,2])**2
    return ratio*BR_BXCE, ratio*BR_BXCE*(0.93)   # Gamma_B+/Gamma_B0 = 0.93 -> BR_B0 = 0.93*BR_B+
    # (the 0.93 enters through Gamma; BR ratio = Gamma_B+/Gamma_B0)

def pstar(mX, mY, mS):
    if mX < mY+mS: return 0.0
    return math.sqrt(max((mX**2-(mY+mS)**2)*(mX**2-(mY-mS)**2), 0.0))/(2.0*mX)

def M_KpiH(mH, tanb, neutral=False):
    """Coherent K -> pi H amplitude, Eq.(Matrix_4FH), absolute value in GeV (see notes)."""
    c, s = mixing(tanb)
    mK = M_K0L if neutral else M_KP
    mpi = M_PI0 if neutral else M_PIP
    xiW = c
    xi = xi_ds(mH, tanb)
    fourq = math.sqrt(GF)*2.0**0.25*xiW*(7.0*LAMBDA_K*(mK**2+mpi**2-mH**2)/18.0
                                         - 7.0*A_K*mK**2/9.0)
    fcnc = xi/(2.0*V_VEV)*MS_QUARK*(mK**2-mpi**2)/(MS_QUARK-MD_QUARK)*F0_KPI
    return fourq+fcnc, fourq, fcnc

def br_K_to_piH(mH, tanb, which='K+'):
    if which == 'K+':
        mK, mpi, GK = M_KP, M_PIP, GAM_KP
    elif which == 'K0L':
        mK, mpi, GK = M_K0L, M_PI0, GAM_K0L
    else:
        raise KeyError(which)
    p = pstar(mK, mpi, mH)
    if p <= 0.0: return 0.0
    M = M_KpiH(mH, tanb, neutral=(which == 'K0L'))[0]
    return p*abs(M)**2/(8.0*math.pi*mK**2*GK)

def br_semileptonic(mH, tanb, which='pi'):
    c, s = mixing(tanb)
    if which == 'pi':
        mX, Br = M_PIP, BR_PI_MUNU
    elif which == 'K':
        mX, Br = M_KP, BR_K_MUNU
    else:
        raise KeyError(which)
    x = mH**2/mX**2
    if x >= 1.0: return 0.0
    return (math.sqrt(2.0)*GF*mX**4*c**2
            /(96.0*math.pi**2*M_MU**2*(1.0-M_MU**2/mX**2)**2))*Br*f_phase(x)*(7.0/9.0)**2

# ---------------------------------------------------------------- driver
if __name__ == '__main__':
    print('CKM |V| ='); print(np.abs(VCKM).round(6))
    for mH in (0.021, 1.0, 2.0):
        for tanb in (5.0, 100.0, 1000.0):
            xi, terms = xi_bs(mH, tanb, return_terms=True)
            print(f'mH={mH:7.3f} tanb={tanb:8.1f} |xi_bs|={abs(xi):.6e} '
                  f'Br(B->XsH)={br_B_to_XsH(mH,tanb)[0]:.4e} |xi_ds|={abs(xi_ds(mH,tanb)):.6e}')
    hs = higgs_sector(1.0, 100.0)
    print('higgs sector at mH=1, tanb=100:', json.dumps(
        {k: v for k, v in hs.items() if k != 'sigma'}, indent=1, default=str))
    print('sigma:', hs['sigma'], 'total', hs['sigma_total'])
