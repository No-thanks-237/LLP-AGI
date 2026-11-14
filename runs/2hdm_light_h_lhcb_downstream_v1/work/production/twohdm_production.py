# -*- coding: utf-8 -*-
"""Type-I 2HDM light CP-even H: production module (DESIGN.md section 3).

Pure model physics, no detector, no MC.  Implements

  3.1  the effective FCNC coupling xi_H^{ij} (paper eq. xi_H_2HDM + App. D),
  3.2  B -> X_s H  (inclusive FASER-paper formula; exclusive mixture only as shape),
  3.3  K -> pi H   (coherent four-quark + FCNC amplitude, paper Eq. Matrix_4FH),
  3.4  pi/K -> H e nu (paper formula after n_h,n_l),
  3.5  h couplings, h -> HH and the 2HDM-scaled SM Higgs cross sections.

Reference: Li, Song, Su, Wang, arXiv:2212.06186 (inputs/references/FASER_light.tex).
Boiarska et al. arXiv:1904.10447 for the exclusive B -> X_s channel table/shape only.

Conventions (fixed by DESIGN.md section 1):
  m_h = 125, v = 246, m_A = m_Hpm = 600, lambda v^2 = 0,
  c = cos(beta-alpha) = 1/tanb, s = sin(beta-alpha) = +sqrt(1-c^2),
  xi_H^V = c, xi_H^f = c - s/tanb (= 1/(tanb^3 (1+s))), xi_h^V = s, xi_h^f = s + c/tanb,
  lambda_HH+H- = -(2 m_Hpm^2 - m_H^2) c / v.

Sign conventions of the xi_H sum (paper eq. xi_H_2HDM, lower row = H):
    xi_H^{ij} = -(4 sqrt2 G_F)/(16 pi^2) sum_k V*_ki m_k^2
                [ g1 c - g2 s - g0 (2v/m_W^2) lambda_HH+H- ] V_kj
(the g0 term is *minus* g0 times the trilinear; lambda_HH+H- < 0 for c > 0, so the
g0 contribution flips sign relative to g1 c at large tanb.  Do not drop it.)

Public API required by response/scan (fixed by the Lead):
  br_B_to_XsH(mH, tanb, meson='B+', mode=None) -> {'total', 'channels', 'total_exclusive', 'mode'}
        mode = 'foresee' (DEFAULT: B-meson phase space, paper's FORESEE 2HDM_H treatment)
             | 'paper_mb' (literal Sec. 2.2 formula, closes at m_b=4.18)
             | 'exclusive' (Boiarska exclusive sum); DEFAULT_B_MODE is module-level.
  br_B_to_XsH_exclusive(mH, tanb, meson='B+') -> {'total', 'channels'} (absolute per channel)
  channel_masses()                    -> {label: m_X' [GeV]}
  br_K_to_piH(mH, tanb, which='K+')   -> float   (K0L: per-K0L BR, NO 0.5 projection)
  br_K_to_Henu(mH, tanb, which='K+')  -> float   (K0L: 0.0)
  br_pi_to_Henu(mH, tanb)             -> float
  higgs_sector(mH, tanb)              -> dict (sigma_process, sigma_total, BR_hHH, ...)
  evaluate(mH, tanb)                  -> flat dict with everything
  evaluate_grid(points)               -> list of flat dicts (batch helper)
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import b_channels  # noqa: E402

# ======================================================================================
# model version / fixed model parameters (DESIGN.md section 1 / 3)
# ======================================================================================
MODEL_TAG = '2HDM-lightH-production-1.0.1'  # 1.0.1: path refactor only (work/common -> work/analysis/grids), physics unchanged
__version__ = '1.0.0'
GF = 1.1663787e-5          # GeV^-2
V_VEV = 246.0              # GeV
M_W = 80.385               # GeV
M_H_SM = 125.0             # m_h [GeV]
M_HPM = 600.0              # m_{H^pm} [GeV]
M_A_HEAVY = 600.0          # m_A [GeV] (only enters through lambda_HH+H- via m_Hpm)
LAMBDA_V2 = 0.0            # lambda v^2 [GeV^2]

# loop-quark masses (MSbar / PDG choices, as fixed in DESIGN.md 3.1)
M_QUARK = {'u': 2.2e-3, 'c': 1.28, 't': 173.1}      # GeV
# quark masses entering the meson rates
M_B_QUARK = 4.18           # b->s inclusive formula
M_C_QUARK = 1.28           # f(m_c^2/m_b^2)
M_K0_NEUTRAL = 0.497611    # PDG K0 (lightest recoil of a neutral b -> s transition)
M_S_QUARK = 0.096          # K->pi H FCNC term (DESIGN.md 3.3; the decay program's
M_D_QUARK = 0.0047         #   chiral Lagrangian uses 0.0924/0.00467 - recorded)

# CKM (PDG 2022).  Standard parametrisation with (s12,s23,s13) from |Vus|,|Vcb|,|Vub|
# and delta = gamma = 1.144 rad; this set is exactly unitary.
CKM_S12, CKM_S23, CKM_S13 = 0.2243, 0.0408, 0.00382
CKM_GAMMA = 1.144          # rad
CKM_MAGNITUDES_PDG = {     # DESIGN.md 3.1 nominal list (rounded, for documentation)
    (0, 0): 0.97373, (0, 1): 0.2243, (0, 2): 0.00382,
    (1, 0): 0.221,   (1, 1): 0.975,  (1, 2): 0.0408,
    (2, 0): 0.0086,  (2, 1): 0.0415, (2, 2): 0.999,
}

# K -> pi H constants (DESIGN.md 3.3)
LAMBDA_QCD = 3.1e-7
A_K = 0.0
F0_KPI = 0.96
M_KP = 0.493677
M_K0L = 0.497611
M_PI_P = 0.13957039
M_PI0 = 0.1349768
GAM_KP = 5.314e-17          # GeV
GAM_K0L = 1.287e-17         # GeV
BR_PI_MUNU = 0.999877
BR_K_MUNU = 0.6356
M_MU = 0.1056583755

# B -> X_s H
BR_B_XC_ENU = 0.104
M_BP = 5.27934            # B+ mass
M_B0 = 5.27965            # B0 mass
# Rate mode: 'foresee' uses the B-meson phase space (1-mH^2/m_B^2)^2, i.e. the paper's
# own FORESEE 2HDM_H model treatment (B endpoint ~5.28 GeV); 'paper_mb' is the literal
# Sec. 2.2 formula with m_b (closes at 4.18); 'exclusive' is the Boiarska exclusive sum.
DEFAULT_B_MODE = 'foresee'
B_MODES = ('foresee', 'paper_mb', 'exclusive')
# Gamma(B->XsH) is meson independent while the total widths differ:
#   BR(B0->XsH) = BR(B+->XsH) * Gamma_B+/Gamma_B0
# with tau_B+ = 1.638 ps, tau_B0 = 1.519 ps (same convention as Boiarska's caption and
# FORESEE's mB/Gamma_B normalisation); the module's B widths give 0.92735.
B0_LIFETIME_SCALE = b_channels.GAM_BP / b_channels.GAM_B0

# h -> HH / Higgs sector (DESIGN.md 3.5)
GAM_H_SM = 4.07e-3          # GeV
BR_H_SM = {                 # HXSWG-like SM branching ratios at m_h = 125 GeV
    'bb': 0.5824, 'WW': 0.2137, 'gg': 0.0818, 'tautau': 0.0627, 'cc': 0.0288,
    'ZZ': 0.0262, 'gamma gamma': 0.00227, 'Zgamma': 0.00153, 'ss': 0.0003,
    'mumu': 0.00022,
}
BR_H_SM_NORM = sum(BR_H_SM.values())

SIGMA_H = {                 # 13 TeV SM Higgs production cross sections [pb]
    'ggF': 48.58, 'VBF': 3.782, 'WH': 1.369, 'ZH': 0.883, 'bbH': 0.2384,
    'ttH': 0.5065, 'tHj': 0.0740, 'tHW': 0.0152, 'tHb': 0.0046,
}
SIGMA_H_TOTAL_SM = sum(SIGMA_H.values())

HBAR_C = 1.973269804e-16    # GeV m


# ======================================================================================
# CKM
# ======================================================================================
def ckm_matrix(s12=CKM_S12, s23=CKM_S23, s13=CKM_S13, delta=CKM_GAMMA):
    """Complex CKM matrix, standard (PDG) parametrisation; exactly unitary.

    With |Vus|=0.2243, |Vcb|=0.0408, |Vub|=0.00382 and gamma=1.144 rad this gives
        |Vud|=0.97374 |Vus|=0.22430 |Vub|=0.00382
        |Vcd|=0.22104 |Vcs|=0.97486 |Vcb|=0.04080
        |Vtd|=0.008333|Vts|=0.040082|Vtb|=0.999167
    The DESIGN.md nominal list quotes |Vtd|=0.0086, |Vts|=0.0415 (B-mixing
    measurements); those are not exactly unitary together with |Vus|,|Vcb|,|Vub|.
    See logs/production_validation.md: the |Vtd Vts| difference is +6.5 % in
    amplitude (+13 % in BR(K->pi H)); ckm_matrix_from_magnitudes() exposes it.
    """
    c12, c23, c13 = math.sqrt(1 - s12 ** 2), math.sqrt(1 - s23 ** 2), math.sqrt(1 - s13 ** 2)
    e = np.exp(1j * delta)
    return np.array([
        [c12 * c13, s12 * c13, s13 / e],
        [-s12 * c23 - c12 * s23 * s13 * e, c12 * c23 - s12 * s23 * s13 * e, s23 * c13],
        [s12 * s23 - c12 * c23 * s13 * e, -c12 * s23 - s12 * c23 * s13 * e, c23 * c13],
    ], dtype=complex)


def ckm_matrix_from_magnitudes(mag=None, delta=CKM_GAMMA):
    """The DESIGN.md nominal |V| list with the standard parametrisation phases.

    Not exactly unitary; kept only for the sensitivity cross-check in validation.
    """
    mag = mag or CKM_MAGNITUDES_PDG
    V = ckm_matrix()
    out = np.zeros((3, 3), dtype=complex)
    for (i, j), m in mag.items():
        out[i, j] = m * V[i, j] / abs(V[i, j])
    return out


V_CKM = ckm_matrix()
CKM_MAG = np.abs(V_CKM)
_QFLAV_DOWN = {'d': 0, 's': 1, 'b': 2}


def ckm_sums(V=None):
    """S_ij = sum_k V*_ki m_k^2 V_kj for i,j in {d,s,b}  (complex, GeV^2)."""
    V = V_CKM if V is None else V
    m2 = np.array([M_QUARK['u'] ** 2, M_QUARK['c'] ** 2, M_QUARK['t'] ** 2])
    S = np.zeros((3, 3), dtype=complex)
    for i in ('d', 's', 'b'):
        for j in ('d', 's', 'b'):
            S[_QFLAV_DOWN[i], _QFLAV_DOWN[j]] = np.sum(np.conj(V[:, _QFLAV_DOWN[i]]) * m2
                                                       * V[:, _QFLAV_DOWN[j]])
    return S


_S_CACHE = {}


def _S():
    if 'default' not in _S_CACHE:
        _S_CACHE['default'] = ckm_sums(V_CKM)
    return _S_CACHE['default']


# ======================================================================================
# 3.1  xi_H^{ij}
# ======================================================================================
_G_CACHE = {}


def g_aux_functions(tanb):
    """App. D Type-I functions g0,g1,g2 for k=u,c,t at fixed tanb.

    g0 = -cot^2 b (3xH^2-4xH x+x^2-2x(2xH-x)ln(xH/x)) / (16 (xH-x)^3)
    g1 = -3/4 + cot^2 b x[5xH^2-8xH x+3x^2-2xH(2xH-x)ln(xH/x)] / (4 (xH-x)^3)
    g2 = cot b X1 + cot^3 b X2
    X1 = -1/4{ xH/(xH-x) - 6/(x-1)^2 + 3 - xH(3xH-2x)/(xH-x)^2 ln xH
               + [xH(3xH-2x)/(xH-x)^2 + 3(x+1)/(x-1)^3] ln x }
    X2 = x(5xH-3x)/(4(xH-x)^2) - xH x(2xH-x)/(2(xH-x)^3) ln(xH/x)
    with x_k = m_k^2/m_W^2, x_H = m_Hpm^2/m_W^2.
    """
    key = round(tanb, 10)
    if key in _G_CACHE:
        return _G_CACHE[key]
    cotb = 1.0 / tanb
    cotb2 = cotb * cotb
    cotb3 = cotb2 * cotb
    xH = M_HPM ** 2 / M_W ** 2
    lnxH = math.log(xH)
    out = {}
    for name, m in M_QUARK.items():
        x = m * m / M_W ** 2
        d = xH - x
        lnr = math.log(xH / x)
        # g0
        g0 = -cotb2 * (3.0 * xH * xH - 4.0 * xH * x + x * x
                       - 2.0 * x * (2.0 * xH - x) * lnr) / (16.0 * d ** 3)
        # g1
        g1 = -0.75 + cotb2 * x * (5.0 * xH * xH - 8.0 * xH * x + 3.0 * x * x
                                  - 2.0 * xH * (2.0 * xH - x) * lnr) / (4.0 * d ** 3)
        # X1
        A = xH * (3.0 * xH - 2.0 * x) / d ** 2
        X1 = -0.25 * ((xH / d - 6.0 / (x - 1.0) ** 2 + 3.0)
                      - A * lnxH
                      + (A + 3.0 * (x + 1.0) / (x - 1.0) ** 3) * math.log(x))
        # X2
        X2 = (x * (5.0 * xH - 3.0 * x) / (4.0 * d ** 2)
              - xH * x * (2.0 * xH - x) / (2.0 * d ** 3) * lnr)
        g2 = cotb * X1 + cotb3 * X2
        out[name] = {'g0': g0, 'g1': g1, 'g2': g2, 'X1': X1, 'X2': X2, 'x': x}
    _G_CACHE[key] = out
    return out


def lambda_HHpmHpm(mH, tanb):
    """lambda_{HH+H-} = -(2 m_Hpm^2 - m_H^2) c / v   (DESIGN.md section 1)."""
    c = 1.0 / tanb
    return -(2.0 * M_HPM ** 2 - mH * mH) * c / V_VEV


def xi_H_terms(mH, tanb):
    """The three bracket terms per loop quark: (g1 c, -g2 s, -g0 2v/mW^2 lambda_HH)."""
    c = 1.0 / tanb
    s = math.sqrt(max(1.0 - c * c, 0.0))
    lam = lambda_HHpmHpm(mH, tanb)
    k = 2.0 * V_VEV / (M_W * M_W) * lam
    out = {}
    G = g_aux_functions(tanb)
    for q in ('u', 'c', 't'):
        g = G[q]
        out[q] = {'g1_c': g['g1'] * c,
                  'minus_g2_s': -g['g2'] * s,
                  'minus_g0_k': -g['g0'] * k,
                  'bracket': g['g1'] * c - g['g2'] * s - g['g0'] * k}
    return {'c': c, 's': s, 'lambda_HH': lam, 'k': k, 'per_quark': out}


_TERM_CACHE = {}


def _terms_cached(mH, tanb):
    key = (round(mH, 10), round(tanb, 10))
    if key not in _TERM_CACHE:
        _TERM_CACHE[key] = xi_H_terms(mH, tanb)
    return _TERM_CACHE[key]


def xi_H(mH, tanb, i='s', j='b'):
    """Complex xi_H^{ij}, i,j in {d,s,b}: sum_k V*_ki m_k^2 [g1 c - g2 s - g0 k] V_kj."""
    c = 1.0 / tanb
    s = math.sqrt(max(1.0 - c * c, 0.0))
    lam = lambda_HHpmHpm(mH, tanb)
    k = 2.0 * V_VEV / (M_W * M_W) * lam
    G = g_aux_functions(tanb)
    m2 = np.array([M_QUARK['u'] ** 2, M_QUARK['c'] ** 2, M_QUARK['t'] ** 2])
    V = V_CKM
    ii, jj = _QFLAV_DOWN[i], _QFLAV_DOWN[j]
    tot = 0.0 + 0.0j
    for iq, q in enumerate(('u', 'c', 't')):
        g = G[q]
        br = g['g1'] * c - g['g2'] * s - g['g0'] * k
        tot += np.conj(V[iq, ii]) * m2[iq] * br * V[iq, jj]
    return -(4.0 * math.sqrt(2.0) * GF / (16.0 * math.pi ** 2)) * tot


def xi_H_details(mH, tanb, i='s', j='b', V=None):
    """Per-term breakdown of xi_H^{ij} for the hand check (validation item a)."""
    c = 1.0 / tanb
    s = math.sqrt(max(1.0 - c * c, 0.0))
    lam = lambda_HHpmHpm(mH, tanb)
    k = 2.0 * V_VEV / (M_W * M_W) * lam
    G = g_aux_functions(tanb)
    V = V_CKM if V is None else V
    m2 = {'u': M_QUARK['u'] ** 2, 'c': M_QUARK['c'] ** 2, 't': M_QUARK['t'] ** 2}
    ii, jj = _QFLAV_DOWN[i], _QFLAV_DOWN[j]
    pref = -(4.0 * math.sqrt(2.0) * GF / (16.0 * math.pi ** 2))
    rows = {}
    tot = 0.0 + 0.0j
    for iq, q in enumerate(('u', 'c', 't')):
        g = G[q]
        t1, t2, t3 = g['g1'] * c, -g['g2'] * s, -g['g0'] * k
        br = t1 + t2 + t3
        ckm = np.conj(V[iq, ii]) * V[iq, jj]
        contrib = m2[q] * ckm * br
        rows[q] = {'m': M_QUARK[q], 'm2': m2[q], 'x_k': g['x'], 'g0': g['g0'],
                   'g1': g['g1'], 'g2': g['g2'], 'X1': g['X1'], 'X2': g['X2'],
                   'g1_c': t1, 'minus_g2_s': t2, 'minus_g0_k': t3, 'bracket': br,
                   'Vki_conj': np.conj(V[iq, ii]), 'Vkj': V[iq, jj],
                   'ckm_product': ckm, 'm2_ckm_bracket': contrib}
        tot += contrib
    return {'prefactor': pref, 'sum': tot, 'xi_H': pref * tot, 'lambda_HH': lam,
            'k': k, 'c': c, 's': s, 'rows': rows}


# ======================================================================================
# 3.2  B -> X_s H
# ======================================================================================
def _f_phase_space(x):
    """f(x) = (1-8x+x^2)(1-x^2) - 12 x^2 ln x  (paper phase-space factor)."""
    if x <= 0.0:
        return 1.0
    return (1.0 - 8.0 * x + x * x) * (1.0 - x * x) - 12.0 * x * x * math.log(x)


F_PHASE_B = _f_phase_space((M_C_QUARK / M_B_QUARK) ** 2)


def _b_ratio_and_total(mH, tanb, meson='B+', mode=None):
    """Br(B->XsH)/Br(B->Xc e nu) and Br(B->XsH) for one rate mode (see DEFAULT_B_MODE)."""
    if meson not in ('B+', 'B0'):
        raise KeyError(meson)
    mode = DEFAULT_B_MODE if mode is None else mode
    if mode not in B_MODES:
        raise ValueError('mode must be one of %s' % (B_MODES,))
    mB = M_BP if meson == 'B+' else M_B0
    # lightest strangeness-carrying recoil: K+- for B+-, K0 for B0 (differs by 4 MeV!)
    m_recoil = M_KP if meson == 'B+' else M_K0_NEUTRAL
    scale = 1.0 if meson == 'B+' else B0_LIFETIME_SCALE   # Gamma_B+/Gamma_B0
    if mode == 'exclusive':
        # the exclusive sum already carries the meson's own width mX/Gamma_X AND its recoil
        # masses (K+- for B+-, K0 for B0), so evaluate it directly for this meson: the old
        # "B+ basis x lifetime scale" is only equivalent while both mesons share the same
        # channel phase space (review F2; differs at the endpoint where B0 -> K0 H closes)
        total, _ch = br_B_to_XsH_exclusive_parts(mH, tanb, meson)
        return total / BR_B_XC_ENU, total
    m_ps = mB if mode == 'foresee' else M_B_QUARK
    if mode == 'foresee' and mH >= mB - m_recoil:
        # mesonic endpoint: no B -> X_s H channel is open above mB - m_recoil
        # (B+- : 4.785663 GeV with K+-, B0 : 4.782039 GeV with K0).  The FORESEE mB
        # phase space would otherwise extrapolate unphysically up to mB.
        return 0.0, 0.0
    if mH >= m_ps:
        return 0.0, 0.0
    xib = abs(xi_H(mH, tanb, 's', 'b'))
    vcb = abs(V_CKM[1, 2])
    ratio = (12.0 * math.pi ** 2 * (V_VEV ** 2 / M_B_QUARK ** 2)
             * (1.0 - mH * mH / m_ps ** 2) ** 2 / F_PHASE_B
             * (xib / vcb) ** 2)
    return ratio * scale, ratio * scale * BR_B_XC_ENU


def br_B_to_XsH(mH, tanb, meson='B+', mode=None):
    """Absolute BR(B -> X_s H), inclusive FASER-paper formula; exclusive channels
    renormalised to the *same total* so that they only carry the kinematic shape.

    Br(B->XsH)/Br(B->Xc e nu) = 12 pi^2 (v^2/m_b^2)(1-mH^2/m_b^2)^2 / f(mc^2/mb^2)
                                * |xi_H^{bs}/Vcb|^2

    Gamma(B->XsH) is meson independent while the total widths differ, so in ALL modes
        BR(B0->XsH) = BR(B+->XsH) * Gamma_B+/Gamma_B0 = 0.92735 * BR(B+->XsH)
    (tau_B+ = 1.638 ps, tau_B0 = 1.519 ps; same convention as Boiarska's caption and
    FORESEE's mB/Gamma_B).  The normalised 'channels' fractions are identical for both
    mesons and sum(channels) == total for each.  B_s is NOT included.

    mode='foresee' returns total = 0 for mH >= m_B - m_K (mesonic endpoint, since K is
    the lightest strangeness-carrying recoil); 'beyond_mesonic_endpoint' flags it.
    Returns {'total', 'channels', 'sum_channels', 'total_exclusive', 'meson', 'mode',
             'ratio_to_BXcenu', 'beyond_mesonic_endpoint', 'channel_shape_extrapolated'}.
    """
    mode = DEFAULT_B_MODE if mode is None else mode
    ratio, total = _b_ratio_and_total(mH, tanb, meson, mode)
    w = b_channels.channel_weights(mH, meson)
    if total > 0.0 and sum(w.values()) <= 0.0:
        # mH in (m_B - m_K, m_B): the two-body b->s channels are all closed, the true
        # recoil X_s is multibody and the H is soft in the B rest frame.  Use the K
        # template (softest available) for the kinematic shape only; flagged here.
        w = {k: (1.0 if k == 'K' else 0.0) for k in w}
    chans = {k: v * total for k, v in w.items()}
    # make sum(channels) == total to machine precision
    if chans and total > 0.0:
        key = max(chans, key=lambda k: chans[k])
        chans[key] += total - sum(chans.values())
    excl_total, _excl_ch = br_B_to_XsH_exclusive_parts(mH, tanb, meson)
    mB = M_BP if meson == 'B+' else M_B0
    m_recoil = M_KP if meson == 'B+' else M_K0_NEUTRAL
    return {'total': total, 'channels': chans, 'sum_channels': sum(chans.values()),
            'total_exclusive': excl_total, 'meson': meson, 'mode': mode,
            'ratio_to_BXcenu': ratio, 'branching_ratio_over_BXcenu': ratio,
            'beyond_mesonic_endpoint': bool(mH >= mB - m_recoil),
            'channel_shape_extrapolated': bool(
                total > 0.0
                and sum(b_channels.channel_weights(mH, meson).values()) <= 0.0)}


def br_B_to_XsH_exclusive_parts(mH, tanb, meson='B+'):
    """Boiarska exclusive b->s channel BRs with the 2HDM xi_H^{bs} (VARIANT rate).

    theta^2-normalised Boiarska channel BRs (inputs/bc5_reuse/bc5_production.py
    br_B_to_S, only the b->s channels: the b->d 'pi' entry is dropped) rescaled by
    (xi_H^{bs}/xi_Boiarska_sb)^2.  Same xi_H^{bs} and same B+/B0 lifetime as the
    nominal inclusive formula.  Extends to the kinematic closing masses of each
    channel (B->K S closes at m_B - m_K = 4.786 GeV) instead of the free-quark
    m_b = 4.18 threshold.  Returns (total, {label: BR})."""
    if meson not in ('B+', 'B0'):
        raise KeyError(meson)
    xib2 = abs(xi_H(mH, tanb, 's', 'b')) ** 2
    if meson == 'B+':
        shape = b_channels.shape_channels_cached(mH)
    else:
        shape = b_channels.shape_channels(mH, 'B0')
    ch = {k: v * xib2 for k, v in shape.items()}
    return sum(ch.values()), ch


def br_B_to_XsH_exclusive(mH, tanb, meson='B+'):
    """SYSTEMATIC VARIANT of the B -> X_s H rate (see br_B_to_XsH_exclusive_parts).

    Returns {'total', 'channels', 'meson'}; sum(channels) == total by construction.
    The nominal/paper prescription stays br_B_to_XsH()['total']; this variant is the
    Boiarska-exclusive estimate and differs by the known factor ~2 model ambiguity
    (Boiarska eq:exclusive-to-inclusive)."""
    total, ch = br_B_to_XsH_exclusive_parts(mH, tanb, meson)
    return {'total': total, 'channels': ch, 'meson': meson}


def channel_masses(meson='B+'):
    """{exclusive label: m_X' [GeV]} for the b->s channels only (K0 for a neutral B)."""
    return b_channels.channel_masses(meson)


def b_channel_weights(mH, meson='B+'):
    """Normalised kinematic mixture (sum 1) over the open b->s exclusive channels."""
    return b_channels.channel_weights(mH, meson)


# --------------------------------------------------------------------------------------
# validity of the effective/inclusive description: Br(B->XsH) must stay << 1
# --------------------------------------------------------------------------------------
def _br_total(mH, tanb):
    """Br(B->XsH) in the DEFAULT_B_MODE (used by the validity boundaries)."""
    return _b_ratio_and_total(mH, tanb, 'B+', DEFAULT_B_MODE)[1]


def tanb_valid(mH, br_target=0.1, tanb_lo=5.0, tanb_hi=1e5):
    """tanb at which Br(B->XsH) = br_target (BR is monotonically decreasing in tanb).

    The region tanb < tanb_valid is outside the validity of the inclusive/effective
    description (the branching ratio is no longer a small correction).  Returns None
    when there is no crossing inside [tanb_lo, tanb_hi] (BR already below the target at
    tanb_lo, or mH above the m_b threshold)."""
    f = lambda tb: _br_total(mH, tb) - br_target   # noqa: E731
    flo, fhi = f(tanb_lo), f(tanb_hi)
    if flo <= 0.0:
        return None
    if fhi > 0.0:
        return None
    lo, hi = math.log(tanb_lo), math.log(tanb_hi)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(math.exp(mid)) > 0.0:
            lo = mid
        else:
            hi = mid
    return math.exp(0.5 * (lo + hi))


def validity_boundaries(mH):
    """{'tanb_BR_gt_0.1': t, 'tanb_BR_gt_1': t} (None if no crossing)."""
    return {'tanb_BR_gt_0.1': tanb_valid(mH, 0.1),
            'tanb_BR_gt_1': tanb_valid(mH, 1.0)}


# ======================================================================================
# 3.3  K -> pi H
# ======================================================================================
def _k_meson(which):
    if which in ('K+', 'Kp', 'K^+'):
        return M_KP, M_PI_P, GAM_KP
    if which in ('K0L', 'KL', 'K_L'):
        return M_K0L, M_PI0, GAM_K0L
    raise KeyError(which)


def k_to_pi_amplitude(mH, tanb, which='K+'):
    """Complex amplitude M(K -> pi H).

    M = G_F^{1/2} 2^{1/4} xi_H^W [ 7 lambda (mK^2+m_pi^2-mH^2)/18 - 7 A_K mK^2/9 ]
        + (xi_H^{ds}/(2v)) m_s (mK^2-m_pi^2)/(m_s-m_d) f_0^{K pi}(q^2)

    A_K = 0.  xi_H^W = c is the W-loop four-quark term; the FCNC term carries the
    complex xi_H^{ds}, so the two amplitudes *coherently* interfere through
    Re(xi_H^{ds}).  Returns (M, dict of components, p_star, mK).
    """
    mK, mPi, _gam = _k_meson(which)
    c = 1.0 / tanb
    fourq = (math.sqrt(GF) * 2.0 ** 0.25 * c
             * (7.0 * LAMBDA_QCD * (mK ** 2 + mPi ** 2 - mH ** 2) / 18.0
                - 7.0 * A_K * mK ** 2 / 9.0))
    xids = xi_H(mH, tanb, 'd', 's')
    fcnc_coef = (M_S_QUARK * (mK ** 2 - mPi ** 2) / (M_S_QUARK - M_D_QUARK) * F0_KPI
                 / (2.0 * V_VEV))
    fcnc = xids * fcnc_coef
    p_star = 0.0
    if mH < mK - mPi:
        p_star = math.sqrt(max((mK * mK - (mH + mPi) ** 2) * (mK * mK - (mH - mPi) ** 2), 0.0)) / (2.0 * mK)
    return {'M': fourq + fcnc, 'four_quark': fourq, 'fcnc': fcnc,
            'fcnc_abs': abs(fcnc), 'four_quark_abs': abs(fourq),
            'interference': 2.0 * (fourq.real * fcnc.real + fourq.imag * fcnc.imag),
            'p_star': p_star, 'mK': mK}


def br_K_to_piH(mH, tanb, which='K+'):
    """Absolute BR(K -> pi H).

    Br = (1/Gamma_K) (2 p*/mK) |M|^2/(16 pi mK).
    K+ -> pi+ H and K0L -> pi0 H use the same amplitude with the neutral masses (the
    paper gives no separate SU(2) coefficient for K0; the K+-normalised isospin
    relation, factor 1, is used - documented).  K0S -> pi0 H is CP suppressed -> 0.
    The returned K0L number is the per-K0L branching ratio: the response layer applies
    the 1/2 PDG-311 flavour projection itself.
    """
    if which in ('K0S', 'KS', 'K_S'):
        return 0.0
    a = k_to_pi_amplitude(mH, tanb, which)
    if a['p_star'] <= 0.0:
        return 0.0
    mK = a['mK']
    _, _, gam = _k_meson(which)
    return (1.0 / gam) * (2.0 * a['p_star'] / mK) * abs(a['M']) ** 2 / (16.0 * math.pi * mK)


def br_K_to_piH_components(mH, tanb, which='K+'):
    """Diagnostic breakdown of the K->pi H rate (for validation/reports)."""
    a = k_to_pi_amplitude(mH, tanb, which)
    br = br_K_to_piH(mH, tanb, which)
    return {'br': br, 'br_fcnc_only': None, 'br_fourq_only': None, 'components': a}


# ======================================================================================
# 3.4  semileptonic pi/K -> H e nu
# ======================================================================================
def br_pi_to_Henu(mH, tanb):
    """Br(pi+ -> H e+ nu) (charge conjugate included), paper formula after n_h,n_l."""
    mX = M_PI_P
    if mH >= mX:
        return 0.0
    x = (M_MU / mX) ** 2
    c = 1.0 / tanb
    pref = (math.sqrt(2.0) * GF * mX ** 4 * c * c
            / (96.0 * math.pi ** 2 * M_MU ** 2 * (1.0 - x) ** 2))
    return pref * BR_PI_MUNU * _f_phase_space((mH / mX) ** 2) * (7.0 / 9.0) ** 2


def br_K_to_Henu(mH, tanb, which='K+'):
    """Br(K+ -> H e+ nu) (charge conjugate included).  K0L -> H e nu is not a
    physical two-body leptonic channel (charge) -> 0.0."""
    if which in ('K0L', 'KL', 'K_L', 'K0S', 'KS', 'K_S'):
        return 0.0
    if which not in ('K+', 'Kp', 'K^+'):
        raise KeyError(which)
    mX = M_KP
    if mH >= mX:
        return 0.0
    x = (M_MU / mX) ** 2
    c = 1.0 / tanb
    pref = (math.sqrt(2.0) * GF * mX ** 4 * c * c
            / (96.0 * math.pi ** 2 * M_MU ** 2 * (1.0 - x) ** 2))
    return pref * BR_K_MUNU * _f_phase_space((mH / mX) ** 2) * (7.0 / 9.0) ** 2


# ======================================================================================
# 3.5  Higgs sector: h couplings, h -> HH, 2HDM Higgs production
# ======================================================================================
def _f_loop(tau):
    """Paper App. A f(tau) (complex for tau > 1)."""
    if tau <= 1.0:
        return math.asin(math.sqrt(max(tau, 0.0))) ** 2
    r = math.sqrt(1.0 - 1.0 / tau)
    return -0.25 * (math.log((1.0 + r) / (1.0 - r)) - 1j * math.pi) ** 2


def A12(tau):
    """A_{1/2}(tau) = [tau + (tau-1) f(tau)]/tau^2  (paper App. A)."""
    f = _f_loop(tau)
    return (tau + (tau - 1.0) * f) / tau ** 2


def A1(tau):
    """A_1(tau) = -1/2 [2 tau^2 + 3 tau + 3(2 tau - 1) f(tau)]/tau^2."""
    f = _f_loop(tau)
    return -0.5 * (2.0 * tau ** 2 + 3.0 * tau + 3.0 * (2.0 * tau - 1.0) * f) / tau ** 2


def g_hHH(mH, tanb):
    """g_C = s/(2v)[(mH^2-3 lambda v^2-mh^2)(2 cot2b s c - c^2 + s^2) + (lambda v^2-mH^2)].

    lambda v^2 = 0.  Paper App. C eq. ghHH.
    """
    c = 1.0 / tanb
    s = math.sqrt(max(1.0 - c * c, 0.0))
    cot2b = (1.0 - tanb ** 2) / (2.0 * tanb)
    bracket = ((mH * mH - 3.0 * LAMBDA_V2 - M_H_SM ** 2)
               * (2.0 * cot2b * s * c - c * c + s * s) + (LAMBDA_V2 - mH * mH))
    return s / (2.0 * V_VEV) * bracket


def kappa_g(mH, tanb):
    """kappa_g = [xi_h^t A12(t) + xi_h^b A12(b)] / [A12(t) + A12(b)]  (App. A A12).

    In Type-I all up- and down-type fermions share xi_h^f, so the ratio is exactly
    xi_h^f; the explicit loop ratio is kept (and unit-tested to equal xi_h^f to 1e-15)
    so the formula matches DESIGN.md 3.5 and generalises to other 2HDM types.
    """
    c = 1.0 / tanb
    s = math.sqrt(max(1.0 - c * c, 0.0))
    xi_f = s + c / tanb
    tau_t = M_H_SM ** 2 / (4.0 * M_QUARK['t'] ** 2)
    tau_b = M_H_SM ** 2 / (4.0 * M_B_QUARK ** 2)
    num = xi_f * A12(tau_t) + xi_f * A12(tau_b)
    den = A12(tau_t) + A12(tau_b)
    return abs(num) / abs(den)


def kappa_gamma(mH, tanb):
    """kappa_gamma from the full one-loop h gamma gamma amplitude (W + t,b,c,tau)."""
    c = 1.0 / tanb
    s = math.sqrt(max(1.0 - c * c, 0.0))
    xi_f = s + c / tanb        # xi_h^f (Type-I: same for all fermions)
    xi_V = s
    tauW = M_H_SM ** 2 / (4.0 * M_W ** 2)
    amps = [('W', 1.0, xi_V, A1(tauW), 1)]
    for (name, m, NcQ2) in (('t', M_QUARK['t'], 3.0 * (2.0 / 3.0) ** 2),
                            ('b', M_B_QUARK, 3.0 * (-1.0 / 3.0) ** 2),
                            ('c', M_QUARK['c'], 3.0 * (2.0 / 3.0) ** 2),
                            ('tau', 1.77686, 1.0)):
        amps.append((name, m, xi_f, A12(M_H_SM ** 2 / (4.0 * m ** 2)), NcQ2))
    num = 1j * 0.0
    den = 1j * 0.0
    for (_n, _m, kf, amp, nq) in amps:
        if _n == 'W':
            num += kf * amp
            den += amp
        else:
            num += nq * kf * amp
            den += nq * amp
    return abs(num) / abs(den)


def higgs_sector(mH, tanb):
    """h couplings, h -> HH branching ratio and 2HDM-scaled SM Higgs cross sections.

    Returns dict with keys: sigma_process, sigma_total, BR_hHH, g_C, g3,
    Gamma_h_other, Gamma_hHH, kappa, kappa_g, kappa_gamma, sigma_tHW_range.
    """
    c = 1.0 / tanb
    s = math.sqrt(max(1.0 - c * c, 0.0))
    xi_h_f = s + c / tanb
    xi_h_V = s

    gC = g_hHH(mH, tanb)
    g3 = -2.0 * gC
    gam_HH = 0.0
    if mH < 0.5 * M_H_SM:
        beta = math.sqrt(max(1.0 - 4.0 * mH * mH / M_H_SM ** 2, 0.0))
        gam_HH = abs(g3) ** 2 / (32.0 * math.pi * M_H_SM) * beta

    kg = kappa_g(mH, tanb)
    kgam = kappa_gamma(mH, tanb)
    # kappa_Zgamma: Z gamma is W-loop dominated; use kappa_V (documented approximation,
    # 0.15 % of the SM width -> negligible everywhere)
    kZgam = xi_h_V

    kap = {'bb': xi_h_f, 'WW': xi_h_V, 'gg': kg, 'tautau': xi_h_f, 'cc': xi_h_f,
           'ZZ': xi_h_V, 'gamma gamma': kgam, 'Zgamma': kZgam, 'ss': xi_h_f,
           'mumu': xi_h_f}
    gam_other = GAM_H_SM * sum(BR_H_SM[k] / BR_H_SM_NORM * kap[k] ** 2 for k in BR_H_SM)
    br_hHH = gam_HH / (gam_other + gam_HH) if gam_HH > 0.0 else 0.0

    # per-process kappa^2 (no global kappa: each process scales with its own coupling)
    sig = {}
    sig['ggF'] = SIGMA_H['ggF'] * kg ** 2
    for p in ('VBF', 'WH', 'ZH'):
        sig[p] = SIGMA_H[p] * xi_h_V ** 2
    for p in ('ttH', 'tHj', 'tHb'):
        sig[p] = SIGMA_H[p] * xi_h_f ** 2
    sig['bbH'] = SIGMA_H['bbH'] * xi_h_f ** 2
    # tHW: interference-aware bracket with cos(phi) = 1 baseline; phi in [0,pi] brackets
    kt, kV = xi_h_f, xi_h_V
    sig_tHW_base = SIGMA_H['tHW'] * (kt ** 2 + kV ** 2 + 2.0 * kt * kV) / 4.0
    sig_tHW_lo = SIGMA_H['tHW'] * (kt ** 2 + kV ** 2 - 2.0 * kt * kV) / 4.0
    sig_tHW_hi = sig_tHW_base
    sig['tHW'] = sig_tHW_base
    sigma_total = sum(sig.values())

    return {
        'sigma_process': sig,
        'sigma_total': sigma_total,
        'BR_hHH': br_hHH,
        'g_C': gC,
        'g3': g3,
        'Gamma_h_other': gam_other,
        'Gamma_hHH': gam_HH,
        'kappa': kap,
        'kappa_g': kg,
        'kappa_gamma': kgam,
        'kappa_V': xi_h_V,
        'kappa_f': xi_h_f,
        'sigma_tHW_range': (sig_tHW_lo, sig_tHW_hi),
        'c': c, 's': s,
    }


# ======================================================================================
# couplings + batch driver
# ======================================================================================
def couplings(tanb):
    """All tree-level 2HDM mixing/coupling modifiers at this tanb."""
    c = 1.0 / tanb
    s = math.sqrt(max(1.0 - c * c, 0.0))
    return {'tanb': tanb, 'c': c, 's': s,
            'xi_H_V': c, 'xi_H_f': c - s / tanb,
            'xi_H_f_stable': 1.0 / (tanb ** 3 * (1.0 + s)),
            'xi_h_V': s, 'xi_h_f': s + c / tanb}


def production_rates(mH, tanb, meson='B+'):
    """Production branching ratios used by the scan (B, K, pi) at one point."""
    b = br_B_to_XsH(mH, tanb, meson)
    kp = br_K_to_piH(mH, tanb, 'K+')
    kl = br_K_to_piH(mH, tanb, 'K0L')
    kenu = br_K_to_Henu(mH, tanb, 'K+')
    penu = br_pi_to_Henu(mH, tanb)
    return {'br_B_to_XsH': b['total'],
            'br_B_to_XsH_mode': b['mode'],
            'br_B_to_XsH_paper_mb': _b_ratio_and_total(mH, tanb, meson, 'paper_mb')[1],
            'br_B_to_XsH_exclusive': b['total_exclusive'],
            'br_B_to_XsH_channels': b['channels'],
            'br_K_to_piH_K+': kp,
            'br_K_to_piH_K0L': kl,
            'br_K_to_Henu_K+': kenu,
            'br_pi_to_Henu': penu,
            'br_B_ratio_to_BXcenu': b['ratio_to_BXcenu']}


def evaluate(mH, tanb, meson='B+'):
    """Everything a scan point needs from this module, in one flat dict."""
    cpl = couplings(tanb)
    rates = production_rates(mH, tanb, meson)
    hs = higgs_sector(mH, tanb)
    xib = xi_H(mH, tanb, 's', 'b')
    xids = xi_H(mH, tanb, 'd', 's')
    xidb = xi_H(mH, tanb, 'd', 'b')
    out = {'mH': float(mH), 'tanb': float(tanb)}
    out.update(cpl)
    out.update(rates)
    out['xi_H_bs'] = xib
    out['xi_H_ds'] = xids
    out['xi_H_db'] = xidb
    out['xi_H_bs_abs'] = abs(xib)
    out['xi_H_ds_abs'] = abs(xids)
    out['lambda_HH'] = lambda_HHpmHpm(mH, tanb)
    out['br_B_to_XsH_validity'] = _b_validity_flag(rates['br_B_to_XsH'])
    out['outside_validity_B'] = rates['br_B_to_XsH'] > 0.1
    out.update(hs)
    return out


def _b_validity_flag(br):
    """Model-validity flag: the inclusive b->s H rate must stay a small correction."""
    if br >= 1.0:
        return 'UNPHYSICAL_BR_GE_1'
    if br > 0.1:
        return 'SATURATED_BR_GT_0.1'
    return 'OK'


def evaluate_arrays(mH_list, tanb_list, meson='B+'):
    """Vectorised-over-tanb evaluation for a fixed mass list? Not needed; kept simple:
    returns {(mH,tanb): dict} for the full product grid (used by the scan cache)."""
    return {(mH, tb): evaluate(mH, tb, meson) for mH in mH_list for tb in tanb_list}


def evaluate_grid(points, meson='B+'):
    """Batch helper: points = iterable of (mH, tanb).  Returns a list of flat dicts."""
    out = []
    for p in points:
        if isinstance(p, dict):
            mH, tb = p['mH'], p['tanb']
        else:
            mH, tb = p[0], p[1]
        out.append(evaluate(mH, tb, meson))
    return out


def mass_grid():
    """The run's fixed mass grid (work/analysis/grids/mass_grid.txt)."""
    path = os.path.abspath(os.path.join(_HERE, '..', 'analysis', 'grids', 'mass_grid.txt'))
    with open(path, 'r', encoding='utf-8') as f:
        return [float(l.split()[0]) for l in f if l.strip() and not l.startswith('#')]


if __name__ == '__main__':
    for (mH, tb) in ((1.0, 100.0), (0.1, 100.0), (1.0, 1e5)):
        e = evaluate(mH, tb)
        print(f"mH={mH} tanb={tb:g}: xi_H^bs={e['xi_H_bs']:.6e} "
              f"BR(B->XsH)={e['br_B_to_XsH']:.6e} BR(K->piH)={e['br_K_to_piH_K+']:.4e} "
              f"BR(pi->He nu)={e['br_pi_to_Henu']:.4e} BR(h->HH)={e['BR_hHH']:.4e} "
              f"sigma_tot={e['sigma_total']:.4f} pb")
