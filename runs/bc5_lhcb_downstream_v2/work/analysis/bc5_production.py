# -*- coding: utf-8 -*-
"""
Meson -> S and Higgs -> S production for BC5, following Boiarska et al. arXiv:1904.10447.

Master formula (their Eq. productionBRmesons):
    BR(X -> X' S) = theta^2 * |xi_ij|^2 m_Qi^2 |M_XX'(m_S^2)|^2 / (32 pi v^2 Gamma_X) * |p_S|/m_X^2
with |p_S| the scalar momentum in the X rest frame.  We reproduce their tab:BR to a few
percent at m_S = 0 (see validate_production.py).
"""
from __future__ import annotations
import math, os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'widths')))
import scalar_widths as sw

V = sw.V_HIGGS
PI = math.pi

# meson masses / lifetimes
M_BP, M_B0, M_BS = 5.27934, 5.27965, 5.36688
M_KP, M_K0L = 0.493677, 0.497611
HBAR_GEV_S = 6.582119569e-25           # GeV s
GAM_BP = HBAR_GEV_S / 1.638e-12
GAM_B0 = HBAR_GEV_S / 1.519e-12
GAM_BS = HBAR_GEV_S / 1.520e-12
GAM_KP = HBAR_GEV_S / 1.2380e-8
GAM_KL = HBAR_GEV_S / 5.116e-8

XI_SB, XI_DS, XI_DB = 3.6e-4, 3.3e-6, 7.9e-5   # b->s, d->s, b->d FCNC loop factors
M_B_QUARK, M_S_QUARK, M_D_QUARK = 4.18, 0.0924, 0.00467

M_KST892, M_KST1410, M_KST1680 = 0.89555, 1.414, 1.718
M_K0ST700, M_K0ST1430 = 0.824, 1.425
M_K1_1270, M_K1_1400 = 1.253, 1.408
M_K2ST1430 = 1.4324


def pS(mX, mXp, mS):
    if mX < mXp + mS:
        return 0.0
    return math.sqrt(max((mX * mX - (mXp + mS) ** 2) * (mX * mX - (mXp - mS) ** 2), 0.0)) / (2.0 * mX)


# --------------------------------------------------------------------------------------
# form factors (Boiarska Appendix C)
# --------------------------------------------------------------------------------------
def ff_ps(q2, F0, m_fit):
    return F0 if m_fit is None else F0 / (1.0 - q2 / m_fit ** 2)


def ff_sca(q2, F0, a, b, mB=5.3):
    return F0 / (1.0 - a * q2 / mB ** 2 + b * (q2 / mB ** 2) ** 2)


def ff_vec892(q2):
    r1, r2, mR, mfit = 1.364, -0.99, 5.279, math.sqrt(36.8)
    return r1 / (1.0 - q2 / mR ** 2) + r2 / (1.0 - q2 / mfit ** 2)


def ff_vec_hi(q2, mV, xi_perp, xi_par, mB=5.3):
    xp = xi_perp / (1.0 - q2 / mB ** 2)
    xl = xi_par / (1.0 - q2 / mB ** 2)
    return (1.0 - 2.0 * mV ** 2 / (mB ** 2 + mV ** 2 - q2)) * xl + (mV / mB) * xp


def ff_psv(q2, mV, which, mB=5.3):
    F0A, aA, bA = 0.22, 2.4, 1.78
    F0B, aB, bB = -0.45, 1.34, 0.69
    th = math.radians(-34.0)
    mA, mBv = 1.31, 1.34
    VA = F0A / (1.0 - aA * q2 / mB ** 2 + bA * (q2 / mB ** 2) ** 2)
    VB = F0B / (1.0 - aB * q2 / mB ** 2 + bB * (q2 / mB ** 2) ** 2)
    if which == 1270:
        return (math.sin(th) * mA * VA + math.cos(th) * mBv * VB) / mV
    return (math.cos(th) * mA * VA - math.sin(th) * mBv * VB) / mV


def ff_ten(q2, mX, mT, F0=0.23, aT=1.23, bT=0.76):
    return F0 / ((1.0 - q2 / mX ** 2) * (1.0 - aT * q2 / mX ** 2 + bT * (q2 / mX ** 2) ** 2))


def matrix_element(channel, mS, mX, mXp):
    """|M_XX'(q^2 = m_S^2)|."""
    q2 = mS ** 2
    p = pS(mX, mXp, mS)
    if channel == 'K':
        return abs((mX ** 2 - mXp ** 2) / (M_S_QUARK - M_B_QUARK)) * ff_ps(q2, 0.33, 6.16)
    if channel == 'pi':
        return abs((mX ** 2 - mXp ** 2) / (M_D_QUARK - M_B_QUARK)) * ff_ps(q2, 0.258, 6.16)
    if channel == 'K*892':
        return abs(-2.0 * mX * p / (M_B_QUARK + M_S_QUARK)) * ff_vec892(q2)
    if channel == 'K*1410':
        return abs(-2.0 * mX * p / (M_B_QUARK + M_S_QUARK)) * ff_vec_hi(q2, 1.414, 0.28, 0.22)
    if channel == 'K*1680':
        return abs(-2.0 * mX * p / (M_B_QUARK + M_S_QUARK)) * ff_vec_hi(q2, 1.718, 0.24, 0.18)
    if channel == 'K0*700':
        return abs((mX ** 2 - mXp ** 2 - q2) / (M_B_QUARK + M_S_QUARK)) * ff_sca(q2, 0.46, 1.6, 1.35)
    if channel == 'K0*1430':
        return abs((mX ** 2 - mXp ** 2 - q2) / (M_B_QUARK + M_S_QUARK)) * ff_sca(q2, 0.17, 4.4, 6.4)
    if channel == 'K1_1270':
        return abs(2.0 * mX * p / (M_S_QUARK - M_B_QUARK)) * abs(ff_psv(q2, 1.253, 1270))
    if channel == 'K1_1400':
        return abs(2.0 * mX * p / (M_S_QUARK - M_B_QUARK)) * abs(ff_psv(q2, 1.408, 1400))
    if channel == 'K2*1430':
        if p <= 0:
            return 0.0
        return abs(-1.0 / (M_B_QUARK + M_S_QUARK) * math.sqrt(2.0 / 3.0) * (mX * p * p / 1.4324) * 2.0
                   * ff_ten(q2, mX, 1.4324))
    raise KeyError(channel)


# B -> X_s S channels: (label, X' mass, Gamma_X, m_Q, xi)
B_CHANNELS = [
    ('K',        M_KP,       GAM_BP, M_B_QUARK, XI_SB),
    ('pi',       0.13957039, GAM_BP, M_B_QUARK, XI_DB),   # b -> d transition
    ('K*892',    M_KST892,   GAM_BP, M_B_QUARK, XI_SB),
    ('K*1410',   M_KST1410,  GAM_BP, M_B_QUARK, XI_SB),
    ('K*1680',   M_KST1680,  GAM_BP, M_B_QUARK, XI_SB),
    ('K0*700',   M_K0ST700,  GAM_BP, M_B_QUARK, XI_SB),
    ('K0*1430',  M_K0ST1430, GAM_BP, M_B_QUARK, XI_SB),
    ('K1_1270',  M_K1_1270,  GAM_BP, M_B_QUARK, XI_SB),
    ('K1_1400',  M_K1_1400,  GAM_BP, M_B_QUARK, XI_SB),
    ('K2*1430',  M_K2ST1430, GAM_BP, M_B_QUARK, XI_SB),
]
# quark-level mass entering |M| for each channel (m_Qi of the decaying meson)
B_MQI = {'K': M_B_QUARK, 'pi': M_B_QUARK, 'K*892': M_B_QUARK, 'K*1410': M_B_QUARK,
         'K*1680': M_B_QUARK, 'K0*700': M_B_QUARK, 'K0*1430': M_B_QUARK,
         'K1_1270': M_B_QUARK, 'K1_1400': M_B_QUARK, 'K2*1430': M_B_QUARK}


def br_B_to_S(ms, meson='B+'):
    """theta^2-normalised BRs of B -> X_s S channels.

    The partial width is the same for B+, B0 and Bs (same quark transition), so the branching
    ratio is partial/Gamma_total and the only lifetime dependence is the GX in the denominator
    below.  An earlier version also multiplied by Gamma_B+/GX, which squared the lifetime
    factor and made the B0 and Bs mixing branching ratios 0.927 times too small.  Note that
    Gamma_B+/Gamma_B0 = 0.927 is exactly the 0.93 of the caption of tab:b-production-quartic,
    so that factor must be applied ONCE, and it is already contained in GX here.
    """
    if meson == 'B+':
        mX, GX = M_BP, GAM_BP
    elif meson == 'B0':
        mX, GX = M_B0, GAM_B0
    elif meson == 'Bs':
        mX, GX = M_BS, GAM_BS
    else:
        raise KeyError(meson)
    out = {}
    for (name, mXp, Gref, mq, xi) in B_CHANNELS:
        if ms + mXp >= mX:
            out[name] = 0.0
            continue
        p = pS(mX, mXp, ms)
        M = matrix_element(name, ms, mX, mXp)
        # BR = partial width / total width; GX already carries the whole lifetime dependence
        br = xi ** 2 * mq ** 2 * M ** 2 / (32.0 * PI * V ** 2 * GX) * p / mX ** 2
        out[name] = br
    return out


def br_K_to_S(ms, which='K+'):
    """theta^2-normalised BR(K -> pi S).  K_S^0 -> pi S is assumed to vanish (CP suppression)."""
    if which == 'K+':
        mX, GX, mP = M_KP, GAM_KP, 0.13957039
    elif which == 'K0L':
        mX, GX, mP = M_K0L, GAM_KL, 0.1349768
    else:
        return 0.0
    if ms + mP >= mX:
        return 0.0
    p = pS(mX, mP, ms)
    M = abs((mX ** 2 - mP ** 2) / (M_D_QUARK - M_S_QUARK)) * 0.96
    return XI_DS ** 2 * M_S_QUARK ** 2 * M ** 2 / (32.0 * PI * V ** 2 * GX) * p / mX ** 2


# --------------------------------------------------------------------------------------
# Higgs portal quartic coupling alpha fixed by BR(h -> SS) = 0.01
# --------------------------------------------------------------------------------------
B_H = 0.01
GAM_H_SM = 4.07e-3
M_H = 125.25


def alpha_from_bh(ms, b_h=B_H, gam_h_sm=GAM_H_SM):
    """alpha such that BR(h -> SS) = b_h, from Gamma(h->SS) = alpha^2 beta/(32 pi m_h)."""
    if ms >= M_H / 2.0:
        return 0.0
    beta = math.sqrt(1.0 - 4.0 * ms ** 2 / M_H ** 2)
    gam_ss = b_h / (1.0 - b_h) * gam_h_sm
    return math.sqrt(32.0 * PI * M_H * gam_ss / beta)


# --------------------------------------------------------------------------------------
# Quartic-coupling production: X -> SS and X -> X' SS  (theta-independent)
#
# Table 4 of Boiarska et al. (label tab:b-production-quartic) gives, at alpha = 1 GeV,
#   BR(m_S = 0)                -- column 2
#   the "closing mass" m_close -- column 3, the mass where BR has dropped by a factor 10
# and the decay width formulas are Eqs. (eq:bstossbranching) and (eq:br-quartic-3-body)
# in Appendix sec:quartic-coupling:
#
#   2-body:  Gamma(X -> SS)        ~ sqrt(1 - 4 m_S^2 / m_X^2)                      (exact)
#   3-body:  Gamma(X -> X' SS)     ~ int_{4m_S^2}^{(m_X-m_X')^2} dq^2
#                                        sqrt(E2*^2 - m_S^2) sqrt(E3*^2 - m_X'^2)
#                                    with E2* = sqrt(q^2)/2 ,
#                                         E3* = (m_X^2 - q^2 - m_X'^2)/(2 sqrt(q^2)).
#
# The previous code replaced BOTH shapes by an exponential, BR(0) 10^(-m_S/m_close), and used
# a single kaon mass as the spectator for every 3-body mode.  The exact 2-body form reproduces
# the paper's closing masses (0.4975 m_X = 2.670 for Bs, 2.627 for B0 against 2.624 in the
# table), which confirms that the exponential was the wrong shape.
#
# For the 3-body modes the q^2 dependence of |M_XX'|^2 (the form-factor models of the paper's
# appendices) is not tabulated here, so |M_XX'|^2 is taken constant.  The phase-space shape is
# then raised to a power fixed by the paper's own closing mass, so that the curve passes
# through BR(0) at m_S = 0 AND through 0.1 BR(0) at m_close, while keeping the exact threshold
# (m_X - m_X')/2 where the phase space vanishes.
#
# The paper's caption: "in the case of decays of B^0 mesons, all the given branching ratios
# should be multiplied by a factor 0.93".  Those B0 counterparts of the B+- channels were
# missing entirely.
# --------------------------------------------------------------------------------------
M_BS_MESON = 5.36688
M_B0_MESON = 5.27963
B0_SCALE = 0.93               # paper caption of tab:b-production-quartic
# Table 4: (meson, kind, final-state meson mass, BR(0) at alpha=1, closing mass [GeV])
QUARTIC = [
    ('Bs', '2b', 0.0,          4.0e-10, 2.670),
    ('B+', '3b', 0.493677,     1.4e-10, 1.998),   # B+ -> SS K+
    ('B+', '3b', 0.700,        1.2e-10, 1.621),   # B+ -> SS K*_0(700)
    ('B+', '3b', 1.270,        1.2e-10, 1.478),   # B+ -> SS K_1(1270)
    ('B+', '3b', 0.892,        9.1e-11, 1.701),   # B+ -> SS K*(892)
    ('B+', '3b', 1.430,        3.5e-11, 1.621),   # B+ -> SS K*_0(1430)
    ('B+', '3b', 1.410,        1.9e-11, 1.358),   # B+ -> SS K*(1410)
    ('B+', '3b', 1.430,        2.5e-11, 1.499),   # B+ -> SS K*_2(1430)
    ('B0', '2b', 0.0,          1.1e-11, 2.624),
    ('B+', '3b', 1.680,        9.9e-12, 1.240),   # B+ -> SS K*(1680)
    ('B+', '3b', 0.1395704,    4.7e-12, 2.149),   # B+ -> SS pi+
    ('B+', '3b', 1.400,        7.3e-13, 1.545),   # B+ -> SS K_1(1400)
]


def _ps2_ratio(ms, mX):
    """Exact two-body phase-space ratio sqrt(1 - 4 m_S^2/m_X^2), Eq. (eq:bstossbranching)."""
    if 2.0 * ms >= mX:
        return 0.0
    return math.sqrt(max(1.0 - 4.0 * ms * ms / (mX * mX), 0.0))


def _ps3_integral(ms, mX, mXp, n=120):
    """Phase-space integral of Eq. (eq:br-quartic-3-body) with |M_XX'|^2 taken constant."""
    lo = 4.0 * ms * ms
    hi = (mX - mXp) ** 2
    if hi <= lo:
        return 0.0
    dq = (hi - lo) / n
    tot = 0.0
    for i in range(n):
        q2 = lo + (i + 0.5) * dq
        sq = math.sqrt(q2)
        e2 = 0.5 * sq
        e3 = (mX * mX - q2 - mXp * mXp) / (2.0 * sq)
        a2 = e2 * e2 - ms * ms
        b2 = e3 * e3 - mXp * mXp
        if a2 > 0.0 and b2 > 0.0:
            tot += math.sqrt(a2) * math.sqrt(b2)
    return tot * dq


def _quartic_ratio(ms, kind, mX, mXp, mclose):
    """BR(m_S)/BR(0) for one quartic channel, anchored on the paper's BR(0) and closing mass."""
    if kind == '2b':
        return _ps2_ratio(ms, mX)
    i0 = _ps3_integral(0.0, mX, mXp)
    if i0 <= 0.0:
        return 0.0
    r_close = _ps3_integral(mclose, mX, mXp) / i0
    r = _ps3_integral(ms, mX, mXp) / i0
    if not (1e-6 < r_close < 0.999999) or r <= 0.0:
        return 10.0 ** (-ms / mclose)          # fall back to the old shape
    return r ** (math.log(0.1) / math.log(r_close))


# Kaon quartic modes, Boiarska tab:K->SS.  Same BR(m_S = 0) at alpha = 1 GeV and closing-mass
# prescription as the B modes above; K0S is 560 times smaller than K0L in the 2-body mode and is
# only listed for completeness (the minimum-bias sample carries K0L, not K0S).
M_K0_NEUTRAL = 0.497611
M_PI0 = 0.1349768
M_PI_CHARGED = 0.13957039
M_K_CHARGED = 0.493677
QUARTIC_KAON = [
    # (kaon, mode, BR(0) at alpha = 1, closing mass)
    ('K0L', 'SS',    4.4e-13, 0.252),
    ('K0L', 'SSpi0', 6.6e-15, 0.140),
    ('K+',  'SSpip', 1.4e-15, 0.136),
    ('K0S', 'SS',    7.8e-16, 0.252),
    ('K0S', 'SSpi0', 1.2e-17, 0.140),
]


def quartic_kaon_channels(ms, kaon=None):
    """(kaon, mode, BR) for the open K -> SS and K -> SS pi modes at this mass.

    The BR already includes alpha(m_S)^2 and the closing-mass factor; every entry produces two
    scalars.
    """
    a2 = alpha_from_bh(ms) ** 2
    mK = {'K0L': M_K0_NEUTRAL, 'K0S': M_K0_NEUTRAL, 'K+': M_K_CHARGED}  # noqa: F821
    out = []
    for (k, mode, br0, mclose) in QUARTIC_KAON:
        if kaon is not None and k != kaon:
            continue
        if mode == 'SS':
            if 2.0 * ms >= mK[k]:
                continue
        else:
            spec = M_PI0 if mode == 'SSpi0' else M_PI_CHARGED
            if 2.0 * ms + spec >= mK[k]:
                continue
        out.append((k, mode, br0 * a2 * 10.0 ** (-ms / mclose)))
    return out


def quartic_channels(ms):
    """List of (meson, kind, final_state_mass, BR) for X -> SS ('2b') or X -> X'SS ('3b').

    The BR already includes alpha(m_S)^2 and the exact mass dependence of the paper; every
    entry produces two scalars.  Each B+- three-body channel is accompanied by its B0
    counterpart at 0.93 of the branching ratio, as the caption of tab:b-production-quartic
    requires.
    """
    a2 = alpha_from_bh(ms) ** 2
    mX = {'B+': M_BP, 'B0': M_B0, 'Bs': M_BS_MESON}
    out = []
    for (mes, kind, mfin, br0, mclose) in QUARTIC:
        mx = mX[mes]
        if ms * 2.0 + mfin >= mx:
            continue
        r = _quartic_ratio(ms, kind, mx, mfin, mclose)
        if r <= 1e-30:
            continue
        out.append((mes, kind, mfin, br0 * a2 * r))
        if mes == 'B+' and kind == '3b':
            # B0 counterpart: same hadronic transition scaled by 0.93, own parent mass
            r0 = _quartic_ratio(ms, kind, mX['B0'], mfin, mclose)
            if r0 > 1e-30:
                out.append(('B0', kind, mfin, br0 * B0_SCALE * a2 * r0))
    return out


# 13 TeV SM Higgs production cross sections in pb (LHC HXSWG YR4 / PDG values)
SIGMA_H = {
    'ggF': 48.58, 'VBF': 3.782, 'WH': 1.369, 'ZH': 0.883,
    'bbH': 0.2384, 'ttH': 0.5065, 'tHj': 0.0740, 'tHW': 0.0152, 'tHb': 0.0046,
}
SIGMA_H_TOTAL = sum(SIGMA_H.values())

# minimum-bias inelastic pp cross section at 13 TeV (Pythia 8.316, SoftQCD:inelastic)
SIGMA_INEL_MB = 79.5

# cross-section normalisation references (to be filled from validated measurements)
REFS = {
    'sigma_inel_13TeV_mb': SIGMA_INEL_MB,
    'sigma_h_total_13TeV_pb': SIGMA_H_TOTAL,
    'sigma_h_sources': SIGMA_H,
}
