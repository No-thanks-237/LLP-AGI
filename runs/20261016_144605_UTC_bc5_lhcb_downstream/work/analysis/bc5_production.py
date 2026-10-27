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
    """theta^2-normalised BRs of B -> X_s S channels (B0/Bs scaled by their lifetimes)."""
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
        br = xi ** 2 * mq ** 2 * M ** 2 / (32.0 * PI * V ** 2 * GX) * p / mX ** 2
        # scale from the B+ reference width to the actual meson width
        if GX != GAM_BP:
            br *= GAM_BP / GX
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
# Boiarska tab:b-production-quartic gives BR(m_S = 0) at alpha = 1 GeV and the "closing mass"
# at which the branching ratio has dropped by a factor of 10.  We interpolate with
#   BR(m_S) = BR(0) * alpha(m_S)^2 * 10^(-m_S/m_close)
# which is the paper's own description of the closing mass, and enforce the kinematic
# thresholds explicitly.  All 3-body modes are treated with a kaon spectator.
# --------------------------------------------------------------------------------------
M_BS_MESON = 5.36688
M_K_SPECT = 0.493677          # kaon spectator used for every 3-body quartic mode
QUARTIC = [
    # (meson, spectator mass or None for X -> SS, BR(0) at alpha=1, closing mass)
    ('Bs', None,      4.0e-10, 2.670),
    ('B+', M_K_SPECT, 1.4e-10, 1.998),
    ('B+', M_K_SPECT, 1.2e-10, 1.621),
    ('B+', M_K_SPECT, 1.2e-10, 1.478),
    ('B+', M_K_SPECT, 9.1e-11, 1.701),
    ('B+', M_K_SPECT, 3.5e-11, 1.621),
    ('B+', M_K_SPECT, 1.9e-11, 1.358),
    ('B+', M_K_SPECT, 2.5e-11, 1.499),
    ('B0', None,      1.1e-11, 2.624),
    ('B+', M_K_SPECT, 9.9e-12, 1.240),
    ('B+', M_K_SPECT, 4.7e-12, 2.149),
    ('B+', M_K_SPECT, 7.3e-13, 1.545),
]


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
    """List of (meson, kind, spectator_mass, BR) for X -> SS (kind '2b') or X -> X'SS ('3b').

    The returned BR already includes the alpha(m_S)^2 factor; each entry produces two scalars.
    """
    a2 = alpha_from_bh(ms) ** 2
    mX = {'B+': M_BP, 'B0': M_B0, 'Bs': M_BS_MESON}
    out = []
    for (mes, spec, br0, mclose) in QUARTIC:
        if spec is None:
            if ms * 2.0 >= mX[mes]:
                continue
            kind = '2b'
        else:
            if ms * 2.0 + spec >= mX[mes]:
                continue
            kind = '3b'
        br = br0 * a2 * 10.0 ** (-ms / mclose)
        out.append((mes, kind, spec, br))
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
