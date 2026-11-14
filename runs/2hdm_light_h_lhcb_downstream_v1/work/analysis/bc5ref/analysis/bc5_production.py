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

XI_DS, XI_DB = 3.3e-6, 7.9e-5   # d->s, b->d FCNC loop factors; xi_sb is set from
                             # Winkler's eq:effcoupling in the B -> X_s S block below
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


# ======================================================================================
# B -> X_s S : Boiarska's exclusive channels (arXiv:1904.10447 tab:BR)
# ======================================================================================
# The bc5 benchmark's production prescription is Boiarska's channel-by-channel table (prompt.md
# section 3.1: "按 Boiarska 文献逐道纳入 ... 排他道与对应包容率不重复相加"), and the reference tool
# SensCalc uses exactly that table by default (`SelectedBrChoice="Exclusive"` ->
# `BtoXsStotal.dat`; the free-quark inclusive formula is only the optional alternative).
#
# Master formula (their eq:productionBRmesons):
#   BR(X_{Q_i} -> X'_{Q_j} S) = theta^2 |xi_ij|^2 m_{Q_i}^2 |M_XX'(m_S^2)|^2 /(32 pi v^2 Gamma_X)
#                               * |p_S| / m_X^2
# with xi_ij from their Table 2 (eq:effcoupling): xi_ds = 3.3e-6, xi_db = 7.9e-5,
# xi_sb = 3.6e-4, and the |M_XX'|^2 form factors of their appendix C.
#
# B_s is NOT a mixing channel here: Boiarska state explicitly (main.tex L682) that the scalar
# production from D_s, J/psi, B_s, B_c and Upsilon has no advantage over D^0, D^+, B^-, B^0
# "while their amount at any experiment is significantly lower. Therefore, we will discuss below
# only production from [the] latter mesons."  B_s therefore contributes only through the quartic
# B_s -> S S entry of tab:b-production-quartic (see QUARTIC below); SensCalc does the same
# (its "B-mixing" term uses fbtob0 + fbtobplus only, with a separate "Bs-quartic" term).
XI_SB = 3.6e-4          # Boiarska Table 2 (their eq:effcoupling, evaluated with their inputs)

# (label, X' mass, xi).  The parent lifetime is not a per-channel constant: B+ and B0 differ
# only through Gamma_X (Gamma_B0 = 1/0.93 x Gamma_B+), which br_B_to_S applies once.
B_CHANNELS = [
    ('K',        M_KP,       XI_SB),
    ('pi',       0.13957039, XI_DB),      # b -> d transition
    ('K*892',    M_KST892,   XI_SB),
    ('K*1410',   M_KST1410,  XI_SB),
    ('K*1680',   M_KST1680,  XI_SB),
    ('K0*700',   M_K0ST700,  XI_SB),
    ('K0*1430',  M_K0ST1430, XI_SB),
    ('K1_1270',  M_K1_1270,  XI_SB),
    ('K1_1400',  M_K1_1400,  XI_SB),
    ('K2*1430',  M_K2ST1430, XI_SB),
]
# recoil mass used by the two-body generator of run_scan.sample_meson_decay_S
B_CHANNEL_MASS = {n[0]: n[1] for n in B_CHANNELS}


def br_B_to_S(ms, meson='B+'):
    """theta^2-normalised BRs of B -> X_s S exclusive channels (Boiarska tab:BR).

    The partial width is the same for B+ and B0 (same quark transition); the branching ratio
    differs only through the total width Gamma_X in the denominator, which already contains the
    Gamma_B+/Gamma_B0 = 0.93 of the table caption exactly once.

    B_s returns {} by construction, see the block comment above.
    """
    if meson == 'Bs':
        return {}
    if meson == 'B+':
        mX, GX = M_BP, GAM_BP
    elif meson == 'B0':
        mX, GX = M_B0, GAM_B0
    else:
        raise KeyError(meson)
    out = {}
    for (name, mXp, xi) in B_CHANNELS:
        if ms + mXp >= mX:
            out[name] = 0.0
            continue
        p = pS(mX, mXp, ms)
        M = matrix_element(name, ms, mX, mXp)
        out[name] = xi ** 2 * M_B_QUARK ** 2 * M ** 2 / (32.0 * PI * V ** 2 * GX) * p / mX ** 2
    return out


# The free-quark (inclusive) estimate of Boiarska eq:inclusive / Winkler eq:inclusive is NOT
# used: their own caption and table are the exclusive ones, and prompt.md forbids adding an
# inclusive estimate on top of the exclusive sum ("排他道与对应包容率不重复相加").


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
# Table 4: (meson, kind, final-state meson mass, BR(0) at alpha=1, closing mass [GeV], channel)
# The last column names the hadronic matrix element M_XX'(q^2) that enters the exact three-body
# integral; it is the same form-factor model as the corresponding mixing channel of App. C.  The
# final-state masses are the effective values that reproduce the paper's BR(0) and closing masses
# (validated to a few percent, see logs/production_audit_v3.md F11/F12) and are deliberately NOT
# replaced by the PDG resonance masses of the matrix_element() channels.
QUARTIC = [
    ('Bs', '2b', 0.0,          4.0e-10, 2.670, 'Bs'),
    ('B+', '3b', 0.493677,     1.4e-10, 1.998, 'K'),        # B+ -> SS K+
    ('B+', '3b', 0.700,        1.2e-10, 1.621, 'K0*700'),   # B+ -> SS K*_0(700)
    ('B+', '3b', 1.270,        1.2e-10, 1.478, 'K1_1270'),  # B+ -> SS K_1(1270)
    ('B+', '3b', 0.892,        9.1e-11, 1.701, 'K*892'),    # B+ -> SS K*(892)
    ('B+', '3b', 1.430,        3.5e-11, 1.621, 'K0*1430'),  # B+ -> SS K*_0(1430)
    ('B+', '3b', 1.410,        1.9e-11, 1.358, 'K*1410'),   # B+ -> SS K*(1410)
    ('B+', '3b', 1.430,        2.5e-11, 1.499, 'K2*1430'),  # B+ -> SS K*_2(1430)
    ('B0', '2b', 0.0,          1.1e-11, 2.624, 'B0'),
    ('B+', '3b', 1.680,        9.9e-12, 1.240, 'K*1680'),   # B+ -> SS K*(1680)
    ('B+', '3b', 0.1395704,    4.7e-12, 2.149, 'pi'),       # B+ -> SS pi+
    ('B+', '3b', 1.400,        7.3e-13, 1.545, 'K1_1400'),  # B+ -> SS K_1(1400)
]


def _ps2_ratio(ms, mX):
    """Exact two-body phase-space ratio sqrt(1 - 4 m_S^2/m_X^2), Eq. (eq:bstossbranching)."""
    if 2.0 * ms >= mX:
        return 0.0
    return math.sqrt(max(1.0 - 4.0 * ms * ms / (mX * mX), 0.0))


def _ps3_integral_exact(ms, mX, mXp, channel=None, n=400):
    """The paper's three-body phase-space integral with the real q^2 dependence.

        I(m_S) = int_{4 m_S^2}^{(m_X-m_X')^2} |M_XX'(q^2)|^2
                 sqrt(E2*^2 - m_S^2) sqrt(E3*^2 - m_X'^2) dq^2 ,
        E2* = sqrt(q^2)/2 ,  E3* = (m_X^2 - q^2 - m_X'^2)/(2 sqrt(q^2))

    (App. D.1, Eq. eq:br-quartic-3-body).  channel=None means |M_XX'|^2 = 1, which is exact for
    the kaon modes because their f0 form factor 0.96 carries no q^2 dependence.  The substitution
    q^2 = lo + (hi-lo) sin^2(pi t/2) removes both square-root endpoints.
    """
    lo = 4.0 * ms * ms
    hi = (mX - mXp) ** 2
    if hi <= lo:
        return 0.0
    tot = 0.0
    for i in range(n):
        t = (i + 0.5) / n
        s2 = math.sin(0.5 * math.pi * t) ** 2
        q2 = lo + (hi - lo) * s2
        dq2 = (hi - lo) * 0.5 * math.pi * math.sin(math.pi * t) / n
        sq = math.sqrt(q2)
        e2 = 0.5 * sq
        e3 = (mX * mX - q2 - mXp * mXp) / (2.0 * sq)
        a2 = e2 * e2 - ms * ms
        b2 = e3 * e3 - mXp * mXp
        if a2 > 0.0 and b2 > 0.0:
            if channel is None:
                mod2 = 1.0
            else:
                m = matrix_element(channel, sq, mX, mXp)
                mod2 = m * m
            tot += mod2 * math.sqrt(a2) * math.sqrt(b2) * dq2
    return tot


_I0_CACHE = {}


def _ps3_ratio(ms, channel, mX, mXp):
    """I(m_S)/I(0) for one three-body channel; I(0) depends only on the channel and is cached."""
    key = (channel, round(mX, 9), round(mXp, 9))
    i0 = _I0_CACHE.get(key)
    if i0 is None:
        i0 = _ps3_integral_exact(0.0, mX, mXp, channel)
        _I0_CACHE[key] = i0
    if i0 <= 0.0:
        return 0.0
    return _ps3_integral_exact(ms, mX, mXp, channel) / i0


def _quartic_ratio(ms, kind, mX, mXp, mclose=None, channel=None):
    """BR(m_S)/BR(0) for one quartic channel from the paper's exact phase space.

    kind '2b': Eq. (eq:bstossbranching), sqrt(1 - 4 m_S^2/m_X^2), already exact.
    kind '3b': Eq. (eq:br-quartic-3-body) with the channel's own |M_XX'(q^2)|^2.  The earlier
    version replaced this by a constant-|M|^2 integral raised to a power fitted through the
    paper's closing mass; that heuristic was up to 2x off for K1(1400)/K2*(1430) individually
    (see logs/production_audit_v3.md section 4.1).  The absolute normalisation still comes from
    tab:b-production-quartic BR(0), which is validated to 3-4 %.
    """
    if kind == '2b':
        return _ps2_ratio(ms, mX)
    return _ps3_ratio(ms, channel, mX, mXp)


# --------------------------------------------------------------------------------------
# differential version of the three-body quartic rate (F07)
# --------------------------------------------------------------------------------------
_Q2_ENV = {}


def _q2_weight(ms, mX, mXp, channel, q2):
    """Integrand of eq:br-quartic-3-body: |M_XX'(q^2)|^2 sqrt(E2*^2-m_S^2) sqrt(E3*^2-m_X'^2)."""
    sq = math.sqrt(q2)
    e2 = 0.5 * sq
    e3 = (mX * mX - q2 - mXp * mXp) / (2.0 * sq)
    a2 = e2 * e2 - ms * ms
    b2 = e3 * e3 - mXp * mXp
    if a2 <= 0.0 or b2 <= 0.0:
        return 0.0
    if channel is None:
        mod2 = 1.0
    else:
        mm = matrix_element(channel, sq, mX, mXp)
        mod2 = mm * mm
    return mod2 * math.sqrt(a2) * math.sqrt(b2)


def sample_dipion_mass(ms, mX, mXp, channel, rng, n_grid=200):
    """Draw m_SS = sqrt(q^2) for X -> X' S S with the density of the rate integral.

    The rate uses the exact integral of |M_XX'(q^2)|^2 over the three-body phase space
    (bc5_production._ps3_integral_exact) but the generator used to draw an unweighted nbody()
    configuration, i.e. a different distribution: |M|^2 varies by factors 2-9 over q^2.  Drawing
    q^2 with this rejection makes the generated events the differential version of the same
    rate.  The angles are isotropic in the S S rest frame, exactly as the rate integral assumes.

    Returns sqrt(q^2), or None if the mode is closed.
    """
    lo = 4.0 * ms * ms
    hi = (mX - mXp) ** 2
    if hi <= lo:
        return None
    key = (round(ms, 6), round(mX, 6), round(mXp, 6), channel)
    wmax = _Q2_ENV.get(key)
    if wmax is None:
        wmax = 0.0
        for k in range(1, n_grid):
            q2 = lo + (hi - lo) * k / n_grid
            wmax = max(wmax, _q2_weight(ms, mX, mXp, channel, q2))
        wmax *= 1.02
        _Q2_ENV[key] = wmax
    if wmax <= 0.0:
        return None
    for _ in range(200):
        q2 = lo + (hi - lo) * rng.random()
        if _q2_weight(ms, mX, mXp, channel, q2) >= rng.random() * wmax:
            return math.sqrt(q2)
    return None


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


_QK_CACHE = {}
_QB_CACHE = {}


def quartic_kaon_channels(ms, kaon=None):
    """(kaon, mode, BR) for the open K -> SS and K -> SS pi modes at this mass.

    The BR already includes alpha(m_S)^2 and the exact phase-space shape of the paper:
    Eq. (eq:bstossbranching) for K -> SS, and Eq. (eq:br-quartic-3-body) with a constant
    |M|^2 = 1 for K -> SS pi (exact there, because the kaon f0 form factor 0.96 has no q^2
    dependence).  The previous 10^(-m_S/m_close) exponential was up to 3.7x too steep at
    0.2 GeV; the K channel is <= 0.08 % of the signal, so this is a correctness fix only.
    Every entry produces two scalars.
    """
    key = (round(ms, 9), kaon)
    if key in _QK_CACHE:
        return _QK_CACHE[key]
    a2 = alpha_from_bh(ms) ** 2
    mK = {'K0L': M_K0_NEUTRAL, 'K0S': M_K0_NEUTRAL, 'K+': M_K_CHARGED}
    out = []
    for (k, mode, br0, mclose) in QUARTIC_KAON:
        if kaon is not None and k != kaon:
            continue
        if mode == 'SS':
            if 2.0 * ms >= mK[k]:
                continue
            r = _ps2_ratio(ms, mK[k])
        else:
            spec = M_PI0 if mode == 'SSpi0' else M_PI_CHARGED
            if 2.0 * ms + spec >= mK[k]:
                continue
            r = _ps3_ratio(ms, None, mK[k], spec)
        if r <= 0.0:
            continue
        out.append((k, mode, br0 * a2 * r))
    _QK_CACHE[key] = out
    return out


def quartic_channels(ms):
    """List of (meson, kind, final_state_mass, BR) for X -> SS ('2b') or X -> X'SS ('3b').

    The BR already includes alpha(m_S)^2 and the exact mass dependence of the paper; every
    entry produces two scalars.  Each B+- three-body channel is accompanied by its B0
    counterpart at 0.93 of the branching ratio, as the caption of tab:b-production-quartic
    requires.
    """
    key = round(ms, 9)
    if key in _QB_CACHE:
        return _QB_CACHE[key]
    a2 = alpha_from_bh(ms) ** 2
    mX = {'B+': M_BP, 'B0': M_B0, 'Bs': M_BS_MESON}
    out = []
    for (mes, kind, mfin, br0, mclose, chan) in QUARTIC:
        mx = mX[mes]
        if ms * 2.0 + mfin >= mx:
            continue
        r = _quartic_ratio(ms, kind, mx, mfin, mclose, chan)
        if r <= 1e-30:
            continue
        # The channel label (the hadronic matrix element |M_XX'(q^2)|^2 that the rate integral
        # used) is part of the tuple so that the KINEMATIC generator can sample q^2 with the same
        # weight (F07: the rate carried |M(q^2)|^2 while the generator used plain nbody).
        out.append((mes, kind, mfin, br0 * a2 * r, chan))
        if mes == 'B+' and kind == '3b':
            # B0 counterpart: same hadronic transition scaled by 0.93, own parent mass
            r0 = _quartic_ratio(ms, kind, mX['B0'], mfin, mclose, chan)
            if r0 > 1e-30:
                out.append(('B0', kind, mfin, br0 * B0_SCALE * a2 * r0, chan))
    _QB_CACHE[key] = out
    return out


# 13 TeV SM Higgs production cross sections in pb (LHC HXSWG YR4 / PDG values)
SIGMA_H = {
    'ggF': 48.58, 'VBF': 3.782, 'WH': 1.369, 'ZH': 0.883,
    'bbH': 0.2384, 'ttH': 0.5065, 'tHj': 0.0740, 'tHW': 0.0152, 'tHb': 0.0046,
}
SIGMA_H_TOTAL = sum(SIGMA_H.values())

# Pythia's "HiggsSM:all = on" does not produce the 13 TeV process mixture that SIGMA_H above
# describes: ggF comes out ~80 % of the events instead of 88 %, and the fermion-annihilation
# channel (code 901) is ~2.7 % instead of the 0.43 % bbH it stands in for.  Multiplying the whole
# sample by SIGMA_H_TOTAL therefore reproduces the total rate but not the per-process weights;
# run_scan.py reweights each group by (target fraction)/(sample fraction).  Codes follow the
# Pythia 8 Higgs-process listing; 903 (8 events) is counted as ggF.
HIGGS_PROC_GROUP = {902: 'ggF', 903: 'ggF', 906: 'VBF', 907: 'VBF', 905: 'WH',
                    904: 'ZH', 901: 'bbH', 908: 'ttH', 909: 'ttH'}

# minimum-bias inelastic pp cross section at 13 TeV (Pythia 8.316, SoftQCD:inelastic)
SIGMA_INEL_MB = 79.5

# cross-section normalisation references (to be filled from validated measurements)
REFS = {
    'sigma_inel_13TeV_mb': SIGMA_INEL_MB,
    'sigma_h_total_13TeV_pb': SIGMA_H_TOTAL,
    'sigma_h_sources': SIGMA_H,
}