# -*- coding: utf-8 -*-
"""B -> X_s H exclusive-channel kinematic mixture for the 2HDM light-H run.

Provenance
----------
Channel list, hadronic matrix elements |M_XX'(q^2)| and meson constants are copied
verbatim (numerically) from inputs/bc5_reuse/bc5_production.py of this run, which in
turn implements Boiarska et al. arXiv:1904.10447 (Table tab:BR and App. C form
factors).  Only the b->s channels are kept (the b->d 'pi' entry is dropped), as
specified by DESIGN.md section 3.2: the exclusive channels describe the *recoil
system* X_s that the inclusive B->X_s H rate integrates over, so they are used ONLY
as a kinematic shape (relative weights renormalised to 1).  They are never added to
the inclusive rate.

Master formula (Boiarska eq:productionBRmesons, theta^2-normalised):
    BR(X_Qi -> X'_Qj S)/theta^2 = |xi_ij|^2 m_Qi^2 |M_XX'(m_S^2)|^2/(32 pi v^2 Gamma_X)
                                  * |p_S| / m_X^2
with the identification theta*xi_ij  <->  xi_H^{bs} of the FASER paper
(arXiv:2212.06186 Eq. eq:Lfceff_phi / eq:xi_H_2HDM), i.e. every extra factor is a
common |xi_H^{bs}|^2 which cancels in the normalised mixture.

Public API (used by twohdm_production.py and by the response layer):
    channel_weights(mH)          -> {label: weight}, sum = 1 (0 for closed channels)
    channel_masses()             -> {label: m_X' in GeV}
    shape_channels(mH, meson)    -> {label: BR/|xi|^2}  (Boiarska units)
    sample_channel(mH, rng)      -> label (weighted by channel_weights)
"""
from __future__ import annotations

import math

PI = math.pi
V_HIGGS = 246.0                    # GeV, same v as everywhere else in the run

# ---- meson / quark constants (identical to inputs/bc5_reuse/bc5_production.py) ----
M_BP, M_B0 = 5.27934, 5.27965
HBAR_GEV_S = 6.582119569e-25       # GeV s
GAM_BP = HBAR_GEV_S / 1.638e-12    # GeV
GAM_B0 = HBAR_GEV_S / 1.519e-12    # GeV
M_B_QUARK, M_S_QUARK, M_D_QUARK = 4.18, 0.0924, 0.00467
XI_SB = 3.6e-4                     # Boiarska Table 2 theta^2-normalised xi_sb

M_KST892, M_KST1410, M_KST1680 = 0.89555, 1.414, 1.718
M_K_NEUTRAL = 0.497611            # PDG K0: lightest s-carrying recoil of a NEUTRAL b -> s
M_K0ST700, M_K0ST1430 = 0.824, 1.425
M_K1_1270, M_K1_1400 = 1.253, 1.408
M_K2ST1430 = 1.4324

# Only b -> s channels (DESIGN.md 3.2: drop the b -> d 'pi' entry of BC5's list).
B_CHANNELS = [
    ('K',       0.493677,   XI_SB),
    ('K*892',   M_KST892,   XI_SB),
    ('K*1410',  M_KST1410,  XI_SB),
    ('K*1680',  M_KST1680,  XI_SB),
    ('K0*700',  M_K0ST700,  XI_SB),
    ('K0*1430', M_K0ST1430, XI_SB),
    ('K1_1270', M_K1_1270,  XI_SB),
    ('K1_1400', M_K1_1400,  XI_SB),
    ('K2*1430', M_K2ST1430, XI_SB),
]
B_CHANNEL_MASS = {n[0]: n[1] for n in B_CHANNELS}


def pS(mX, mXp, mS):
    """|p_S| in the X rest frame (GeV)."""
    if mX < mXp + mS:
        return 0.0
    return math.sqrt(max((mX * mX - (mXp + mS) ** 2) * (mX * mX - (mXp - mS) ** 2), 0.0)) / (2.0 * mX)


# ---- form factors (Boiarska Appendix C; copied from bc5_production.py) ----
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
    """|M_XX'(q^2 = m_S^2)| for the B -> X_s transition."""
    q2 = mS ** 2
    p = pS(mX, mXp, mS)
    if channel == 'K':
        return abs((mX ** 2 - mXp ** 2) / (M_S_QUARK - M_B_QUARK)) * ff_ps(q2, 0.33, 6.16)
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


def channel_masses(meson='B+'):
    """{channel: m_X' [GeV]}; for a NEUTRAL B the lightest recoil is K0, not K+-."""
    out = dict(B_CHANNEL_MASS)
    if meson == 'B0':
        out['K'] = M_K_NEUTRAL
    return out


def shape_channels(mS, meson='B+'):
    """BR/|xi_H^{bs}|^2 per exclusive channel (Boiarska master formula, one power of
    Gamma_X and the open-channel phase space included).  Sum over channels equals the
    exclusive B->X_s S rate divided by |xi|^2; it is *not* the inclusive rate."""
    if meson == 'B+':
        mX, GX = M_BP, GAM_BP
    elif meson == 'B0':
        mX, GX = M_B0, GAM_B0
    else:
        raise KeyError(meson)
    out = {}
    for (name, mXp, _xi) in B_CHANNELS:
        if name == 'K' and meson == 'B0':
            mXp = M_K_NEUTRAL      # B0 -> K0 H, not B0 -> K+- H
        if mS + mXp >= mX:
            out[name] = 0.0
            continue
        p = pS(mX, mXp, mS)
        M = matrix_element(name, mS, mX, mXp)
        out[name] = M_B_QUARK ** 2 * M ** 2 / (32.0 * PI * V_HIGGS ** 2 * GX) * p / mX ** 2
    return out


_SHAPE_CACHE = {}


def shape_channels_cached(mS, meson='B+'):
    """Cached shape_channels(mS, meson) (depends on mH and the B charge, not on tanb)."""
    key = (round(mS, 9), meson)
    if key not in _SHAPE_CACHE:
        _SHAPE_CACHE[key] = shape_channels(mS, meson)
    return _SHAPE_CACHE[key]


def channel_weights(mS, meson='B+'):
    """Normalised (sum = 1) kinematic mixture over the open b->s exclusive channels.

    Independent of tanb (all b->s channels carry the same xi_H^{bs}); depends on mH
    only through the phase space and the q^2-dependent form factors.  Channels that
    are kinematically closed get weight 0.  Returns all-zero weights if every channel
    is closed (mH >= m_B - m_K), in which case the caller has no B rate anyway.
    """
    raw = shape_channels_cached(mS, meson)
    tot = sum(raw.values())
    if tot <= 0.0:
        return {k: 0.0 for k in raw}
    return {k: v / tot for k, v in raw.items()}


def sample_channel(mS, rng, meson='B+'):
    """Draw one exclusive channel label with the normalised mixture (None if closed)."""
    w = channel_weights(mS, meson)
    labels = [k for k, v in w.items() if v > 0.0]
    if not labels:
        return None
    u = rng.random()
    acc = 0.0
    for k in labels:
        acc += w[k]
        if u <= acc:
            return k
    return labels[-1]
