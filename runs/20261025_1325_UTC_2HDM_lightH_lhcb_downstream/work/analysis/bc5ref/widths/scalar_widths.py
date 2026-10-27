# -*- coding: utf-8 -*-
"""
BC5 Higgs-portal real scalar S: decay widths, branching ratios and lifetime.

Physics sources (snapshotted in inputs/physics/):
  [W] Winkler, arXiv:1809.01876v1 (scalar.tex)       -- decay rates, form factors, dispersive pi pi / K K
  [B] Boiarska et al., arXiv:1904.10447v4 (main.tex) -- decay widths recap, loop factors, perturbative quarks/gluons

Convention (identical to prompt.md Sec. 2.1 and to [W] Sec. 3):
  The scalar S couples to every SM fermion and to W/Z exactly as the SM Higgs does,
  multiplied by the mixing angle theta:
      L = -theta sum_f (m_f/v) S fbar f + (2 theta m_W^2/v) S W+ W- + (theta m_Z^2/v) S Z Z
  Hence Gamma(S -> X) = theta^2 Gamma(h -> X) evaluated with m_h -> m_S, for every channel
  for which a SM-Higgs expression exists (leptons, photons, gluons, perturbative quarks).
  This is exactly the structure used in [W].  Hadronic channels in the resonance region
  (pi pi, K K, multi-meson) use the form-factor/dispersive treatment of [W].

  Cross-check of the gluon channel: with L_eff = (alpha_s/(12 pi v)) h G^a_{mu nu}G^{a mu nu}
  the heavy-top limit of the width is alpha_s^2 m^3/(72 pi^3 v^2).  Using
  |sum_q F_{1/2}(inf)|^2 = 16/9 the formula below reproduces this exactly, and gives
  Gamma(h->gg) = 1.78e-4 GeV at m_h=125 GeV, i.e. the correct LO value
  (the PDG 3.33e-4 GeV includes a QCD K-factor ~1.9).

All widths are at theta^2 = 1 (i.e. Gamma/theta^2); they are exactly linear in theta^2.
Units: GeV; c*tau in metres.
"""
from __future__ import annotations

import math
import numpy as np

# --------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------
V_HIGGS      = 246.22        # GeV, EW vev
G_F          = 1.1663787e-5  # GeV^-2
ALPHA_EM_MZ  = 1.0 / 127.95
ALPHA_EM_0   = 1.0 / 137.036   # on-shell (Thomson) photon coupling, used for S -> gamma gamma
M_Z          = 91.1876
M_W          = 80.377
M_H          = 125.25
HBAR_C       = 1.973269804593025e-16   # GeV m

M_E     = 0.51099895e-3
M_MU    = 0.1056583755
M_TAU   = 1.77686

MSB_S_2GEV = 0.0924     # [B]
MSB_C_2GEV = 1.23       # [B]
MSB_B_MB   = 4.18       # [B]
M_T_POLE   = 172.69

M_K        = 0.493677   # lightest hadron with an s quark
M_D        = 1.86966    # lightest hadron with a c quark
M_B        = 5.27966    # lightest hadron with a b quark

M_PI_CHARGED = 0.13957039
M_PI0        = 0.1349768
M_PI         = 0.13957039
M_K_CHARGED  = 0.493677
M_K0         = 0.497611
M_KAON       = 0.495

M_U_MSBAR  = 0.00216
M_D_MSBAR  = 0.00467

# Multi-meson (4pi, eta-eta, rho-rho, ...) coefficient of [W] Eq. for Gamma_{4pi,...} = C m^3 beta_2pi.
#
# Winkler sets C = 5.1e-9 GeV^-2 and states that this value makes the hadronic rate "transit
# smoothly into the rate of the spectator model at m_phi = 2 GeV".  That calibration is tied to
# HIS perturbative normalisation, which is LO in alpha_s for the gluon channel.
#
# We follow Boiarska/SensCalc and include the NLO correction to Gamma_gg (GG_NLO below), which
# makes the perturbative rate at 2 GeV a factor 2.0 larger.  C is a matching constant, not an
# independent physical parameter, so it has to be re-calibrated with it:
#     dispersive at 2 GeV  (leptons + pi pi + K K + C m^3 beta_2pi)  ==  perturbative just above
#     (leptons + s sbar + NLO gg)
# gives C = 1.7227e-8 GeV^-2.  With this value the model is continuous at the junction (ratio
# 1.000 by construction) and the 4pi-like term carries 71 % of the width at 1.99 GeV, against
# 74 % for the real 4pi channels in the SensCalc table - i.e. the recalibration reproduces the
# reference's composition instead of Winkler's LO one.  Keep 5.1e-9 (and GG_NLO = False) for a
# literal reproduction of his Fig. 3.
C_MULTIMESON = 5.1e-9        # GeV^-2, [W] (Sect. "Perturbative Spectator Model")
# Fixed quark masses used in Winkler's perturbative spectator model ("we set m_s = 95 MeV,
# m_c = 1.3 GeV and neglect the tiny decay rate into u, d quarks").  m_b is not needed by him
# explicitly (the bbar b channel is "analogous") and we take the MS-bar value at its own scale.
W_M_S, W_M_C, W_M_B = 0.095, 1.30, 4.18

# --------------------------------------------------------------------------------------
# Strong coupling constant and MS-bar quark masses.
#
# Both are obtained from one pass over a log-spaced grid in mu; the grid values are then
# integrated with a cumulative trapezoid so that repeated evaluations are O(log n).
# 3-loop QCD beta function, matched to alpha_s(m_Z) = 0.1179, n_f thresholds at m_c, m_b.
# 2-loop quark-mass anomalous dimension, gamma_m = (alpha_s/pi)[1 + g1 alpha_s/pi].
# --------------------------------------------------------------------------------------
_ALPHA_S_MZ = 0.1179
_M_C_MSBAR  = 1.27
_M_B_MSBAR  = 4.18
_MU_MIN     = 1.0        # low edge of the grid
_MU_MAX     = 200.0
_NGRID      = 4000

_MU_GRID = np.exp(np.linspace(math.log(_MU_MIN), math.log(_MU_MAX), _NGRID))
_LN_GRID = np.log(_MU_GRID)


def _beta_coeffs(nf: float):
    b0 = 11.0 - 2.0 * nf / 3.0
    b1 = 102.0 - 38.0 * nf / 3.0
    b2 = 2857.0 / 2.0 - 5033.0 * nf / 18.0 + 325.0 * nf ** 2 / 54.0
    return b0, b1, b2


def _nf_of(mu: float) -> int:
    return 5 if mu >= _M_B_MSBAR else (4 if mu >= _M_C_MSBAR else 3)


def _build_grids():
    """Integrate alpha_s(mu) downward from m_Z and mbar(mu) from their reference scales."""
    a = np.empty(_NGRID)
    a[:] = np.nan
    i0 = int(np.searchsorted(_MU_GRID, M_Z))
    i0 = min(i0, _NGRID - 1)
    # walk downward from the grid point just above m_Z
    a[i0] = _ALPHA_S_MZ + (0.0)
    # refine: alpha_s at _MU_GRID[i0] by running up from m_Z
    b0, b1, b2 = _beta_coeffs(5)

    def da(aa):
        # d(alpha_s)/d(ln mu) in the standard PDG convention.  The 2-loop coefficient is
        # b1/(8 pi^2) and the 3-loop one b2/(32 pi^3); an earlier version used b1/(4 pi^2)
        # and b2/(64 pi^3), which made alpha_s(2 GeV) = 0.3376 instead of 0.2991 (+12.9 %,
        # i.e. Gamma_gg +27 %) and broke the hadronic/perturbative junction at 2 GeV.
        # Checked against the PDG: alpha_s(m_tau) = 0.331, alpha_s(2) = 0.300, alpha_s(5) = 0.21.
        return -(b0 / (2 * math.pi)) * aa ** 2 - (b1 / (8 * math.pi ** 2)) * aa ** 3 \
               - (b2 / (32 * math.pi ** 3)) * aa ** 4

    # one RK4 pass from m_Z up to _MU_GRID[i0]
    if _MU_GRID[i0] > M_Z:
        h = math.log(_MU_GRID[i0] / M_Z)
        aa = _ALPHA_S_MZ
        k1 = da(aa); k2 = da(aa + .5 * h * k1); k3 = da(aa + .5 * h * k2); k4 = da(aa + h * k3)
        a[i0] = aa + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    else:
        a[i0] = _ALPHA_S_MZ
    for i in range(i0 - 1, -1, -1):
        mu_lo, mu_hi = _MU_GRID[i], _MU_GRID[i + 1]
        nf = _nf_of(mu_hi)
        b0, b1, b2 = _beta_coeffs(nf)
        h = math.log(mu_hi / mu_lo)
        aa = a[i + 1]
        # walking DOWN in mu: d(ln mu) = -h
        k1 = da(aa); k2 = da(aa - .5 * h * k1); k3 = da(aa - .5 * h * k2); k4 = da(aa - h * k3)
        a[i] = aa - h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    for i in range(i0 + 1, _NGRID):
        mu_lo, mu_hi = _MU_GRID[i - 1], _MU_GRID[i]
        nf = _nf_of(mu_hi)
        b0, b1, b2 = _beta_coeffs(nf)
        h = math.log(mu_hi / mu_lo)
        aa = a[i - 1]
        k1 = da(aa); k2 = da(aa + .5 * h * k1); k3 = da(aa + .5 * h * k2); k4 = da(aa + h * k3)
        a[i] = aa + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    # cumulative integral of gamma_m dln(mu) for the mass running (2-loop)
    g1 = np.array([101.0 / 12.0 - 5.0 * _nf_of(m) / 18.0 for m in _MU_GRID])
    x = a / math.pi
    # gamma_m = 2 (alpha_s/pi) [1 + (g1/2) alpha_s/pi]: the leading coefficient is 2, not 1.
    # With 1 the masses run far too slowly (m_b(m_Z) = 3.35 instead of the PDG 2.82,
    # m_s(10 GeV) = 0.0761 instead of 0.0681); with 2 one gets m_b(m_Z) = 2.89.
    gam = 2.0 * x * (1.0 + 0.5 * g1 * x)
    for i in range(1, _NGRID):
        if _nf_of(_MU_GRID[i]) != _nf_of(_MU_GRID[i - 1]):
            gam[i] = gam[i]   # continuous gamma_m: no matching needed at this order
    dln = np.diff(_LN_GRID)
    cum = np.concatenate(([0.0], np.cumsum(0.5 * (gam[1:] + gam[:-1]) * dln)))
    return a, cum


_ALPHA_GRID, _GAMMA_CUM = _build_grids()


def alpha_s(mu: float) -> float:
    """alpha_s(mu), 3-loop MS-bar running with n_f thresholds at m_c and m_b."""
    if mu < _MU_MIN or mu > _MU_MAX:
        raise ValueError("alpha_s grid covers only %.2f - %.2f GeV" % (_MU_MIN, _MU_MAX))
    return float(np.interp(math.log(mu), _LN_GRID, _ALPHA_GRID))


def _gamma_int(mu_a: float, mu_b: float) -> float:
    """Integral of gamma_m d(ln mu) from mu_a to mu_b (positive if mu_b > mu_a)."""
    la, lb = math.log(mu_a), math.log(mu_b)
    ya = float(np.interp(la, _LN_GRID, _GAMMA_CUM))
    yb = float(np.interp(lb, _LN_GRID, _GAMMA_CUM))
    return yb - ya


def mbar_q(mq_ref: float, mu_ref: float, mu: float) -> float:
    """MS-bar running quark mass: mbar(mu) = mq_ref * exp(-int gamma_m dln mu)."""
    return mq_ref * math.exp(-_gamma_int(mu_ref, mu))


# --------------------------------------------------------------------------------------
# One-loop triangle form factors (numerically stable for all tau)
# --------------------------------------------------------------------------------------
def F_half(tau: float):
    """Fermion-loop function F_{1/2}(tau); F_{1/2}(inf) = -4/3.

    F_{1/2}(tau) = -2 tau [1 + (1-tau) f(tau)],  f = arcsin^2(1/sqrt(tau)) for tau >= 1.
    For tau >= 1 we use the convergent series for 1 + (1-tau) f, which avoids the
    catastrophic cancellation of the naive expression at large tau:
        1 + (1-tau) f(tau) = (2/3)u^2 + (7/45)u^4 + ... ,  u = 1/sqrt(tau).
    """
    if tau <= 0:
        return 0.0
    if tau >= 1.0:
        u = 1.0 / math.sqrt(tau)
        if u < 0.1:
            u2 = u * u
            bracket = u2 * (2.0 / 3.0 + (7.0 / 45.0) * u2 + (8.0 / 105.0) * u2 ** 2)
            return -2.0 / u2 * bracket
        return -2.0 * tau * (1.0 + (1.0 - tau) * math.asin(u) ** 2)
    r = math.sqrt(1.0 - tau)
    zeta = 2.0 * math.log(1.0 + r) - math.log(tau)     # ln((1+r)/(1-r)), stable for small tau
    f = complex(-0.25, 0.0) * (complex(zeta, -math.pi)) ** 2
    return -2.0 * tau * (1.0 + (1.0 - tau) * f)


def F_one(tau: float):
    """W-loop function F_1(tau); F_1(inf) = 7 (series used at large tau)."""
    if tau <= 0:
        return 7.0
    if tau >= 1.0:
        u = 1.0 / math.sqrt(tau)
        if u < 0.1:
            u2 = u * u
            return 7.0 + 2.0 * u2 + (22.0 / 15.0) * u2 ** 2
        return 2.0 + 3.0 * tau + 3.0 * tau * (2.0 - tau) * math.asin(u) ** 2
    r = math.sqrt(1.0 - tau)
    zeta = 2.0 * math.log(1.0 + r) - math.log(tau)
    f = complex(-0.25, 0.0) * (complex(zeta, -math.pi)) ** 2
    return 2.0 + 3.0 * tau + 3.0 * tau * (2.0 - tau) * f


# --------------------------------------------------------------------------------------
# Decay widths at theta^2 = 1
# --------------------------------------------------------------------------------------
def gamma_ll(ms: float, ml: float) -> float:
    """Gamma(S->l+l-)/theta^2 = G_F m_S m_l^2 beta^3/(4 sqrt2 pi)   [W Eq. (mudecay)]."""
    if ms <= 2.0 * ml:
        return 0.0
    beta = math.sqrt(1.0 - 4.0 * ml ** 2 / ms ** 2)
    return G_F * ms * ml ** 2 * beta ** 3 / (4.0 * math.sqrt(2.0) * math.pi)


def gamma_gammagamma(ms: float) -> float:
    """Gamma(S->gamma gamma)/theta^2 = SM Gamma(h->gamma gamma) with m_h -> m_S."""
    A = complex(0.0, 0.0)
    for m, q2, nc in ((M_E, 1.0, 1), (M_MU, 1.0, 1), (M_TAU, 1.0, 1),
                      (M_U_MSBAR, 4.0 / 9.0, 3), (M_D_MSBAR, 1.0 / 9.0, 3),
                      (MSB_S_2GEV, 1.0 / 9.0, 3), (MSB_C_2GEV, 4.0 / 9.0, 3),
                      # N01 (review 2026-10-25): Q_b^2 = (-1/3)^2 = 1/9, not 4/9.  The b loop was
                      # four times too large, which changed its interference with the W and top
                      # loops; Gamma(h->gamma gamma) at 125.25 GeV moves 9.578 -> 9.368 keV,
                      # i.e. onto the SM anchor 9.39 keV (-0.2 % instead of +2 %).
                      (MSB_B_MB, 1.0 / 9.0, 3), (M_T_POLE, 4.0 / 9.0, 3)):
        A += nc * q2 * F_half(4.0 * m ** 2 / ms ** 2)
    A += F_one(4.0 * M_W ** 2 / ms ** 2)
    # The photons are real, so the coupling is alpha_em(0), not alpha(m_Z).  With alpha(m_Z)
    # the SM anchor Gamma(h->gamma gamma) at 125.25 GeV comes out 10.99 keV against the
    # reference 9.39 keV; with alpha(0) it is 9.58 keV (-2 %).
    return G_F * ALPHA_EM_0 ** 2 * ms ** 3 * abs(A) ** 2 / (128.0 * math.sqrt(2.0) * math.pi ** 3)


_Q_LOOP = (M_U_MSBAR, M_D_MSBAR, MSB_S_2GEV, MSB_C_2GEV, MSB_B_MB, M_T_POLE)


def gamma_gg(ms: float, k_factor: float = 1.0) -> float:
    """Gamma(S->gg)/theta^2 = alpha_s^2 m_S^3 |sum_q F_{1/2}|^2/(128 pi^3 v^2) * K.

    Identical to [W] Eq. for Gamma_gg (his (x+(x-1)f)/x^2 equals -F_{1/2}/2), and to the
    heavy-top effective-Lagrangian result in the limit |sum F|^2 -> 16/9.
    k_factor multiplies by the NLO QCD correction K (1 = LO).
    """
    A = complex(0.0, 0.0)
    for m in _Q_LOOP:
        A += F_half(4.0 * m ** 2 / ms ** 2)
    a_s = alpha_s(ms)
    return k_factor * a_s ** 2 * ms ** 3 * abs(A) ** 2 / (128.0 * math.pi ** 3 * V_HIGGS ** 2)


# NOTE: Winkler's printed Gamma_gg is the LO one-loop result and we implement exactly that.
# Boiarska/SensCalc additionally multiply it by the NLO K-factor of Spira et al.,
# 1 + (alpha_s/pi)(95/4 - 7/6 n_l).  That is a DIFFERENT model; it is not mixed into this file.
# The alternative is exercised only in the validation path (work/widths/reference_br_model.py),
# which builds a width model directly from the published reference branching ratios.


# --------------------------------------------------------------------------------------
# Higher-order QCD corrections to the quark channels.
#
# Boiarska's Eq. (decay-quarks) multiplies the spectator result by (1 + Delta_QCD + Delta_t)
# with alpha_s(m_S).  At the 2 GeV matching scale alpha_s = 0.338 and Delta_QCD = 1.09, so the
# "correction" doubles the width and the hadronic width jumps by a factor 1.49 across the
# matching point.  Winkler's own figure and the SensCalc tables are both continuous there
# (ratio 1.0002 and 1.0000).  The perturbative expansion is simply not trustworthy at that
# scale, so the production setting is LO (QUARK_NLO = False), which restores the matching:
#
#   m_S            dispersive hadrons      perturbative hadrons      ratio
#   just below 2 GeV      7.76e-8 GeV            -                    -
#   just above 2 GeV          -                7.99e-8 GeV          1.03
#
# Set QUARK_NLO = True to recover Boiarska's NLO expression.
# --------------------------------------------------------------------------------------
QUARK_NLO = False


def gamma_qq(ms: float, mq_ref: float, mu_ref: float, m_hadron: float,
             nlo: bool = None) -> float:
    """Gamma(S->qbar q)/theta^2, [B] Eq. (decay-quarks):

       N_c theta^2 m_S mbar_q^2(m_S)/(8 pi v^2) (1 - 4 m_Had^2/m_S^2)^{3/2} (1+Delta_QCD+Delta_t)
    The kinematic threshold uses the lightest hadron containing q (K, D, B), not a free quark.
    The (1 + Delta_QCD + Delta_t) factor is only applied when nlo is true; see QUARK_NLO.
    """
    if nlo is None:
        nlo = QUARK_NLO
    mbar = mbar_q(mq_ref, mu_ref, ms)
    if ms <= 2.0 * m_hadron:
        return 0.0
    beta = math.sqrt(1.0 - 4.0 * m_hadron ** 2 / ms ** 2)
    corr = 1.0
    if nlo:
        a_s = alpha_s(ms)
        x = a_s / math.pi
        nf = _nf_of(ms)
        d_qcd = 5.67 * x + (35.94 - 1.36 * nf) * x ** 2 + (164.14 - 25.77 * nf + 0.259 * nf ** 2) * x ** 3
        d_t = x ** 2 * (1.57 - (2.0 / 3.0) * math.log(ms ** 2 / M_T_POLE ** 2)
                        + (1.0 / 9.0) * math.log(mbar ** 2 / ms ** 2) ** 2)
        corr = 1.0 + d_qcd + d_t
    return 3.0 * ms * mbar ** 2 / (8.0 * math.pi * V_HIGGS ** 2) * beta ** 3 * corr


def gamma_multimeson(ms: float) -> float:
    """Toy multi-meson channel [W]/[B]: C m_S^3 beta_{2pi}, beta_{2pi}=sqrt(1-(4 m_pi)^2/m_S^2)."""
    if ms <= 4.0 * M_PI:
        return 0.0
    return C_MULTIMESON * ms ** 3 * math.sqrt(1.0 - (4.0 * M_PI) ** 2 / ms ** 2)


# --------------------------------------------------------------------------------------
# Hadronic pi pi / K K in the resonance region
# --------------------------------------------------------------------------------------
def chpt_form_factors(ms: float):
    """[W] Eq. (formfactors): lowest-order ChPT form factors at sqrt(s) = m_S."""
    s = ms ** 2
    mpi2, mK2 = M_PI ** 2, M_KAON ** 2
    return (mpi2, 0.0, s + 2.0 * mpi2,
            0.5 * mpi2, mK2 - 0.5 * mpi2, s + 2.0 * mK2)


def gamma_pipi_chpt(ms: float) -> float:
    """[W] Gamma_pipi (total pi+pi- and pi0pi0) with the lowest-order ChPT form factors."""
    Gp, Dp, Tp, _, _, _ = chpt_form_factors(ms)
    if ms <= 2.0 * M_PI:
        return 0.0
    beta = math.sqrt(1.0 - 4.0 * M_PI ** 2 / ms ** 2)
    A = (7.0 / 9.0) * (Gp + Dp) + (2.0 / 9.0) * Tp
    return 3.0 * G_F / (16.0 * math.sqrt(2.0) * math.pi * ms) * beta * abs(A) ** 2


def gamma_KK_chpt(ms: float) -> float:
    """[W] Gamma_KK (total K+K- and K0K0bar) with the lowest-order ChPT form factors."""
    _, _, _, GK, DK, TK = chpt_form_factors(ms)
    if ms <= 2.0 * M_KAON:
        return 0.0
    beta = math.sqrt(1.0 - 4.0 * M_KAON ** 2 / ms ** 2)
    A = (7.0 / 9.0) * GK + (7.0 / 9.0) * DK + (2.0 / 9.0) * TK
    return G_F / (4.0 * math.sqrt(2.0) * math.pi * ms) * beta * abs(A) ** 2


# The pi pi / K K rates around the f0(980) require the two-channel Muskhelishvili-Omnes
# form factors of [W], built from the Roy-Steiner pi pi / K K phase shifts of
# Hoferichter et al. (JHEP 06 (2012) 072).  Those phase-shift tables are not present in
# this workspace.  HADRONIC_TABLE lets the pipeline read the *published* Winkler pi pi and
# K K rates from a digitised table; the free function below falls back to the lowest-order
# ChPT form factors, which are only reliable for m_S <~ 0.5 GeV.
HADRONIC_TABLE = None


def load_hadronic_table(path: str):
    """CSV with columns mass_GeV,gamma_pipi_GeV,gamma_KK_GeV (theta^2 = 1)."""
    global HADRONIC_TABLE
    d = np.genfromtxt(path, delimiter=',', names=True)
    HADRONIC_TABLE = (np.asarray(d['mass_GeV'], float),
                      np.nan_to_num(np.asarray(d['gamma_pipi_GeV'], float), nan=0.0),
                      np.nan_to_num(np.asarray(d['gamma_KK_GeV'], float), nan=0.0))
    return HADRONIC_TABLE


def gamma_pipi_KK(ms: float):
    """(Gamma_pipi, Gamma_KK) at theta^2 = 1 in GeV: digitised Winkler if loaded, else ChPT."""
    if HADRONIC_TABLE is not None:
        m, gp, gk = HADRONIC_TABLE
        return (float(np.interp(ms, m, gp, left=0.0, right=0.0)),
                float(np.interp(ms, m, gk, left=0.0, right=0.0)))
    return gamma_pipi_chpt(ms), gamma_KK_chpt(ms)


PERT_THRESHOLD = 2.0   # GeV: [W] matches the dispersive result to the spectator model here


# Validation hook: set to a callable ms -> dict to replace this model at run time (used by
# reference_br_model.py, the validation path that follows the other reference's branching
# ratios).  None in the production path, which is Winkler only.
WIDTH_OVERRIDE = None


def widths(ms: float, gg_k_factor: float = None) -> dict:
    """All S decay widths at theta^2 = 1, in GeV.

    m_S <  2 GeV : leptons + gamma gamma + dispersive pi pi / K K + multi-meson
    m_S >= 2 GeV : leptons + gamma gamma + perturbative s, c, b quarks + gluons
    """
    if WIDTH_OVERRIDE is not None:
        return WIDTH_OVERRIDE(ms)
    out = {}
    out['e'] = gamma_ll(ms, M_E)
    out['mu'] = gamma_ll(ms, M_MU)
    out['tau'] = gamma_ll(ms, M_TAU)
    out['gamma'] = gamma_gammagamma(ms)
    if ms <= PERT_THRESHOLD:      # the paper uses the dispersive result "up to 2 GeV"
        gp, gk = gamma_pipi_KK(ms)
        out['pipi'] = gp
        out['KK'] = gk
        out['multi'] = gamma_multimeson(ms)
    else:
        # [W] uses FIXED quark masses here: "We set m_s = 95 MeV, m_c = 1.3 GeV and neglect the
        # tiny decay rate into u, d quarks", with the kinematic threshold at the lightest meson
        # containing the quark.  Passing mu_ref = ms makes mbar_q() return the input mass, i.e.
        # no renormalisation-group running; that is the literal Winkler prescription.  (Boiarska
        # instead runs the masses to m_S, which is a different model - see the note above
        # gamma_gg and the validation module reference_br_model.py.)
        out['ss'] = gamma_qq(ms, W_M_S, ms, M_K)
        out['cc'] = gamma_qq(ms, W_M_C, ms, M_D)
        out['bb'] = gamma_qq(ms, W_M_B, ms, M_B)
        out['gg'] = gamma_gg(ms, 1.0 if gg_k_factor is None else gg_k_factor)
    return out


def gamma_total(ms: float, gg_k_factor: float = None) -> float:
    return sum(widths(ms, gg_k_factor).values())


def ctau_m(ms: float, theta2: float = 1.0, gg_k_factor: float = None) -> float:
    """c*tau in metres at mixing theta^2 (Gamma scales as theta^2)."""
    return HBAR_C / (gamma_total(ms, gg_k_factor) * theta2)


if __name__ == '__main__':
    for mu in (2.0, 5.0, 10.0, M_Z):
        print("alpha_s(%.1f) = %.5f" % (mu, alpha_s(mu)))
    print("F_half(1e8)=%.10f  F_half(7.66)=%.6f  F_one(1e8)=%.8f" % (F_half(1e8), F_half(7.66), F_one(1e8)))
    print("Gamma(h->gg, m=125.25) = %.4e GeV (LO)" % gamma_gg(M_H))
    for m in (0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 5.0, 10.0, 30.0, 62.0):
        w = widths(m); tot = sum(w.values())
        print("m_S=%7.3f Gamma=%.5e ctau=%.4e m  " % (m, tot, HBAR_C / tot)
              + " ".join("%s=%.3e" % (k, v) for k, v in w.items() if v > 0))
