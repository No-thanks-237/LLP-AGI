# -*- coding: utf-8 -*-
"""
Muskhelishvili-Omnes machinery for the S -> pi pi / K K form factors.

Implements the equations of Winkler arXiv:1809.01876 (section "Dispersive Analysis"):
the unitarity relation, the two-channel T matrix, the once-subtracted Omnes dispersion
relation, and the form-factor assembly that feeds the decay rates.

The T-matrix input (the isoscalar S-wave pi pi / K K phases delta, psi and the pi pi -> K K
amplitude g) is supplied by phases.py.  This module is model independent.
"""
from __future__ import annotations
import math
import numpy as np

M_PI = 0.13957039
M_K = 0.493677
G_F = 1.1663788e-5

S_PIPI = (2.0 * M_PI) ** 2
S_KK = (2.0 * M_K) ** 2


def _beta(s):
    return (math.sqrt(max(1.0 - S_PIPI / s, 0.0)),
            math.sqrt(max(1.0 - S_KK / s, 0.0)))


def t_matrix(s, delta, psi, g):
    """Winkler's two-channel isoscalar S-wave T matrix and the diagonal beta factors."""
    bp, bk = _beta(s)
    op = 1.0 if s > S_PIPI else 0.0
    ok = 1.0 if s > S_KK else 0.0
    eta = math.sqrt(max(1.0 - 4.0 * bp * bk * g * g * ok, 0.0))
    T = np.zeros((2, 2), complex)
    if op > 0.0:
        T[0, 0] = (eta * np.exp(2j * delta) - 1.0) / (2j * bp)
        T[0, 1] = g * np.exp(1j * psi) * op * ok
        T[1, 0] = T[0, 1]
    if ok > 0.0:
        T[1, 1] = (eta * np.exp(2j * (psi - delta)) - 1.0) / (2j * bk)
    return T, np.array([bp, bk])


def _im_matrix(s, delta, psi, g):
    """M(s) with Im Omega = M Re Omega, from Im Omega = T* beta Omega."""
    T, beta = t_matrix(s, delta, psi, g)
    B = T.conj() * beta[None, :]
    return np.linalg.solve(np.eye(2, dtype=complex) - 1j * B, B)


def omnes(delta, psi, g, s_grid, n_iter=200, tol=1e-11):
    """Two-channel Omnes solution -- NOT YET VALIDATED, hence disabled.

    A discretised solution of
        Re Omega(s) = 1 + (s/pi) P int ds' Im Omega(s') / (s' (s' - s)),
        Im Omega(s) = M(s) Re Omega(s),
        M(s) = (1 - i T* beta)^{-1} T* beta
    is implemented below (direct linear solve; a fixed-point iteration is unusable because
    Im Omega = tan(delta) Re Omega diverges at delta = pi/2), but it does not yet reproduce
    the analytic single-channel benchmark
        Omega(s) = exp( (s/pi) P int delta(s')/(s'(s'-s)) ds' ) * exp(i delta(s)),
    so it is not wired into the production path.  Until it is, the pi pi / K K rates come
    from the digitised Winkler figure (hadronic_table.csv).

    What is missing to finish this:
      1. the two-channel T-matrix input: the isoscalar S-wave pi pi / K K phases delta, psi
         and the inelasticity eta from the Roy-Steiner analysis (Hoferichter et al., cited by
         Winkler as the source of his phase shifts).  This is an EXTERNAL data input, exactly
         as it is for Winkler -- it is not a formula in the paper;
      2. a validated discretisation of the principal-value integral (the Cauchy kernel near
         the endpoints and near delta = pi/2 needs a subtraction-aware quadrature, not the
         simple log formula currently used);
      3. regression against the digitised Winkler curves (work/widths/hadronic_table.csv),
         which is the natural validation target once (1) and (2) are in place.
    """
    raise NotImplementedError(
        'the two-channel Omnes solver is not validated yet; see the docstring. '
        'Use hadronic_table.csv through scalar_widths.load_hadronic_table() instead.')

# ---------------------------------------------------------------------------------------
# Form factors and rates  (Winkler eqs. after omnessolution)
# ---------------------------------------------------------------------------------------
P_OMNES = 0.73      # Winkler: p = 1 - 2 m_pi^2 Omega11'(0) - (4 m_K^2/sqrt3) Omega12'(0)
Q_OMNES = 0.52      # Winkler: q = 1 - sqrt3 m_pi^2 Omega21'(0) - 2 m_K^2 Omega22'(0)


def form_factors(s, Om, p=P_OMNES, q=Q_OMNES):
    """(Gamma_pi, Delta_pi, Theta_pi, Gamma_K, Delta_K, Theta_K) at s = m_S^2."""
    o11, o12, o21, o22 = Om[0, 0], Om[0, 1], Om[1, 0], Om[1, 1]
    mpi2, mK2 = M_PI ** 2, M_K ** 2
    G_pi = mpi2 * (o11 + o12 / math.sqrt(3.0))
    D_pi = (2.0 / math.sqrt(3.0)) * (mK2 - mpi2 / 2.0) * o12
    T_pi = (2.0 * mpi2 + p * s) * o11 + (2.0 / math.sqrt(3.0)) * (2.0 * mK2 + q * s) * o12
    G_K = (mpi2 / 2.0) * (math.sqrt(3.0) * o21 + o22)
    D_K = (mK2 - mpi2 / 2.0) * o22
    T_K = (math.sqrt(3.0) / 2.0) * (2.0 * mpi2 + p * s) * o21 + (2.0 * mK2 + q * s) * o22
    return G_pi, D_pi, T_pi, G_K, D_K, T_K


def rates(ms, Om):
    """(Gamma_pipi, Gamma_KK) at theta^2 = 1, in GeV, from Winkler's rate equations."""
    s = ms ** 2
    G_pi, D_pi, T_pi, G_K, D_K, T_K = form_factors(s, Om)
    bp, bk = _beta(s)
    A_pi = (7.0 / 9.0) * G_pi + (7.0 / 9.0) * D_pi + (2.0 / 9.0) * T_pi
    A_K = (7.0 / 9.0) * G_K + (7.0 / 9.0) * D_K + (2.0 / 9.0) * T_K
    g_pipi = 3.0 * G_F / (16.0 * math.sqrt(2.0) * math.pi * ms) * bp * abs(A_pi) ** 2
    g_KK = G_F / (4.0 * math.sqrt(2.0) * math.pi * ms) * bk * abs(A_K) ** 2
    return g_pipi, g_KK


# ---------------------------------------------------------------------------------------
def _selftest():
    """Unit test: with g = 0 and psi = delta the K K channel decouples, so Omega is a pure
    single-channel Omnès function and Watson's theorem requires arg Omega = delta."""
    s_grid = np.linspace(S_PIPI + 1e-6, 4.0, 1200)

    def delta(s):
        x = min(max((s - S_PIPI) / (1.69 - S_PIPI), 0.0), 1.0)
        return math.pi * (x * x * (3.0 - 2.0 * x))       # smooth 0 -> pi, saturates at pi

    Om = omnes(delta, delta, lambda s: 0.0, s_grid)
    ph = np.unwrap(np.angle(Om[:, 0, 0]))
    tgt = np.array([delta(s) for s in s_grid])
    d = np.abs(ph - tgt)
    sel = s_grid < 1.6
    print('single-channel Omnès check: max |arg Omega - delta| = %.3e rad  (mean %.3e, n=%d)'
          % (d[sel].max(), d[sel].mean(), int(sel.sum())))
    return d[sel].max()


if __name__ == '__main__':
    _selftest()

# ---------------------------------------------------------------------------------------
# Single-channel Omnes function (elastic region, s < 4 M_K^2)
#
# Below the K K threshold only one channel is open and Watson's theorem fixes the phase of
# Omega to be the pi pi isoscalar S-wave phase shift, so the once-subtracted dispersion
# relation has the closed-form solution
#
#     Omega(s) = exp{ (s/pi) P int ds' delta(s') / (s' (s' - s)) } * exp(i delta(s))
#
# (Winkler arXiv:1809.01876 eq. (2.33)-equivalent for the single channel).  Above the
# parametrisation range of the phase the integrand is continued smoothly to 2 pi, which makes
# it integrable; the sensitivity to that continuation is ~10% on |Omega|.
# ---------------------------------------------------------------------------------------
def omnes_single(s, delta, s_threshold, s_to_2pi=25.0, s_upper=400.0):
    """Complex single-channel Omnes function at s, from the phase shift delta(t) (radians).

    s_threshold : threshold of the channel in GeV^2
    s_to_2pi    : scale above which the phase is held at 2 pi (must exceed the range in which
                  delta is known)
    """
    import cmath
    from scipy.integrate import quad
    if s <= s_threshold:
        return complex(1.0, 0.0)
    d = lambda t: delta(t, s_to_2pi) if delta.__code__.co_argcount > 1 else delta(t)
    f = lambda t: d(t) / (t * (t - s))
    I = 0.0
    I += quad(f, s_threshold + 1e-12, max(s - 1e-9, s_threshold + 1e-12), limit=6000)[0]
    I += quad(f, s + 1e-9, s_upper, limit=6000)[0]
    return cmath.exp((s / math.pi) * I + 1j * d(s))
