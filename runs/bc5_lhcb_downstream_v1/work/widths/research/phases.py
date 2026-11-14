# -*- coding: utf-8 -*-
"""
Isoscalar S-wave pi pi phase shift and inelasticity.

Source: Garcia-Martin, Kaminski, Pelaez, Ruiz de Elvira, Yndurain,
        Phys. Rev. D 83 (2011) 074004 [arXiv:1102.2183], Appendix A section 1,
        Eqs. (A1)-(A4) and Table V, page 25-26.  The "CFD" column of Table V is used
        (constrained fit: satisfies forward dispersion relations, Roy and GKPY equations).

    delta_0^0(s): Eq. (A1)-(A2) for 4 M_pi^2 <= s <= (0.85 GeV)^2
                  Eq. (A3)     for (0.85 GeV)^2 < s <= (1.42 GeV)^2
    eta_0^0(s):   Eq. (A4), equal to 1 below the K K threshold 4 M_K^2

Validity: 2 M_pi <= sqrt(s) <= 1.42 GeV.  Above 1.42 GeV the phase is not determined by
this parametrisation; guide_to_2pi() continues it smoothly to 2 pi, which is the asymptotic
behaviour Winkler uses (arXiv:1809.01876, section "Dispersive Analysis").

This module replaces nothing by itself: it supplies the input for the Omnes dispersion in
omnes.py.  Original implementation of the parametrisation, transcribed verbatim from the
paper; validated against the paper's own Table XII phases.
"""
from __future__ import annotations
import math

M_PI = 0.13957
M_K = 0.496
M_ETA = 0.54751

S_MATCH = 0.85 ** 2
SQRT_S_MAX = 1.42

# Table V, CFD column
CFD = dict(B0=7.14, B1=-25.3, B2=-33.2, B3=-26.2, z0=M_PI,
           d0=226.5, c=-81.0, B=93.3, C=48.7, D=-88.3,
           e1=4.9, e2=-15.1, e3=4.7, e4=0.32)

_K_M = math.sqrt(M_K ** 2 - S_MATCH / 4.0)          # |k2(s_M)|


def _k2abs(s):
    return math.sqrt(abs(s / 4.0 - M_K ** 2))


def _k3(s):
    x = s / 4.0 - M_ETA ** 2
    return math.sqrt(x) if x > 0 else 0.0


def _w(s, p=CFD):
    s0 = 4 * M_K ** 2
    rs = math.sqrt(s)
    if s <= s0:
        return (rs - math.sqrt(s0 - s)) / (rs + math.sqrt(s0 - s))
    return complex(rs, -math.sqrt(s - s0)) / complex(rs, math.sqrt(s - s0))


def _cot_delta_low(s, p=CFD):
    br = (p['z0'] ** 2 / (M_PI * math.sqrt(s)) + p['B0'] + p['B1'] * _w(s, p)
          + p['B2'] * _w(s, p) ** 2 + p['B3'] * _w(s, p) ** 3)
    return (math.sqrt(s) / (2 * math.sqrt(s / 4.0 - M_PI ** 2))) * \
           (M_PI ** 2 / (s - 0.5 * p['z0'] ** 2)) * br


def _delta_low(s, p=CFD):
    return math.atan2(1.0, _cot_delta_low(s, p))


_h = 1e-6
_D_M_DEG = math.degrees(_delta_low(S_MATCH))
_D_MP_DEG = math.degrees((_delta_low(S_MATCH + _h) - _delta_low(S_MATCH - _h)) / (2 * _h))


def delta00(s, p=CFD):
    """delta_0^0(s) in radians, 4 M_pi^2 <= s <= (1.42 GeV)^2.

    Eq. (A3) is written in degrees in the paper, so every term of the intermediate branch is
    evaluated in degrees and converted once at the end.
    """
    if s <= S_MATCH:
        return _delta_low(s, p)
    a = _k2abs(s)
    if s < 4 * M_K ** 2:
        return math.radians(p['d0'] * (1 - a / _K_M) ** 2
                            + _D_M_DEG * (a / _K_M) * (2 - a / _K_M)
                            + a * (_K_M - a) * (8 * _D_MP_DEG + p['c'] * (_K_M - a) / M_K ** 3))
    return math.radians(p['d0'] + p['B'] * (a ** 2 / M_K ** 2) + p['C'] * (a ** 4 / M_K ** 4)
                        + p['D'] * (1.0 if s > 4 * M_ETA ** 2 else 0.0) * (a ** 2 / M_ETA ** 2))


def eta00(s, p=CFD):
    """eta_0^0(s); 1 below the K K threshold."""
    if s <= 4 * M_K ** 2:
        return 1.0
    a, k3, rs = _k2abs(s), _k3(s), math.sqrt(s)
    ex = -(a / rs) * (p['e1'] + p['e2'] * a / rs + p['e3'] * a ** 2 / s) ** 2 \
         - p['e4'] * (1.0 if s > 4 * M_ETA ** 2 else 0.0) * k3 ** 3 / rs
    return math.exp(ex)


def delta_guided(s, s_to_2pi=4.0 ** 2):
    """delta_0^0 continued to 2 pi above the parametrisation range (Winkler's asymptotics).

    Between 1.42 GeV and sqrt(s_to_2pi) the phase is interpolated linearly in ln(s) from
    the last CFD value to 2 pi and then held constant.
    """
    s_max = SQRT_S_MAX ** 2
    if s <= s_max:
        return delta00(s)
    d_end = delta00(s_max)
    if s >= s_to_2pi or s_to_2pi <= s_max:
        return 2 * math.pi
    t = math.log(s / s_max) / math.log(s_to_2pi / s_max)
    return d_end + (2 * math.pi - d_end) * t
