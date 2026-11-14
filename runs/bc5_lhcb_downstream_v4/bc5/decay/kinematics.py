# -*- coding: utf-8 -*-
"""Relativistic phase space and fragmentation for the BC5 scalar decays.

Exact uniform n-body sampler (n = 2, 3, 4) built from binary-cluster factorisation: every node
contributes lambda^{1/2}(s; s_L, s_R)/s <= 1, so plain rejection is exact.  The previous
sequential generator was missing the sub-system phase-space volume, which biased the three-body
<s23> by 7 % and every four-body energy sharing (see docs/model.md, F03).

Units: GeV, c = 1.
"""
from __future__ import annotations

import math

from bc5.model.widths import (M_E, M_MU, M_TAU, M_PI_CHARGED, M_PI0, M_K_CHARGED, M_K0)

PI = math.pi
M_PI = M_PI_CHARGED
M_K = M_K_CHARGED
M_K0L = M_K0
M_ETA = 0.547862

# ======================================================================================
# relativistic kinematics
# ======================================================================================
def boost(parts, e_b, p_b):
    """Boost 5-tuples (m,e,px,py,pz) by the four-momentum (e_b, p_b)."""
    pb2 = p_b[0] ** 2 + p_b[1] ** 2 + p_b[2] ** 2
    if pb2 < 1e-24:
        return list(parts)
    b = [x / e_b for x in p_b]                 # boost velocity beta
    b2 = sum(x * x for x in b)
    g = 1.0 / math.sqrt(max(1e-15, 1.0 - min(b2, 1 - 1e-15)))
    out = []
    for (m, e, px, py, pz) in parts:
        bp = b[0] * px + b[1] * py + b[2] * pz
        ep = g * (e + bp)
        k = ((g - 1.0) * bp / b2) + g * e          # multiplies beta, not p_b
        out.append((m, ep, px + k * b[0], py + k * b[1], pz + k * b[2]))
    return out


def _isotropic_dir(rng):
    ct = 2.0 * rng.random() - 1.0
    st = math.sqrt(max(0.0, 1.0 - ct * ct))
    ph = 2.0 * PI * rng.random()
    return (st * math.cos(ph), st * math.sin(ph), ct)


def two_body(ms, m1, m2, rng):
    if ms < m1 + m2:
        return None
    p = math.sqrt(max((ms * ms - (m1 + m2) ** 2) * (ms * ms - (m1 - m2) ** 2), 0.0)) / (2.0 * ms)
    u = _isotropic_dir(rng)
    return [(m1, math.sqrt(m1 * m1 + p * p), p * u[0], p * u[1], p * u[2]),
            (m2, math.sqrt(m2 * m2 + p * p), -p * u[0], -p * u[1], -p * u[2])]


# ======================================================================================
# exact uniform n-body phase space  (F03, validated in work/review_check/numerics)
# ======================================================================================
def _lam(s, a2, b2):
    """Kallen function lambda(s, a^2, b^2) with a, b the masses."""
    sa, sb = math.sqrt(max(a2, 0.0)), math.sqrt(max(b2, 0.0))
    return (s - (sa + sb) ** 2) * (s - (sa - sb) ** 2)


_G3_CACHE = {}


def _g3_envelope(ms, m1, m2, m3):
    """(smin, smax, gmax) for the exact three-body rejection; cached per mass/masses."""
    key = (round(ms, 9), round(m1, 9), round(m2, 9), round(m3, 9))
    got = _G3_CACHE.get(key)
    if got is not None:
        return got
    smin, smax = (m2 + m3) ** 2, (ms - m1) ** 2
    gmax = 0.0
    for k in range(1, 201):
        s = smin + (smax - smin) * k / 200.0
        if s > 0.0:
            g = math.sqrt(max(_lam(ms * ms, m1 * m1, s), 0.0)) * \
                math.sqrt(max(_lam(s, m2 * m2, m3 * m3), 0.0)) / s
            gmax = max(gmax, g)
    out = (smin, smax, gmax * 1.02)
    _G3_CACHE[key] = out
    return out


def _three_body(ms, masses, rng, max_try=1000):
    m1, m2, m3 = masses
    smin, smax, gmax = _g3_envelope(ms, m1, m2, m3)
    if smax <= smin or gmax <= 0.0:
        return None
    while max_try > 0:
        max_try -= 1
        s = smin + (smax - smin) * rng.random()
        g = math.sqrt(max(_lam(ms * ms, m1 * m1, s), 0.0)) * \
            math.sqrt(max(_lam(s, m2 * m2, m3 * m3), 0.0)) / s
        if g >= rng.random() * gmax:
            break
    else:
        return None
    p = math.sqrt(max(_lam(ms * ms, m1 * m1, s), 0.0)) / (2.0 * ms)
    e1 = math.sqrt(m1 * m1 + p * p)
    e2 = math.sqrt(s + p * p)
    u = _isotropic_dir(rng)
    sub = two_body(math.sqrt(s), m2, m3, rng)
    if sub is None:
        return None
    sub = boost(sub, e2, (-p * u[0], -p * u[1], -p * u[2]))
    return [(m1, e1, p * u[0], p * u[1], p * u[2])] + sub


_G4_CACHE = {}


def _g4_envelope(ms, masses, n_grid=40):
    """(lo12, hi12, lo34, hi34, wmax) for the exact four-body rejection; cached.

    N04 (review 2026-10-25): the node weight w = prod lambda^{1/2}/s is much smaller than 1 close
    to threshold, so rejecting against a flat envelope of 1 made the sampler fail up to 78-100 %
    of the time at m_S ~ 0.56 GeV (four pions) and the caller counted those failures as "no
    charged daughters".  The acceptance is now w/wmax with wmax computed once per (mass, masses)
    on a 40x40 grid, which makes the failure probability negligible and the sampler exact.
    """
    key = (round(ms, 9),) + tuple(round(x, 9) for x in masses)
    got = _G4_CACHE.get(key)
    if got is not None:
        return got
    m1, m2, m3, m4 = masses
    lo12, hi12 = (m1 + m2) ** 2, (ms - m3 - m4) ** 2
    lo34, hi34 = (m3 + m4) ** 2, (ms - m1 - m2) ** 2
    wmax = 0.0
    if hi12 > lo12 and hi34 > lo34:
        ms2 = ms * ms
        for a in range(n_grid):
            s12 = lo12 + (hi12 - lo12) * (a + 0.5) / n_grid
            l12 = _lam(s12, m1 * m1, m2 * m2)
            if l12 <= 0.0:
                continue
            for b in range(n_grid):
                s34 = lo34 + (hi34 - lo34) * (b + 0.5) / n_grid
                l34 = _lam(s34, m3 * m3, m4 * m4)
                if l34 <= 0.0:
                    continue
                lout = _lam(ms2, s12, s34)
                if lout <= 0.0:
                    continue
                w = (math.sqrt(l12) / s12) * (math.sqrt(l34) / s34) * (math.sqrt(lout) / ms2)
                if w > wmax:
                    wmax = w
    out = (lo12, hi12, lo34, hi34, wmax * 1.02)
    _G4_CACHE[key] = out
    return out


def _four_body(ms, masses, rng, max_try=10000):
    """Binary tree (12)|(34): exact, because every node factor lambda^{1/2}(s;sL,sR)/s <= 1."""
    lo12, hi12, lo34, hi34, wmax = _g4_envelope(ms, masses)
    if wmax <= 0.0 or hi12 <= lo12 or hi34 <= lo34:
        return None
    m1, m2, m3, m4 = masses
    ms2 = ms * ms
    while max_try > 0:
        max_try -= 1
        s12 = lo12 + (hi12 - lo12) * rng.random()
        s34 = lo34 + (hi34 - lo34) * rng.random()
        l12 = _lam(s12, m1 * m1, m2 * m2)
        if l12 <= 0.0:
            continue
        l34 = _lam(s34, m3 * m3, m4 * m4)
        if l34 <= 0.0:
            continue
        lout = _lam(ms2, s12, s34)
        if lout <= 0.0:
            continue
        w = (math.sqrt(l12) / s12) * (math.sqrt(l34) / s34) * (math.sqrt(lout) / ms2)
        if w >= rng.random() * wmax:
            break
    else:
        raise RuntimeError('four-body rejection failed after %d tries (m_S=%.6f, masses=%r); '
                           'this is a generator failure, not a physics outcome' % (max_try, ms, masses))
    pout = math.sqrt(lout) / (2.0 * ms)
    e12 = math.sqrt(s12 + pout * pout)
    e34 = math.sqrt(s34 + pout * pout)
    u = _isotropic_dir(rng)
    pq = (pout * u[0], pout * u[1], pout * u[2])
    d12 = two_body(math.sqrt(s12), m1, m2, rng)
    d34 = two_body(math.sqrt(s34), m3, m4, rng)
    if d12 is None or d34 is None:
        return None
    out12 = boost(d12, e12, pq)
    out34 = boost(d34, e34, (-pq[0], -pq[1], -pq[2]))
    return out12 + out34


def nbody(ms, masses, rng):
    """Exact uniform n-body phase space for n = 1..4, unweighted.

    The previous implementation sampled the recursive sub-system mass uniformly and applied a
    rejection weight only for the FIRST two-body factor |p*|/M, while the sub-generator was
    normalised to unit probability.  The missing factor is the (n-1)-body phase-space VOLUME
    Phi_{n-1}(q), which is constant only for a massless two-body sub-system - exactly the case
    the old "massless three-body" test covered, which is why the bias went unnoticed.  It is
    large for a massive three-body system (3 pi, M = 1 GeV: <s23>/M^2 = 0.327 instead of
    0.3528) and for every four-body system (massless four-body: <E/M> = 0.333,0.222,0.222,0.222
    instead of 0.25 for all four).

    The replacement uses the standard cluster factorisation
        dPhi_n(P) = ds_A/(2pi) ds_B/(2pi) dPhi_2(P;q_A,q_B) dPhi_A(q_A) dPhi_B(q_B),
    where for a binary tree the joint density of the internal masses is the product of the node
    weights lambda^{1/2}(s; s_L, s_R)/s (each <= 1, so plain rejection is exact).
    n = 3 uses one node against the (2,3) pair; n = 4 uses the tree (12)|(34).
    Validated: massless four-body <E_i/M> = 0.2500 (all), massive three-body <s23> =
    0.35278 +- 0.00026 vs the exact 0.352813, four-momentum conservation 1e-12.
    """
    n = len(masses)
    if n == 1:
        return [(masses[0], masses[0], 0.0, 0.0, 0.0)]
    if ms < sum(masses):
        return None
    if n == 2:
        return two_body(ms, masses[0], masses[1], rng)
    if n == 3:
        return _three_body(ms, masses, rng)
    if n == 4:
        return _four_body(ms, masses, rng)
    raise NotImplementedError('nbody: n = %d is not implemented' % n)


def _poisson(nmean, rng):
    L = math.exp(-nmean); k = 0; p = 1.0
    while True:
        k += 1; p *= rng.random()
        if p <= L or k > 60:
            return k - 1


def jet_fragment(e_jet, rng, mult_factor=1.0):
    """Charged hadrons from a jet of energy e_jet, in the jet rest frame (jet along +z).

    Longitudinal phase space dN/dz ~ z(1-z) with a Gaussian transverse momentum of
    <pT> = 0.35 GeV; the charged multiplicity follows the measured jet-charge scaling
    <n_ch> = 1.0 + 1.5 ln(E_jet/1 GeV), scaled by mult_factor (gluon jets fragment harder).
    This is an explicit approximation to Pythia/Lund fragmentation.
    """
    if e_jet < 0.35:
        return []
    nmean = max(0.0, (1.0 + 1.5 * math.log(max(e_jet, 1e-3))) * mult_factor)
    n = _poisson(nmean, rng)
    out = []
    e_left = e_jet
    for _ in range(n):
        z = max(rng.random() ** 2, 0.02)
        e = z * e_left
        if e < 0.2:
            break
        pt = abs(rng.gauss(0.0, 0.35))
        ph = 2 * PI * rng.random()
        pz = math.sqrt(max(e * e - M_PI ** 2 - pt * pt, 1e-8))
        out.append((M_PI, e, pt * math.cos(ph), pt * math.sin(ph), pz))
        e_left -= e
        if e_left < 0.2:
            break
    return out


def _rotate_jet(jet, axis, rng):
    """Rotate jet momenta (jet along +z) onto the direction 'axis'."""
    ax = axis
    tmp = (0.0, 0.0, 1.0) if abs(ax[2]) < 0.9 else (1.0, 0.0, 0.0)
    e1 = (ax[1] * tmp[2] - ax[2] * tmp[1], ax[2] * tmp[0] - ax[0] * tmp[2], ax[0] * tmp[1] - ax[1] * tmp[0])
    n1 = math.sqrt(sum(x * x for x in e1)) or 1.0
    e1 = tuple(x / n1 for x in e1)
    e2 = (ax[1] * e1[2] - ax[2] * e1[1], ax[2] * e1[0] - ax[0] * e1[2], ax[0] * e1[1] - ax[1] * e1[0])
    out = []
    for (m, e, px, py, pz) in jet:
        vx = px * e1[0] + py * e2[0] + pz * ax[0]
        vy = px * e1[1] + py * e2[1] + pz * ax[1]
        vz = px * e1[2] + py * e2[2] + pz * ax[2]
        out.append((m, e, vx, vy, vz))
    return out
