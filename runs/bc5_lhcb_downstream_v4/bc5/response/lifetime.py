# -*- coding: utf-8 -*-
"""LHCb downstream detector response and the exact decay-in-window integral.

A(z) is the per-decay probability that the S daughters pass the downstream selection; it is
averaged over sampled decays and interpolated piecewise linearly in z.  The decay-in-window
probability is then integrated exactly in the exponential (F10), which is what makes the
short-lifetime edge unbiased (the old fixed-grid trapezoid was off by up to 1.5x).
"""
from __future__ import annotations

import math

import numpy as np

from bc5.decay import channels as C
from bc5.decay.kinematics import boost
from . import kernel as dk

Z_DEC_MIN, Z_DEC_MAX = 1.0, 2.5
ETA_S_MIN, ETA_S_MAX = 2.0, 5.0
RECO_EFF = 0.4
ZGRID = np.linspace(Z_DEC_MIN, Z_DEC_MAX, 13)

def _passes_pair(daughters_lab, x, y, z):
    """1.0 if at least one + and one - track of this decay survive the downstream selection.

    Each daughter may carry its own birth offset from the S decay vertex (7th element, metres),
    which is how the charged daughters of a flying neutral kaon are handled (F04): the transport
    kernel is given the true production point and rejects it itself when it is not
    reconstructable (BornAfterFront, outside the SciFi envelope, beam hole, ...).
    """
    good = []
    for a in daughters_lab:
        m, e, px, py, pz, ch = a[0], a[1], a[2], a[3], a[4], a[5]
        off = a[6] if len(a) > 6 else None
        if not dk.passes_birth_kinematics(px, py, pz, e, dk.BIRTH_CUTS['minE'], dk.BIRTH_CUTS['minPT']):
            continue
        if off is None:
            xb, yb, zb = x, y, z
        else:
            xb, yb, zb = x + off[0], y + off[1], z + off[2]
        acc, status, front, end = dk.transport_track(xb, yb, zb, px, py, pz, e, ch)
        if acc:
            good.append(ch)
    if len(good) < 2:
        return 0.0
    return 1.0 if (any(c > 0 for c in good) and any(c < 0 for c in good)) else 0.0


def acceptance_profile(ms, p4, v0_mm, rng, ndecay=12, zgrid=ZGRID):
    """A(z_dec) for one scalar, averaged over ndecay sampled decays."""
    e, px, py, pz = p4
    if pz <= 0:
        return np.zeros(len(zgrid))
    pmod = math.sqrt(px * px + py * py + pz * pz)
    if pmod <= abs(pz):
        return np.zeros(len(zgrid))
    eta = 0.5 * math.log((pmod + pz) / max(pmod - pz, 1e-300))    # pseudorapidity
    if not (ETA_S_MIN < eta < ETA_S_MAX):
        return np.zeros(len(zgrid))
    br, _ = C.channel_branching(ms)
    names = list(br.keys())
    cum = np.cumsum([br[n] for n in names])
    tot = np.zeros(len(zgrid))
    x0, y0, z0 = v0_mm[0] * 1e-3, v0_mm[1] * 1e-3, v0_mm[2] * 1e-3
    dirx, diry = px / pz, py / pz
    for _ in range(ndecay):
        u = rng.random()
        idx = int(np.searchsorted(cum, u))
        idx = min(idx, len(names) - 1)
        rest = C.s_decay_products(ms, names[idx], rng, s_p4=(e, px, py, pz))
        if not rest:
            # decay into no reconstructable product: it contributes zero response, but it MUST
            # stay in the denominator.  Dividing by the number of non-empty decays instead
            # would give a conditional efficiency and overestimate the visible response.
            continue
        lab5 = boost([(a[0], a[1], a[2], a[3], a[4]) for a in rest], e, (px, py, pz))
        lab = [(b[0], b[1], b[2], b[3], b[4], a[5], a[6])
               for b, a in zip(lab5, rest)]
        for k, zd in enumerate(zgrid):
            tot[k] += _passes_pair(lab, x0 + (zd - z0) * dirx, y0 + (zd - z0) * diry, zd)
    if ndecay <= 0:
        return np.zeros(len(zgrid))
    return tot / ndecay


def _p_window(alpha, z0, A, zgrid):
    """Decay-in-window probability, integrated EXACTLY in the exponential (F10).

        p = int dz alpha exp(-alpha (z - z0)) A(z) Theta(z - z0),   alpha = m_S /(c tau p_z)

    The previous code trapezoided the product on the fixed 13-node grid ZGRID.  For
    alpha*h >> 1 (short lifetime) the first cell contributes 0.5*alpha*h instead of
    1 - exp(-alpha h), an O(1) error per cell: with A = 1 the code gave ratios 1.127
    (decay length 0.1 m) / 0.555 (birth point inside the window) against the analytic answer.
    It also could not represent a birth point z0 between grid nodes.

    A(z) is known only on the grid, so it is taken piecewise LINEAR; then each cell integrates in
    closed form:

        cell_k = f1 * ( A(lo) * c0 + s_k * c1 ),
        f1 = exp(-alpha (lo - z0)),  lo = max(z_k, z0),  d = z_{k+1} - lo,
        c0 = 1 - exp(-alpha d),                       ( = alpha * int_0^d exp(-alpha u) du )
        c1 = (1 - (1 + x) exp(-x))/x * d,  x = alpha d ( = alpha * int_0^d u exp(-alpha u) du )
        s_k = (A_{k+1} - A_k)/h_k.

    This is exact for piecewise-linear A for any node count, reduces to the analytic answer for
    A = 1 to 2e-16, and is O(h^2) in the A interpolation.  Vectorised over rows.

    alpha, z0 : (N,)   A : (N, K)   zgrid : (K,)  ->  (N,)
    """
    z = np.asarray(zgrid, dtype=float)
    alpha = np.asarray(alpha, dtype=float).reshape(-1, 1)
    z0 = np.asarray(z0, dtype=float).reshape(-1, 1)
    A = np.asarray(A, dtype=float)
    if A.ndim == 1:
        A = A.reshape(1, -1)
    z1, z2 = z[:-1], z[1:]
    h = z2 - z1
    lo = np.maximum(z1[None, :], z0)
    d = np.clip(z2[None, :] - lo, 0.0, None)
    f1 = np.exp(-alpha * np.clip(lo - z0, 0.0, None))
    s_k = (A[:, 1:] - A[:, :-1]) / h[None, :]
    A_lo = A[:, :-1] + s_k * (lo - z1[None, :])
    x = alpha * d
    c0 = -np.expm1(-x)
    # N05 (review 2026-10-25): the slope term c1 = d (1-(1+x)e^{-x})/x has an O(x^2) numerator
    # formed by subtracting two numbers both close to one, so the direct evaluation loses digits
    # as x -> 0.  The old cut-over to a two-term series was at x < 1e-8, which left a relative
    # error of ~1.9 at x ~ 1e-8 (a long lifetime with a non-constant A(z)).  Use a four-term
    # series below x = 1e-3 (where it is accurate to 1e-16) and the expm1-based form above:
    #     1-(1+x)e^{-x} = x^2/2 - x^3/3 + x^4/8 - x^5/30 + O(x^6).
    c1_series = d * (x * x * (0.5 - x / 3.0 + x * x / 8.0 - x ** 3 / 30.0)) / np.where(x > 0.0, x, 1.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        # (1-(1+x)e^{-x})/x = c0/x - e^{-x}, with c0 = -expm1(-x) evaluated stably
        c1_stable = (c0 / np.where(x > 0.0, x, 1.0) - np.exp(-x)) * d
    c1 = np.where(x < 1e-3, c1_series, c1_stable)
    cell = f1 * (A_lo * c0 + s_k * c1)
    return np.clip(np.sum(cell, axis=1), 0.0, 1.0)

def p_window(ms, pz, z0, A, ctau_m, zgrid=ZGRID):
    """Decay-in-window probability for a VECTOR of scalars (F10).

    alpha = m_S /(c tau p_z), z0 = birth point in metres, A = A(z) on zgrid, ctau_m = c tau.
    See _p_window for the exact cell integration.
    """
    pz = np.asarray(pz, dtype=float)
    alpha = np.where(pz > 0.0, ms / np.maximum(ctau_m * np.maximum(pz, 1e-9), 1e-30), 0.0)
    if not np.any(alpha > 0.0):
        return np.zeros(len(pz))
    return _p_window(alpha, np.asarray(z0, dtype=float), A, zgrid)


def lifetime_probability(ms, p4, z0_m, a_z, ctau_m, zgrid=ZGRID):
    """p = int dz (m/(c tau pz)) exp(-(z-z0) m/(c tau pz)) A(z),  z > z0.

    Derived from p = int dl/lambda exp(-l/lambda) A(l) with l = (z-z0) p/pz and
    lambda = (p/m) c tau, which gives dl/lambda = dz m/(c tau pz).  See _p_window for the exact
    cell integration (piecewise-linear A, exact exponential).
    """
    e, px, py, pz = p4
    if pz <= 0 or ctau_m <= 0:
        return 0.0
    alpha = ms / (ctau_m * pz)
    return float(_p_window(np.array([alpha]), np.array([z0_m]),
                           np.asarray(a_z, dtype=float).reshape(1, -1), zgrid)[0])
