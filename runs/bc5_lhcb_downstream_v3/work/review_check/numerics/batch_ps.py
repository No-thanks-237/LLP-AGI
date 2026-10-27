# -*- coding: utf-8 -*-
"""Vectorised (numpy) version of the exact uniform n-body sampler (BC5 F03 audit).

Same exact algorithm as correct_ps.py (documented there), but samples a whole batch at
once: 1e7 decays cost a few seconds instead of 1e7 Python calls.  Each batch is exact
(unweighted) - the rejection uses only the phase-space node factors lambda^{1/2}/s <= 1.
"""
import math
import numpy as _np
from correct_ps import _g3max, _lam

PI = math.pi


def _lam_np(s, a2, b2):
    return (s - (_np.sqrt(a2) + _np.sqrt(b2)) ** 2) * (s - (_np.sqrt(a2) - _np.sqrt(b2)) ** 2)


def _dirs_np(n, rng):
    ct = 2.0 * rng.random(n) - 1.0
    st = _np.sqrt(_np.maximum(0.0, 1.0 - ct * ct))
    ph = 2.0 * PI * rng.random(n)
    return st * _np.cos(ph), st * _np.sin(ph), ct


def _boost_np(e, px, py, pz, eb, bx, by, bz):
    """Boost a batch of 4-vectors from the rest frame of (eb,bx,by,bz) into the lab."""
    b2 = bx * bx + by * by + bz * bz
    g = eb / _np.sqrt(_np.maximum(eb * eb - b2, 1e-30))
    bp = bx * px + by * py + bz * pz
    ep = g * (e + bp / eb)
    k = ((g - 1.0) * bp / _np.maximum(b2, 1e-300)) + g * e / eb
    return ep, px + k * bx, py + k * by, pz + k * bz


_BLOCK_MAX = 2_000_000


def three_body_batch(ms, masses, n, rng):
    """n events of exact uniform 3-body phase space; returns (E,PX,PY,PZ), each (n,3)."""
    m1, m2, m3 = masses
    if ms < m1 + m2 + m3:
        return None
    smin, smax, gmax = _g3max(ms, m1, m2, m3)
    if smax <= smin or gmax <= 0.0:
        return None
    s = _np.empty(0)
    est = 0.5
    while len(s) < n:
        need = n - len(s)
        k = min(max(1024, int(need / est * 1.05) + 256), _BLOCK_MAX)
        cand = smin + (smax - smin) * rng.random(k)
        g = _np.sqrt(_np.maximum(_lam_np(ms * ms, m1 * m1, cand), 0.0)) * \
            _np.sqrt(_np.maximum(_lam_np(cand, m2 * m2, m3 * m3), 0.0)) / cand
        good = g >= rng.random(k) * gmax
        s = _np.concatenate([s, cand[good]])
        if good.any():
            est = min(0.99, max(0.02, 0.5 * est + 0.5 * good.mean()))
    s = s[:n]
    p = _np.sqrt(_np.maximum(_lam_np(ms * ms, m1 * m1, s), 0.0)) / (2.0 * ms)
    e1 = _np.sqrt(m1 * m1 + p * p); e2 = _np.sqrt(s + p * p)
    ux, uy, uz = _dirs_np(n, rng)
    sq = _np.sqrt(s)
    p12 = _np.sqrt(_np.maximum((s - (m2 + m3) ** 2) * (s - (m2 - m3) ** 2), 0.0)) / (2.0 * sq)
    ux2, uy2, uz2 = _dirs_np(n, rng)
    em2 = _np.sqrt(m2 * m2 + p12 * p12); em3 = _np.sqrt(m3 * m3 + p12 * p12)
    # subsystem momentum in the P frame is -p*u
    bx, by, bz = -p * ux, -p * uy, -p * uz
    e2a, x2, y2, z2 = _boost_np(em2, p12 * ux2, p12 * uy2, p12 * uz2, e2, bx, by, bz)
    e3a, x3, y3, z3 = _boost_np(em3, -p12 * ux2, -p12 * uy2, -p12 * uz2, e2, bx, by, bz)
    return (_np.stack([e1, e2a, e3a], axis=1), _np.stack([p * ux, x2, x3], axis=1),
            _np.stack([p * uy, y2, y3], axis=1), _np.stack([p * uz, z2, z3], axis=1))


def four_body_batch(ms, masses, n, rng):
    """n events of exact uniform 4-body phase space; returns (E,PX,PY,PZ), each (n,4)."""
    m1, m2, m3, m4 = masses
    if ms < m1 + m2 + m3 + m4:
        return None
    lo12 = (m1 + m2) ** 2; hi12 = (ms - m3 - m4) ** 2
    lo34 = (m3 + m4) ** 2; hi34 = (ms - m1 - m2) ** 2
    if hi12 <= lo12 or hi34 <= lo34:
        return None
    ms2 = ms * ms
    s12 = _np.empty(0); s34 = _np.empty(0)
    est = 0.05
    while len(s12) < n:
        need = n - len(s12)
        k = min(max(1024, int(need / est * 1.05) + 256), _BLOCK_MAX)
        a = lo12 + (hi12 - lo12) * rng.random(k)
        b = lo34 + (hi34 - lo34) * rng.random(k)
        l12 = _lam_np(a, m1 * m1, m2 * m2)
        l34 = _lam_np(b, m3 * m3, m4 * m4)
        lout = _lam_np(ms2, a, b)
        good = (l12 > 0) & (l34 > 0) & (lout > 0)
        w = _np.zeros(k)
        w[good] = (_np.sqrt(l12[good]) / a[good]) * (_np.sqrt(l34[good]) / b[good]) * \
                  (_np.sqrt(lout[good]) / ms2)
        acc = w >= rng.random(k)
        s12 = _np.concatenate([s12, a[acc]]); s34 = _np.concatenate([s34, b[acc]])
        if acc.any():
            est = min(0.99, max(0.005, 0.5 * est + 0.5 * acc.mean()))
    s12 = s12[:n]; s34 = s34[:n]
    lout = _lam_np(ms2, s12, s34)
    pout = _np.sqrt(_np.maximum(lout, 0.0)) / (2.0 * ms)
    e12 = _np.sqrt(s12 + pout * pout); e34 = _np.sqrt(s34 + pout * pout)
    Ux, Uy, Uz = _dirs_np(n, rng)
    bx, by, bz = pout * Ux, pout * Uy, pout * Uz
    q12 = _np.sqrt(s12)
    p12 = _np.sqrt(_np.maximum(_lam_np(s12, m1 * m1, m2 * m2), 0.0)) / (2.0 * q12)
    ux, uy, uz = _dirs_np(n, rng)
    a1 = _np.sqrt(m1 * m1 + p12 * p12); a2 = _np.sqrt(m2 * m2 + p12 * p12)
    e1, x1, y1, z1 = _boost_np(a1, p12 * ux, p12 * uy, p12 * uz, e12, bx, by, bz)
    e2, x2, y2, z2 = _boost_np(a2, -p12 * ux, -p12 * uy, -p12 * uz, e12, bx, by, bz)
    q34 = _np.sqrt(s34)
    p34 = _np.sqrt(_np.maximum(_lam_np(s34, m3 * m3, m4 * m4), 0.0)) / (2.0 * q34)
    vx, vy, vz = _dirs_np(n, rng)
    a3 = _np.sqrt(m3 * m3 + p34 * p34); a4 = _np.sqrt(m4 * m4 + p34 * p34)
    e3, x3, y3, z3 = _boost_np(a3, p34 * vx, p34 * vy, p34 * vz, e34, -bx, -by, -bz)
    e4, x4, y4, z4 = _boost_np(a4, -p34 * vx, -p34 * vy, -p34 * vz, e34, -bx, -by, -bz)
    return (_np.stack([e1, e2, e3, e4], axis=1), _np.stack([x1, x2, x3, x4], axis=1),
            _np.stack([y1, y2, y3, y4], axis=1), _np.stack([z1, z2, z3, z4], axis=1))


def nbody_batch(ms, masses, n, rng):
    k = len(masses)
    if k == 2:
        if ms < masses[0] + masses[1]:
            return None
        p = _np.sqrt(_np.maximum((ms * ms - (masses[0] + masses[1]) ** 2) *
                                 (ms * ms - (masses[0] - masses[1]) ** 2), 0.0)) / (2.0 * ms)
        ux, uy, uz = _dirs_np(n, rng)
        e1 = _np.sqrt(masses[0] ** 2 + p * p); e2 = _np.sqrt(masses[1] ** 2 + p * p)
        return (_np.stack([e1, e2], axis=1), _np.stack([p * ux, -p * ux], axis=1),
                _np.stack([p * uy, -p * uy], axis=1), _np.stack([p * uz, -p * uz], axis=1))
    if k == 3:
        return three_body_batch(ms, masses, n, rng)
    if k == 4:
        return four_body_batch(ms, masses, n, rng)
    raise NotImplementedError
