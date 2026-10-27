# -*- coding: utf-8 -*-
"""Exact uniform n-body phase-space sampler for n = 2, 3, 4 (audit of BC5 F03).

Algorithm (exact, unweighted, no matrix element):
  * n = 2 : two-body, isotropic (same as production).
  * n = 3 : draw s23 in [(m2+m3)^2, (M-m1)^2] with density
                propto  sqrt(lam(M^2,m1^2,s)) * sqrt(lam(s,m2^2,m3^2)) / s
            by rejection inside a constant envelope (cached); then two-body decays.
  * n = 4 : binary tree (12)|(34).  Draw s12 and s34 uniformly inside their allowed
            boxes and accept with
                sqrt(lam(s12,m1^2,m2^2))/s12 * sqrt(lam(s34,m3^2,m4^2))/s34
              * sqrt(lam(M^2,s12,s34))/M^2   (each factor <= 1, so valid envelope)
            then two two-body decays and two boosts.
  This is the standard exact cluster factorisation
      dPhi_n(P) = dsA/(2pi) dsB/(2pi) dPhi_2(P;qA,qB) dPhi_A(qA) dPhi_B(qB).
  For a binary tree the joint density of the internal masses is prod_nodes
  lam^{1/2}(s_node; s_L, s_R)/s_node and every node factor is <= 1, so the
  rejection acceptance probability is the product itself.  Unweighted, exact.
"""
import math, random

PI = math.pi
_sqrt = math.sqrt

def _lam(s, a2, b2):
    sa = _sqrt(a2); sb = _sqrt(b2)
    return (s - (sa + sb) ** 2) * (s - (sa - sb) ** 2)

def _dir(rng):
    ct = 2.0 * rng.random() - 1.0
    st = _sqrt(max(0.0, 1.0 - ct * ct))
    ph = 2.0 * PI * rng.random()
    return (st * math.cos(ph), st * math.sin(ph), ct)

def _boost(p4, e_b, pb):
    """Boost 5-tuple (m,e,px,py,pz) from the rest frame of (e_b,pb) to the lab."""
    m, e, px, py, pz = p4
    b2 = pb[0] * pb[0] + pb[1] * pb[1] + pb[2] * pb[2]
    if b2 < 1e-30:
        return p4
    g = e_b / _sqrt(max(e_b * e_b - b2, 1e-30))
    bp = pb[0] * px + pb[1] * py + pb[2] * pz
    ep = g * (e + bp / e_b)
    k = ((g - 1.0) * bp / b2) + g * e / e_b
    return (m, ep, px + k * pb[0], py + k * pb[1], pz + k * pb[2])

def two_body(ms, m1, m2, rng):
    if ms < m1 + m2:
        return None
    p = _sqrt(max((ms * ms - (m1 + m2) ** 2) * (ms * ms - (m1 - m2) ** 2), 0.0)) / (2.0 * ms)
    u = _dir(rng)
    e1 = _sqrt(m1 * m1 + p * p); e2 = _sqrt(m2 * m2 + p * p)
    return [(m1, e1, p * u[0], p * u[1], p * u[2]),
            (m2, e2, -p * u[0], -p * u[1], -p * u[2])]

# ---------------------------------------------------------------- n = 3
_G3 = {}

def _g3max(ms, m1, m2, m3):
    key = (ms, m1, m2, m3)
    v = _G3.get(key)
    if v is not None:
        return v
    smin = (m2 + m3) ** 2; smax = (ms - m1) ** 2
    gmax = 0.0
    for k in range(1, 201):
        s = smin + (smax - smin) * k / 200.0
        if s <= 0.0:
            continue
        g = _sqrt(max(_lam(ms * ms, m1 * m1, s), 0.0)) * _sqrt(max(_lam(s, m2 * m2, m3 * m3), 0.0)) / s
        if g > gmax:
            gmax = g
    gmax *= 1.02
    _G3[key] = (smin, smax, gmax)
    return _G3[key]

def three_body(ms, masses, rng, max_try=1000):
    m1, m2, m3 = masses
    if ms < m1 + m2 + m3:
        return None
    smin, smax, gmax = _g3max(ms, m1, m2, m3)
    if smax <= smin or gmax <= 0.0:
        return None
    while max_try > 0:
        max_try -= 1
        s = smin + (smax - smin) * rng.random()
        g = _sqrt(max(_lam(ms * ms, m1 * m1, s), 0.0)) * _sqrt(max(_lam(s, m2 * m2, m3 * m3), 0.0)) / s
        if g >= rng.random() * gmax:
            break
    else:
        return None
    # P rest frame: particle 1 against the (2,3) system of invariant mass sqrt(s)
    p = _sqrt(max(_lam(ms * ms, m1 * m1, s), 0.0)) / (2.0 * ms)
    e1 = _sqrt(m1 * m1 + p * p)
    e2 = _sqrt(s + p * p)
    u = _dir(rng)
    sub = two_body(_sqrt(s), m2, m3, rng)
    if sub is None:
        return None
    q = (-p * u[0], -p * u[1], -p * u[2])
    sub = [_boost(a, e2, q) for a in sub]
    return [(m1, e1, p * u[0], p * u[1], p * u[2])] + sub

# ---------------------------------------------------------------- n = 4
def four_body(ms, masses, rng, max_try=500):
    m1, m2, m3, m4 = masses
    if ms < m1 + m2 + m3 + m4:
        return None
    lo12 = (m1 + m2) ** 2; hi12 = (ms - m3 - m4) ** 2
    lo34 = (m3 + m4) ** 2; hi34 = (ms - m1 - m2) ** 2
    if hi12 <= lo12 or hi34 <= lo34:
        return None
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
        w = (_sqrt(l12) / s12) * (_sqrt(l34) / s34) * (_sqrt(lout) / ms2)
        if w >= rng.random():
            break
    else:
        return None
    # root: P -> q12 + q34
    pout = _sqrt(lout) / (2.0 * ms)
    e12 = _sqrt(s12 + pout * pout); e34 = _sqrt(s34 + pout * pout)
    U = _dir(rng)
    pq = (pout * U[0], pout * U[1], pout * U[2])
    # cluster 12
    d12 = two_body(_sqrt(s12), m1, m2, rng)
    # cluster 34
    d34 = two_body(_sqrt(s34), m3, m4, rng)
    if d12 is None or d34 is None:
        return None
    out12 = [_boost(a, e12, pq) for a in d12]
    out34 = [_boost(a, e34, (-pq[0], -pq[1], -pq[2])) for a in d34]
    return out12 + out34

def nbody(ms, masses, rng):
    n = len(masses)
    if n == 1:
        return [(masses[0], masses[0], 0.0, 0.0, 0.0)]
    if n == 2:
        return two_body(ms, masses[0], masses[1], rng)
    if n == 3:
        return three_body(ms, masses, rng)
    if n == 4:
        return four_body(ms, masses, rng)
    raise NotImplementedError("n=%d" % n)
