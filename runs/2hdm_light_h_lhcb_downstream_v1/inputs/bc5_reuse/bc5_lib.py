# -*- coding: utf-8 -*-
"""
BC5 scalar: decay kinematics, LHCb downstream response, lifetime integral.
Units: GeV, metres, c = 1.  Input sample vertices are mm and are converted to m.

Physics fixes of 2026-10-24 (response to logs/physics_review_20261024.md):

  F01  the Pythia hadronisation library stores every S decay with the parton/ tau pair
       along +-z.  A spin-0 S must emit its products isotropically with respect to the
       LABORATORY axes, so the whole library event is now rotated by one random SO(3)
       matrix per decay (never particle by particle).  Without this the jet axis is
       always the lab z axis and the high-mass acceptance is wrong by up to 60 %.
  F02  HadronLib.sample() now requires the mass node to match the requested mass; the
       library must be generated on the scan grid.  If a node is missing (within
       NODE_MTOL) the four-momenta are dilated by m/m_node, which keeps four-momentum
       conservation and the total invariant mass exact and only rescales the daughter
       masses by the (small) mismatch; a warning is printed once per channel/mass.
  F03  nbody() is the exact uniform n-body phase-space sampler for n = 2,3,4 (binary
       cluster factorisation with the lambb^{1/2}(s;sL,sR)/s node weights), replacing
       the previous sequential generator that was missing the sub-system volume.
  F04  the neutral kaons of S -> K0 K0bar are given a real flight: one decay length per
       kaon is sampled in the lab, the charged daughters are boosted from the true decay
       point and the transport kernel decides whether that point is reconstructable.
  F05  the merged multi-meson width is represented by an isoscalar four-pion system,
       decomposed into its charge states by the isospin weights (4/9 pi+pi-pi+pi-,
       4/9 pi+pi-pi0pi0, 1/9 four pi0 -> invisible) instead of an arbitrary 3-way split.
  F12  the photon channel is an explicit channel (no charged daughters) so that the
       branching-ratio vector sums to one.

Every changed piece is documented at the function that implements it.
"""
from __future__ import annotations
import math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WIDTHS = os.path.abspath(os.path.join(HERE, '..', 'widths'))
RESP = os.path.abspath(os.path.join(HERE, '..', 'response'))
sys.path.insert(0, WIDTHS)
sys.path.insert(0, RESP)
import scalar_widths as sw
import downstream_ctypes as dk

# The pi pi / K K rates in the resonance region come from the digitised Winkler figure.
_HADRONIC = os.path.join(WIDTHS, 'hadronic_table.csv')
if os.path.isfile(_HADRONIC):
    sw.load_hadronic_table(_HADRONIC)

PI = math.pi
M_E, M_MU, M_TAU = sw.M_E, sw.M_MU, sw.M_TAU
M_PI, M_PI0, M_K = sw.M_PI_CHARGED, sw.M_PI0, sw.M_K_CHARGED
M_K0L = sw.M_K0
M_ETA = 0.547862

Z_DEC_MIN, Z_DEC_MAX = 1.0, 2.5
ETA_S_MIN, ETA_S_MAX = 2.0, 5.0
RECO_EFF = 0.4
ZGRID = np.linspace(Z_DEC_MIN, Z_DEC_MAX, 13)


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


# ======================================================================================
# S decay channels
# ======================================================================================
# F12: 'gamma' is an explicit channel.  widths() always contained it, but the channel list did
# not, so sum(BR over CHANNELS) < 1 and the remainder was silently handed to the LAST channel
# by idx = min(...).  It contributes no charged daughter but it must be sampled as such.
CHANNELS = ['e', 'mu', 'tau', 'pipi', 'KK', 'multi', 'ss', 'cc', 'bb', 'gg', 'gamma']

# F04: neutral-kaon lifetimes, used to give K_S / K_L a real flight from the S decay point.
# c*tau = c * tau directly (tau in s, c in m/s); the first version of this line divided by 1/tau,
# which gives 1.8e-26 m instead of 2.7 cm and silently switched the whole kaon-flight fix off -
# caught by comparing against the independent implementation in work/review_check/kaon.
C_LIGHT = 2.99792458e8                        # m/s
CTAU_KS = C_LIGHT * 0.8954e-10                # 2.684 cm
CTAU_KL = C_LIGHT * 5.116e-8                  # 15.34 m
# Mass node tolerance of the hadronisation library (F02).  The library is generated on the
# scan grid, so a real mismatch means the grid changed and the library must be regenerated.
NODE_MTOL = 2e-3

_TAU_MIX = [(0.35, 'pi'), (0.30, 'l'), (0.15, 'rho'), (0.20, '3pi')]


def _tau_charged(m_tau, rng):
    """One or three charged tracks from a tau decay, in the tau rest frame."""
    x = rng.random(); c = 0.0; mode = 'pi'
    for w, m in _TAU_MIX:
        c += w
        if x <= c:
            mode = m; break
    if mode == '3pi':
        d = nbody(m_tau, [M_PI, M_PI, M_PI], rng)
        return [(M_PI, a[1], a[2], a[3], a[4]) for a in d]
    if mode == 'pi':
        d = two_body(m_tau, M_PI, 0.0, rng)
        return [(M_PI, d[0][1], d[0][2], d[0][3], d[0][4])]
    if mode == 'rho':
        d = two_body(m_tau, 0.775, 0.0, rng)
        return [(M_PI, d[0][1], d[0][2], d[0][3], d[0][4])]
    # three-body l nu nu : E_l spectrum ~ (m^2/2)(1-2E/m)^5 style; use ~ (1-2E/m)^3
    u = rng.random()
    e = 0.5 * m_tau * (1.0 - u ** 0.25)
    p = math.sqrt(max(e * e - M_PI ** 2, 0.0))
    ct = 2 * rng.random() - 1; st = math.sqrt(max(0.0, 1 - ct * ct)); ph = 2 * PI * rng.random()
    return [(M_PI, e, p * st * math.cos(ph), p * st * math.sin(ph), p * ct)]


def _K0L_charged(m_k, rng):
    d = nbody(m_k, [M_PI, 0.0, 0.0], rng)
    return [(M_PI, d[0][1], d[0][2], d[0][3], d[0][4])]


def _K0_charged(m_k, is_short, rng):
    """Charged daughters of ONE neutral kaon, in the kaon rest frame.

    K_S : pi+ pi- (69.2 %) or pi0 pi0 (30.8 %) -> two charged pions in the first case.
    K_L : pi+ pi- pi0 (12.5 %) -> two charged; pi l nu (68.1 %) -> one charged lepton plus one
          charged pion; the rest (3 pi0) is invisible.
    Returns (m, e, px, py, pz, charge) tuples.  The flight of the kaon itself is handled by the
    caller (F04), which knows the lab frame.
    """
    if is_short:
        if rng.random() >= 0.692:
            return []
        d = two_body(m_k, M_PI, M_PI, rng)
        if d is None:
            return []
        return [(M_PI, d[0][1], d[0][2], d[0][3], d[0][4], +1),
                (M_PI, d[1][1], d[1][2], d[1][3], d[1][4], -1)]
    u = rng.random()
    if u < 0.125:
        d = nbody(m_k, [M_PI, M_PI, M_PI0], rng)
        if d is None:
            return []
        return [(M_PI, d[0][1], d[0][2], d[0][3], d[0][4], +1),
                (M_PI, d[1][1], d[1][2], d[1][3], d[1][4], -1)]
    if u < 0.806:
        d = nbody(m_k, [M_PI, M_E, 0.0], rng)
        if d is None:
            return []
        return [(M_PI, d[0][1], d[0][2], d[0][3], d[0][4], +1),
                (M_E, d[1][1], d[1][2], d[1][3], d[1][4], -1)]
    return []


_LIB_MAP = {'gg': 'gg', 'ss': 'ss', 'cc': 'cc', 'bb': 'bb', 'tau': 'tautau'}


def s_decay_products(ms, channel, rng, s_p4=None):
    """Charged daughters of S -> channel in the S rest frame.

    Returns (m, e, px, py, pz, charge, offset) 7-tuples.  'offset' is the LAB displacement of
    the production point from the S decay vertex in metres (None if the particle is produced at
    the vertex).  It is only non-zero for the charged daughters of a neutral kaon, whose flight
    cannot be neglected (F04).

    s_p4 = (E, px, py, pz) of the S in the lab is needed to sample that flight; without it the
    kaons are treated as prompt (used by the standalone validation scripts only).

    For gg, ss, cc, bb and tau/tau the Pythia library is used when it has been loaded; the
    analytic models below remain as a fallback and for the elementary channels
    (e, mu, pi pi, K K, multi-meson).
    """
    if HADRON_LIB is not None and channel in _LIB_MAP:
        got = HADRON_LIB.sample(ms, _LIB_MAP[channel], rng)
        if got is not None:
            # F01: the library stores every decay with its parton/tau pair along +-z, because the
            # generator injects them that way.  The decay axis of a spin-0 particle must be
            # isotropic with respect to the LABORATORY axes, so rotate the WHOLE event by one
            # random SO(3) matrix (rotating particle by particle would destroy the internal
            # correlations).  The rotation preserves every invariant mass and the total
            # four-momentum.  Writing the library events already rotated would be equivalent but
            # would leave the requirement invisible in the code.
            got = _rotate_all(got, rng)
            return [(a[0], a[1], a[2], a[3], a[4], a[5], None) for a in got]
    out = []
    if channel in ('e', 'mu'):
        m = M_E if channel == 'e' else M_MU
        d = two_body(ms, m, m, rng)
        if d is None:
            return []
        return [(m, d[0][1], d[0][2], d[0][3], d[0][4], +1, None),
                (m, d[1][1], d[1][2], d[1][3], d[1][4], -1, None)]
    if channel == 'tau':
        d = two_body(ms, M_TAU, M_TAU, rng)
        if d is None:
            return []
        for i, sg in ((0, +1), (1, -1)):
            sub = _tau_charged(M_TAU, rng)
            for (m, e, px, py, pz) in boost(sub, d[i][1], d[i][2:5]):
                out.append((m, e, px, py, pz, sg, None))
        return out
    if channel == 'pipi':
        if rng.random() < 1.0 / 3.0:
            return []
        d = two_body(ms, M_PI, M_PI, rng)
        if d is None:
            return []
        return [(M_PI, d[0][1], d[0][2], d[0][3], d[0][4], +1, None),
                (M_PI, d[1][1], d[1][2], d[1][3], d[1][4], -1, None)]
    if channel == 'KK':
        if rng.random() < 0.5:
            d = two_body(ms, M_K, M_K, rng)
            if d is None:
                return []
            return [(M_K, d[0][1], d[0][2], d[0][3], d[0][4], +1, None),
                    (M_K, d[1][1], d[1][2], d[1][3], d[1][4], -1, None)]
        # Neutral K0 K0bar pair: the S is CP even, so the CP-even superposition gives EQUAL
        # proportions of K_S K_S and K_L K_L (prompt.md line 106).  Each kaon then decays with
        # its real charged branching ratios (see _K0_charged).
        #
        # F04: the kaon is given a real FLIGHT.  c*tau(K_S) = 2.68 cm and c*tau(K_L) = 15.3 m,
        # so at these energies K_L decays hundreds of metres downstream (its charged daughters
        # are never reconstructable) while K_S typically flies 0.7-1.3 m for a B-produced S and
        # 9-12 m for a Higgs-produced one, i.e. often PAST the SciFi front face.  The old code
        # put every daughter at the S decay vertex and always counted it.
        # One decay length per kaon is sampled in the LAB frame; the daughters are boosted by the
        # kaon's lab four-momentum and tagged with their true birth offset from the S vertex.  The
        # transport kernel then decides whether that birth point is reconstructable (it rejects
        # BornAfterFront, outside the SciFi envelope, beam hole, ...), so no special case is
        # needed for a kaon that decays inside the tracking volume.
        d = two_body(ms, M_K0L, M_K0L, rng)
        if d is None:
            return []
        is_short = rng.random() < 0.5
        ctau = CTAU_KS if is_short else CTAU_KL
        for i in (0, 1):
            # N02 (review 2026-10-25): ONE decay length per kaon, drawn ONCE, before the loop over
            # its daughters.  A kaon decays once, so all its direct daughters are produced at the
            # same space-time point; the first version drew a length inside the daughter loop,
            # giving each pion of the same K_S/K_L a different birth point (hence the wrong joint
            # reachability probability P instead of P^2).
            off = None
            if s_p4 is not None:
                kaon_lab = boost([(d[i][0], d[i][1], d[i][2], d[i][3], d[i][4])],
                                 s_p4[0], (s_p4[1], s_p4[2], s_p4[3]))[0]
                pmod = math.sqrt(kaon_lab[2] ** 2 + kaon_lab[3] ** 2 + kaon_lab[4] ** 2)
                if pmod > 0.0:
                    fly = -math.log(max(rng.random(), 1e-12)) * (pmod / M_K0L) * ctau
                    off = (fly * kaon_lab[2] / pmod, fly * kaon_lab[3] / pmod,
                           fly * kaon_lab[4] / pmod)
            for (m, e, px, py, pz, ch) in _K0_charged(M_K0L, is_short, rng):
                b = boost([(m, e, px, py, pz)], d[i][1], d[i][2:5])[0]
                out.append((b[0], b[1], b[2], b[3], b[4], ch, off))
        return out
    if channel == 'multi':
        # F05: the merged multi-meson width C m^3 beta_2pi of Winkler stands for
        # "4pi, eta eta, rho rho, ..." and the paper gives NO exclusive composition, so the
        # charge states of the four-pion system are fixed by isospin: an isoscalar pi-pi pair is
        # 2/3 (pi+pi-) and 1/3 (pi0pi0) - the same rule that already gives the 1/3 pi0pi0
        # fraction of the pi pi channel - and applying it to the two pairs of the 4pi state gives
        #     4/9  (pi+pi-)(pi+pi-)   -> four charged pions
        #     4/9  (pi+pi-)(pi0pi0)   -> pi+ pi- pi0 pi0
        #     1/9  (pi0pi0)(pi0pi0)   -> invisible
        # prompt.md section 3.1 asks for four pions with uniform four-body phase space; the
        # total rate is kept exactly as Winkler's, only the charge composition is made explicit.
        x = rng.random()
        if x < 4.0 / 9.0:
            masses, charges = [M_PI] * 4, [+1, +1, -1, -1]
        elif x < 8.0 / 9.0:
            masses, charges = [M_PI, M_PI, M_PI0, M_PI0], [+1, -1, 0, 0]
        else:
            return []            # four pi0: no charged daughter
        if ms < sum(masses):
            return []            # below 4 m_pi the toy width is zero anyway
        d = nbody(ms, masses, rng)
        if d is None:
            return []
        return [(d[i][0], d[i][1], d[i][2], d[i][3], d[i][4], charges[i], None)
                for i in range(len(d))]
    if channel in ('gg', 'ss', 'cc', 'bb'):
        e_jet = ms / 2.0
        mf = 2.25 if channel == 'gg' else 1.0
        u = _isotropic_dir(rng)
        for sgn in (+1, -1):
            ax = (sgn * u[0], sgn * u[1], sgn * u[2])
            jet = _rotate_jet(jet_fragment(e_jet, rng, mf), ax, rng)
            for (m, e, px, py, pz) in jet:
                out.append((m, e, px, py, pz, +1 if rng.random() < 0.5 else -1, None))
        return out
    return []       # gamma: no charged daughters


def _rotate_all(parts, rng):
    """Rotate every 3-momentum of one decay by the SAME uniform random SO(3) matrix (F01)."""
    q = [rng.gauss(0.0, 1.0) for _ in range(4)]
    n = math.sqrt(sum(x * x for x in q))
    if n <= 0.0:
        return parts
    w, x, y, z = (c / n for c in q)
    r = [[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
         [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
         [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]]
    out = []
    for a in parts:
        px, py, pz = a[2], a[3], a[4]
        out.append((a[0], a[1],
                    r[0][0] * px + r[0][1] * py + r[0][2] * pz,
                    r[1][0] * px + r[1][1] * py + r[1][2] * pz,
                    r[2][0] * px + r[2][1] * py + r[2][2] * pz,
                    a[5]))
    return out


# --------------------------------------------------------------------------------------
# Pythia 8.316 hadronisation / decay library
#
# work/magnus/hadron_lib/hadron_lib_<channel>.txt.gz holds, for every channel and every mass
# on the scan grid, 1500 full Pythia hadronisations of the isolated colour-singlet S decay,
# with the stable charged particles given in the S rest frame.  When present it replaces the
# analytic jet model and the leading-track tau model entirely.
# --------------------------------------------------------------------------------------
class HadronLib:
    """Random access to a Pythia-generated S-decay library."""

    def __init__(self, paths):
        self.chan = {}
        for name, path in paths.items():
            d = np.loadtxt(path, comments='#')
            if d.ndim == 1:
                d = d.reshape(1, -1)
            mass, samp, pid, chg = d[:, 0], d[:, 1].astype(np.int64), d[:, 2].astype(int), d[:, 3].astype(int)
            e, px, py, pz = d[:, 4], d[:, 5], d[:, 6], d[:, 7]
            m = np.sqrt(np.maximum(e ** 2 - px ** 2 - py ** 2 - pz ** 2, 0.0))
            # group rows into (mass node, sample) blocks; the file is written mass- then sample-major
            key = mass * 1e6 + samp
            _, start, counts = np.unique(key, return_index=True, return_counts=True)
            order = np.argsort(start)
            self.chan[name] = {
                'mass': mass, 'key': key, 'start': start[order], 'counts': counts[order],
                'e': e, 'px': px, 'py': py, 'pz': pz, 'm': m, 'chg': chg,
                'mass_nodes': np.unique(mass),
                # number of generation ATTEMPTS per node.  Groups whose decay produced no
                # charged particle are absent from the file, so they must be sampled as
                # explicit empty outcomes instead of being normalised away.
                'n_attempt': int(samp.max()) + 1 if len(samp) else 0,
            }

    def sample(self, ms, channel, rng):
        """One random S -> channel decay as [(m, e, px, py, pz, charge), ...] in the S rest frame."""
        c = self.chan.get(channel)
        if c is None:
            return None
        i = int(np.argmin(np.abs(c['mass_nodes'] - ms)))
        node = float(c['mass_nodes'][i])
        rel = abs(node - ms) / max(ms, node, 1e-30)
        if rel > 0.20:
            return None            # channel closed at this mass: the library has no entry
        if rel > NODE_MTOL:
            # N06 (review 2026-10-25): the earlier fallback scaled (m,E,p) of every stable
            # daughter by m_S/m_node.  That matches the parent four-momentum but changes the
            # PHYSICAL masses of pions, kaons and protons, which is not a legitimate on-shell
            # mapping.  The library is now generated on the exact scan grid, so a mismatch means
            # the grid changed: refuse instead of silently producing wrong final states.
            raise RuntimeError(
                'hadron library has no node for %s at m_S = %.6f (nearest %.6f, %.2f%% away); '
                'regenerate the library on this mass grid (work/magnus/bundles_stage_hadronlib3)'
                % (channel, ms, node, 100.0 * rel))
        # Keys are mass*1e6 + sample_index with sample_index < NSAMPLE_MAX.  The upper bound
        # must be the NEXT MASS NODE, not node*1e6 + 1e6: in key units +1e6 is +1 GeV, which
        # swept in every node within 1 GeV above and mixed their spectra into this one.
        nodes = c['mass_nodes']
        nxt = float(nodes[i + 1]) if i + 1 < len(nodes) else node + 1.0
        # Draw one generation ATTEMPT uniformly over the node's attempts.  Attempts whose decay
        # produced no charged particle have no rows in the file, so they must come back as an
        # explicit empty outcome; drawing only among the non-empty groups would renormalise the
        # visible fraction upwards (review P2: 198 of 1500 attempts missing for the gg node).
        n_att = int(c.get('n_attempt', 0))
        if n_att <= 0:
            return None
        tgt = node * 1e6 + rng.randrange(n_att)
        row = int(np.searchsorted(c['key'], tgt - 0.5, side='left'))
        if row >= len(c['key']) or abs(float(c['key'][row]) - tgt) > 0.4:
            return []
        # row is an index into the ROW array; start/counts are indexed by GROUP, so convert
        g = int(np.searchsorted(c['start'], row, side='right')) - 1
        if g < 0 or g >= len(c['start']):
            return []
        s, n = int(c['start'][g]), int(c['counts'][g])
        return [(float(c['m'][k]), float(c['e'][k]), float(c['px'][k]),
                 float(c['py'][k]), float(c['pz'][k]), int(c['chg'][k]))
                for k in range(s, s + n)]


HADRON_LIB = None


def load_hadron_lib(directory):
    """Load hadron_lib_<channel>.txt(.gz) from a directory; returns True on success."""
    global HADRON_LIB
    paths = {}
    for ch in ('gg', 'ss', 'cc', 'bb', 'tautau'):
        for ext in ('.txt.gz', '.txt'):
            f = os.path.join(directory, 'hadron_lib_%s%s' % (ch, ext))
            if os.path.isfile(f):
                paths[ch] = f
                break
    if not paths:
        return False
    HADRON_LIB = HadronLib(paths)
    return True


_CH_CACHE = {}


def channel_branching(ms):
    key = round(ms, 6)
    if key in _CH_CACHE:
        return _CH_CACHE[key]
    w = sw.widths(ms)
    tot = sum(w.values())
    br = {c: w.get(c, 0.0) / tot for c in CHANNELS}
    _CH_CACHE[key] = (br, tot)
    return _CH_CACHE[key]


# ======================================================================================
# detector response
# ======================================================================================
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
    br, _ = channel_branching(ms)
    names = list(br.keys())
    cum = np.cumsum([br[n] for n in names])
    tot = np.zeros(len(zgrid))
    x0, y0, z0 = v0_mm[0] * 1e-3, v0_mm[1] * 1e-3, v0_mm[2] * 1e-3
    dirx, diry = px / pz, py / pz
    for _ in range(ndecay):
        u = rng.random()
        idx = int(np.searchsorted(cum, u))
        idx = min(idx, len(names) - 1)
        rest = s_decay_products(ms, names[idx], rng, s_p4=(e, px, py, pz))
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
