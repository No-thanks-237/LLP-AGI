# -*- coding: utf-8 -*-
"""S decay channels: composition, charged daughters, hadronisation-library routing.

The channel list must reproduce the branching-ratio vector exactly (the photon channel is
explicit, F12).  Neutral kaons from S -> K0 K0bar are given a real flight (F04), and the merged
multi-meson width is decomposed into its isospin charge states (F05).
"""
from __future__ import annotations

import math

from bc5.model.widths import widths as _widths
from . import hadronlib
from .kinematics import (PI, M_E, M_MU, M_TAU, M_PI, M_PI0, M_K, M_K0L,
                         boost, two_body, nbody, jet_fragment, _isotropic_dir, _rotate_jet)

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
    if hadronlib.HADRON_LIB is not None and channel in _LIB_MAP:
        got = hadronlib.HADRON_LIB.sample(ms, _LIB_MAP[channel], rng)
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

_CH_CACHE = {}


def channel_branching(ms):
    key = round(ms, 6)
    if key in _CH_CACHE:
        return _CH_CACHE[key]
    w = _widths(ms)
    tot = sum(w.values())
    br = {c: w.get(c, 0.0) / tot for c in CHANNELS}
    _CH_CACHE[key] = (br, tot)
    return _CH_CACHE[key]
