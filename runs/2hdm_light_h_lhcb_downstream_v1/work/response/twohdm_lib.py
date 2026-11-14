# -*- coding: utf-8 -*-
"""2HDM light-H: per-DECAY-CHANNEL kinematics and LHCb-downstream response.

Provenance
----------
Adapted from inputs/bc5_reuse/bc5_lib.py (BC5 v3.6, 2026-10-25).  Everything that
describes the H decay and its transport through the downstream detector is reused
VERBATIM from that file (same functions, same constants):

    boost, two_body, exact n-body phase space (_lam/_g3/_three_body/_g4/_four_body/nbody),
    jet_fragment (analytic fallback), _tau_charged, _K0_charged (K_S/K_L treatment),
    _rotate_all (F01 isotropic decay), HadronLib / load_hadron_lib (Pythia library),
    _passes_pair semantics, the exact piecewise-linear lifetime integral (_p_window/p_window).

The ONLY new layer (DESIGN.md section 4) is that the response is evaluated PER DECAY
CHANNEL for every H candidate, instead of sampling one channel with branching-ratio
weights (which is impossible when the BRs depend on tanb):

    A_ch(mH; i, z_j) = probability that H candidate i, decaying into channel ch at
                       z_j (13-node grid 1.0..2.5 m), yields >=1 opposite-charge
                       daughter pair passing BirthKinematics + DownstreamTransport
                       with the exact BC5 geometry (By=1 T 3..7 m, SciFi 7.7..9.4 m,
                       halfX 3.24, halfY 2.415, hole 0.09) and birth cuts
                       E>5 GeV, pT>0.5 GeV.

Channel names (DESIGN.md section 2): gaga, ee, mumu, tautau, pipi, KK, 4pi, cc, ss, gg.
BC5 aliases are accepted too (e, mu, tau, multi, gamma, bb) so one library serves both
the 2HDM engine and the BC5 v3.6 regression adapter.

Units: GeV, metres, c = 1.  Sample vertices are mm and are converted to m by the caller.
"""
from __future__ import annotations
import math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import downstream_ctypes as dk

PI = math.pi

# --- constants, identical to inputs/bc5_reuse/bc5_lib.py + BC5 work/widths/scalar_widths.py
M_E = 0.51099895e-3
M_MU = 0.1056583755
M_TAU = 1.77686
M_PI = 0.13957039
M_PI0 = 0.1349768
M_K = 0.493677
M_K0L = 0.497611
M_ETA = 0.547862
HBAR_C = 1.973269804593025e-16          # GeV m
C_LIGHT = 2.99792458e8                  # m/s
CTAU_KS = C_LIGHT * 0.8954e-10          # 2.684 cm
CTAU_KL = C_LIGHT * 5.116e-8            # 15.34 m
NODE_MTOL = 2e-3                        # hadron-library node tolerance (F02)

Z_DEC_MIN, Z_DEC_MAX = 1.0, 2.5
ETA_S_MIN, ETA_S_MAX = 2.0, 5.0
RECO_EFF = 0.4
ZGRID = np.linspace(Z_DEC_MIN, Z_DEC_MAX, 13)

#: 2HDM channel names fixed by DESIGN.md section 2 (the response must cover all of them).
CHANNELS_2HDM = ['gaga', 'ee', 'mumu', 'tautau', 'pipi', 'KK', '4pi', 'cc', 'ss', 'gg']

#: canonical channel name for every accepted alias (BC5 model adapter uses the aliases).
_ALIAS = {'gaga': 'gaga', 'gamma': 'gaga',
          'ee': 'ee', 'e': 'ee',
          'mumu': 'mumu', 'mu': 'mumu',
          'tautau': 'tautau', 'tau': 'tautau',
          'pipi': 'pipi',
          'KK': 'KK',
          '4pi': '4pi', 'multi': '4pi',
          'cc': 'cc', 'ss': 'ss', 'gg': 'gg', 'bb': 'bb'}

#: canonical channel -> Pythia hadron-library key (BC5: gg, ss, cc, bb, tau tau).
_LIB_OF = {'gg': 'gg', 'ss': 'ss', 'cc': 'cc', 'bb': 'bb', 'tautau': 'tautau'}

#: Mass windows in which the AUTHOR's Light_scalar_decay program gives the channel a
#: non-zero width, read off src/CPevenDecay.cpp (the program itself uses mpion = 0.13497):
#:   Gamma_pipi, Gamma_4pi, Gamma_KK are forced to 1e-40 for coe > 2 and below their
#:   thresholds; Gamma_gg and Gamma_ss are forced to 1e-40 below 2 GeV; Gamma_cc below
#:   2 m_D = 3.729686; Gamma_tautau below 2 m_tau = 3.5536; gaga is always open.
#: This is only the DEFAULT channel mask used when no decay table is available; the
#: production engine passes the exact mask from work/decay (BR > 0).  Channels outside
#: their window have an identically zero detector response (no decay exists).
OPEN_2HDM = {
    'gaga': (0.0, None),
    'ee': (2 * M_E, None),
    'mumu': (2 * M_MU, None),
    'tautau': (2 * 1.7768, None),
    'pipi': (2 * 0.13497, 2.0),
    'KK': (2 * 0.498, 2.0),
    '4pi': (4 * 0.13497, 2.0),
    'cc': (2 * 1.864843, None),
    'ss': (2.0, None),
    'gg': (2.0, None),
}

#: Default decays sampled per channel per H candidate.  The hadronisation channels get
#: fewer because one of their decays generates many daughters and their acceptance is
#: already averaged over many tracks; the MC error of A_ch enters the total signal only
#: through the average over thousands of candidates (see logs/response_validation.md).
DEFAULT_NDECAY = {'*': 10, 'cc': 6, 'ss': 6, 'gg': 6, 'tautau': 6}


def open_channels_2hdm(ms, channels=None):
    """Channels with a non-zero width at ms according to the author program's windows."""
    chans = list(channels) if channels is not None else list(CHANNELS_2HDM)
    out = []
    for ch in chans:
        c = _ALIAS.get(ch, ch)
        lo, hi = OPEN_2HDM.get(c, (0.0, None))
        if ms > lo and (hi is None or ms < hi):
            out.append(ch)
    return out


# ======================================================================================
# relativistic kinematics (verbatim from bc5_lib.py)
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
# exact uniform n-body phase space  (F03, verbatim from bc5_lib.py)
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
    """(lo12, hi12, lo34, hi34, wmax) for the exact four-body rejection; cached (N04)."""
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
    """Exact uniform n-body phase space for n = 1..4, unweighted (F03, verbatim)."""
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
    """Analytic fallback for gg/cc/ss/bb when the Pythia library is absent (verbatim)."""
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
    """Rotate jet momenta (jet along +z) onto the direction 'axis' (verbatim)."""
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
# tau and neutral-kaon decays (verbatim from bc5_lib.py)
# ======================================================================================
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


def _K0_charged(m_k, is_short, rng):
    """Charged daughters of ONE neutral kaon, in the kaon rest frame (F04, verbatim).

    K_S : pi+ pi- (69.2 %) or pi0 pi0 (30.8 %).
    K_L : pi+ pi- pi0 (12.5 %); pi l nu (68.1 %); the rest invisible.
    Returns (m, e, px, py, pz, charge) tuples; the kaon FLIGHT is handled by the caller.
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


# ======================================================================================
# Pythia 8.316 hadronisation / decay library (verbatim from bc5_lib.py)
# ======================================================================================
class HadronLib:
    """Random access to a Pythia-generated H-decay library (one file per channel)."""

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
                # number of generation ATTEMPTS per node.  Groups without charged daughters
                # are absent from the file and must be sampled as explicit empty outcomes.
                'n_attempt': int(samp.max()) + 1 if len(samp) else 0,
            }

    def has(self, channel):
        return channel in self.chan

    def sample(self, ms, channel, rng):
        """One random H -> channel decay as [(m, e, px, py, pz, charge), ...] in the H rest frame.

        Returns None if the channel has no library entry at this mass (closed below threshold);
        raises if a node exists but is further than NODE_MTOL away (grid changed -> regenerate).
        """
        c = self.chan.get(channel)
        if c is None:
            return None
        i = int(np.argmin(np.abs(c['mass_nodes'] - ms)))
        node = float(c['mass_nodes'][i])
        rel = abs(node - ms) / max(ms, node, 1e-30)
        if rel > 0.20:
            return None            # channel closed at this mass: the library has no entry
        if rel > NODE_MTOL:
            raise RuntimeError(
                'hadron library has no node for %s at m_H = %.6f (nearest %.6f, %.2f%% away); '
                'regenerate the library on this mass grid' % (channel, ms, node, 100.0 * rel))
        # Keys are mass*1e6 + sample_index with sample_index < n_attempt.
        nodes = c['mass_nodes']
        nxt = float(nodes[i + 1]) if i + 1 < len(nodes) else node + 1.0
        n_att = int(c.get('n_attempt', 0))
        if n_att <= 0:
            return None
        tgt = node * 1e6 + rng.randrange(n_att)
        row = int(np.searchsorted(c['key'], tgt - 0.5, side='left'))
        if row >= len(c['key']) or abs(float(c['key'][row]) - tgt) > 0.4:
            return []
        g = int(np.searchsorted(c['start'], row, side='right')) - 1
        if g < 0 or g >= len(c['start']):
            return []
        s, n = int(c['start'][g]), int(c['counts'][g])
        return [(float(c['m'][k]), float(c['e'][k]), float(c['px'][k]),
                 float(c['py'][k]), float(c['pz'][k]), int(c['chg'][k]))
                for k in range(s, s + n)]


HADRON_LIB = None
HADRON_LIB_DIR = None


def load_hadron_lib(directory, channels=('gg', 'ss', 'cc', 'bb', 'tautau')):
    """Load hadron_lib_<channel>.txt(.gz) from a directory; returns True on success."""
    global HADRON_LIB, HADRON_LIB_DIR
    paths = {}
    for ch in channels:
        for ext in ('.txt.gz', '.txt'):
            f = os.path.join(directory, 'hadron_lib_%s%s' % (ch, ext))
            if os.path.isfile(f):
                paths[ch] = f
                break
    if not paths:
        return False
    HADRON_LIB = HadronLib(paths)
    HADRON_LIB_DIR = os.path.abspath(directory)
    return True


# ======================================================================================
# H decay products, PER CHANNEL (the new layer)
# ======================================================================================
def decay_products(ms, channel, rng, s_p4=None):
    """Charged daughters of H -> channel in the H rest frame.

    Returns (m, e, px, py, pz, charge, offset) 7-tuples.  'offset' is the LAB displacement
    of the production point from the H decay vertex in metres (None if the particle is
    produced at the vertex); non-zero only for the charged daughters of a flying neutral
    kaon (F04).

    s_p4 = (E, px, py, pz) of the H in the lab, needed for the neutral-kaon flight.

    Channel treatment is EXACTLY the BC5 one (bc5_lib.s_decay_products):
      gaga   -> no charged daughters (A = 0)
      ee     -> e+ e-          (two-body)
      mumu   -> mu+ mu-        (two-body)
      tautau -> Pythia tau-decay library when loaded, else the BC5 analytic tau model
      pipi   -> 1/3 pi0 pi0 (invisible) + 2/3 pi+ pi-
      KK     -> 1/2 K+K- + 1/2 K0K0bar (CP-even), K_S/K_L with a real flight
      4pi    -> isospin 4/9 (pi+pi-)(pi+pi-), 4/9 (pi+pi-)(pi0pi0), 1/9 invisible
      cc,ss,gg -> Pythia hadron library when loaded, else the BC5 analytic jet fallback
    """
    ch = _ALIAS.get(channel)
    if ch is None:
        raise KeyError('unknown decay channel %r' % (channel,))
    if ch == 'gaga':
        return []
    if ch in ('ee', 'mumu'):
        m = M_E if ch == 'ee' else M_MU
        d = two_body(ms, m, m, rng)
        if d is None:
            return []
        return [(m, d[0][1], d[0][2], d[0][3], d[0][4], +1, None),
                (m, d[1][1], d[1][2], d[1][3], d[1][4], -1, None)]
    if ch == 'pipi':
        if rng.random() < 1.0 / 3.0:
            return []
        d = two_body(ms, M_PI, M_PI, rng)
        if d is None:
            return []
        return [(M_PI, d[0][1], d[0][2], d[0][3], d[0][4], +1, None),
                (M_PI, d[1][1], d[1][2], d[1][3], d[1][4], -1, None)]
    if ch == 'KK':
        if rng.random() < 0.5:
            d = two_body(ms, M_K, M_K, rng)
            if d is None:
                return []
            return [(M_K, d[0][1], d[0][2], d[0][3], d[0][4], +1, None),
                    (M_K, d[1][1], d[1][2], d[1][3], d[1][4], -1, None)]
        # Neutral K0 K0bar pair: the H is CP even, so equal K_S K_S and K_L K_L; each kaon
        # flies (F04) and decays with its real charged branching ratios (_K0_charged).
        d = two_body(ms, M_K0L, M_K0L, rng)
        if d is None:
            return []
        is_short = rng.random() < 0.5
        ctau = CTAU_KS if is_short else CTAU_KL
        out = []
        for i in (0, 1):
            # N02: ONE decay length per kaon, drawn once, so all its daughters share the
            # same birth point.
            off = None
            if s_p4 is not None:
                kaon_lab = boost([(d[i][0], d[i][1], d[i][2], d[i][3], d[i][4])],
                                 s_p4[0], (s_p4[1], s_p4[2], s_p4[3]))[0]
                pmod = math.sqrt(kaon_lab[2] ** 2 + kaon_lab[3] ** 2 + kaon_lab[4] ** 2)
                if pmod > 0.0:
                    fly = -math.log(max(rng.random(), 1e-12)) * (pmod / M_K0L) * ctau
                    off = (fly * kaon_lab[2] / pmod, fly * kaon_lab[3] / pmod,
                           fly * kaon_lab[4] / pmod)
            for (m, e, px, py, pz, chg) in _K0_charged(M_K0L, is_short, rng):
                b = boost([(m, e, px, py, pz)], d[i][1], d[i][2:5])[0]
                out.append((b[0], b[1], b[2], b[3], b[4], chg, off))
        return out
    if ch == '4pi':
        # F05: isoscalar four-pion system, isospin weights 4/9, 4/9, 1/9.
        x = rng.random()
        if x < 4.0 / 9.0:
            masses, charges = [M_PI] * 4, [+1, +1, -1, -1]
        elif x < 8.0 / 9.0:
            masses, charges = [M_PI, M_PI, M_PI0, M_PI0], [+1, -1, 0, 0]
        else:
            return []            # four pi0: no charged daughter
        if ms < sum(masses):
            return []
        d = nbody(ms, masses, rng)
        if d is None:
            return []
        return [(d[i][0], d[i][1], d[i][2], d[i][3], d[i][4], charges[i], None)
                for i in range(len(d))]
    if ch == 'tautau':
        if HADRON_LIB is not None and HADRON_LIB.has('tautau'):
            got = HADRON_LIB.sample(ms, 'tautau', rng)
            if got is None:
                return []                    # below 2 m_tau: channel closed
            got = _rotate_all(got, rng)      # F01
            return [(a[0], a[1], a[2], a[3], a[4], a[5], None) for a in got]
        d = two_body(ms, M_TAU, M_TAU, rng)
        if d is None:
            return []
        out = []
        for i, sg in ((0, +1), (1, -1)):
            sub = _tau_charged(M_TAU, rng)
            for (m, e, px, py, pz) in boost(sub, d[i][1], d[i][2:5]):
                out.append((m, e, px, py, pz, sg, None))
        return out
    # cc, ss, gg, bb: Pythia hadron library (DESIGN section 5); analytic fallback only
    # when the library has not been loaded at all.
    libkey = _LIB_OF[ch]
    if HADRON_LIB is not None and HADRON_LIB.has(libkey):
        got = HADRON_LIB.sample(ms, libkey, rng)
        if got is None:
            return []                        # closed at this mass (author's program: BR = 0)
        got = _rotate_all(got, rng)          # F01
        return [(a[0], a[1], a[2], a[3], a[4], a[5], None) for a in got]
    e_jet = ms / 2.0
    mf = 2.25 if ch == 'gg' else 1.0
    u = _isotropic_dir(rng)
    out = []
    for sgn in (+1, -1):
        ax = (sgn * u[0], sgn * u[1], sgn * u[2])
        jet = _rotate_jet(jet_fragment(e_jet, rng, mf), ax, rng)
        for (m, e, px, py, pz) in jet:
            out.append((m, e, px, py, pz, +1 if rng.random() < 0.5 else -1, None))
    return out


# ======================================================================================
# downstream transport: same kernel, same cuts (fast path over the ctypes wrapper)
# ======================================================================================
_CFG = dk.DEFAULT_CONFIG
_BY = float(_CFG['by'])
_FZMIN, _FZMAX = float(_CFG['fieldZMin']), float(_CFG['fieldZMax'])
_FRONT, _BACK = float(_CFG['frontZ']), float(_CFG['backZ'])
_HALFX, _HALFY = float(_CFG['halfX']), float(_CFG['halfY'])
_HOLE = float(_CFG['holeRadius'])
_FULL = 1 if _CFG['fullEnvelope'] else 0
_MIN_E = float(dk.BIRTH_CUTS['minE'])
_MIN_PT = float(dk.BIRTH_CUTS['minPT'])
_MIN_PT2 = _MIN_PT * _MIN_PT

try:                                            # reuse the wrapper's ctypes binding
    _LIB = dk._LIB
    _BUF = dk._OUT_ARRAY()
    _TRANSPORT = _LIB.ds_transport
    _FAST = True
except AttributeError:                          # defensive: fall back to the public API
    _FAST = False


def passes_pair(daughters_lab, x, y, z):
    """1.0 if >=1 + and >=1 - track of this decay survive the downstream selection.

    Each daughter may carry its own birth offset from the H decay vertex (7th element),
    which is how the charged daughters of a flying neutral kaon are handled (F04).
    """
    has_pos = False
    has_neg = False
    for a in daughters_lab:
        e = a[1]
        px = a[2]; py = a[3]; pz = a[4]
        if not (e > _MIN_E and (px * px + py * py) > _MIN_PT2):
            continue
        ch = a[5]
        off = a[6] if len(a) > 6 else None
        if off is None:
            xb, yb, zb = x, y, z
        else:
            xb, yb, zb = x + off[0], y + off[1], z + off[2]
        if _FAST:
            st = _TRANSPORT(xb, yb, zb, px, py, pz, e, ch,
                            _BY, _FZMIN, _FZMAX, _FRONT, _BACK, _HALFX, _HALFY, _HOLE,
                            _FULL, _BUF)
            ok = (st == 0)
        else:
            ok = dk.transport_track(xb, yb, zb, px, py, pz, e, ch)[0]
        if ok:
            if ch > 0:
                has_pos = True
            elif ch < 0:
                has_neg = True
            if has_pos and has_neg:
                return 1.0
    return 0.0


def acceptance_profiles(ms, p4, v0_mm, rng, channels=None, ndecay=DEFAULT_NDECAY, zgrid=ZGRID):
    """Per-channel A_ch(z_dec) for ONE H candidate: array (n_channels, n_z).

    This is DESIGN.md section 4: for each channel sampled independently, ndecay decays of
    H -> ch are generated with exactly the BC5 kinematics, boosted to the lab from the
    candidate's production point, and the fraction that yields an accepted opposite-charge
    pair is recorded at every z node.  Channels with no charged daughters (gaga) or closed
    at this mass contribute a zero row.  Empty decays stay in the DENOMINATOR (the
    efficiency is unconditional, as in BC5).

    Parameters
    ----------
    ms      : H mass [GeV]
    p4      : (E, px, py, pz) of the H in the lab [GeV]
    v0_mm   : production vertex (x, y, z) in MILLIMETRES (sample convention)
    rng     : random.Random
    channels: channel names; default = open_channels_2hdm(ms) (closed channels are 0)
    ndecay  : decays per channel: int or {channel: n, '*' : n} (default DEFAULT_NDECAY)
    """
    chans = list(channels) if channels is not None else open_channels_2hdm(ms)
    nz = len(zgrid)
    out = np.zeros((len(chans), nz), dtype=float)
    e, px, py, pz = p4
    if pz <= 0:
        return out
    pmod = math.sqrt(px * px + py * py + pz * pz)
    if pmod <= abs(pz):
        return out
    eta = 0.5 * math.log((pmod + pz) / max(pmod - pz, 1e-300))    # pseudorapidity
    if not (ETA_S_MIN < eta < ETA_S_MAX):
        return out
    x0, y0, z0 = v0_mm[0] * 1e-3, v0_mm[1] * 1e-3, v0_mm[2] * 1e-3
    dirx, diry = px / pz, py / pz
    zs = np.asarray(zgrid, dtype=float)
    # A candidate born at or after the last node can never decay inside the window:
    # A_ch is identically zero (the lifetime integral of bc5_lib also returns 0).
    if z0 >= float(zs[-1]):
        return out
    # Only the nodes from the one just below z0 onwards can contribute: the piecewise
    # linear A(z) enters the exact lifetime integral only for z > z0 (all earlier cells
    # have zero width, and bc5_lib._p_window integrates them to exactly zero).
    k_start = int(np.searchsorted(zs, z0, side='right')) - 1
    if k_start < 0:
        k_start = 0
    bx = x0 + (zs - z0) * dirx
    by = y0 + (zs - z0) * diry
    bp = [(float(bx[k]), float(by[k]), float(zs[k])) for k in range(k_start, nz)]
    if not isinstance(ndecay, dict):
        nd = {c: int(ndecay) for c in chans}
    else:
        nd = {c: int(ndecay.get(c, ndecay.get('*', 10))) for c in chans}
    s_p4 = (e, px, py, pz)
    for ic, chname in enumerate(chans):
        ch = _ALIAS.get(chname)
        n = nd[chname]
        if n <= 0 or ch == 'gaga':
            continue
        row = np.zeros(nz)
        for _ in range(n):
            rest = decay_products(ms, ch, rng, s_p4=s_p4)
            if not rest:
                continue
            lab5 = boost([(a[0], a[1], a[2], a[3], a[4]) for a in rest], e, (px, py, pz))
            lab = [(b[0], b[1], b[2], b[3], b[4], a[5], a[6])
                   for b, a in zip(lab5, rest)]
            for j, (xb, yb, zb) in enumerate(bp):
                row[k_start + j] += passes_pair(lab, xb, yb, zb)
        out[ic] = row / n
    return out


# ======================================================================================
# decay-in-window probability: exact piecewise-linear integration (F10, verbatim)
# ======================================================================================
def _p_window(alpha, z0, A, zgrid):
    """Decay-in-window probability, integrated EXACTLY in the exponential (F10).

        p = int dz alpha exp(-alpha (z - z0)) A(z) Theta(z - z0),   alpha = m_H /(c tau p_z)

    A(z) is piecewise linear on zgrid, so each cell integrates in closed form.
    Vectorised over rows: alpha, z0 : (N,)   A : (N, K)   zgrid : (K,)  ->  (N,).
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
    c1_series = d * (x * x * (0.5 - x / 3.0 + x * x / 8.0 - x ** 3 / 30.0)) / np.where(x > 0.0, x, 1.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        c1_stable = (c0 / np.where(x > 0.0, x, 1.0) - np.exp(-x)) * d
    c1 = np.where(x < 1e-3, c1_series, c1_stable)
    cell = f1 * (A_lo * c0 + s_k * c1)
    return np.clip(np.sum(cell, axis=1), 0.0, 1.0)


def p_window(ms, pz, z0, A, ctau_m, zgrid=ZGRID):
    """Decay-in-window probability for a VECTOR of H candidates (F10).

    alpha = m_H /(c tau p_z), z0 = birth point [m], A = A(z) on zgrid, ctau_m = c tau [m].
    """
    pz = np.asarray(pz, dtype=float)
    alpha = np.where(pz > 0.0, ms / np.maximum(np.asarray(ctau_m, dtype=float) * np.maximum(pz, 1e-9), 1e-30), 0.0)
    if not np.any(alpha > 0.0):
        return np.zeros(len(pz))
    return _p_window(alpha, np.asarray(z0, dtype=float), A, zgrid)


def lifetime_probability(ms, p4, z0_m, a_z, ctau_m, zgrid=ZGRID):
    """p = int dz (m/(c tau pz)) exp(-(z-z0) m/(c tau pz)) A(z), z > z0 (scalar version)."""
    e, px, py, pz = p4
    if pz <= 0 or ctau_m <= 0:
        return 0.0
    alpha = ms / (ctau_m * pz)
    return float(_p_window(np.array([alpha]), np.array([z0_m]),
                           np.asarray(a_z, dtype=float).reshape(1, -1), zgrid)[0])

# ======================================================================================
# parent-decay-position integrated window probability (Rao-Blackwellised p_window)
# ======================================================================================
def _cell_int_vertex(beta, A, z1, z2, h):
    """int_{z1}^{z2} (A1 + s (z - z1)) exp(-beta z) dz      (beta (N,1), A (N,K))."""
    x = beta * h                                  # (N, K-1)
    e = np.exp(-x)
    small = x < 1e-3
    c0 = np.where(small, h * (1.0 - x / 2.0 + x * x / 6.0 - x ** 3 / 24.0),
                  -np.expm1(-x) / beta)
    c1 = np.where(small, h * h * (0.5 - x / 3.0 + x * x / 8.0 - x ** 3 / 30.0),
                  (1.0 - e * (1.0 + x)) / (beta * beta))
    s = (A[:, 1:] - A[:, :-1]) / h
    return np.exp(-beta * z1) * (A[:, :-1] * c0 + s * c1)


def _cell_int_vertex_z(beta, A, z1, z2, h):
    """int_{z1}^{z2} (A1 + s (z - z1)) * z * exp(-beta z) dz   (beta (N,1), A (N,K))."""
    x = beta * h
    e = np.exp(-x)
    small = x < 1e-3
    J0 = np.where(small, h * (1.0 - x / 2.0 + x * x / 6.0 - x ** 3 / 24.0),
                  -np.expm1(-x) / beta)
    J1 = np.where(small, h * h * (0.5 - x / 3.0 + x * x / 8.0 - x ** 3 / 30.0),
                  (1.0 - e * (1.0 + x)) / (beta * beta))
    J2 = np.where(small, h ** 3 * (1.0 / 3.0 - x / 4.0 + x * x / 10.0 - x ** 3 / 36.0),
                  (2.0 - e * (2.0 + 2.0 * x + x * x)) / (beta ** 3))
    s = (A[:, 1:] - A[:, :-1]) / h
    return np.exp(-beta * z1) * (A[:, :-1] * z1 * J0
                                 + (A[:, :-1] + s * z1) * J1 + s * J2)


def _p_window_vertex(alpha, lam, A, zgrid):
    """p = int dz A(z) W(z),  W(z) = C [exp(-z/lam) - exp(-alpha z)],  C = alpha/(alpha*lam-1).

    This is the parent-decay-position average of the H decay-in-window probability: the
    parent is produced at z=0, flies a distance z0 ~ Exp(mean = lam) (lam = beta_z gamma
    c*tau of the parent), and the H then decays with rate alpha = m_H/(c*tau_H pz_H).
    Integrating z0 analytically removes the exponential sensitivity to the sampled parent
    decay vertex (the dominant variance source of the short-lived crossings).

    lam -> 0 : prompt parent, W -> alpha exp(-alpha z) (the old p_window at z0 = 0).
    lam <= 0 : parent flying backwards -> no forward H -> p = 0.
    alpha*lam -> 1 : removable singularity, W -> alpha^2 z exp(-alpha z) (exact limit
    branch, no parameter perturbation).
    Vectorised over rows: alpha (N,), lam (N,), A (N, K).
    """
    z = np.asarray(zgrid, dtype=float)
    alpha = np.asarray(alpha, dtype=float).reshape(-1, 1)
    lam = np.asarray(lam, dtype=float).reshape(-1, 1)
    A = np.asarray(A, dtype=float)
    z1, z2 = z[:-1], z[1:]
    h = z2 - z1
    n = alpha.shape[0]
    out = np.zeros(n)
    ok = lam[:, 0] > 0.0
    if not ok.any():
        return out
    lm = np.where(ok[:, None], lam, 1.0)
    # alpha*lam -> 1 is a REMOVABLE singularity of W: there W -> alpha^2 z exp(-alpha z).
    # Use the exact limit within |alpha*lam-1| < 1e-5 (the relative error of the limit
    # versus the full expression is then O(1e-5*alpha*z) < 1e-3) and the two-exponential
    # closed form outside (its cancellation error is ~1e-16/1e-5 = 1e-11).
    rel = alpha * lm - 1.0
    bad = np.abs(rel) < 1e-5
    lms = np.where(bad, 1.0, lm)
    with np.errstate(divide='ignore', invalid='ignore'):
        C = alpha / (alpha * lms - 1.0)
        I1 = _cell_int_vertex(1.0 / lms, A, z1, z2, h)
        I2 = _cell_int_vertex(alpha, A, z1, z2, h)
        p = C[:, 0] * (I1 - I2).sum(axis=1)
    if bad.any():
        idx = np.flatnonzero(bad[:, 0])
        al = alpha[idx]
        lim = (al[:, 0] ** 2) * _cell_int_vertex_z(al, A[idx], z1, z2, h).sum(axis=1)
        p[idx] = lim
    p = np.where(np.isfinite(p), p, 0.0)
    out = np.where(ok, p, 0.0)
    return np.clip(out, 0.0, 1.0)


def p_window_vertex(ms, pz, A, ctau_m, lam_z, zgrid=ZGRID):
    """p for a VECTOR of H candidates, with the parent decay position integrated out.

    ms      : H mass [GeV]
    pz      : H longitudinal momentum [GeV] (vector)
    A       : acceptance profile (N, K) on zgrid [1,2.5] m (0 outside the decay window)
    ctau_m  : c*tau_H [m]
    lam_z   : (pz_parent/m_parent) * c*tau_parent [m] (vector; <=0 -> p = 0)
    """
    pz = np.asarray(pz, dtype=float)
    alpha = np.where(pz > 0.0,
                     ms / np.maximum(np.asarray(ctau_m, dtype=float) * np.maximum(pz, 1e-9), 1e-30),
                     0.0)
    lam = np.asarray(lam_z, dtype=float)
    p_fwd = _p_window_vertex(alpha, np.maximum(lam, 0.0), A, zgrid)
    p_prompt = _p_window_vertex(alpha, np.full_like(lam, 1e-12), A, zgrid)
    # backward-flying parent: birth point at z = -u, u ~ Exp(|lam|), so the H has an extra
    # exp(-alpha u) suppression -> p = p_prompt / (1 + alpha |lam|)
    p_back = p_prompt / (1.0 + alpha * np.abs(lam))
    return np.where(lam > 0.0, p_fwd, np.where(lam < 0.0, p_back, p_prompt))
