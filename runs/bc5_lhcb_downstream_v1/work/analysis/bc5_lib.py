# -*- coding: utf-8 -*-
"""
BC5 scalar: decay kinematics, LHCb downstream response, lifetime integral.
Units: GeV, metres, c = 1.  Input sample vertices are mm and are converted to m.
"""
from __future__ import annotations
import math, os, sys, random
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


def nbody(ms, masses, rng):
    """Uniform n-body phase space via the sequential M-generator.

    The base case is n = 2 (a genuine two-body decay); the n = 1 case is only reached when
    the drawn sub-system invariant mass equals that single mass, and is handled by the
    n = 2 branch of the caller.
    """
    n = len(masses)
    if n == 1:
        return [(masses[0], masses[0], 0.0, 0.0, 0.0)]
    if n == 2:
        return two_body(ms, masses[0], masses[1], rng)
    if ms < sum(masses):
        return None
    m_rest = sum(masses[1:])
    # Uniform n-body phase space factorises as dPhi_n ~ dPhi_2 x dPhi_(n-1) x d(m2^2), but the
    # two-body factor dPhi_2 ~ |p*|/M still depends on the drawn subsystem mass m2.  Sampling
    # m2^2 uniformly and NOT weighting by |p*| biases the distribution: for massless three-body
    # the marginal of x = s23/M^2 must be 2(1-x) (<x^2> = 1/6), whereas an unweighted m2^2 draw
    # gives a flat marginal (<x^2> = 1/3).  Importance-sample m2^2 against |p*| by rejection.
    lo2, hi2 = m_rest * m_rest, (ms - masses[0]) ** 2
    if hi2 <= lo2:
        return None

    def _pstar(mm):
        lam = (ms * ms - (masses[0] + mm) ** 2) * (ms * ms - (masses[0] - mm) ** 2)
        return math.sqrt(max(lam, 0.0)) / (2.0 * ms)

    ps_max = 0.0
    for _t in range(17):
        ps_max = max(ps_max, _pstar(math.sqrt(lo2 + (hi2 - lo2) * _t / 16.0)))
    if ps_max <= 0.0:
        return None
    m2 = None
    for _try in range(64):
        cand = math.sqrt(rng.uniform(lo2, hi2))
        if _pstar(cand) >= rng.random() * ps_max:
            m2 = cand
            break
    if m2 is None:
        return None
    p = math.sqrt(max((ms * ms - (masses[0] + m2) ** 2) * (ms * ms - (masses[0] - m2) ** 2), 0.0)) / (2.0 * ms)
    u = _isotropic_dir(rng)
    e1 = math.sqrt(masses[0] ** 2 + p * p)
    e2 = math.sqrt(m2 * m2 + p * p)
    sub = nbody(m2, masses[1:], rng)
    if sub is None:
        return None
    sub = boost(sub, e2, (-p * u[0], -p * u[1], -p * u[2]))
    return [(masses[0], e1, p * u[0], p * u[1], p * u[2])] + sub


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
CHANNELS = ['e', 'mu', 'tau', 'pipi', 'KK', 'multi', 'ss', 'cc', 'bb', 'gg']

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


_LIB_MAP = {'gg': 'gg', 'ss': 'ss', 'cc': 'cc', 'bb': 'bb', 'tau': 'tautau'}


def s_decay_products(ms, channel, rng):
    """Charged daughters of S -> channel in the S rest frame as (m, e, px, py, pz, charge).

    For gg, ss, cc, bb and tau/tau the Pythia library is used when it has been loaded; the
    analytic models below remain as a fallback and for the elementary channels
    (e, mu, pi pi, K K, 4 pi).
    """
    if HADRON_LIB is not None and channel in _LIB_MAP:
        got = HADRON_LIB.sample(ms, _LIB_MAP[channel], rng)
        if got is not None:
            return got
    out = []
    if channel in ('e', 'mu'):
        m = M_E if channel == 'e' else M_MU
        d = two_body(ms, m, m, rng)
        if d is None:
            return []
        return [(m, d[0][1], d[0][2], d[0][3], d[0][4], +1), (m, d[1][1], d[1][2], d[1][3], d[1][4], -1)]
    if channel == 'tau':
        d = two_body(ms, M_TAU, M_TAU, rng)
        if d is None:
            return []
        for i, sg in ((0, +1), (1, -1)):
            sub = _tau_charged(M_TAU, rng)
            for (m, e, px, py, pz) in boost(sub, d[i][1], d[i][2:5]):
                out.append((m, e, px, py, pz, sg))
        return out
    if channel == 'pipi':
        if rng.random() < 1.0 / 3.0:
            return []
        d = two_body(ms, M_PI, M_PI, rng)
        if d is None:
            return []
        return [(M_PI, d[0][1], d[0][2], d[0][3], d[0][4], +1),
                (M_PI, d[1][1], d[1][2], d[1][3], d[1][4], -1)]
    if channel == 'KK':
        if rng.random() < 0.5:
            d = two_body(ms, M_K, M_K, rng)
            if d is None:
                return []
            return [(M_K, d[0][1], d[0][2], d[0][3], d[0][4], +1),
                    (M_K, d[1][1], d[1][2], d[1][3], d[1][4], -1)]
        # Neutral K0 K0bar pair: the S is CP even, so the CP-even superposition gives EQUAL
        # proportions of K_S K_S and K_L K_L (prompt.md line 106).  K_S -> pi+ pi- is
        # essentially 100 %; K_L is invisible most of the time, and its visible fraction keeps
        # the charged-pion model used before.  Both are emitted from the parent vertex.
        d = two_body(ms, M_K0L, M_K0L, rng)
        if d is None:
            return []
        out = []
        kaon_is_short = rng.random() < 0.5
        for i, sg in ((0, +1), (1, -1)):
            if kaon_is_short or rng.random() < 0.80:
                sub = _K0L_charged(M_K0L, rng)
                for (m, e, px, py, pz) in boost(sub, d[i][1], d[i][2:5]):
                    out.append((m, e, px, py, pz, sg))
        return out
    if channel == 'multi':
        x = rng.random()
        if x < 1.0 / 3.0:
            masses, charges = [M_PI] * 4, [+1, +1, -1, -1]
        elif x < 2.0 / 3.0:
            masses, charges = [M_PI, M_PI, M_PI0, M_PI0], [+1, -1, 0, 0]
        else:
            masses, charges = [M_ETA, M_ETA, M_PI, M_PI], [0, 0, +1, -1]
        if ms < sum(masses):
            masses, charges = [M_PI] * 4, [+1, +1, -1, -1]
            if ms < sum(masses):
                return []
        d = nbody(ms, masses, rng)
        if d is None:
            return []
        return [(d[i][0], d[i][1], d[i][2], d[i][3], d[i][4], charges[i]) for i in range(len(d))]
    if channel in ('gg', 'ss', 'cc', 'bb'):
        e_jet = ms / 2.0
        mf = 2.25 if channel == 'gg' else 1.0
        u = _isotropic_dir(rng)
        for sgn in (+1, -1):
            ax = (sgn * u[0], sgn * u[1], sgn * u[2])
            jet = _rotate_jet(jet_fragment(e_jet, rng, mf), ax, rng)
            for (m, e, px, py, pz) in jet:
                out.append((m, e, px, py, pz, +1 if rng.random() < 0.5 else -1))
        return out
    return []       # gamma gamma : no charged daughters


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
            }

    def sample(self, ms, channel, rng):
        """One random S -> channel decay as [(m, e, px, py, pz, charge), ...] in the S rest frame."""
        c = self.chan.get(channel)
        if c is None:
            return None
        i = int(np.argmin(np.abs(c['mass_nodes'] - ms)))
        node = float(c['mass_nodes'][i])
        if abs(node - ms) > 0.20 * max(ms, node):
            return None            # channel closed at this mass: the library has no entry
        # Keys are mass*1e6 + sample_index with sample_index < NSAMPLE_MAX.  The upper bound
        # must be the NEXT MASS NODE, not node*1e6 + 1e6: in key units +1e6 is +1 GeV, which
        # swept in every node within 1 GeV above and mixed their spectra into this one.
        nodes = c['mass_nodes']
        nxt = float(nodes[i + 1]) if i + 1 < len(nodes) else node + 1.0
        lo = int(np.searchsorted(c['key'], node * 1e6 - 0.5, side='left'))
        hi = int(np.searchsorted(c['key'], nxt * 1e6 - 0.5, side='left'))
        # start/counts are indexed by (mass, sample) group, not by row: convert the row range
        glo = int(np.searchsorted(c['start'], lo, side='left'))
        ghi = int(np.searchsorted(c['start'], hi, side='left'))
        if ghi <= glo:
            return None
        j = rng.randrange(glo, ghi)
        s, n = int(c['start'][j]), int(c['counts'][j])
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
    good = []
    for (m, e, px, py, pz, ch) in daughters_lab:
        if not dk.passes_birth_kinematics(px, py, pz, e, dk.BIRTH_CUTS['minE'], dk.BIRTH_CUTS['minPT']):
            continue
        acc, status, front, end = dk.transport_track(x, y, z, px, py, pz, e, ch)
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
        rest = s_decay_products(ms, names[idx], rng)
        if not rest:
            # decay into no reconstructable product: it contributes zero response, but it MUST
            # stay in the denominator.  Dividing by the number of non-empty decays instead
            # would give a conditional efficiency and overestimate the visible response.
            continue
        chg = [a[5] for a in rest]
        lab5 = boost([(a[0], a[1], a[2], a[3], a[4]) for a in rest], e, (px, py, pz))
        lab = [(b[0], b[1], b[2], b[3], b[4], c) for b, c in zip(lab5, chg)]
        for k, zd in enumerate(zgrid):
            tot[k] += _passes_pair(lab, x0 + (zd - z0) * dirx, y0 + (zd - z0) * diry, zd)
    if ndecay <= 0:
        return np.zeros(len(zgrid))
    return tot / ndecay


def lifetime_probability(ms, p4, z0_m, a_z, ctau_m, zgrid=ZGRID):
    """p = int_{zgrid} dz (m/(c tau pz)) exp(-(z-z0) m/(c tau pz)) A(z),  z > z0.

    Derived from p = int dl/lambda exp(-l/lambda) A(l) with l = (z-z0) p/pz and
    lambda = (p/m) c tau, which gives dl/lambda = dz m/(c tau pz).
    """
    e, px, py, pz = p4
    if pz <= 0 or ctau_m <= 0:
        return 0.0
    alpha = ms / (ctau_m * pz)
    w = np.zeros(len(zgrid))
    m = zgrid > z0_m
    w[m] = alpha * np.exp(-(zgrid[m] - z0_m) * alpha)
    return float(np.trapz(w * a_z, zgrid))
