# -*- coding: utf-8 -*-
"""
Main driver:  s(m_S, theta^2) for BC5 at LHCb downstream.

Production channels
  * mesons (B+, B0, Bs through FCNC mixing; K+, K0L through FCNC mixing) taken from the
    Pythia8 SoftQCD:inelastic minimum-bias parent sample; the event OR is evaluated inside
    each simulated pp collision.
  * h -> SS with sigma(pp->h) from the LHC HXSWG and alpha fixed by BR(h->SS) = 0.01.

Signal definition (prompt.md): an event counts once if at least one S in it satisfies
  1 < z_dec < 2.5 m,  2 < eta_S < 5,  and has an opposite-charge daughter pair with
  E > 5 GeV, pT > 0.5 GeV that passes the magnetic-field + SciFi envelope selection.
An overall event-level reconstruction efficiency of 0.4 is applied once.
"""
from __future__ import annotations
import json, math, os, random, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', 'widths')))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', 'response')))
import bc5_lib as L
import bc5_production as P
import scalar_widths as sw

RUN = os.path.abspath(os.path.join(HERE, '..', '..'))
SAMPLES = os.environ.get('BC5_SAMPLES', os.path.join(RUN, 'work', 'magnus', 'samples'))
OUT = os.path.join(RUN, 'outputs', 'data')
LOGDIR = os.path.join(RUN, 'logs')
os.makedirs(OUT, exist_ok=True)

MASS_MIN, MASS_MAX = 0.021, 62.49
THETA2_MIN, THETA2_MAX = 1e-20, 1e-3
# Physical threshold / inflection masses that a purely log-spaced grid misses.  Measured with
# extra real scans at 28 off-grid masses: the 52-point log grid put the drawn long-lived edge
# up to 2.93x too high at 0.2118 GeV (2 m_mu), up to 1.96x at 0.9884 GeV (2 m_K) and up to 1.79x
# at 1.0 GeV, and the short-lived edge up to 4.6x too high in the 0.95-1.0 GeV window.  These
# nodes are merged into every log grid so that the drawn curve turns the corner where the
# physics does; see logs/v3_plan.md section 4.
THRESHOLD_MASSES = [
    0.2118,    # 2 m_mu   : S -> mu mu opens
    0.2802,    # 2 m_pi   : S -> pi pi opens
    0.3541,    # K+  -> pi+ S closes
    0.3626,    # K0L -> pi0 S closes
    0.9500,    # onset of the K K turn-on
    0.9884,    # 2 m_K    : S -> K K opens
    1.0000,    # requested 1 GeV landmark
    1.1000,    # f0(980) region
    1.2000,    # f0(980) region
    2.0000,    # hadronic <-> perturbative switch (PERT_THRESHOLD)
    2.2000,
    2.6834,    # B_s -> SS two-body closes
    3.5567,    # 2 m_tau  : S -> tau tau opens
    3.7579,    # 2 m_D    : S -> c cbar opens
    7.5000,
    10.5834,   # 2 m_b    : S -> b bbar opens
]
_NODE_TOL = 0.006     # ln(m) separation below which two nodes are the same point


def merge_mass_grid(grid, extra=None, tol=_NODE_TOL):
    """Union of a log grid and the physical threshold nodes.

    A threshold node wins over a nearby arbitrary grid point (0.3541 replaces 0.3533), but two
    nodes that are distinct physical features (0.9884 = 2 m_K and 1.0 GeV) both survive because
    they are 1.2 % apart, i.e. more than tol.
    """
    extra = THRESHOLD_MASSES if extra is None else extra
    keep = []
    for m in sorted(grid):
        if any(abs(math.log(m) - math.log(t)) < tol for t in extra):
            continue                       # a physical node sits on top of it
        keep.append(m)
    nodes = []
    for t in sorted(extra):
        if MASS_MIN <= t <= MASS_MAX and not any(abs(math.log(t) - math.log(k)) < tol for k in nodes):
            nodes.append(t)
    return np.array(sorted(keep + nodes))
LUMI = [25.0, 300.0]                     # fb^-1
PB_PER_FB = 1e3                          # 1 fb^-1 = 1e3 pb^-1
# number of generated parent events (set from the Magnus job summary before the scan)
_np_rng = np.random.default_rng(int(os.environ.get('BC5_NPRNG_SEED', '20261016')))
SAMPLE_B = os.environ.get('BC5_SAMPLE_B', 'bbar_heavy.txt')
SAMPLE_K = os.environ.get('BC5_SAMPLE_K', 'minbias_light.txt')
SAMPLE_H = os.environ.get('BC5_SAMPLE_H', 'higgs.txt')
# The number of generated events per sample must match the files on disk, otherwise the MC
# normalisation is wrong.  samples/normalisation.json records it (n_gen = max event id + 1),
# so the defaults are correct for this run; the BC5_NGEN_* variables still override.
_NORM_FILE = os.path.join(SAMPLES, 'normalisation.json')
try:
    with open(_NORM_FILE) as _f:
        _NORM = json.load(_f)
except (OSError, ValueError):
    _NORM = {}

N_GEN_MINBIAS = int(os.environ.get('BC5_NGEN_MINBIAS', _NORM.get('K', {}).get('n_gen', 200000)))
N_GEN_BBAR = int(os.environ.get('BC5_NGEN_BBAR', _NORM.get('B', {}).get('n_gen', 373659)))
N_GEN_HIGGS = int(os.environ.get('BC5_NGEN_HIGGS', _NORM.get('H', {}).get('n_gen', 200000)))
# pp -> b bbar total cross section.  v3 uses 560 ub: LHCb measured sigma(pp -> H_b X) in
# 2 < eta < 5 at 13 TeV as 144 +- 1 +- 21 ub (PRL 118, 052002); our own parent sample has 25.41 %
# of the b hadrons in that window, so sigma_bb(total) = 144/0.2541 = 567 +- 83 ub, and the same
# procedure at 7 TeV reproduces the known ~290 ub.  The previous 500 ub was 12 % low; the
# chi_bb = 8.6e-3 of Boiarska's tab:meson-amounts implies 684 ub, which is 21 % high.
# See logs/v3_plan.md section 2.
SIGMA_BB_PB = float(os.environ.get('BC5_SIGMA_BB_PB',
                              _NORM.get('B', {}).get('sigma_pb', 5.6e8)))   # 560 ub, pp -> b bbar
SIGMA_INEL_PB = float(os.environ.get('BC5_SIGMA_INEL_PB',
                                     _NORM.get('K', {}).get('sigma_pb', 7.95e10)))  # 79.5 mb
# Parent-sample yield normalisation, DERIVED FROM THE SAMPLE FILE ITSELF (no hard-coded fudge).
#
# The B estimator counts the b mesons that the parent file actually contains, so the implied
# B-meson production rate is sigma_bb x (rows per generated event) / 2.  The file must therefore
# carry every B+-/B0/B_s that the generator produced.  A dedicated Pythia 8.316 job with these
# exact generation settings (Magnus 3a9a4a80bb8a162c) shows that the generator produces
#
#     B+-/B0/B_s per pp -> b bbar event : 1.918   (status 91-92 + 83-84 + 81/82/85/86)
#     of which the old parents.cc filter (status 91/92 only) wrote : 1.318
#
# and that the dropped subset has the same <E> (110.8 vs 110.4 GeV), <p> (110.0 vs 109.5 GeV)
# and <|eta|> (2.711 vs 2.717).  The missing hadrons are therefore real and kinematically
# indistinguishable, which is the 27 % b-hadron deficit of the production audit (F1).
#
# The correction is computed at run time as REF_B_MESON_YIELD / (rows per generated event), so it
# is automatically 1.0 for a sample that contains everything (Magnus job f2d79005260dd1b2 writes
# the 83/84 lines) and is printed loudly either way.  BC5_B_YIELD_CORR overrides it for tests.
REF_B_MESON_YIELD = 1.918
B_YIELD_CORR = float(os.environ.get('BC5_B_YIELD_CORR', '0'))     # 0 = derive from the sample
# Bump this whenever the production model changes; a scan written by a different model is never
# merged into (see --only-masses).
# v3.4: the physics-review fixes of 2026-10-24 - isotropic hadronisation library, exact n-body
# phase space, neutral-kaon flight, isospin 4pi multi-meson composition, photon channel, exact
# lifetime integral, |M(q^2)|^2 differential for the quartic three-body modes, Winkler's
# inclusive B -> X_s S prescription, Higgs pre-decay copy, K0 flavour projection (F01-F13).
# v3.5: second review round, 2026-10-25 - Q_b^2 = 1/9 in the photon loop (N01), one decay vertex
# per neutral kaon (N02), exact four-body rejection envelope (N04), numerically stable lifetime
# slope term (N05), no mass-dilating library fallback (N06), per-process Higgs weighting, input
# fingerprints.
# v3.6: B -> X_s S back to BOIARSKA's exclusive channels (xi_sb = 3.6e-4, tab:BR) and no B_s
# mixing channel, matching both the paper (B_s is excluded there) and the reference tool
# SensCalc (SelectedBrChoice="Exclusive").  The Winkler inclusive formula is gone.
MODEL_TAG = 'v3.6-2026-10-25'


# --------------------------------------------------------------------------------------
# Input fingerprint (review 2026-10-25, item 4.8): a scan must be reproducible from the files it
# names, so every scan records the sha256 of the samples, the decay library and the four source
# files.  A partial re-run refuses to merge into a scan whose fingerprint differs.
_FINGERPRINT_TARGETS = ['work/analysis/bc5_lib.py', 'work/analysis/bc5_production.py',
                        'work/analysis/run_scan.py', 'work/widths/scalar_widths.py']


def _sha256(path, nbytes=1 << 20):
    import hashlib
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        while True:
            b = fh.read(nbytes)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def inputs_fingerprint():
    """{relative path: sha256} for the samples, the decay library and the production sources."""
    import hashlib
    out = {}
    for rel in _FINGERPRINT_TARGETS:
        p = os.path.join(RUN, rel)
        if os.path.isfile(p):
            out[rel] = _sha256(p)
    for name in (SAMPLE_B, SAMPLE_K, SAMPLE_H):
        for cand in (name, name + '.gz'):
            p = os.path.join(SAMPLES, cand)
            if os.path.isfile(p):
                out['samples/' + cand] = _sha256(p)
                break
    libdir = os.environ.get('BC5_HADRON_LIB', os.path.join(RUN, 'work', 'magnus', 'hadron_lib'))
    for ch in ('gg', 'ss', 'cc', 'bb', 'tautau'):
        for ext in ('.txt.gz', '.txt'):
            p = os.path.join(libdir, 'hadron_lib_%s%s' % (ch, ext))
            if os.path.isfile(p):
                out['hadron_lib/hadron_lib_%s%s' % (ch, ext)] = _sha256(p)
                break
    return out


_SAMPLE_CACHE = {}


def load_sample(name):
    """Read a parents .txt(.gz) file -> dict of arrays (cached)."""
    if name in _SAMPLE_CACHE:
        return _SAMPLE_CACHE[name]
    base = os.path.join(SAMPLES, name)
    path = None
    for cand in (base, base + '.gz'):
        if os.path.exists(cand):
            path = cand
            break
    if path is None:
        _SAMPLE_CACHE[name] = None
        return None
    ev, pdg, px, py, pz, e, x, y, z, t, code = ([] for _ in range(11))
    has_status = False
    opener = open
    if path.endswith('.gz'):
        import gzip
        opener = gzip.open
    with opener(path, 'rt') as f:
        for line in f:
            if not line or line[0] == '#':
                continue
            v = line.split()
            if len(v) < 10:
                continue
            ev.append(int(v[0])); pdg.append(int(v[1]))
            px.append(float(v[2])); py.append(float(v[3])); pz.append(float(v[4])); e.append(float(v[5]))
            x.append(float(v[6])); y.append(float(v[7])); z.append(float(v[8])); t.append(float(v[9]))
            # optional 11th column: the Pythia STATUS of the record (needed to pick the
            # pre-decay copy of a resonance, see the Higgs block of run_mass)
            if len(v) >= 11:
                code.append(int(v[10])); has_status = True
            else:
                code.append(0)
    out = {'ev': np.asarray(ev, dtype=np.int64), 'pdg': np.asarray(pdg, dtype=np.int32),
           'px': np.asarray(px), 'py': np.asarray(py), 'pz': np.asarray(pz), 'e': np.asarray(e),
           'x': np.asarray(x), 'y': np.asarray(y), 'z': np.asarray(z), 't': np.asarray(t)}
    if has_status:
        out['status'] = np.asarray(code, dtype=np.int32)
    _SAMPLE_CACHE[name] = out
    return out


def sample_meson_decay_S(ms, mes, rng):
    """Return a list of (E,px,py,pz,x_mm,y_mm,z_mm) S candidates from one parent meson."""
    e, px, py, pz = mes['e'], mes['px'], mes['py'], mes['pz']
    m = math.sqrt(max(e * e - px * px - py * py - pz * pz, 1e-9))
    p3 = math.sqrt(px * px + py * py + pz * pz)
    if p3 <= 0:
        return []
    out = []
    if mes['kind'] == 'B':
        brs = P.br_B_to_S(ms, mes['bmeson'])
        names = [n for n, b in brs.items() if b > 0]
        if not names:
            return []
        w = np.array([brs[n] for n in names]); w = w / w.sum()
        idx = int(np.searchsorted(np.cumsum(w), rng.random()))
        idx = min(idx, len(names) - 1)
        mXp = P.B_CHANNEL_MASS[names[idx]]
        d = L.two_body(m, mXp, ms, rng)
        if d is None:
            return []
        s4 = d[1]
        ctau = sw.HBAR_C / P.GAM_BP if mes['bmeson'] == 'B+' else (
            sw.HBAR_C / P.GAM_B0 if mes['bmeson'] == 'B0' else sw.HBAR_C / P.GAM_BS)
    else:
        which = mes['kaon']
        br = P.br_K_to_S(ms, which)
        if br <= 0:
            return []
        mP = 0.13957039 if which == 'K+' else 0.1349768
        d = L.two_body(m, mP, ms, rng)
        if d is None:
            return []
        s4 = d[1]
        ctau = sw.HBAR_C / (P.GAM_KP if which == 'K+' else P.GAM_KL)
    # boost S to the lab
    s_lab = L.boost([s4], e, (px, py, pz))[0]
    # parent decay vertex: displaced by beta*gamma*c*tau with t ~ Exp(1)
    u = -math.log(max(rng.random(), 1e-12))
    # ctau is in METRES (HBAR_C in GeV m) while the sample vertices are in MILLIMETRES, and
    # the caller converts the whole vertex by 1e-3.  The displacement must therefore be
    # converted to mm here, otherwise it is a factor 1000 too small and disagrees with the
    # lifetime integral, which uses the vertex directly.
    disp = (p3 / m) * ctau * u * 1e3
    xd = mes['x'] + disp * px / p3
    yd = mes['y'] + disp * py / p3
    zd = mes['z'] + disp * pz / p3
    return [(s_lab[1], s_lab[2], s_lab[3], s_lab[4], xd, yd, zd)]


def sample_quartic_decay_S(ms, mes, kind, spec, rng, channel=None):
    """Two scalars from X -> SS (kind '2b') or X -> X' SS (kind '3b'), boosted to the lab.

    F07: the three-body case now draws the S S invariant mass with the |M_XX'(q^2)|^2 weight of
    the rate integral (bc5_production.sample_dipion_mass) instead of an unweighted nbody(), so
    the generated events are the differential version of the same rate.  The S S rest frame is
    isotropic, which is exactly the assumption of the rate integral.
    """
    e, px, py, pz = mes['e'], mes['px'], mes['py'], mes['pz']
    m = math.sqrt(max(e * e - px * px - py * py - pz * pz, 1e-9))
    if kind == '2b':
        d = L.two_body(m, ms, ms, rng)
        if d is None:
            return []
    else:
        if m < 2 * ms + spec:
            return []
        q = P.sample_dipion_mass(ms, m, spec, channel, rng)
        if q is None:
            return []
        dss = L.two_body(q, ms, ms, rng)
        dpar = L.two_body(m, q, spec, rng)
        if dss is None or dpar is None:
            return []
        # the dipion system is particle 1 of the X rest frame; boost the two scalars into it
        d = [(a[0], a[1], a[2], a[3], a[4]) for a in
             L.boost([(a[0], a[1], a[2], a[3], a[4]) for a in dss], dpar[0][1], dpar[0][2:5])]
    # The parent meson flies before it decays, exactly as in the mixing channels.  Without this
    # the quartic scalars would all start at the production point and their displacement would
    # be lost (review P2).
    pmod = math.sqrt(px * px + py * py + pz * pz)
    ctau = None
    if mes.get('kaon') is not None:
        ctau = sw.HBAR_C / (P.GAM_KP if mes['kaon'] == 'K+' else P.GAM_KL)
    elif mes.get('bmeson') is not None:
        ctau = sw.HBAR_C / {'B+': P.GAM_BP, 'B0': P.GAM_B0, 'Bs': P.GAM_BS}[mes['bmeson']]
    if ctau is None or pmod <= 0:
        xd, yd, zd = mes['x'], mes['y'], mes['z']
    else:
        disp = (pmod / m) * ctau * (-math.log(max(rng.random(), 1e-12))) * 1e3
        xd = mes['x'] + disp * px / pmod
        yd = mes['y'] + disp * py / pmod
        zd = mes['z'] + disp * pz / pmod
    out = []
    for s4 in d:
        sl = L.boost([s4], e, (px, py, pz))[0]
        out.append((sl[1], sl[2], sl[3], sl[4], xd, yd, zd))
    return out


def decay_vertex_and_selection(ms, s, rng, ndecay=10):
    """Return (p4, v0_mm, A_z) for one S, or None if eta_S is outside the selection."""
    e, px, py, pz, x, y, z = s
    if pz <= 0:
        return None
    pmod = math.sqrt(px * px + py * py + pz * pz)
    if pmod <= abs(pz):
        return None
    eta = 0.5 * math.log((pmod + pz) / max(pmod - pz, 1e-300))    # pseudorapidity
    if not (L.ETA_S_MIN < eta < L.ETA_S_MAX):
        return None
    a = L.acceptance_profile(ms, (e, px, py, pz), (x, y, z), rng, ndecay=ndecay)
    return (e, px, py, pz), (x, y, z), a


# --------------------------------------------------------------------------------------
def run_mass(ms, nsamp_per_channel, rng, ndecay=10, verbose=True):
    """Build the response table for one scalar mass."""
    t0 = time.time()
    rows = {'chan': [], 'ev': [], 'e': [], 'px': [], 'py': [], 'pz': [],
            'x': [], 'y': [], 'z': [], 'A': [], 'w': [], 'par': [], 'hid': [], 'kq': []}
    # For every channel we record how many parent candidates the sample offers and how many of
    # them were actually evaluated.  The efficiency must be normalised by the number of PARENTS
    # SAMPLED and then scaled to the full sample, not divided by the number of generated events:
    # otherwise subsampling the parents silently rescales the whole signal.
    counts = {'B': [0, 0], 'K': [0, 0], 'H': [0, 0]}
    # Collision ids of every SAMPLED parent.  Parents whose scalars all fail the selection leave
    # no row in the table, but their collision still contributes a zero response; it must stay in
    # the sample so that the Monte-Carlo variance is estimated over the full sample.
    samp_ev = {'B': [], 'K': [], 'H': []}

    def add(channel, ev, s, w=1.0, par=0, hid=0, kq=0.0):
        """One scalar row of one decay history of one parent particle.

        Bookkeeping needed to get the probability algebra right (audit items 3.1 and 3.3):
          par  -- index of the parent particle, unique inside the channel; shares lifetime/vertex
          hid  -- index of the decay history of that parent.  Histories of the same parent are
                  mutually EXCLUSIVE (a parent decays once), so they are added, not OR-ed.
          kq   -- 1.0 if this history's rate carries a factor theta^2 (mixing-induced
                  production), 0.0 if it does not (quartic-induced production, or h -> SS).
        The weight w is the branching ratio of the history; one channel is sampled per parent
        with the normalised branching ratios, so the expectation over that sampling equals
        BR_a p_a / sum BR_a and the weight sum BR_a restores the absolute rate.
        """
        r = decay_vertex_and_selection(ms, s, rng, ndecay=ndecay)
        if r is None:
            return
        p4, v0, a = r
        rows['w'].append(w)
        rows['par'].append(par); rows['hid'].append(hid); rows['kq'].append(kq)
        rows['chan'].append(channel); rows['ev'].append(ev)
        rows['e'].append(p4[0]); rows['px'].append(p4[1]); rows['py'].append(p4[2]); rows['pz'].append(p4[3])
        # decay_vertex_and_selection returns the vertex in the sample's MILLIMETRES.
        # per_scalar_p integrates over z in METRES, so convert once here, at the interface.
        rows['x'].append(v0[0] * 1e-3); rows['y'].append(v0[1] * 1e-3)
        rows['z'].append(v0[2] * 1e-3)
        rows['A'].append(a)

    # --- B mesons (from the minimum-bias sample) --------------------------------------
    smp = load_sample(SAMPLE_B)
    if smp is not None:
        sel = np.isin(smp['pdg'], [521, -521, 511, -511, 531, -531])
        idxs = np.where(sel)[0]
        counts['B'][0] = len(idxs)
        if len(idxs) > nsamp_per_channel:
            idxs = _np_rng.choice(idxs, size=nsamp_per_channel, replace=False)
        counts['B'][1] = len(idxs)
        samp_ev['B'] = smp['ev'][idxs]
        for ip, i in enumerate(idxs):
            pdg = smp['pdg'][i]
            bm = {521: 'B+', -521: 'B+', 511: 'B0', -511: 'B0', 531: 'Bs', -531: 'Bs'}[pdg]
            mes = {'kind': 'B', 'bmeson': bm, 'x': smp['x'][i], 'y': smp['y'][i], 'z': smp['z'][i],
                   'e': smp['e'][i], 'px': smp['px'][i], 'py': smp['py'][i], 'pz': smp['pz'][i]}
            wB = sum(P.br_B_to_S(ms, bm).values())
            for s in sample_meson_decay_S(ms, mes, rng) or []:
                add('B', smp['ev'][i], s, wB, par=ip, hid=0, kq=1.0)
            # quartic-coupling channels X -> SS / X -> X' SS (two scalars per parent).
            # One mode is sampled with its normalised branching ratio and the summed branching
            # ratio is used as the weight, exactly as for the mixing channels.
            qmodes = [q for q in P.quartic_channels(ms) if q[0] == bm]
            if qmodes:
                wq = sum(q[3] for q in qmodes)
                u = rng.random(); c = 0.0; sel = qmodes[-1]
                for q in qmodes:
                    c += q[3] / wq
                    if u <= c:
                        sel = q; break
                qs = sample_quartic_decay_S(ms, mes, sel[1], sel[2], rng, channel=sel[4])
                for s in qs or []:
                    # quartic production is proportional to alpha^2, not to theta^2
                    add('B', smp['ev'][i], s, wq, par=ip, hid=1, kq=0.0)

    # --- kaons ------------------------------------------------------------------------
    smp = load_sample(SAMPLE_K)
    if smp is not None:
        k_par_off = 0
        for pdg, which in ((321, 'K+'), (-321, 'K+'), (311, 'K0L'), (-311, 'K0L')):
            idxs = np.where(smp['pdg'] == pdg)[0]
            counts['K'][0] += len(idxs)
            if len(idxs) > nsamp_per_channel:
                idxs = _np_rng.choice(idxs, size=nsamp_per_channel, replace=False)
            counts['K'][1] += len(idxs)
            samp_ev['K'].extend(smp['ev'][idxs])
            # F09-2: a PDG 311 / -311 is a strong-interaction FLAVOUR state, not a mass
            # eigenstate.  The non-coherent superposition is 1/2 K_L + 1/2 K_S; only the K_L
            # component can decay to pi0 S (K_S -> pi S is CP suppressed, see br_K_to_S), so the
            # K_L rate must be weighted by |<K_L|K^0>|^2 = 1/2 and the K_S component contributes
            # nothing.  Applying the full K_L rate to every 311 double-counts this source.
            # (Same projection for the quartic modes: only the K0L entries are used, and the
            # K0S half, which is 560 times smaller, is dropped - documented, negligible.)
            k0_proj = 0.5 if which == 'K0L' else 1.0
            for ip, i in enumerate(idxs):
                ip += k_par_off
                mes = {'kind': 'K', 'kaon': which, 'x': smp['x'][i], 'y': smp['y'][i], 'z': smp['z'][i],
                       'e': smp['e'][i], 'px': smp['px'][i], 'py': smp['py'][i], 'pz': smp['pz'][i]}
                wK = P.br_K_to_S(ms, which) * k0_proj
                for s in sample_meson_decay_S(ms, mes, rng) or []:
                    add('K', smp['ev'][i], s, wK, par=ip, hid=0, kq=1.0)
                # quartic modes K -> SS and K -> SS pi (two scalars per parent).  They are
                # MUTUALLY EXCLUSIVE histories of the same parent, so each mode needs its own
                # hid: they are then added, not OR-ed (F09-3).  Sharing hid = 1 over-estimated
                # this (already tiny) contribution by up to 1.78x.
                for imode, (_k, mode, wq) in enumerate(P.quartic_kaon_channels(ms, which)):
                    kind = '2b' if mode == 'SS' else '3b'
                    spec = 0.0 if mode == 'SS' else (P.M_PI0 if mode == 'SSpi0' else P.M_PI_CHARGED)
                    for s in sample_quartic_decay_S(ms, mes, kind, spec, rng) or []:
                        add('K', smp['ev'][i], s, wq * k0_proj, par=ip, hid=1 + imode, kq=0.0)
            k_par_off += len(idxs)

    # --- Higgs ------------------------------------------------------------------------
    smp = load_sample(SAMPLE_H)
    nHiggsPair = 0
    if smp is not None and ms < P.M_H / 2.0:
        # F06: Pythia's event record contains several copies of the same Higgs.  Their order is
        # -22 (hard process, pT = 0 because it recoils against the whole event) -> -44 (parton
        # shower copies) -> -62 (the last pre-decay copy, i.e. the physical Higgs that decays).
        # The old code kept the FIRST record of every event, which is an early shower copy:
        # <pT> = 55.4 against 60.2 GeV, and 3-5 % on the Higgs-channel efficiency.  If the sample
        # carries the Pythia status code (11th column) the -62 copy is selected explicitly;
        # otherwise the LAST copy of each event is used, which is the -62 copy in the record
        # order (verified against a Pythia 8.316 listing, see work/review_check/model).
        if 'status' in smp:
            cand = np.where(np.abs(smp['status']) == 62)[0]
            if len(cand) == 0:
                raise SystemExit('Higgs sample has a status column and no pre-decay copy')
            ev_c = smp['ev'][cand]
            order = np.argsort(ev_c, kind='stable')
            cand, ev_c = cand[order], ev_c[order]
            last = np.ones(len(cand), dtype=bool)
            last[:-1] = ev_c[1:] != ev_c[:-1]
            idxs = np.sort(cand[last])
        else:
            _, first_idx, cnt = np.unique(smp['ev'], return_index=True, return_counts=True)
            idxs = np.sort(first_idx + cnt - 1)
        counts['H'][0] = len(idxs)
        if len(idxs) > nsamp_per_channel:
            idxs = _np_rng.choice(idxs, size=nsamp_per_channel, replace=False)
        counts['H'][1] = len(idxs)
        # the H rows are keyed by the sample ROW index (one Higgs per row), so the collision
        # key must be the row index too, not smp['ev']
        samp_ev['H'] = np.asarray(idxs, dtype=np.int64)
        # Per-process reweighting of the Higgs sample (review 2026-10-25, section 3 of
        # efficiency_pythia_followup): the Pythia mixture is ggF-heavy (79.8 % against 87.6 % in
        # the target cross sections) and its fermion-annihilation component stands in for bbH at
        # 2.7 % instead of 0.43 %.  Each process group is reweighted by
        # (target fraction)/(sample fraction), so the sample reproduces the SIGMA_H mixture and
        # not only its total.  tHj/tHW/tHb (0.17 % of the total) have no generator code here, so
        # the mean weight is 0.9983; that deficit is documented rather than renormalised away.
        wH = np.ones(len(idxs))
        grp = np.array([P.HIGGS_PROC_GROUP.get(int(c), 'other') for c in smp['pdg'][idxs]])
        for g_ in sorted(set(grp.tolist())):
            if g_ == 'other' or g_ not in P.SIGMA_H:
                continue
            m_ = grp == g_
            frac = float(m_.sum()) / float(len(idxs))
            if frac > 0.0:
                wH[m_] = P.SIGMA_H[g_] / (P.SIGMA_H_TOTAL * frac)
        for ip, i in enumerate(idxs):
            # the S pair is produced with the Higgs' own invariant mass (the sample carries a
            # Breit-Wigner of width 4.1 MeV), not the fixed nominal P.M_H
            mH = math.sqrt(max(smp['e'][i] * smp['e'][i] - smp['px'][i] * smp['px'][i]
                               - smp['py'][i] * smp['py'][i] - smp['pz'][i] * smp['pz'][i], 1e-12))
            d = L.two_body(mH, ms, ms, rng)
            if d is None:
                continue
            pair = L.boost(d, smp['e'][i], (smp['px'][i], smp['py'][i], smp['pz'][i]))
            # the two scalars of one Higgs form ONE history: their OR is taken, and the
            # history is identified by the parent index, not by adjacency in the row array
            for s4 in pair:
                add('H', i, (s4[1], s4[2], s4[3], s4[4], smp['x'][i], smp['y'][i], smp['z'][i]),
                    w=float(wH[ip]), par=ip, hid=0, kq=0.0)
            nHiggsPair += 1

    if verbose:
        print('  m_S=%.4f  N_S=%d  (B=%d K=%d H=%d)  %.1f s' % (
            ms, len(rows['chan']), rows['chan'].count('B'), rows['chan'].count('K'),
            rows['chan'].count('H'), time.time() - t0))
    out = {k: (np.array(v) if k != 'A' else np.array(v)) for k, v in rows.items()}
    out['nHiggsPair'] = nHiggsPair
    out['counts'] = counts
    out['samp_ev'] = {k: np.asarray(v, dtype=np.int64) for k, v in samp_ev.items()}
    return out


# --------------------------------------------------------------------------------------
def per_scalar_p(ms, tab, theta2, pz, z0, A):
    """p_k = decay-in-window probability x daughter acceptance, for a vector of scalars.

    F10: the exponential is integrated exactly and A(z) piecewise linearly on ZGRID (bc5_lib
    .p_window).  The old fixed-grid trapezoid was biased by up to a factor 1.5 when
    alpha * dz >> 1 and could not represent a birth point between grid nodes.
    """
    if len(pz) == 0:
        return np.zeros(0)
    ctau = sw.HBAR_C / (sw.gamma_total(ms) * theta2)          # m
    return L.p_window(ms, pz, z0, A, ctau)


def channel_efficiency(tab, p, chan, n_gen, theta2):
    """Per-pp-collision efficiency, with the probability algebra of audit items 3.1 and 3.3.

    Three levels, in this order:
      1. scalar row -> decay history.  A history's expected response is
              o_h = 1 - prod_rows (1 - w p theta^(2 kq)),
         so the branching ratio is applied ONCE per history, after the OR over its scalars.
         kq = 1 for mixing-induced histories (rate proportional to theta^2) and 0 for
         quartic-induced ones and for h -> SS (rate independent of theta).
      2. history -> parent.  The histories of one parent are mutually EXCLUSIVE, because a
         parent particle decays once, so their expected responses ADD:
              p_par = sum_h o_h .
      3. parent -> collision: o_ev = 1 - prod_parents (1 - p_par).

    Normalisation: the sample offers counts[chan][0] parents of which counts[chan][1] were
    evaluated, and spans n_gen generated collisions, so
        eps = 0.4 * sum_events o_ev * (n_tot/n_smp) / n_gen ,
        se  = 0.4 * (n_tot/(n_smp n_gen)) * sqrt(n_coll * Var(o_ev)) .
    """
    sel = tab['chan'] == chan
    if not np.any(sel) or n_gen <= 0:
        return 0.0, 0.0
    c = tab.get('counts', {}).get(chan, [0, 0])
    n_tot, n_smp = int(c[0]), int(c[1])
    if n_smp <= 0:
        return 0.0, 0.0
    ev = tab['ev'][sel].astype(np.int64)
    par = tab['par'][sel].astype(np.int64)
    hid = tab['hid'][sel].astype(np.int64)
    kq = tab['kq'][sel]
    br = tab['w'][sel]                       # branching ratio of the history (constant inside it)
    # per-scalar probability, WITHOUT the branching ratio: the ratio must be applied once per
    # history, after the OR over that history's scalars (review P1-b).  Mixing histories carry
    # one factor theta^2; quartic histories and h -> SS carry none.
    pe = p[sel] * np.where(kq > 0.5, theta2, 1.0)
    log1mp = np.log1p(-np.clip(pe, 0.0, 1.0 - 1e-15))

    # 1. rows -> histories (group by ev, par, hid)
    order = np.lexsort((hid, par, ev))
    ev_s, par_s, hid_s, log_s = ev[order], par[order], hid[order], log1mp[order]
    br_s = br[order]
    new_h = np.ones(len(ev_s), dtype=bool)
    new_h[1:] = ((ev_s[1:] != ev_s[:-1]) | (par_s[1:] != par_s[:-1])
                 | (hid_s[1:] != hid_s[:-1]))
    start_h = np.flatnonzero(new_h)
    # -expm1 is numerically stable for very small probabilities
    o_inner = -np.expm1(np.add.reduceat(log_s, start_h))
    o_hist = br_s[start_h] * o_inner        # branch ratio applied ONCE, after the OR
    ev_h, par_h = ev_s[start_h], par_s[start_h]

    # 2. histories -> parents (exclusive histories ADD)
    new_p = np.ones(len(o_hist), dtype=bool)
    new_p[1:] = (ev_h[1:] != ev_h[:-1]) | (par_h[1:] != par_h[:-1])
    start_p = np.flatnonzero(new_p)
    p_par = np.add.reduceat(o_hist, start_p)
    ev_p = ev_h[start_p]

    # 3. parents -> collisions
    new_e = np.ones(len(p_par), dtype=bool)
    new_e[1:] = ev_p[1:] != ev_p[:-1]
    start_e = np.flatnonzero(new_e)
    o_found = -np.expm1(np.add.reduceat(np.log1p(-np.clip(p_par, 0.0, 1.0 - 1e-15)), start_e))

    # 4. restore the sampled collisions that produced no surviving scalar: their response is
    #    exactly zero, but they must stay in the sample or the variance is underestimated
    #    (review P1-c).  A collision is identified by its sample event id.
    samp = tab.get('samp_ev', {}).get(chan)
    if samp is not None and len(samp) > 0:
        all_ev = np.unique(samp)
        have = ev_p[start_e]
        o_all = np.zeros(len(all_ev))
        if len(have):
            pos = np.searchsorted(all_ev, have)
            o_all[pos] = o_found
        o_ev = o_all
    else:
        o_ev = o_found
    n_coll = len(o_ev)
    norm = (n_tot / (n_smp * n_gen)) * L.RECO_EFF
    eps = float(np.sum(o_ev)) * norm
    var_o = float(np.var(o_ev)) if n_coll > 1 else 0.0
    se = norm * math.sqrt(max(var_o, 0.0) * n_coll)
    return eps, se


def signal_per_lumi(ms, tab, theta2, lumi_fb, sigma_bb_pb, sigma_inel_pb,
                    n_gen_bbar, n_gen_minbias, n_gen_higgs):
    """Expected signal events s = 1e3 L sum_c sigma_c eps_c (prompt.md Eq. 4.2).

    Every channel is normalised by its own measured/generated cross section:
      B  : pp -> b bbar total cross section (the bbar sample is unbiased in eta)
      K  : minimum-bias inelastic cross section
      H  : sum of the SM Higgs production cross sections times BR(h -> SS)
    """
    if len(tab['pz']) == 0:
        return {'B': 0.0, 'K': 0.0, 'H': 0.0, 'eps_B': 0.0, 'eps_K': 0.0, 'eps_H': 0.0,
                'relerr_B': 0.0, 'relerr_K': 0.0, 'relerr_H': 0.0}
    p = per_scalar_p(ms, tab, theta2, tab['pz'], tab['z'], tab['A'])
    out = {}
    # The theta^2 dependence now lives inside channel_efficiency, per decay history: mixing
    # histories carry one factor of theta^2, quartic histories and h -> SS carry none.
    epsB, seB = channel_efficiency(tab, p, 'B', n_gen_bbar, theta2)
    epsB *= B_YIELD_CORR
    seB *= B_YIELD_CORR
    epsK, seK = channel_efficiency(tab, p, 'K', n_gen_minbias, theta2)
    epsH, seH = channel_efficiency(tab, p, 'H', n_gen_higgs, theta2)
    out['B'] = PB_PER_FB * lumi_fb * sigma_bb_pb * epsB
    out['K'] = PB_PER_FB * lumi_fb * sigma_inel_pb * epsK
    out['H'] = PB_PER_FB * lumi_fb * P.SIGMA_H_TOTAL * P.B_H * epsH
    out['eps_B'] = epsB
    out['eps_K'] = epsK
    out['eps_H'] = epsH
    out['eps_geom_B'] = epsB
    out['eps_geom_K'] = epsK
    out['relerr_B'] = (seB / epsB) if epsB > 0 else 0.0
    out['relerr_K'] = (seK / epsK) if epsK > 0 else 0.0
    out['relerr_H'] = (seH / epsH) if epsH > 0 else 0.0
    return out


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--nmass', type=int, default=40)
    ap.add_argument('--nsamp', type=int, default=1500)
    ap.add_argument('--ndecay', type=int, default=10)
    ap.add_argument('--seed', type=int, default=20261016)
    ap.add_argument('--mmin', type=float, default=0.3)
    ap.add_argument('--mmax', type=float, default=MASS_MAX)
    ap.add_argument('--masses', type=str, default=None,
                    help='explicit comma-separated mass list in GeV; overrides mmin/mmax/nmass')
    ap.add_argument('--no-threshold-nodes', action='store_true',
                    help='use the bare log grid without the physical threshold nodes')
    ap.add_argument('--nsamp-map', type=str, default=None,
                    help='JSON {"mass_GeV": nsamp} giving a per-mass parent sample size; '
                         'masses not listed use --nsamp.  This is what makes the sampling '
                         'adaptive: only the masses whose crossing is too noisy are re-run.')
    ap.add_argument('--force-mix', action='store_true',
                    help='allow merging into a scan produced by a different model (dangerous)')
    ap.add_argument('--only-masses', type=str, default=None,
                    help='comma-separated masses to (re-)evaluate; every other mass is copied '
                         'from the existing output file, so a re-run of a few points merges into '
                         'the production scan instead of replacing it')
    args = ap.parse_args()

    # ---- decay-width model ------------------------------------------------------------
    # Production path = Winkler (arXiv:1809.01876) only.  BC5_WIDTH_MODEL=reference_br swaps in
    # the published reference branching-ratio table for the validation scan.
    _wmodel = os.environ.get('BC5_WIDTH_MODEL', 'winkler').strip().lower()
    if _wmodel == 'reference_br':
        import reference_br_model
        sw.WIDTH_OVERRIDE = reference_br_model.widths
        print('  WIDTH MODEL: reference branching-ratio table (validation only)')
    else:
        sw.WIDTH_OVERRIDE = None
        print('  WIDTH MODEL: Winkler arXiv:1809.01876')

    # Pythia 8.316 hadronisation / decay library for S -> gg, ss, cc, bb, tau tau
    libdir = os.environ.get('BC5_HADRON_LIB', os.path.join(RUN, 'work', 'magnus', 'hadron_lib'))
    if L.load_hadron_lib(libdir):
        print('  hadronisation library loaded from %s' % libdir)
    else:
        print('  WARNING: no hadronisation library found; using the analytic fallback models')

    rng = random.Random(args.seed)
    if args.masses:
        masses = np.array(sorted(float(x) for x in args.masses.split(',') if x.strip()))
    else:
        masses = np.exp(np.linspace(math.log(args.mmin), math.log(args.mmax), args.nmass))
        if not args.no_threshold_nodes:
            n_log = len(masses)
            masses = merge_mass_grid(masses)
            print('  mass grid: %d log-spaced points + physical threshold nodes -> %d masses'
                  % (n_log, len(masses)))
    thetas = np.exp(np.linspace(math.log(THETA2_MIN), math.log(THETA2_MAX), 240))

    # --- parent-sample yield normalisation (see REF_B_MESON_YIELD) -----------------------
    global B_YIELD_CORR
    if B_YIELD_CORR <= 0.0:
        smpB = load_sample(SAMPLE_B)
        if smpB is not None and N_GEN_BBAR > 0:
            sel = np.isin(smpB['pdg'], [521, -521, 511, -511, 531, -531])
            per_event = float(sel.sum()) / float(N_GEN_BBAR)
            B_YIELD_CORR = REF_B_MESON_YIELD / per_event if per_event > 0 else 1.0
            print('  b mesons in the parent sample: %.4f per generated event '
                  '(reference %.3f) -> B-channel yield correction %.4f'
                  % (per_event, REF_B_MESON_YIELD, B_YIELD_CORR))
            if abs(B_YIELD_CORR - 1.0) < 0.02:
                print('  the sample carries the full b-meson yield; correction is 1.00')
        else:
            B_YIELD_CORR = 1.0
    else:
        print('  B-channel yield correction fixed by BC5_B_YIELD_CORR = %.4f' % B_YIELD_CORR)

    nsamp_map = {}
    if args.nsamp_map:
        with open(args.nsamp_map) as fh:
            nsamp_map = {float(k): int(v) for k, v in json.load(fh).items()}
        print('  per-mass nsamp map: %d entries' % len(nsamp_map))
    only = None
    if args.only_masses:
        only = np.array(sorted(float(x) for x in args.only_masses.split(',') if x.strip()))

    scan_name = os.environ.get('BC5_SCAN_OUT', 'signal_scan.json')
    scan_path = os.path.join(OUT, scan_name)
    res = {}
    nsamp_by_mass = {}
    if only is not None and os.path.exists(scan_path):
        with open(scan_path) as fh:
            prev = json.load(fh)
        # A partial re-run must never mix results produced by a different physics model: the
        # sigma_bb normalisation, the B_s treatment, the quartic phase space and the parent
        # sample all changed between v2 and v3, and a curve built from a mixture would be wrong
        # in a way that is invisible afterwards.
        problems = []
        if prev.get('model_tag') != MODEL_TAG:
            problems.append('model_tag %r != %r' % (prev.get('model_tag'), MODEL_TAG))
        if abs(float(prev.get('sigma_bb_pb', -1.0)) / SIGMA_BB_PB - 1.0) > 1e-3:
            problems.append('sigma_bb_pb %r != %r' % (prev.get('sigma_bb_pb'), SIGMA_BB_PB))
        if abs(float(prev.get('b_yield_corr', -1.0)) / B_YIELD_CORR - 1.0) > 1e-3:
            problems.append('b_yield_corr %r != %r' % (prev.get('b_yield_corr'), B_YIELD_CORR))
        fp_old = prev.get('inputs') or {}
        fp_new = inputs_fingerprint()
        diff = sorted(k for k in set(fp_old) | set(fp_new) if fp_old.get(k) != fp_new.get(k))
        if diff:
            problems.append('input files changed (%d): %s' % (len(diff), ', '.join(diff[:4])))
        if problems and not args.force_mix:
            raise SystemExit('refusing to merge into %s, which was produced by a different '
                             'model:\n  %s\nRe-run the full grid, or pass --force-mix.'
                             % (scan_name, '\n  '.join(problems)))
        res = {round(float(k), 6): v for k, v in prev['signal_per_fb_per_channel'].items()}
        nsamp_by_mass = {round(float(k), 6): float(v)
                         for k, v in prev.get('nsamp_by_mass', {}).items()}
        print('  merging into %s: %d masses kept, %d re-evaluated (provenance ok)'
              % (scan_name, len(res), len(only)))

    run_list = masses if only is None else only
    sigma_inel_pb = SIGMA_INEL_PB
    for ms in run_list:
        ms = round(float(ms), 6)          # canonical key, so a re-run REPLACES the old point
        n = int(nsamp_map.get(ms, args.nsamp))
        tab = run_mass(ms, n, rng, ndecay=args.ndecay)
        row = []
        for t2 in thetas:
            d = signal_per_lumi(ms, tab, float(t2), 1.0, SIGMA_BB_PB, sigma_inel_pb,
                                N_GEN_BBAR, N_GEN_MINBIAS, N_GEN_HIGGS)
            row.append(d)
        res[ms] = row
        nsamp_by_mass[ms] = n
        del tab

    out_masses = sorted(res.keys())
    with open(scan_path, 'w') as f:
        json.dump({'masses': out_masses, 'theta2': [float(t) for t in thetas],
                   'channels': ['mesons', 'H'],
                   'n_gen_minbias': N_GEN_MINBIAS, 'n_gen_higgs': N_GEN_HIGGS,
                   'nsamp_by_mass': {str(m): nsamp_by_mass.get(m) for m in out_masses},
                   'sigma_bb_pb': SIGMA_BB_PB,
                   'b_yield_corr': B_YIELD_CORR,
                   'model_tag': MODEL_TAG,
                   'inputs': inputs_fingerprint(),
                   'signal_per_fb_per_channel': {str(m): res[m] for m in out_masses}}, f)
    print('wrote %s (%d masses, model %s)' % (scan_name, len(out_masses), MODEL_TAG))


if __name__ == '__main__':
    main()