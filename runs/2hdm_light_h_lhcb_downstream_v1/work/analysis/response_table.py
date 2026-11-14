# -*- coding: utf-8 -*-
"""Parent-sample reader, per-mass response-table builder and versioned npz cache.

A response table contains ONLY the parts of the calculation that depend on the H mass
(production kinematics + detector response).  The model-dependent numbers (BRs, lifetime,
production branching ratios, cross sections) are supplied by work/analysis/models.py at
evaluate time, so ONE table serves all tanb (2HDM) or all theta^2 (BC5 regression).

Table schema (all numpy arrays; see TABLE_VERSION)
    channels     list[str]                     decay channels covered
    mH           float
    pz, z0       float64 (n_scalar,)           H lab momentum z / birth point [m]
    A            float32 (n_scalar, n_ch, 13)  per-channel A_ch(z) on twohdm_lib.ZGRID
    scal_hist    int32   (n_scalar,)           scalar -> history
    hist_event   int64   (n_hist,)             history -> collision (globally unique ids)
    hist_parent  int64   (n_hist,)             history -> parent particle
    hist_group   int8    (n_hist,)             0 = B, 1 = inelastic pi/K, 2 = h -> HH
    hist_wkey    list[str]                     history -> model weight key
    samp_events  dict {group: int64 array}     sampled collisions, including zero-response
    n_gen        dict {group: int}             generated collisions of the parent sample
    yield_corr   dict {group: float}           sample-yield correction (1.0 except BC5 B)
    h_frac       dict {process: fraction}      Pythia Higgs process mixture of the sample
    nsamp        dict {group: int}             parents sampled per channel group
    meta         dict                          provenance (sha256, config, seeds, timing, ...)
"""
from __future__ import annotations
import gzip, hashlib, json, os, sys, time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(RUN, 'work', 'response'))

import twohdm_lib as L
import models as M

TABLE_VERSION = 3   # v3: per-history parent flight length lam (vertex-integrated decay window)
CACHE_DIR = os.path.join(HERE, 'cache')
SAMPLES_DEFAULT = os.path.join(RUN, 'work', 'magnus', 'samples')
SAMPLE_FILES = {'B': 'bbar_heavy.txt.gz',
                'inel': 'minbias_pik.txt.gz',
                'h': 'higgs.txt.gz'}
DEFAULT_NSAMP = {'B': 1500, 'inel': 1500, 'h': 1500}
DEFAULT_SEED = 20261025
# Event-id / parent-id offsets per group: the three parent samples are independent, so
# collision ids and parent indices must be made globally unique BEFORE any sorting, or the
# event-level OR would merge unrelated histories (review BUG B/C).
_GROUP_OFFSET = {0: 0, 1: 10 ** 12, 2: 2 * 10 ** 12}
_PARENT_OFFSET = {0: 0, 1: 10 ** 9, 2: 2 * 10 ** 9}


def sha256_file(path, nbytes=1 << 20):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        while True:
            b = fh.read(nbytes)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def sha256_code(paths):
    h = hashlib.sha256()
    for p in paths:
        if os.path.isfile(p):
            h.update(sha256_file(p).encode())
    return h.hexdigest()


class Samples:
    """Lazy reader for the three parent samples; arrays cached for the process lifetime.

    Reads ANY whitespace-separated 10- or 11-column parents file in the BC5 format
        event pdg px py pz E x y z t [status]
    and filters by pdg afterwards, so it is agnostic to which particle species the file
    carries (B sample, pi/K minbias sample, Higgs sample with process codes in column 2).
    """

    def __init__(self, directory=SAMPLES_DEFAULT, files=None, verbose=True,
                 inel_ev_sub=100000, inel_seed=20261025, inel_n_gen_full=None):
        self.dir = os.path.abspath(directory)
        self.files = dict(SAMPLE_FILES if files is None else files)
        self.verbose = verbose
        self._cache = {}
        self.n_gen = {}
        self._norm = None
        self.paths = {}
        # The minbias pi/K sample is ~46 rows/event (520 MB gz): keep ALL rows of a random
        # EVENT subsample (default 50k of 5e5 collisions) and normalise by n_gen_sub.
        self.inel_ev_sub = inel_ev_sub
        self.inel_seed = inel_seed
        norm0 = os.path.join(self.dir, 'normalisation.json')
        self.inel_n_gen_full = inel_n_gen_full
        if self.inel_n_gen_full is None and os.path.isfile(norm0):
            try:
                with open(norm0) as fh:
                    blk = json.load(fh).get('inel', {})
                self.inel_n_gen_full = int(blk['n_gen']) if 'n_gen' in blk else None
            except (OSError, ValueError, KeyError):
                self.inel_n_gen_full = None
        norm = os.path.join(self.dir, 'normalisation.json')
        if os.path.isfile(norm):
            try:
                with open(norm) as fh:
                    self._norm = json.load(fh)
            except (OSError, ValueError):
                self._norm = None

    def path(self, key):
        if key in self.paths:
            return self.paths[key]
        base = os.path.join(self.dir, self.files[key])
        for cand in (base, base + '.gz'):
            if os.path.isfile(cand):
                self.paths[key] = cand
                return cand
        self.paths[key] = None
        return None

    def available(self):
        return {k: self.path(k) for k in self.files}

    def _norm_n_gen(self, key, path):
        if not self._norm:
            return None
        base = os.path.basename(path)
        for blk_key, grp in (('B', 'B'), ('K', 'inel'), ('H', 'h'), ('inel', 'inel'), ('h', 'h')):
            blk = self._norm.get(blk_key)
            if isinstance(blk, dict) and 'n_gen' in blk:
                fn = str(blk.get('file', ''))
                if fn and base.startswith(fn.split('.')[0]):
                    return int(blk['n_gen'])
                if grp == key and not fn:
                    return int(blk['n_gen'])
        return None

    def _load_full(self, path, key):
        opener = gzip.open if path.endswith('.gz') else open
        with opener(path, 'rt') as fh:
            data = np.loadtxt(fh, comments='#', dtype=np.float64, ndmin=2)
        if data.shape[1] < 10:
            raise RuntimeError('%s: expected >= 10 columns, got %d' % (path, data.shape[1]))
        # float32 keeps the memory per process small enough for a many-core parallel scan
        # (positions mm and momenta GeV both keep 7 significant digits; the kernel is fed
        # float64 afterwards and the precision loss is ~1e-7 relative).
        out = {
            'ev': data[:, 0].astype(np.int64),
            'pdg': data[:, 1].astype(np.int32),
            'px': data[:, 2].astype(np.float32), 'py': data[:, 3].astype(np.float32),
            'pz': data[:, 4].astype(np.float32), 'e': data[:, 5].astype(np.float32),
            'x': data[:, 6].astype(np.float32), 'y': data[:, 7].astype(np.float32),
            'z': data[:, 8].astype(np.float32),
        }
        if data.shape[1] >= 11:
            out['status'] = data[:, 10].astype(np.int32)
        return out

    def _load_event_subsample(self, path, key):
        """Stream an oversized parents file and keep ALL rows of a random EVENT subsample.

        The minbias pi/K sample has ~46 rows per event (22.8 M rows, 520 MB gz), so
        materialising it is wasteful.  A uniform random subset of COLLISIONS is statistically
        equivalent for an event-level analysis: the event OR is exact inside the subsample and
        the efficiency is normalised by the number of SUB-sampled collisions, so it stays
        unbiased.  Rows of the chosen collisions are all kept (no eta pre-cut beyond the
        sample's own content; the BC5 N03 lesson: cut on the signal scalar, never the parent).
        """
        n_full = self.inel_n_gen_full or self._norm_n_gen(key, path)
        if not n_full:
            raise RuntimeError('cannot subsample %s: unknown number of generated events' % path)
        n_sub = min(int(self.inel_ev_sub), int(n_full))
        rng = np.random.default_rng(self.inel_seed)
        chosen = np.sort(rng.choice(int(n_full), size=n_sub, replace=False))
        mask = np.zeros(int(n_full) + 2, dtype=bool)
        mask[chosen] = True
        import os.path as _op
        tag = '%s_%d' % (_op.splitext(_op.basename(path))[0].replace('.txt', ''),
                         os.path.getsize(path))
        cache = os.path.join(CACHE_DIR, 'sample_%s_evsub%d_seed%d.npz'
                             % (tag, n_sub, int(self.inel_seed)))
        if os.path.isfile(cache):
            d = np.load(cache, allow_pickle=False)
            out = {k: d[k] for k in d.files}
            out['ev'] = out['ev'].astype(np.int64)
            out['pdg'] = out['pdg'].astype(np.int32)
            if self.verbose:
                print('  loaded inel from subsample cache %s rows=%d n_gen=%d'
                      % (os.path.basename(cache), len(out['ev']), n_sub))
            return out, n_sub
        t0 = time.time()
        opener = gzip.open if path.endswith('.gz') else open
        chunks = {k: [] for k in ('ev', 'pdg', 'px', 'py', 'pz', 'e', 'x', 'y', 'z')}
        block = 400000
        with opener(path, 'rt') as fh:
            while True:
                lines = fh.readlines(block)
                if not lines:
                    break
                buf = ''.join(l for l in lines if l and l[0] != '#')
                if not buf.strip():
                    continue
                d = np.loadtxt(buf.splitlines(), dtype=np.float64, ndmin=2,
                               comments=None)
                if d.size == 0:
                    continue
                ev = d[:, 0].astype(np.int64)
                keep = mask[np.minimum(ev, int(n_full) + 1)]
                if not keep.any():
                    continue
                d = d[keep]
                chunks['ev'].append(d[:, 0].astype(np.int64))
                chunks['pdg'].append(d[:, 1].astype(np.int32))
                for i, k in enumerate(('px', 'py', 'pz', 'e', 'x', 'y', 'z')):
                    chunks[k].append(d[:, 2 + i].astype(np.float32))
        out = {k: np.concatenate(v) if v else np.zeros(0) for k, v in chunks.items()}
        os.makedirs(CACHE_DIR, exist_ok=True)
        np.savez_compressed(cache, **out)
        if self.verbose:
            print('  streamed inel %-20s rows=%-9d n_gen(sub)=%-7d %.1f s -> cache %s'
                  % (os.path.basename(path), len(out['ev']), n_sub, time.time() - t0,
                     os.path.basename(cache)))
        return out, n_sub

    def get(self, key):
        if key in self._cache:
            return self._cache[key]
        path = self.path(key)
        if path is None:
            self._cache[key] = None
            return None
        t0 = time.time()
        n_norm = self._norm_n_gen(key, path)
        if key == 'inel' and self.inel_ev_sub:
            out, n_gen = self._load_event_subsample(path, key)
            self.n_gen[key] = n_gen
            self._cache[key] = out
            if self.verbose:
                print('  inel: event subsample of %d collisions (n_gen used for normalisation)'
                      % n_gen)
            return out
        out = self._load_full(path, key)
        self._cache[key] = out
        n_max = int(out['ev'].max()) + 1
        self.n_gen[key] = n_norm if n_norm is not None else n_max
        if self.verbose:
            print('  loaded %-5s %-22s rows=%-9d n_gen=%-8d (max_ev+1=%-9d) %.1f s'
                  % (key, os.path.basename(path), len(out['ev']), self.n_gen[key], n_max,
                     time.time() - t0))
        return out


# ======================================================================================
# table construction
# ======================================================================================
def _canonical_order(cand):
    """Sort histories by (event, parent, original index) and remap the scalar->history map."""
    hev, hpar = cand['hist_event'], cand['hist_parent']
    order = np.lexsort((np.arange(len(hev)), hpar, hev))
    inv = np.empty_like(order)
    inv[order] = np.arange(len(order))
    out = dict(cand)
    for k in ('hist_event', 'hist_parent', 'hist_group', 'hist_scale', 'hist_lam'):
        out[k] = cand[k][order]
    out['hist_wkey'] = [cand['hist_wkey'][i] for i in order]
    out['scal_hist'] = inv[cand['scal_hist']]
    for g in M.GROUP_NAMES:
        ev = cand['samp_events'].get(g)
        out['samp_events'][g] = np.asarray(ev if ev is not None else [], dtype=np.int64)
    return out


def _hadron_lib_sha_now():
    """sha256 of the hadron-library manifest (or None when no library is loaded).

    The detector response of the hadronisation channels (ss, gg, cc, bb, tautau) depends on
    the Pythia hadron library; a table built WITHOUT the library (analytic jet fallback) must
    never be reused when the library is loaded.  The manifest lists the sha256 of all four
    library files, so hashing it identifies the library state.
    """
    d = getattr(L, 'HADRON_LIB_DIR', None)
    if not d:
        return None
    man = os.path.join(d, 'manifest.json')
    return sha256_file(man) if os.path.isfile(man) else 'no-manifest:' + os.path.basename(d)


def _meta(model, mH, samples, nsamp, ndecay, seed):
    import platform
    dll = os.path.join(RUN, 'work', 'response', 'downstream_kernel.dll')
    return {
        'dll_sha256': sha256_file(dll) if os.path.isfile(dll) else None,
        # single source of truth: the same fingerprint the cache/skip checks compare against
        'code_sha256': _code_sha_now(),
        'config': {k: getattr(L, k) for k in
                   ('Z_DEC_MIN', 'Z_DEC_MAX', 'ETA_S_MIN', 'ETA_S_MAX', 'RECO_EFF', 'NODE_MTOL')},
        'detector_config': dict(L.dk.DEFAULT_CONFIG),
        'birth_cuts': dict(L.dk.BIRTH_CUTS),
        'zgrid': [float(z) for z in L.ZGRID],
        'hadron_lib_dir': getattr(L, 'HADRON_LIB_DIR', None),
        'ndecay': ndecay,
        'nsamp': nsamp,
        'seed': seed,
        'hadron_lib': getattr(L, 'HADRON_LIB_DIR', None),
        'hadron_lib_sha': _hadron_lib_sha_now(),
        'python': platform.python_version(),
        'samples': {k: (sha256_file(p) if p else None) for k, p in samples.available().items()},
        'samples_stat': {k: [os.path.getsize(p), int(os.stat(p).st_mtime_ns)]
                         for k, p in samples.available().items() if p},
        'samples_n_gen': dict(samples.n_gen),
    }


class _NoSamples:
    """Empty sample store, used by the unit tests (no parent samples needed)."""
    n_gen = {}

    def get(self, key):
        return None

    def available(self):
        return {}


def build_table(model, mH, samples, seed=DEFAULT_SEED, nsamp=None, ndecay=None, verbose=True):
    """Build the response table for one H mass (no caching; use get_table for the cache)."""
    t0 = time.time()
    if samples is None:
        samples = _NoSamples()
    mH = float(mH)
    nsamp = dict(DEFAULT_NSAMP if nsamp is None else nsamp)
    if ndecay is None:
        ndecay = dict(L.DEFAULT_NDECAY)
    chans = list(model.channels(mH))
    import random as _random
    s = (seed + int(round(mH * 1e6))) % 1000003
    rng = _random.Random(s)
    np_rng = np.random.default_rng(s)
    cand = model.build_candidates(samples, mH, rng, np_rng, nsamp)
    # ---- review BUG B/C: make collision ids and parent indices globally unique FIRST,
    # then sort histories by (event, parent, original index). ---------------------------
    hgroup = np.asarray(cand['hist_group'], dtype=np.int8)
    ev_off = np.zeros(len(hgroup), dtype=np.int64)
    par_off = np.zeros(len(hgroup), dtype=np.int64)
    for gi, o in _GROUP_OFFSET.items():
        ev_off[hgroup == gi] = o
    for gi, o in _PARENT_OFFSET.items():
        par_off[hgroup == gi] = o
    cand['hist_event'] = np.asarray(cand['hist_event'], dtype=np.int64) + ev_off
    cand['hist_parent'] = np.asarray(cand['hist_parent'], dtype=np.int64) + par_off
    cand = _canonical_order(cand)
    # ---- review BUG A: the scalars must be ordered by their (remapped) history index,
    # otherwise the run-length grouping in event_engine is invalid. -------------------
    perm = np.argsort(cand['scal_hist'], kind='stable')
    cand['scal_hist'] = cand['scal_hist'][perm]
    cand['scal_p4'] = cand['scal_p4'][perm]
    cand['scal_v0'] = cand['scal_v0'][perm]
    if len(cand['scal_hist']):
        assert np.all(np.diff(cand['scal_hist']) >= 0), 'scalars not sorted by history'
    # ---- review BUG D: sampled collisions sorted and unique (a collision appears once) --
    samp_events = {}
    for g, name in enumerate(M.GROUP_NAMES):
        ev = cand['samp_events'].get(name)
        if ev is None or len(ev) == 0:
            samp_events[name] = np.zeros(0, dtype=np.int64)
        else:
            samp_events[name] = np.unique(np.asarray(ev, dtype=np.int64) + _GROUP_OFFSET[g])
    n_s = len(cand['scal_hist'])
    n_z = len(L.ZGRID)
    A = np.zeros((n_s, len(chans), n_z), dtype=np.float32)
    p4 = cand['scal_p4']
    v0 = cand['scal_v0']
    tA = time.time()
    n_used = 0
    for i in range(n_s):
        Ai = L.acceptance_profiles(mH, (p4[i, 0], p4[i, 1], p4[i, 2], p4[i, 3]),
                                   (v0[i, 0], v0[i, 1], v0[i, 2]), rng,
                                   channels=chans, ndecay=ndecay)
        if Ai.any():
            n_used += 1
        A[i] = Ai
    ycorr = model.yield_corr(samples)
    m = _meta(model, mH, samples, nsamp, ndecay, seed)
    m.update({'t_build_s': time.time() - t0, 't_response_s': time.time() - tA,
              'n_scalar': int(n_s), 'n_scalar_nonzero_response': int(n_used),
              'n_hist': int(len(cand['hist_event'])), 'channels': chans,
              'stratification': getattr(model, 'stratification', None),
              'b_yield_check': getattr(model, 'b_yield_check', None)})
    table = {
        'version': TABLE_VERSION,
        'model': model.name,
        'mH': mH,
        'channels': chans,
        'pz': p4[:, 3].astype(np.float64),
        'z0': (v0[:, 2] * 1e-3).astype(np.float64),
        'A': A,
        'scal_hist': cand['scal_hist'].astype(np.int32),
        'hist_event': cand['hist_event'].astype(np.int64),
        'hist_parent': cand['hist_parent'].astype(np.int64),
        'hist_group': cand['hist_group'].astype(np.int8),
        'hist_wkey': list(cand['hist_wkey']),
        'hist_scale': np.asarray(cand['hist_scale'], dtype=np.float64),
        'hist_lam': np.asarray(cand['hist_lam'], dtype=np.float64),
        'samp_events': samp_events,
        'n_gen': {g: int(samples.n_gen[g]) for g in M.GROUP_NAMES if g in samples.n_gen},
        'yield_corr': {g: float(ycorr[g]) for g in M.GROUP_NAMES},
        'h_frac': dict(cand.get('h_frac', {})),
        'counts': {g: [int(v[0]), int(v[1])] for g, v in cand.get('counts', {}).items()},
        'nsamp': nsamp,
        'meta': m,
    }
    if verbose:
        print('  [table %s mH=%.6g] scalars=%d (nonzero response %d) histories=%d '
              'response %.1f s, total %.1f s'
              % (model.name, mH, n_s, n_used, len(table['hist_event']),
                 m['t_response_s'], m['t_build_s']))
    return table


# ======================================================================================
# cache
# ======================================================================================
def cache_path(model, mH):
    return os.path.join(CACHE_DIR, 'response_%s_m%.6f_v%d.npz' % (model.name, float(mH), TABLE_VERSION))


def save_table(table, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    meta = {k: table[k] for k in ('version', 'model', 'mH', 'channels', 'hist_wkey', 'n_gen',
                                  'yield_corr', 'h_frac', 'counts', 'nsamp', 'meta')}
    arrays = {k: table[k] for k in ('pz', 'z0', 'A', 'scal_hist', 'hist_event',
                                    'hist_parent', 'hist_group', 'hist_scale', 'hist_lam')}
    for g in M.GROUP_NAMES:
        arrays['samp_events_' + g] = table['samp_events'][g]
    np.savez_compressed(path, meta_json=np.array(json.dumps(meta)), **arrays)
    return path


def load_table(path):
    d = np.load(path, allow_pickle=False)
    table = json.loads(str(d['meta_json']))
    for k in ('pz', 'z0', 'A', 'scal_hist', 'hist_event', 'hist_parent', 'hist_group'):
        table[k] = d[k]
    table['hist_scale'] = (d['hist_scale'] if 'hist_scale' in d.files
                           else np.ones(len(table['hist_event']), dtype=np.float64))
    table['hist_lam'] = (d['hist_lam'] if 'hist_lam' in d.files
                         else np.zeros(len(table['hist_event']), dtype=np.float64))
    table['samp_events'] = {g: d['samp_events_' + g] for g in M.GROUP_NAMES}
    return table


def _code_sha_now():
    # b_channels.py / twohdm_production.py feed channel_masses() and the kinematic mixture,
    # i.e. the table content itself, so they are part of the response-code fingerprint
    return sha256_code([os.path.join(RUN, 'work', 'response', 'twohdm_lib.py'),
                        os.path.join(HERE, 'models.py'),
                        os.path.join(HERE, 'response_table.py'),
                        os.path.join(RUN, 'work', 'production', 'b_channels.py'),
                        os.path.join(RUN, 'work', 'production', 'twohdm_production.py')])


def _SAMPLES_AVAILABLE():
    """[(group, path)] for the default sample directory (cheap; no file hashing)."""
    out = []
    for k, fn in SAMPLE_FILES.items():
        base = os.path.join(SAMPLES_DEFAULT, fn)
        p = None
        for cand in (base, base + '.gz'):
            if os.path.isfile(cand):
                p = cand
                break
        out.append((k, p))
    return out


def _cache_matches(table, nsamp, ndecay, seed):
    m = table.get('meta', {})
    if table.get('version') != TABLE_VERSION:
        return 'table version'
    if m.get('code_sha256') != _code_sha_now():
        return 'response code changed'
    if int(m.get('seed', -1)) != int(seed):
        return 'seed'
    # the meta records the detector configuration, the DLL and the parent-sample fingerprint;
    # compare them, otherwise a changed detector/sample would silently reuse the old response
    dll = os.path.join(RUN, 'work', 'response', 'downstream_kernel.dll')
    if m.get('dll_sha256') != (sha256_file(dll) if os.path.isfile(dll) else None):
        return 'detector kernel dll'
    for key, now in (('detector_config', dict(L.dk.DEFAULT_CONFIG)),
                     ('birth_cuts', dict(L.dk.BIRTH_CUTS)),
                     ('config', {k: getattr(L, k) for k in
                                 ('Z_DEC_MIN', 'Z_DEC_MAX', 'ETA_S_MIN', 'ETA_S_MAX',
                                  'RECO_EFF', 'NODE_MTOL')}),
                     ('zgrid', [float(z) for z in L.ZGRID])):
        if m.get(key) != now:
            return 'detector/birth configuration (%s)' % key
    if 'samples_stat' in m:
        for k, p in _SAMPLES_AVAILABLE():
            st = m['samples_stat'].get(k)
            cur = [os.path.getsize(p), int(os.stat(p).st_mtime_ns)] if p else None
            if st != cur:
                return 'parent sample changed (%s)' % k
    # v2 caches built before this field existed are trusted; newly built ones must match the
    # currently loaded hadron library exactly (review finding: cache could otherwise poison
    # the hadronisation response with the analytic fallback).
    if 'hadron_lib' in m:
        if m.get('hadron_lib') != getattr(L, 'HADRON_LIB_DIR', None):
            return 'hadron library'
        if m.get('hadron_lib_sha') != _hadron_lib_sha_now():
            return 'hadron library content'
    if dict(m.get('nsamp', {})) != dict(nsamp):
        return 'nsamp'
    if dict(m.get('ndecay', {})) != dict(ndecay):
        return 'ndecay'
    return None


def get_table(model, mH, samples=None, seed=DEFAULT_SEED, nsamp=None, ndecay=None,
              cache=True, force=False, verbose=True):
    """Build (or load from cache) the response table for one mass."""
    nsamp = dict(DEFAULT_NSAMP if nsamp is None else nsamp)
    if ndecay is None:
        ndecay = dict(L.DEFAULT_NDECAY)
    path = cache_path(model, mH)
    if cache and os.path.isfile(path) and not force:
        try:
            tab = load_table(path)
        except Exception as exc:                      # corrupt cache: rebuild
            if verbose:
                print('  cache %s unreadable (%s): rebuilding' % (os.path.basename(path), exc))
            tab = None
        if tab is not None:
            why = _cache_matches(tab, nsamp, ndecay, seed)
            if why is None:
                if verbose:
                    print('  [cache hit %s mH=%.6g] scalars=%d'
                          % (model.name, float(mH), len(tab['scal_hist'])))
                return tab
            if verbose:
                print('  cache %s stale (%s): rebuilding' % (os.path.basename(path), why))
    if samples is None:
        samples = Samples(verbose=verbose)
    tab = build_table(model, mH, samples, seed=seed, nsamp=nsamp, ndecay=ndecay, verbose=verbose)
    if cache:
        save_table(tab, path)
        if verbose:
            print('  cached -> %s (%.1f MB)' % (os.path.basename(path),
                                                os.path.getsize(path) / 1e6))
    return tab
