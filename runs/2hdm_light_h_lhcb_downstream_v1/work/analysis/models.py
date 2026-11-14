# -*- coding: utf-8 -*-
"""Model adapters for the response/event engine.

Two models share ONE engine (work/analysis/event_engine.py):

  TwoHDMModel  -- the physics of THIS run: production rates from
                  work/production/twohdm_production.py, decay widths/BRs from
                  work/decay/decay_widths.py (author Light_scalar_decay).  Scan
                  parameter is tanb.

  BC5Model     -- the BC5 v3.6 reference model (scalar_widths + BC5 production, mixing
                  angle theta^2), used ONLY for the end-to-end regression against the
                  published BC5 25 fb^-1 long-lived boundary.  It lives on the frozen
                  copies in work/analysis/bc5ref/ (sha256 recorded in the validation log).

A model must provide:
    channels(mH)                 -> list of decay-channel names stored in the response table
    open_channels(mH, param)     -> list of channels with non-zero width
    decay_br(mH, param)          -> {channel: BR}   (sums to 1; keys in table['channels'])
    ctau_m(mH, param)            -> c * tau_H [m]
    build_candidates(samples, mH, rng, np_rng, nsamp, mass_grid)
                                 -> candidate/history records (production kinematics)
    history_weights(mH, param, table) -> np.ndarray over table['hist_*']
    sigma_pb(mH, param)          -> {group: pb}
    yield_corr                   -> {group: float}   (1.0 except the BC5 B correction)

Parent-sample normalisation (DESIGN.md section 4):
    s = 1e3 * L * sum_c sigma_c * eps_c ,  eps_c = yield_corr_c * 0.4/N_gen_c * sum_events o
"""
from __future__ import annotations
import math, os, sys

import numpy as np

RUN = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
RESP_DIR = os.path.join(RUN, 'work', 'response')
PROD_DIR = os.path.join(RUN, 'work', 'production')
DECAY_DIR = os.path.join(RUN, 'work', 'decay')
BC5REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bc5ref')
for _p in (RESP_DIR, PROD_DIR, DECAY_DIR,
           os.path.join(BC5REF, 'analysis'), os.path.join(BC5REF, 'widths'),
           os.path.join(BC5REF, 'response')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import twohdm_lib as L

HBAR_C = 1.973269804593025e-16          # GeV m
HBAR_GEV_S = 6.582119569e-25           # GeV s
# Electron mass used in the three-body semileptonic phase space X -> H e nu.  The paper's
# Br(X->H e nu) formula uses m_mu only through the NORMALISATION Br(X->mu nu); the final
# state is e + nu with m_e -> 0.  Using m_mu in the phase space would distort the H energy
# spectrum badly at low mH, so m_e is used here (m_e = 0.511 MeV, numerically negligible).
M_E_SEMILEPTONIC = L.M_E

# Parent-meson properties (identical to the BC5 v3.6 values, inputs/bc5_reuse/bc5_production.py)
M_BP, M_B0, M_BS = 5.27934, 5.27965, 5.36688
M_KP, M_K0L = 0.493677, 0.497611
M_PI_CHARGED = L.M_PI
GAM_BP = HBAR_GEV_S / 1.638e-12
GAM_B0 = HBAR_GEV_S / 1.519e-12
GAM_BS = HBAR_GEV_S / 1.520e-12
GAM_KP = HBAR_GEV_S / 1.2380e-8
GAM_KL = HBAR_GEV_S / 5.116e-8
GAM_PI = HBAR_GEV_S / 2.6033e-8
CTAU_BP = HBAR_C / GAM_BP
CTAU_B0 = HBAR_C / GAM_B0
CTAU_BS = HBAR_C / GAM_BS
CTAU_KP = HBAR_C / GAM_KP
CTAU_KL = HBAR_C / GAM_KL
CTAU_PI = HBAR_C / GAM_PI

#: mass of the on-shell SM Higgs in the higgs.txt sample (Pythia 8.316, HiggsSM:all, mh=125.25)
M_HIGGS_SAMPLE = 125.25

# Pythia 8 Higgs-process codes of the higgs.txt sample -> target process name
# (identical to BC5 work/analysis/bc5_production.py HIGGS_PROC_GROUP)
HIGGS_PROC_GROUP = {902: 'ggF', 903: 'ggF', 906: 'VBF', 907: 'VBF', 905: 'WH',
                    904: 'ZH', 901: 'bbH', 908: 'ttH', 909: 'ttH'}

# sigma(pp -> b bbar) [pb]: 560 ub, BC5 v3.6 value (samples are unbiased in eta)
SIGMA_BB_PB = 5.6e8
# inelastic pp cross section [pb]: 79.5 mb (measured normalisation of the minbias sample)
SIGMA_INEL_PB = 7.95e10


# ======================================================================================
# helpers shared by both models
# ======================================================================================
def _p4_and_mass(e, px, py, pz):
    m = math.sqrt(max(e * e - px * px - py * py - pz * pz, 1e-12))
    return m


def _parent_decay_vertex(mes, ctau_m, rng, mass=None):
    """Displace the parent decay point by beta*gamma*c*tau*Exp(1) [mm sample units].

    mass: explicit parent mass.  For the 22.8 M-row inelastic sample the %8g text quantises
    E strongly above ~100 GeV, so sqrt(E^2-p^2) is noisy for a 0.14-0.50 GeV parent; the
    PDG mass is used there instead (samples-owner note).
    """
    e, px, py, pz = mes['e'], mes['px'], mes['py'], mes['pz']
    m = float(mass) if mass is not None else _p4_and_mass(e, px, py, pz)
    p3 = math.sqrt(px * px + py * py + pz * pz)
    if p3 <= 0.0 or ctau_m <= 0.0:
        return mes['x'], mes['y'], mes['z']
    # BC5 convention: c*tau is in METRES, the sample vertex in MILLIMETRES -> x1e3
    disp = (p3 / m) * ctau_m * (-math.log(max(rng.random(), 1e-12))) * 1e3
    return (mes['x'] + disp * px / p3, mes['y'] + disp * py / p3, mes['z'] + disp * pz / p3)


def _boost_to_lab(p4_rest, mes):
    lab = L.boost([p4_rest], mes['e'], (mes['px'], mes['py'], mes['pz']))[0]
    return (lab[1], lab[2], lab[3], lab[4])


def _boost_on_shell(p4_rest, mes, mass):
    """Boost a rest-frame particle into the lab with a mass-consistent parent momentum.

    The parent samples are stored as float32 (memory), so the stored E carries a ~1e-7
    relative error.  In the lab boost beta = p/E and gamma = 1/sqrt(1-beta^2) amplify that
    error by gamma^2: for the TeV-scale parents of the high-boost tail this makes the
    boosted four-momenta (and even the daughter mass shell) wrong by O(1).  Rebuilding
    E = sqrt(p^2 + m^2) from the 3-momentum and the KNOWN parent mass keeps the boost
    mass-shell consistent, and the boosted particle is re-set on shell as well.
    """
    px, py, pz = float(mes['px']), float(mes['py']), float(mes['pz'])
    e = math.sqrt(px * px + py * py + pz * pz + float(mass) * float(mass))
    lab = L.boost([p4_rest], e, (px, py, pz))[0]
    m = float(p4_rest[0])
    p2 = lab[2] * lab[2] + lab[3] * lab[3] + lab[4] * lab[4]
    return (math.sqrt(p2 + m * m), lab[2], lab[3], lab[4])


def _sample_categorical(weights, rng):
    """Index sampled proportionally to a dict {key: weight} (all weights > 0)."""
    items = [(k, float(w)) for k, w in weights.items() if w > 0.0]
    if not items:
        return None
    tot = sum(w for _k, w in items)
    u = rng.random() * tot
    c = 0.0
    for k, w in items:
        c += w
        if u <= c:
            return k
    return items[-1][0]


class HistBuilder:
    """Accumulates histories and their scalars in a compact, sortable form."""

    def __init__(self):
        self.h_event = []
        self.h_parent = []
        self.h_group = []
        self.h_wkey = []
        self.h_scale = []
        self.h_lam = []
        self.s_hist = []
        self.s_p4 = []
        self.s_v0 = []
        self.samp_events = {'B': [], 'inel': [], 'h': []}
        self.h_frac = {}
        # (n_tot, n_smp) of the parent particles available / evaluated, per group.  The
        # signal is normalised by the number of parents SAMPLED and then scaled to the full
        # sample (BC5 channel_efficiency); without this factor a subsample would silently
        # rescale the whole signal.
        self.counts = {'B': [0, 0], 'inel': [0, 0], 'h': [0, 0]}

    def add(self, event, parent, group, wkey, scalars, scale=1.0, lam=0.0):
        hidx = len(self.h_event)
        self.h_event.append(int(event))
        self.h_parent.append(int(parent))
        self.h_group.append(group)
        self.h_wkey.append(wkey)
        # stratified/importance weight of this parent (Horvitz-Thompson); 1.0 = plain uniform
        self.h_scale.append(float(scale))
        # (pz_parent/m_parent) * c*tau_parent [m]: the parent decay position is integrated
        # out analytically in the engine (0 -> prompt parent)
        self.h_lam.append(float(lam))
        for (p4, v0) in scalars:
            self.s_hist.append(hidx)
            self.s_p4.append(p4)
            self.s_v0.append(v0)
        return hidx

    def finalize(self):
        return {
            'hist_event': np.asarray(self.h_event, dtype=np.int64),
            'hist_parent': np.asarray(self.h_parent, dtype=np.int64),
            'hist_group': np.asarray(self.h_group, dtype=np.int8),
            'hist_wkey': list(self.h_wkey),
            'hist_scale': np.asarray(self.h_scale, dtype=np.float64),
            'hist_lam': np.asarray(self.h_lam, dtype=np.float64),
            'scal_hist': np.asarray(self.s_hist, dtype=np.int32),
            'scal_p4': np.asarray(self.s_p4, dtype=np.float64).reshape(-1, 4),
            'scal_v0': np.asarray(self.s_v0, dtype=np.float64).reshape(-1, 3),
            'samp_events': {k: np.asarray(v, dtype=np.int64) for k, v in self.samp_events.items()},
            'h_frac': dict(self.h_frac),
            'counts': {k: [int(v[0]), int(v[1])] for k, v in self.counts.items()},
        }


GROUP_NAMES = ('B', 'inel', 'h')
_GROUP_ID = {'B': 0, 'inel': 1, 'h': 2}

#: Parent-momentum strata used by the stratified (Horvitz-Thompson) subsampling.
#: The short-lived crossings are driven by the extreme high-pz tail of the parent
#: sample: c*tau*beta*gamma ~ 1 m needs pz of a few TeV, so a uniform subsample
#: represents that tail with only a handful of events (65-99 % MC error at the short
#: crossings).  The strata below census the tail (up to the sample population) and keep
#: a uniform bulk sample for the long-lived edges; every sampled parent carries its
#: exact inverse-inclusion-probability weight, so the estimator remains unbiased.
B_STRAT_EDGES = (0.0, 200.0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 1e30)
B_STRAT_ALLOC = (5000, 3500, 3000, 2500, 1800, 1014, 337, 142)
H_STRAT_EDGES = (0.0, 200.0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 1e30)
H_STRAT_ALLOC = (800, 500, 500, 600, 700, 1061, 346, 113)


def _scaled_alloc(alloc, n_target):
    """Rescale a stratification allocation so that it sums to n_target (shape preserved).

    The module-level B_STRAT_ALLOC / H_STRAT_ALLOC define the *relative* allocation (tail
    census); the --nsamp value is the requested total number of parents.  This makes the
    recorded nsamp actually control the sample, instead of being a label only.
    """
    a = np.asarray(alloc, dtype=float)
    try:
        n = int(n_target)
    except (TypeError, ValueError):
        return tuple(int(x) for x in a)
    if n <= 0 or a.sum() <= 0:
        return tuple(int(x) for x in a)
    f = float(n) / float(a.sum())
    return tuple(max(0, int(round(x * f))) for x in a)


def _stratified_rows(pz, alloc, edges, np_rng):
    """Stratified subsample of row positions by |pz| with Horvitz-Thompson weights.

    Uniform sampling inside each pz bin, WITHOUT replacement; a bin whose target exceeds
    its population is censused.  Returns (positions, weights, n_sampled), where the weight
    already includes the n_smp/N_tot factor that pairs with the engine's n_tot/n_smp
    normalisation, so that sum_sampled weight*o == sum_all o in expectation.
    """
    pz = np.abs(np.asarray(pz, dtype=float))
    edges = np.asarray(edges, dtype=float)
    bin_of = np.digitize(pz, edges[1:-1], right=False)
    n_bin = len(alloc)
    N_b = np.bincount(bin_of, minlength=n_bin).astype(np.int64)
    n_b = np.minimum(np.asarray(alloc, dtype=np.int64), N_b)
    pos = []
    for b in range(n_bin):
        idx = np.flatnonzero(bin_of == b)
        if n_b[b] <= 0 or len(idx) == 0:
            continue
        if n_b[b] >= len(idx):
            pos.append(idx)
        else:
            pos.append(np_rng.choice(idx, size=int(n_b[b]), replace=False))
    pos = np.sort(np.concatenate(pos)) if pos else np.zeros(0, dtype=np.int64)
    n_smp = int(len(pos))
    N_tot = int(len(pz))
    b_of = bin_of[pos]
    w = (N_b[b_of] / np.maximum(n_b[b_of], 1)).astype(float) * (float(n_smp) / float(N_tot))
    return pos, w, n_smp


# ======================================================================================
# 2HDM light-H model
# ======================================================================================
class TwoHDMModel:
    """Type-I 2HDM light CP-even H (Li et al. 2212.06186), LHCb-downstream production."""

    name = '2hdm'
    param_name = 'tanb'
    # larger tanb -> xi_H^f ~ tanb^-3 -> Gamma_H ~ tanb^-6 -> c*tau_H ~ tanb^6: LONGER lifetime
    param_lifetime_sign = +1

    def __init__(self, verbose=True, decay_dir=None, b_mode=None):
        import importlib
        import twohdm_production as P
        if decay_dir is not None and decay_dir not in sys.path:
            sys.path.insert(0, decay_dir)
        D = importlib.import_module('decay_widths')
        self.P = P
        self.D = D
        self.verbose = verbose
        self._br_cache = {}
        self._mix_cache = {}
        self._open_check = {}
        # B rate prescription: production default ('foresee' nominal) unless overridden.
        self.b_mode = b_mode if b_mode is not None else getattr(P, 'DEFAULT_B_MODE', None)
        self._b_mode_how = 'version'
        # The B kinematic mixture is tanb-independent (production contract); verify once.
        self.b_mix_ref_tanb = 100.0

    # ---- channels and decay ------------------------------------------------------
    def channels(self, mH):
        """Channel list stored in the response table = author-program mass windows."""
        return L.open_channels_2hdm(mH)

    def open_channels(self, mH, tanb):
        return self.channels(mH)

    def decay_br(self, mH, tanb):
        key = (round(float(mH), 9), round(float(tanb), 9))
        got = self._br_cache.get(key)
        if got is not None:
            return got
        widths, br = self.D.widths_and_br(mH, tanb)
        tot = float(sum(widths.values()))
        br = {k: float(v) for k, v in br.items()}
        # Consistency check: the decay module must agree with the author program's channel
        # windows (twohdm_lib.OPEN_2HDM) that define the response-table channel set.  The
        # program uses 1e-40 as its "closed" placeholder, so a channel whose WIDTH exceeds
        # 1e-35 GeV is genuinely open and must be covered by the table.
        chans = set(self.channels(mH))
        for ch, w in widths.items():
            if float(w) > 1e-35 and ch not in chans and ch != 'gaga':
                raise RuntimeError(
                    'decay module gives Gamma(%s)=%.3e GeV at mH=%.6g, tanb=%.6g but the '
                    'author-program window says the channel is closed; response table wrong'
                    % (ch, float(w), mH, tanb))
        got = (br, tot)
        self._br_cache[key] = got
        return got

    def ctau_m(self, mH, tanb):
        _br, tot = self.decay_br(mH, tanb)
        if tot <= 0:
            return float('inf')
        return HBAR_C / tot

    # ---- B production kinematics ------------------------------------------------
    def b_kin_mixture(self, mH, meson='B+'):
        """Normalised kinematic mixture over the b->s exclusive channels (tanb-independent).

        meson matters near the endpoint: the lightest recoil is K+- for B+- but K0 for B0,
        so the K channel closes 3.6 MeV earlier for the neutral B.
        """
        key = (round(float(mH), 9), meson)
        got = self._mix_cache.get(key)
        if got is not None:
            return got
        if hasattr(self.P, 'b_channel_weights'):
            w = self.P.b_channel_weights(mH, meson)
            got = {k: float(v) for k, v in w.items() if v > 0.0}
        else:
            res = self.P.br_B_to_XsH(mH, self.b_mix_ref_tanb, meson)
            ch = res['channels']
            tot = float(sum(ch.values()))
            got = {k: float(v) / tot for k, v in ch.items() if v > 0.0} if tot > 0 else {}
        self._mix_cache[key] = got
        return got

    def _b_total(self, mH, tanb, meson):
        """Absolute Br(B -> X_s H) for the active B rate prescription (b_mode)."""
        P = self.P
        fn = P.br_B_to_XsH
        if self.b_mode is not None:
            import inspect
            params = inspect.signature(fn).parameters
            for kw in ('mode', 'b_mode'):
                if kw in params:
                    self._b_mode_how = 'kwarg:' + kw
                    return float(fn(mH, tanb, meson, **{kw: self.b_mode})['total'])
            for attr in ('B_MODE', 'b_mode', 'B_MODE_DEFAULT'):
                if hasattr(P, attr):
                    old = getattr(P, attr)
                    setattr(P, attr, self.b_mode)
                    try:
                        self._b_mode_how = 'module-attr:' + attr
                        return float(fn(mH, tanb, meson)['total'])
                    finally:
                        setattr(P, attr, old)
        self._b_mode_how = 'default'
        return float(fn(mH, tanb, meson)['total'])

    def b_parent_types(self):
        return ('B+', 'B0')

    def open_modes(self, mH, ref_tanb=100.0):
        """Which production modes can contribute at all at this mass (mass thresholds only).

        Every production branching ratio in this model is proportional to a xi_H^{ij}^2 and
        vanishes identically outside the channel's mass threshold, so evaluating it at one
        reference tanb is enough to decide whether a history can EVER carry weight.  Modes
        that are closed are not sampled at all: this saves most of the CPU above 0.5 GeV
        (the pi/K channel) and above 4.786 GeV (the B channel).
        """
        key = round(float(mH), 9)
        got = getattr(self, '_open_cache', {}).get(key)
        if got is not None:
            return got
        P = self.P
        got = {
            'B:mix': self._b_total(mH, ref_tanb, 'B+') > 0.0,
            'K+:piH': P.br_K_to_piH(mH, ref_tanb, 'K+') > 0.0,
            'K+:enu': P.br_K_to_Henu(mH, ref_tanb, 'K+') > 0.0,
            'K0L:piH': P.br_K_to_piH(mH, ref_tanb, 'K0L') > 0.0,
            'K0L:enu': P.br_K_to_Henu(mH, ref_tanb, 'K0L') > 0.0,
            'pi+:enu': P.br_pi_to_Henu(mH, ref_tanb) > 0.0,
        }
        if not hasattr(self, '_open_cache'):
            self._open_cache = {}
        self._open_cache[key] = got
        return got

    def build_candidates(self, samples, mH, rng, np_rng, nsamp, mass_grid=None):
        hb = HistBuilder()
        # per-meson kinematic mixture and recoil masses (K+- for B+-, K0 for B0)
        mix_by = {bm: self.b_kin_mixture(mH, bm) for bm in ('B+', 'B0')}
        masses_by = {bm: self.P.channel_masses(bm) for bm in ('B+', 'B0')}
        mix = mix_by['B+'] or mix_by['B0']
        modes = self.open_modes(mH)
        self._modes_used = modes

        # ---------------- B mesons (bbar_heavy sample) --------------------------
        smp = samples.get('B')
        if smp is not None and mix and modes['B:mix']:
            sel = np.isin(smp['pdg'], [521, -521, 511, -511])
            idxs_all = np.where(sel)[0]
            n_tot = len(idxs_all)
            b_alloc = _scaled_alloc(B_STRAT_ALLOC, nsamp.get('B'))
            pick, w_strat, n_smp = _stratified_rows(
                smp['pz'][idxs_all], b_alloc, B_STRAT_EDGES, np_rng)
            idxs = idxs_all[pick]
            hb.counts['B'] = [n_tot, n_smp]
            hb.samp_events['B'] = smp['ev'][idxs].copy()
            self.stratification = {'B': {'edges': list(B_STRAT_EDGES), 'alloc': list(b_alloc),
                                         'n_tot': n_tot, 'n_smp': n_smp}}
            for ip, i in enumerate(idxs):
                pdg = int(smp['pdg'][i])
                bm = 'B+' if abs(pdg) == 521 else 'B0'
                mes = {'e': smp['e'][i], 'px': smp['px'][i], 'py': smp['py'][i], 'pz': smp['pz'][i],
                       'x': smp['x'][i], 'y': smp['y'][i], 'z': smp['z'][i]}
                got = self._b_mixing_scalars(mH, mes, bm, mix_by[bm], masses_by[bm], rng)
                if got is not None:
                    p4, v0, lam = got
                    hb.add(smp['ev'][i], ip, _GROUP_ID['B'], 'B:mix:' + bm, [(p4, v0)],
                           scale=float(w_strat[ip]), lam=lam)
        elif smp is not None:
            hb.samp_events['B'] = np.zeros(0, dtype=np.int64)

        # ---------------- pi / K mesons (minbias sample) ------------------------
        smp = samples.get('inel')
        if smp is not None:
            base = 0
            for pdg, which in ((321, 'K+'), (-321, 'K+'), (311, 'K0L'), (-311, 'K0L'),
                               (211, 'pi+'), (-211, 'pi+')):
                if which == 'K+' and not (modes['K+:piH'] or modes['K+:enu']):
                    continue
                if which == 'K0L' and not (modes['K0L:piH'] or modes['K0L:enu']):
                    continue
                if which == 'pi+' and not modes['pi+:enu']:
                    continue
                idxs = np.where(smp['pdg'] == pdg)[0]
                n_tot = len(idxs)
                n_smp = min(int(nsamp['inel']), n_tot)
                if n_smp > 0:
                    idxs = np_rng.choice(idxs, size=n_smp, replace=False)
                hb.counts['inel'][0] += n_tot
                hb.counts['inel'][1] += n_smp
                hb.samp_events['inel'].extend(smp['ev'][idxs].tolist())
                for ip, i in enumerate(idxs):
                    mes = {'e': smp['e'][i], 'px': smp['px'][i], 'py': smp['py'][i], 'pz': smp['pz'][i],
                           'x': smp['x'][i], 'y': smp['y'][i], 'z': smp['z'][i]}
                    for wkey, scal, lam in self._inel_histories(mH, mes, which, rng):
                        hb.add(smp['ev'][i], base + ip, _GROUP_ID['inel'], wkey, [scal],
                               lam=lam)
                base += n_smp
            hb.samp_events['inel'] = np.asarray(hb.samp_events['inel'], dtype=np.int64)

        # ---------------- Higgs (higgs sample, h -> HH) -------------------------
        smp = samples.get('h')
        if smp is not None:
            idxs_all, proc_all = self._higgs_predecay_rows(smp)
            n_tot = len(idxs_all)
            # The process mixture used for the reweighting must be the FULL-sample one
            # (the pz-stratified subsample is not process-blind: ggF and VBF have very
            # different pz spectra).  The stratification weights then restore the unbiased
            # sum over events.
            groups_all = np.array([HIGGS_PROC_GROUP.get(int(c), 'other') for c in proc_all])
            h_alloc = _scaled_alloc(H_STRAT_ALLOC, nsamp.get('h'))
            pick, w_strat, n_smp = _stratified_rows(
                smp['pz'][idxs_all], h_alloc, H_STRAT_EDGES, np_rng)
            idxs = idxs_all[pick]
            groups = groups_all[pick]
            hb.counts['h'] = [n_tot, n_smp]
            hb.samp_events['h'] = smp['ev'][idxs].copy()
            self.stratification = dict(getattr(self, 'stratification', {}),
                                       h={'edges': list(H_STRAT_EDGES), 'alloc': list(h_alloc),
                                          'n_tot': n_tot, 'n_smp': n_smp})
            for ip, i in enumerate(idxs):
                scal = self._higgs_scalars(mH, smp, i, rng)
                if scal is not None:
                    grp = groups[ip]
                    hb.add(smp['ev'][i], ip, _GROUP_ID['h'], 'h:' + str(grp), scal,
                           scale=float(w_strat[ip]), lam=0.0)
            hb.h_frac = {g: float(np.mean(groups_all == g)) for g in sorted(set(groups_all.tolist()))}
        else:
            hb.h_frac = {}
        return hb.finalize()

    # ---- candidate builders ------------------------------------------------------
    def _b_mixing_scalars(self, mH, mes, bm, mix, masses, rng):
        label = _sample_categorical(mix, rng)
        if label is None:
            return None
        mXp = masses[label]
        # PDG B mass: the parent rest-frame energy scale must be exact (the %.8g sample text
        # makes sqrt(E^2-p^2) noisy in the high-E tail; the boost still uses the measured
        # E,px,py,pz, whose direction/velocity is accurate).
        m = M_BP if bm == 'B+' else M_B0
        d = L.two_body(m, mXp, mH, rng)
        if d is None:
            return None
        p4 = _boost_on_shell(d[1], mes, m)
        ctau = CTAU_BP if bm == 'B+' else CTAU_B0
        # The B decay position is integrated out analytically in the engine: store the
        # production vertex and the z-flight scale lam_z = (pz_B/m_B) c*tau_B [m].
        pz = float(mes['pz'])
        lam = (pz / m) * ctau          # signed: >0 forward flight, <0 backward flight
        v0 = (mes['x'], mes['y'], mes['z'])
        return (p4, v0, lam)

    def _inel_histories(self, mH, mes, which, rng):
        """One or two mutually exclusive histories of one K / pi parent.

        The parent mass is the PDG mass (not sqrt(E^2-p^2)): the minbias file stores E with
        8 significant digits, which makes E^2-p^2 unacceptably noisy for a 0.14-0.50 GeV
        meson at E > 100 GeV (samples-owner note (a)).
        """
        out = []
        if which == 'K+':
            ctau = CTAU_KP
            mpar = M_KP
            specs = (('K+:piH', 'piH', M_PI_CHARGED), ('K+:enu', 'enu', None))
        elif which == 'K0L':
            ctau = CTAU_KL
            mpar = M_K0L
            # no K0L -> H e nu: the contract sets br_K_to_Henu('K0L') = 0 (there is no
            # physical two-body K0L leptonic decay to normalise against), so the history is
            # not created at all rather than carried with zero weight.
            specs = (('K0L:piH', 'piH', L.M_PI0),)
        else:
            ctau = CTAU_PI
            mpar = M_PI_CHARGED
            specs = (('pi+:enu', 'enu', None),)
        m = float(mpar)
        modes = getattr(self, '_modes_used', {})
        for wkey, mode, mrecoil in specs:
            if modes and not modes.get(wkey, True):
                continue
            if mode == 'piH':
                d = L.two_body(m, mrecoil, mH, rng)
                if d is None:
                    continue
                p4 = _boost_on_shell(d[1], mes, m)
            else:
                # three-body X -> H e nu with the exact phase space (m_e ~ 0, massless nu)
                d = L.nbody(m, [mH, M_E_SEMILEPTONIC, 0.0], rng)
                if d is None:
                    continue
                p4 = _boost_on_shell(d[0], mes, m)
            pz = float(mes['pz'])
            lam = (pz / m) * ctau      # signed parent flight length [m]
            v0 = (mes['x'], mes['y'], mes['z'])
            out.append((wkey, (p4, v0), lam))
        return out

    def _higgs_predecay_rows(self, smp):
        if 'status' not in smp:
            raise RuntimeError('higgs sample must carry the Pythia status column (11th)')
        cand = np.where(np.abs(smp['status']) == 62)[0]
        if len(cand) == 0:
            raise RuntimeError('higgs sample has a status column and no |status|=62 pre-decay copy')
        ev_c = smp['ev'][cand]
        order = np.argsort(ev_c, kind='stable')
        cand, ev_c = cand[order], ev_c[order]
        last = np.ones(len(cand), dtype=bool)
        last[:-1] = ev_c[1:] != ev_c[:-1]
        # keep one pre-decay copy per event -> the physical Higgs
        idxs = cand[last]
        proc = smp['pdg'][idxs]
        return np.sort(idxs), proc[np.argsort(idxs)]

    def _higgs_scalars(self, mH, smp, i, rng):
        px, py, pz = float(smp['px'][i]), float(smp['py'][i]), float(smp['pz'][i])
        mhig = M_HIGGS_SAMPLE
        if mhig < 2.0 * mH:
            return None
        d = L.two_body(mhig, mH, mH, rng)
        if d is None:
            return None
        e = math.sqrt(px * px + py * py + pz * pz + mhig * mhig)
        pair = L.boost(d, e, (px, py, pz))
        v0 = (smp['x'][i], smp['y'][i], smp['z'][i])
        out = []
        for a in pair:
            p2 = a[2] * a[2] + a[3] * a[3] + a[4] * a[4]
            out.append(((math.sqrt(p2 + a[0] * a[0]), a[2], a[3], a[4]), v0))
        return out

    # ---- weights ------------------------------------------------------------------
    def hist_weight_keys(self):
        """Possible wkey prefixes (for documentation/tests)."""
        return ['B:mix:B+', 'B:mix:B0', 'K+:piH', 'K+:enu', 'K0L:piH', 'K0L:enu',
                'pi+:enu', 'h:<process>']

    def history_weights(self, mH, tanb, table):
        """Production probability of every history (applied ONCE, outside its scalars' OR)."""
        br, _tot = self.decay_br(mH, tanb)
        sigma = self.sigma_pb(mH, tanb)
        proc = sigma['sigma_process']
        # normalise the process mixture by the SUM of the h production cross sections (WITHOUT
        # BR(h->HH), which lives in sigma['h']): then eps_h is the pure event-OR probability,
        # the same convention as eps_B / eps_inel and as BC5's published eps_H.
        s_tot = sigma['sigma_h_total']
        frac = table.get('h_frac', {})
        out = {}
        for wkey in set(table['hist_wkey']):
            if wkey.startswith('B:mix:'):
                bm = wkey.split(':')[2]
                w = self._b_total(mH, tanb, 'B0' if bm == 'B0' else 'B+')
            elif wkey == 'K+:piH':
                w = self.P.br_K_to_piH(mH, tanb, 'K+')
            elif wkey == 'K+:enu':
                w = self.P.br_K_to_Henu(mH, tanb, 'K+')
            elif wkey == 'K0L:piH':
                # PDG 311 is 50 % K_L: the 1/2 flavour projection belongs HERE (response layer),
                # the production module returns the per-K0L branching ratio.
                w = 0.5 * self.P.br_K_to_piH(mH, tanb, 'K0L')
            elif wkey == 'pi+:enu':
                w = self.P.br_pi_to_Henu(mH, tanb)
            elif wkey.startswith('h:'):
                g = wkey.split(':', 1)[1]
                if g not in proc or frac.get(g, 0.0) <= 0.0:
                    w = 0.0
                else:
                    # reweight the Pythia process mixture to the 2HDM sigma_i.  BR(h->HH)
                    # stays OUT of the history weight and enters s through sigma['h'].
                    w = proc[g] / (s_tot * frac[g])
            else:
                raise KeyError('unknown history key %r' % (wkey,))
            out[wkey] = float(w)
        return out

    def sigma_pb(self, mH, tanb):
        hs = self.P.higgs_sector(mH, tanb)
        return {'B': SIGMA_BB_PB, 'inel': SIGMA_INEL_PB,
                'h': float(hs['sigma_total']) * float(hs['BR_hHH']),
                'sigma_h_total': float(hs['sigma_total']),
                'sigma_process': {k: float(v) for k, v in hs['sigma_process'].items()},
                'BR_hHH_detail': float(hs['BR_hHH'])}

    REF_B_MESON_YIELD = 1.918   # generator truth: B+-, B0, Bs per pp -> b bbar event

    def yield_corr(self, samples):
        """B yield correction, derived at run time (Lead ruling 2026-10-25).

        Generator truth: 1.918 B+-, B0, Bs per pp -> b bbar event.  The sample's Bs
        fraction sets the truth B+-,B0 yield: truth = 1.918 * (1 - f_Bs); the correction
        is truth / (N_BpmB0 / N_gen) and is applied ONCE, on eps_B.  (It reduces to
        BC5's 1.918 / (all-b rows per event); because B_s has no production channel in
        the 2HDM, the missing B_s rate is re-injected as B+-,B0.  Recorded in meta.)
        """
        smp = samples.get('B')
        corr = {'B': 1.0, 'inel': 1.0, 'h': 1.0}
        if smp is None:
            return corr
        n_gen = float(samples.n_gen['B'])
        n_bp0 = int(np.isin(smp['pdg'], [521, -521, 511, -511]).sum())
        n_bs = int(np.isin(smp['pdg'], [531, -531]).sum())
        n_all = n_bp0 + n_bs
        if n_all > 0:
            f_bs = n_bs / float(n_all)
            truth_bp0 = self.REF_B_MESON_YIELD * (1.0 - f_bs)
            measured = n_bp0 / n_gen
            corr['B'] = truth_bp0 / measured if measured > 0 else 1.0
        self.b_yield_check = {'n_bp0': n_bp0, 'n_bs': n_bs, 'n_gen': n_gen,
                              'f_bs': (n_bs / float(n_all)) if n_all else 0.0,
                              'bp0_per_event': n_bp0 / n_gen,
                              'all_b_per_event': n_all / n_gen,
                              'reference_all_b_per_event': self.REF_B_MESON_YIELD,
                              'corr': corr['B']}
        if self.verbose:
            print('  [B yield] B+/B0/event = %.4f, f_Bs = %.4f -> truth B+/B0 = %.4f, '
                  'corr = %.4f' % (n_bp0 / n_gen, self.b_yield_check['f_bs'],
                                   self.REF_B_MESON_YIELD * (1 - self.b_yield_check['f_bs']),
                                   corr['B']))
        return corr

    def extra_meta(self, mH, tanb):
        """Provenance / validity payload attached to every scan point (Lead integration contract)."""
        import hashlib, os
        key = (round(float(mH), 9), round(float(tanb), 9))
        cache = getattr(self, '_extra_cache', None)
        if cache is None:
            cache = self._extra_cache = {}
        if key in cache:
            return dict(cache[key])
        out = {'b_mode': self.b_mode, 'b_mode_via': self._b_mode_how}
        try:
            hs = self.P.higgs_sector(mH, tanb)
            out.update({'BR_hHH': float(hs['BR_hHH']),
                        'Gamma_h_other_GeV': float(hs['Gamma_h_other']),
                        'g_C': float(hs['g_C']), 'g3': float(hs['g3']),
                        'sigma_h_total_pb': float(hs['sigma_total'])})
        except Exception as exc:
            out['higgs_sector_error'] = repr(exc)
        try:
            out['validity'] = self.P.validity_boundaries(mH)
        except Exception:
            out['validity'] = None
        out['br_B_to_XsH'] = self._b_total(mH, tanb, 'B+')
        for name, mod in (('production', self.P), ('decay', self.D)):
            f = getattr(mod, '__file__', None)
            if f and os.path.isfile(f):
                with open(f, 'rb') as fh:
                    out[name + '_sha256'] = hashlib.sha256(fh.read()).hexdigest()
        for name, rel in (('decay_binary',
                           os.path.join(RUN, 'work', 'decay', 'build', 'light_scalar_decay.exe')),
                          ('decay_binary_nix',
                           os.path.join(RUN, 'work', 'decay', 'build', 'light_scalar_decay'))):
            if os.path.isfile(rel):
                with open(rel, 'rb') as fh:
                    out[name + '_sha256'] = hashlib.sha256(fh.read()).hexdigest()
                out[name + '_path'] = rel
                break
        cache[key] = dict(out)
        return out


# ======================================================================================
# BC5 v3.6 reference model (regression only)
# ======================================================================================
class BC5Model:
    """BC5 dark scalar: mixing angle theta^2, strict-Winkler widths, Boiarska exclusive B."""

    name = 'bc5'
    param_name = 'theta2'
    # larger theta^2 -> Gamma_S = theta^2 Gamma_S(1) -> c*tau_S ~ theta^-2: SHORTER lifetime
    param_lifetime_sign = -1
    B_H = 0.01                  # BC5: BR(h -> SS) fixed at 0.01
    REF_B_MESON_YIELD = 1.918   # BC5 run_scan.REF_B_MESON_YIELD

    def __init__(self, verbose=True):
        import scalar_widths as sw
        import bc5_production as P
        self.sw = sw
        self.P = P
        self.verbose = verbose
        # the digitised Winkler hadronic table is loaded by bc5_lib; do the same here
        htab = os.path.join(BC5REF, 'widths', 'hadronic_table.csv')
        if os.path.isfile(htab):
            sw.load_hadronic_table(htab)
        self._br_cache = {}

    def channels(self, mS):
        w = self.sw.widths(mS)
        return [c for c, v in w.items() if v > 0.0]

    def open_channels(self, mS, theta2):
        return self.channels(mS)

    def decay_br(self, mS, theta2):
        key = round(float(mS), 9)
        got = self._br_cache.get(key)
        if got is None:
            w = self.sw.widths(mS)
            tot = float(sum(w.values()))
            got = ({c: float(v) / tot for c, v in w.items()}, tot)
            self._br_cache[key] = got
        return got

    def ctau_m(self, mS, theta2):
        _br, tot = self.decay_br(mS, 0.0)
        return HBAR_C / (tot * float(theta2))

    # ---- candidates ---------------------------------------------------------------
    def _b_mixing_scalars(self, mS, mes, bm, rng):
        w = self.P.br_B_to_S(mS, bm)
        tot = sum(w.values())
        if tot <= 0:
            return None
        label = _sample_categorical(w, rng)
        if label is None:
            return None
        mXp = self.P.B_CHANNEL_MASS[label]
        m = _p4_and_mass(mes['e'], mes['px'], mes['py'], mes['pz'])
        d = L.two_body(m, mXp, mS, rng)
        if d is None:
            return None
        p4 = _boost_to_lab(d[1], mes)
        ctau = {'B+': CTAU_BP, 'B0': CTAU_B0, 'Bs': CTAU_BS}[bm]
        v0 = _parent_decay_vertex(mes, ctau, rng)
        return (p4, v0)

    def _sample_quartic(self, mS, mes, kind, spec, channel, rng, wq):
        """Two scalars from one quartic history (polarisation: exactly BC5's generator)."""
        e, px, py, pz = mes['e'], mes['px'], mes['py'], mes['pz']
        m = _p4_and_mass(e, px, py, pz)
        if kind == '2b':
            d = L.two_body(m, mS, mS, rng)
            if d is None:
                return None
        else:
            if m < 2 * mS + spec:
                return None
            q = self.P.sample_dipion_mass(mS, m, spec, channel, rng)
            if q is None:
                return None
            dss = L.two_body(q, mS, mS, rng)
            dpar = L.two_body(m, q, spec, rng)
            if dss is None or dpar is None:
                return None
            d = [(a[0], a[1], a[2], a[3], a[4]) for a in
                 L.boost([(a[0], a[1], a[2], a[3], a[4]) for a in dss],
                         dpar[0][1], dpar[0][2:5])]
        ctau = {'B+': CTAU_BP, 'B0': CTAU_B0, 'Bs': CTAU_BS}[mes['bmeson']]
        v0 = _parent_decay_vertex(mes, ctau, rng)
        out = []
        for s4 in d:
            sl = L.boost([s4], e, (px, py, pz))[0]
            out.append(((sl[1], sl[2], sl[3], sl[4]), v0))
        return out

    def build_candidates(self, samples, mS, rng, np_rng, nsamp, mass_grid=None):
        hb = HistBuilder()
        # ---------------- B mesons: BC5 selects B+-, B0, Bs rows -----------------
        smp = samples.get('B')
        self._n_b_rows_all = 0
        if smp is not None:
            sel = np.isin(smp['pdg'], [521, -521, 511, -511, 531, -531])
            idxs = np.where(sel)[0]
            self._n_b_rows_all = len(idxs)
            n_tot = len(idxs)
            n_smp = min(int(nsamp['B']), n_tot)
            if n_smp > 0:
                idxs = np_rng.choice(idxs, size=n_smp, replace=False)
            hb.counts['B'] = [n_tot, n_smp]
            hb.samp_events['B'] = smp['ev'][idxs].copy()
            for ip, i in enumerate(idxs):
                pdg = int(smp['pdg'][i])
                bm = {521: 'B+', -521: 'B+', 511: 'B0', -511: 'B0', 531: 'Bs', -531: 'Bs'}[pdg]
                mes = {'e': smp['e'][i], 'px': smp['px'][i], 'py': smp['py'][i], 'pz': smp['pz'][i],
                       'x': smp['x'][i], 'y': smp['y'][i], 'z': smp['z'][i], 'bmeson': bm}
                scal = self._b_mixing_scalars(mS, mes, bm, rng)
                if scal is not None:
                    hb.add(smp['ev'][i], ip, _GROUP_ID['B'], 'B:mix:' + bm, [scal])
                qmodes = [q for q in self.P.quartic_channels(mS) if q[0] == bm]
                if qmodes:
                    wq = sum(q[3] for q in qmodes)
                    u = rng.random(); c = 0.0; chosen = qmodes[-1]
                    for q in qmodes:
                        c += q[3] / wq
                        if u <= c:
                            chosen = q; break
                    scal = self._sample_quartic(mS, mes, chosen[1], chosen[2], chosen[4], rng, wq)
                    if scal is not None:
                        hb.add(smp['ev'][i], ip, _GROUP_ID['B'], 'B:q:' + bm, scal)
        # ---------------- kaons (BC5 minbias sample carries K only) --------------
        smp = samples.get('inel')
        if smp is not None:
            base = 0
            for pdg, which in ((321, 'K+'), (-321, 'K+'), (311, 'K0L'), (-311, 'K0L')):
                idxs = np.where(smp['pdg'] == pdg)[0]
                n_tot = len(idxs)
                n_smp = min(int(nsamp['inel']), n_tot)
                if n_smp > 0:
                    idxs = np_rng.choice(idxs, size=n_smp, replace=False)
                hb.counts['inel'][0] += n_tot
                hb.counts['inel'][1] += n_smp
                hb.samp_events['inel'].extend(smp['ev'][idxs].tolist())
                k0_proj = 0.5 if which == 'K0L' else 1.0
                for ip, i in enumerate(idxs):
                    mes = {'e': smp['e'][i], 'px': smp['px'][i], 'py': smp['py'][i], 'pz': smp['pz'][i],
                           'x': smp['x'][i], 'y': smp['y'][i], 'z': smp['z'][i]}
                    m = _p4_and_mass(mes['e'], mes['px'], mes['py'], mes['pz'])
                    br = self.P.br_K_to_S(mS, which)
                    if br > 0:
                        mP = M_PI_CHARGED if which == 'K+' else L.M_PI0
                        d = L.two_body(m, mP, mS, rng)
                        if d is not None:
                            p4 = _boost_to_lab(d[1], mes)
                            ctau = CTAU_KP if which == 'K+' else CTAU_KL
                            v0 = _parent_decay_vertex(mes, ctau, rng)
                            hb.add(smp['ev'][i], base + ip, _GROUP_ID['inel'],
                                   'K:%s:ko' % which, [(p4, v0)])
                    for imode, (k_, mode, wq) in enumerate(self.P.quartic_kaon_channels(mS, which)):
                        kind = '2b' if mode == 'SS' else '3b'
                        spec = 0.0 if mode == 'SS' else (L.M_PI0 if mode == 'SSpi0' else L.M_PI_CHARGED)
                        mes2 = dict(mes, bmeson=None)
                        mes2['kaon'] = which
                        ctau = CTAU_KP if which == 'K+' else CTAU_KL
                        sc = self._quartic_kaon_scalars(mS, mes, kind, spec, ctau, rng)
                        if sc is not None:
                            hb.add(smp['ev'][i], base + ip, _GROUP_ID['inel'],
                                   'K:%s:q%d' % (which, imode), sc)
                base += n_smp
            hb.samp_events['inel'] = np.asarray(hb.samp_events['inel'], dtype=np.int64)
        # ---------------- Higgs -------------------------------------------------
        smp = samples.get('h')
        if smp is not None:
            cand = np.where(np.abs(smp['status']) == 62)[0]
            ev_c = smp['ev'][cand]
            order = np.argsort(ev_c, kind='stable')
            cand, ev_c = cand[order], ev_c[order]
            last = np.ones(len(cand), dtype=bool)
            last[:-1] = ev_c[1:] != ev_c[:-1]
            idxs = np.sort(cand[last])
            n_smp = min(int(nsamp['h']), len(idxs))
            keep = np_rng.choice(len(idxs), size=n_smp, replace=False) if n_smp < len(idxs) else np.arange(len(idxs))
            keep = np.sort(keep)
            idxs = idxs[keep]
            hb.counts['h'] = [len(cand), n_smp]
            hb.samp_events['h'] = smp['ev'][idxs].copy()
            groups = np.array([HIGGS_PROC_GROUP.get(int(c), 'other') for c in smp['pdg'][idxs]])
            for ip, i in enumerate(idxs):
                e, px, py, pz = smp['e'][i], smp['px'][i], smp['py'][i], smp['pz'][i]
                mhig = _p4_and_mass(e, px, py, pz)
                d = L.two_body(mhig, mS, mS, rng)
                if d is None:
                    continue
                pair = L.boost(d, e, (px, py, pz))
                v0 = (smp['x'][i], smp['y'][i], smp['z'][i])
                hb.add(smp['ev'][i], ip, _GROUP_ID['h'], 'h:' + str(groups[ip]),
                       [((a[1], a[2], a[3], a[4]), v0) for a in pair])
            hb.h_frac = {g: float(np.mean(groups == g)) for g in sorted(set(groups.tolist()))}
        else:
            hb.h_frac = {}
        out = hb.finalize()
        return out

    def _quartic_kaon_scalars(self, mS, mes, kind, spec, ctau, rng):
        e, px, py, pz = mes['e'], mes['px'], mes['py'], mes['pz']
        m = _p4_and_mass(e, px, py, pz)
        if kind == '2b':
            d = L.two_body(m, mS, mS, rng)
            if d is None:
                return None
        else:
            if m < 2 * mS + spec:
                return None
            q = self.P.sample_dipion_mass(mS, m, spec, None, rng)
            if q is None:
                return None
            dss = L.two_body(q, mS, mS, rng)
            dpar = L.two_body(m, q, spec, rng)
            if dss is None or dpar is None:
                return None
            d = [(a[0], a[1], a[2], a[3], a[4]) for a in
                 L.boost([(a[0], a[1], a[2], a[3], a[4]) for a in dss],
                         dpar[0][1], dpar[0][2:5])]
        v0 = _parent_decay_vertex(mes, ctau, rng)
        out = []
        for s4 in d:
            sl = L.boost([s4], e, (px, py, pz))[0]
            out.append(((sl[1], sl[2], sl[3], sl[4]), v0))
        return out

    # ---- weights ------------------------------------------------------------------
    def history_weights(self, mS, theta2, table):
        th2 = float(theta2)
        out = {}
        for wkey in set(table['hist_wkey']):
            if wkey.startswith('B:mix:'):
                bm = wkey.split(':')[2]
                w = th2 * sum(self.P.br_B_to_S(mS, bm).values())
            elif wkey.startswith('B:q:'):
                bm = wkey.split(':')[2]
                w = sum(q[3] for q in self.P.quartic_channels(mS) if q[0] == bm)
            elif wkey.startswith('K:') and wkey.endswith(':ko'):
                which = wkey.split(':')[1]
                k0 = 0.5 if which == 'K0L' else 1.0
                w = th2 * self.P.br_K_to_S(mS, which) * k0
            elif wkey.startswith('K:') and ':q' in wkey:
                parts = wkey.split(':')
                which = parts[1]
                imode = int(parts[2][1:])
                k0 = 0.5 if which == 'K0L' else 1.0
                modes = self.P.quartic_kaon_channels(mS, which)
                w = (modes[imode][2] if imode < len(modes) else 0.0) * k0
            elif wkey.startswith('h:'):
                g = wkey.split(':', 1)[1]
                frac = table.get('h_frac', {})
                proc = self.P.SIGMA_H
                # BC5 reweighting: (target fraction)/(sample fraction), target = SIGMA_H
                if g in proc and frac.get(g, 0.0) > 0.0:
                    # BR(h->SS) = B_H enters through sigma_pb['h'], as in BC5 run_scan
                    w = proc[g] / (self.P.SIGMA_H_TOTAL * frac[g])
                else:
                    w = 0.0
            else:
                raise KeyError('unknown history key %r' % (wkey,))
            out[wkey] = float(w)
        return out

    def sigma_pb(self, mS, theta2):
        return {'B': SIGMA_BB_PB, 'inel': SIGMA_INEL_PB,
                'h': self.P.SIGMA_H_TOTAL * self.B_H,
                'sigma_h_total': self.P.SIGMA_H_TOTAL,
                'sigma_process': dict(self.P.SIGMA_H), 'BR_hHH_detail': self.B_H}

    def yield_corr(self, samples):
        """Exactly BC5 run_scan: 1.918 / (B+-, B0, Bs rows per generated event)."""
        smp = samples.get('B')
        corr = {'B': 1.0, 'inel': 1.0, 'h': 1.0}
        if smp is None:
            return corr
        n_gen = samples.n_gen['B']
        sel = np.isin(smp['pdg'], [521, -521, 511, -511, 531, -531])
        per_event = float(sel.sum()) / float(n_gen)
        corr['B'] = self.REF_B_MESON_YIELD / per_event if per_event > 0 else 1.0
        self.b_yield_check = {'per_event_all_b': per_event, 'corr': corr['B']}
        if self.verbose:
            print('  [BC5 B yield] all b mesons per event = %.4f -> B_YIELD_CORR = %.4f'
                  % (per_event, corr['B']))
        return corr

    def extra_meta(self, mS, theta2):
        return {}
