# -*- coding: utf-8 -*-
"""Unit tests for the response library, the probability algebra and the event engine.

Run:  python test_engine.py            (from work/analysis; needs the BC5 reference copy
                                        in work/analysis/bc5ref for the kinematics identity test)

Covers DESIGN.md section 4 cross-checks (i)-{kernel}, (ii) p_window vs analytic A=1,
(iii) event-OR algebra toy test, plus:
  * transport fast path vs the public downstream_ctypes API (bitwise)
  * per-channel decay kinematics IDENTICAL to bc5_lib for the same random stream
    (this is the "decay kinematics per channel must be exactly BC5's" requirement)
  * response-table invariants (sorted histories/scalars, unique sampled collisions)
"""
from __future__ import annotations
import math, os, random, subprocess, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(RUN, 'work', 'response'))

import twohdm_lib as L
import models as M
import response_table as RT
import event_engine as E

RESULTS = []


def check(name, ok, detail=''):
    RESULTS.append((name, bool(ok)))
    print(('PASS: ' if ok else 'FAIL: ') + name + (('  [' + str(detail) + ']') if detail else ''))


# ======================================================================================
# 1. kernel validation suite (copied verbatim from BC5)
# ======================================================================================
def test_kernel_suite():
    path = os.path.join(RUN, 'work', 'response', 'test_kernel.py')
    p = subprocess.run([sys.executable, path], cwd=os.path.dirname(path),
                       capture_output=True, text=True)
    out = p.stdout
    npass = int(out.split('SUMMARY: ')[1].split('/')[0])
    ntot = int(out.split('SUMMARY: ')[1].split('/')[1].split()[0])
    # a2 is the documented spec conflict (axis track -> BeamHole), identical to BC5
    check('kernel suite: %d/%d checks pass (a2 = documented spec conflict)' % (npass, ntot),
          npass == ntot - 1 and 'FAILED: a2' in out)


# ======================================================================================
# 2. transport fast path vs public API
# ======================================================================================
def test_transport_fastpath():
    rng = random.Random(20261025)
    worst = 0
    for _ in range(300):
        x, y, z = rng.uniform(-3, 3), rng.uniform(-2, 2), rng.uniform(0, 2.4)
        px, py, pz = rng.gauss(0, 1), rng.gauss(0, 1), rng.uniform(1, 60)
        e = math.sqrt(px * px + py * py + pz * pz + 0.25)
        ch = rng.choice([-1.0, 1.0])
        a = L.passes_pair([(0.0, e, px, py, pz, ch, None)], x, y, z)
        b = L.passes_pair([(0.0, e, px, py, pz, ch, None)], x, y, z)
        worst = max(worst, abs(a - b))
        if abs(a - b) > 0:
            break
    # compare the fast path with the public transport_track on an opposite-charge pair
    x, y, z = 0.4, 0.2, 0.1
    px, py, pz = 0.7, -0.3, 12.0
    e = math.sqrt(px * px + py * py + pz * pz)
    daughters = [(0.0, e, px, py, pz, 1.0, None), (0.0, e, -px, -py, pz, -1.0, None)]
    fast = L.passes_pair(daughters, x, y, z)
    okp, sp, _f, _e = L.dk.transport_track(x, y, z, px, py, pz, e, 1.0)
    okm, sm, _f, _e = L.dk.transport_track(x, y, z, -px, -py, pz, e, -1.0)
    expect = 1.0 if (okp and okm) else 0.0
    check('transport fast path equals the public API', fast == expect,
          'fast=%r expect=%r status=(%s,%s)' % (fast, expect, sp, sm))


# ======================================================================================
# 3. p_window: analytic A = 1 and piecewise-linear A vs quadrature
# ======================================================================================
def test_p_window():
    z = L.ZGRID
    A1 = np.ones((1, len(z)))
    ms, pz = 1.0, 20.0
    # The lifetime integral runs over the z-node range [1.0, 2.5] m only (A(z) is not defined
    # outside it), so for a birth point z0 < 1 the survival factor exp(-alpha*(1-z0)) multiplies
    # the in-window integral.  Fill in the analytic answer on that domain.
    ctau = 0.5
    alpha = ms / (ctau * pz)
    z0 = 0.0
    expect = math.exp(-alpha * (z[0] - z0)) * (1.0 - math.exp(-alpha * (z[-1] - z[0])))
    got = L.p_window(ms, np.array([pz]), np.array([z0]), A1, ctau)[0]
    check('p_window vs analytic A=1 (z0 before window)', abs(got / expect - 1) < 1e-12,
          'got %.15g expect %.15g' % (got, expect))
    z0 = 1.7
    expect = 1.0 - math.exp(-alpha * (z[-1] - z0))
    got = L.p_window(ms, np.array([pz]), np.array([z0]), A1, ctau)[0]
    check('p_window vs analytic A=1 (z0 inside window)', abs(got / expect - 1) < 1e-12,
          'got %.15g expect %.15g' % (got, expect))
    got = L.p_window(ms, np.array([pz]), np.array([3.0]), A1, ctau)[0]
    check('p_window is zero for a birth point after the window', got == 0.0, 'got %r' % got)
    # short lifetime (alpha*h >> 1): the exact first-cell answer, not 0.5*alpha*h
    ctau = 1e-4
    alpha = ms / (ctau * pz)
    expect = math.exp(-alpha * (z[0] - 0.0)) * (1.0 - math.exp(-alpha * (z[-1] - z[0])))
    got = L.p_window(ms, np.array([pz]), np.array([0.0]), A1, ctau)[0]
    check('p_window exact for a short lifetime (alpha*h >> 1)',
          abs(got / expect - 1) < 1e-10, 'got %.12g expect %.12g' % (got, expect))
    # random piecewise-linear A vs a fine numerical quadrature on the node domain
    rng = np.random.default_rng(7)
    worst = 0.0
    for _ in range(20):
        A = rng.random((1, len(z)))
        ctau = 10.0 ** rng.uniform(-3, 1)
        alpha = ms / (ctau * pz)
        for z0 in (1.0, 1.3, 1.7, 2.2):
            got = L.p_window(ms, np.array([pz]), np.array([z0]), A, ctau)[0]
            zz = np.linspace(z0, z[-1], 200001)
            Av = np.interp(zz, z, A[0])
            num = np.trapezoid(alpha * np.exp(-alpha * (zz - z0)) * Av, zz)
            if num > 1e-12:
                worst = max(worst, abs(got / num - 1))
    check('p_window piecewise-linear A vs numerical quadrature', worst < 2e-4,
          'worst rel diff %.2e' % worst)


# ======================================================================================
# 4. per-channel decay kinematics identical to bc5_lib
# ======================================================================================
_DECAY_PAIRS = [('ee', 'e'), ('mumu', 'mu'), ('tautau', 'tau'), ('pipi', 'pipi'),
                ('KK', 'KK'), ('4pi', 'multi'), ('cc', 'cc'), ('ss', 'ss'), ('gg', 'gg')]


def _compare_decays(ms, channels, ndraw=3):
    """Draw decays with mine and with bc5_lib from the SAME RNG stream and compare exactly."""
    import bc5_lib as B5          # staged in work/analysis/bc5ref/analysis
    bad = []
    for mine, theirs in _DECAY_PAIRS:
        if mine not in channels:
            continue
        for _ in range(ndraw):
            r1 = random.Random(424242)
            r2 = random.Random(424242)
            a = L.decay_products(ms, mine, r1)
            b = B5.s_decay_products(ms, theirs, r2)
            aa = [(round(x[0], 12), round(x[1], 9), round(x[2], 9), round(x[3], 9),
                   round(x[4], 9), x[5], None if x[6] is None else
                   tuple(round(v, 12) for v in x[6])) for x in (a or [])]
            bb = [(round(x[0], 12), round(x[1], 9), round(x[2], 9), round(x[3], 9),
                   round(x[4], 9), x[5], None if x[6] is None else
                   tuple(round(v, 12) for v in x[6])) for x in (b or [])]
            if aa != bb:
                bad.append((mine, len(aa), len(bb)))
                break
    return bad


def test_decay_identity():
    """The per-channel decay kinematics must be EXACTLY bc5_lib's (DESIGN section 4)."""
    import bc5_lib as B5
    saved = L.HADRON_LIB
    libdir = os.path.join(RUN, '..', 'bc5_lhcb_downstream_v3', 'work', 'magnus', 'hadron_lib')
    # (a) analytic fallbacks: every channel is defined analytically
    L.HADRON_LIB = None
    B5.HADRON_LIB = None
    allch = [m for m, _t in _DECAY_PAIRS]
    for ms in (1.0, 3.5):
        bad = _compare_decays(ms, allch)
        check('decay kinematics identical to bc5_lib (analytic fallback, mH=%.2f)' % ms, not bad,
              ('channels differing: %r' % bad) if bad else '')
    # (b) with the Pythia library loaded in BOTH modules, on a mass that is a library node.
    #     Channels whose library has no node at this mass are masked out by the engine
    #     (their BR is zero there); bc5_lib would silently fall back to its analytic model,
    #     which is exactly the behaviour the per-channel design must avoid.
    if os.path.isdir(libdir):
        if L.load_hadron_lib(libdir) and B5.load_hadron_lib(libdir):
            for ms in (3.5567,):
                chans = [m for m, _t in _DECAY_PAIRS if m not in ('cc',)]
                bad = _compare_decays(ms, chans)
                check('decay kinematics identical to bc5_lib (with library, mH=%.4f)' % ms,
                      not bad, ('channels differing: %r' % bad) if bad else '')
        else:
            print('  (BC5 hadron library not available: skipping the with-library test)')
    L.HADRON_LIB = saved
    B5.HADRON_LIB = None


# ======================================================================================
# 5. event-OR algebra toy test
# ======================================================================================
def _toy_table():
    """Two collisions:
       event 10: parent 0 -> two EXCLUSIVE histories (h0: 1 scalar, h1: 1 scalar)
       event 20: parent 1 -> one history h2 (2 scalars) ; parent 2 -> one history h3 (1 scalar)
    """
    z = L.ZGRID
    nz = len(z)
    A = np.zeros((5, 1, nz))
    A[0, 0, :] = 1.0; A[1, 0, :] = 1.0
    A[2, 0, :] = 0.6; A[3, 0, :] = 0.4; A[4, 0, :] = 1.0 * (z > 1.7)
    tab = {
        'channels': ['ch'], 'A': A.astype(np.float32),
        'pz': np.full(5, 20.0), 'z0': np.array([0.5, 0.6, 0.7, 0.8, 0.9]),
        'scal_hist': np.array([0, 1, 2, 2, 3], dtype=np.int32),
        'hist_event': np.array([10, 10, 20, 20], dtype=np.int64),
        'hist_parent': np.array([0, 0, 1, 2], dtype=np.int64),
        'hist_group': np.zeros(4, dtype=np.int8),
        'hist_wkey': ['w0', 'w1', 'w2', 'w3'],
        'samp_events': {'B': np.array([10, 20], dtype=np.int64), 'inel': np.array([], dtype=np.int64),
                        'h': np.array([], dtype=np.int64)},
        'n_gen': {'B': 2, 'inel': 0, 'h': 0},
        'yield_corr': {'B': 1.0, 'inel': 1.0, 'h': 1.0},
        'counts': {'B': [2, 2], 'inel': [0, 0], 'h': [0, 0]},
    }
    w = {'w0': 0.5, 'w1': 0.3, 'w2': 0.9, 'w3': 0.7}
    return tab, w


def test_or_algebra():
    tab, w = _toy_table()
    p = np.array([0.2, 0.4, 0.1, 0.3, 0.25])
    ev, o = E.event_responses(tab, p, w)
    # brute force
    def hist_or(ps):
        r = 1.0
        for q in ps:
            r *= (1.0 - q)
        return 1.0 - r
    ev_brute = {}
    # event 10
    par0 = w['w0'] * hist_or([p[0]]) + w['w1'] * hist_or([p[1]])
    ev_brute[10] = hist_or([par0])
    # event 20: parent 1 (one history, 2 scalars), parent 2 (one history)
    par1 = w['w2'] * hist_or([p[2], p[3]])
    par2 = w['w3'] * hist_or([p[4]])
    ev_brute[20] = hist_or([par1, par2])
    got = dict(zip(ev.tolist(), o.tolist()))
    ok = (set(got) == set(ev_brute)
          and all(abs(got[k] - ev_brute[k]) < 1e-12 for k in ev_brute))
    check('event OR algebra vs brute force', ok,
          'mine %r brute %r' % ({k: round(v, 12) for k, v in got.items()},
                                {k: round(v, 12) for k, v in ev_brute.items()}))
    eff = E.group_efficiency(tab, ev, o, tab['n_gen'], tab['yield_corr'], counts=tab['counts'])
    eps_B = eff['B'][0]
    expect = L.RECO_EFF * (tab['counts']['B'][0] / tab['counts']['B'][1]) / tab['n_gen']['B'] * (
        ev_brute[10] + ev_brute[20])
    check('eps normalisation (0.4 once, n_tot/n_smp, 1/N_gen)', abs(eps_B - expect) < 1e-12,
          'got %.12g expect %.12g' % (eps_B, expect))
    # sensitivity: a parent-id collision across events must NOT merge (review BUG C)
    tab2 = dict(tab)
    # collision event 20's first parent and event 10's parent now share index 0; a grouping
    # on hist_parent alone (the review BUG C) would merge the two collisions' histories.
    tab2['hist_parent'] = np.array([0, 0, 0, 1], dtype=np.int64)
    ev2, o2 = E.event_responses(tab2, p, w)
    got2 = dict(zip(ev2.tolist(), o2.tolist()))
    check('events with equal parent indices are not merged', abs(got2[10] - ev_brute[10]) < 1e-12
          and abs(got2[20] - ev_brute[20]) < 1e-12, 'got %r' % ({k: round(v, 12) for k, v in got2.items()}))


# ======================================================================================
# 6. response-table invariants (uses the 2HDM is not available -> synthetic model)
# ======================================================================================
class _TinyModel:
    """Synthetic model used only to exercise the table builder's ordering invariants."""
    name = 'tiny'
    param_name = 'x'
    param_lifetime_sign = +1

    def channels(self, mH):
        return ['pipi', 'KK']

    def open_channels(self, mH, x):
        return self.channels(mH)

    def decay_br(self, mH, x):
        return {'pipi': 0.6, 'KK': 0.4}, 1.0

    def ctau_m(self, mH, x):
        return 0.5

    def sigma_pb(self, mH, x):
        return {'B': 1.0, 'inel': 1.0, 'h': 1.0}

    def yield_corr(self, samples):
        return {'B': 1.0, 'inel': 1.0, 'h': 1.0}

    def history_weights(self, mH, x, table):
        return {k: 1.0 for k in set(table['hist_wkey'])}

    def build_candidates(self, samples, mH, rng, np_rng, nsamp, mass_grid=None):
        hb = M.HistBuilder()
        # event 5: parent 0 with 2 histories; event 3: parent 1
        for (ev, par, wk, nsc) in ((5, 0, 'a', 1), (5, 0, 'b', 2), (3, 1, 'a', 1)):
            scal = []
            for _ in range(nsc):
                px, py, pz = rng.gauss(0, 0.1), rng.gauss(0, 0.1), rng.uniform(2, 40)
                e = math.sqrt(px * px + py * py + pz * pz + mH * mH)
                scal.append(((e, px, py, pz), (rng.uniform(-.01, .01), rng.uniform(-.01, .01),
                                               rng.uniform(-.01, .01))))
            hb.add(ev, par, 0, wk, scal)
        hb.samp_events['B'] = [5, 3]
        hb.counts['B'] = [100, 2]
        return hb.finalize()


def test_table_invariants():
    model = _TinyModel()
    tab = RT.build_table(model, 1.0, samples=None, verbose=False)
    # histories must be ordered by the PAIR (event, parent); parents are only unique inside a
    # group, so a global ordering of hist_parent is neither required nor produced.
    key = list(zip(tab['hist_event'].tolist(), tab['hist_parent'].tolist()))
    ok = (bool(np.all(np.diff(tab['scal_hist']) >= 0))
          and key == sorted(key)
          and bool(np.array_equal(tab['samp_events']['B'], np.array([3, 5]))))
    check('table invariants: sorted scalars/histories/events, unique sorted samp_events', ok,
          'scal_hist=%r hist_event=%r parents=%r samp=%r' % (
              tab['scal_hist'].tolist(), tab['hist_event'].tolist(),
              tab['hist_parent'].tolist(), tab['samp_events']['B'].tolist()))


# ======================================================================================
def main():
    print('-- 1. kernel suite')
    test_kernel_suite()
    print('-- 2. transport fast path')
    test_transport_fastpath()
    print('-- 3. lifetime integral')
    test_p_window()
    print('-- 4. decay kinematics identity with bc5_lib')
    test_decay_identity()
    print('-- 5. event OR algebra')
    test_or_algebra()
    print('-- 6. table invariants')
    test_table_invariants()
    npass = sum(1 for _n, ok in RESULTS if ok)
    print('\nSUMMARY: %d/%d checks passed' % (npass, len(RESULTS)))
    for n, ok in RESULTS:
        if not ok:
            print('  FAILED: ' + n)
    return 0 if npass == len(RESULTS) else 1


if __name__ == '__main__':
    sys.exit(main())
