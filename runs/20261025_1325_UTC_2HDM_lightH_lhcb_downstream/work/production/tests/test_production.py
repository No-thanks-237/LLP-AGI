# -*- coding: utf-8 -*-
"""Unit tests for work/production/twohdm_production.py (+ b_channels.py).

Run either way:
    python work/production/tests/test_production.py
    python -m pytest work/production/tests/test_production.py
Self-contained (no pytest dependency required).
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROD = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, PROD)
import twohdm_production as P   # noqa: E402
import b_channels              # noqa: E402

WOLFRAM_G = {
    'u': {'g0': -3.36549386575e-07, 'g1': -0.750000000000, 'g2': 0.0351515305077,
          'X1': 3.51515305077, 'X2': -3.19741896996e-10},
    'c': {'g0': -3.36526819507e-07, 'g1': -0.750000005029, 'g2': 0.0350966706407,
          'X1': 3.50966706910, 'X2': -5.02910117637e-05},
    't': {'g0': -2.73828701384e-07, 'g1': -0.750013977299, 'g2': 0.0105130174230,
          'X1': 1.05131571960, 'X2': -0.139772989813},
}


def test_couplings():
    for tb in (5.0, 10.0, 100.0, 1e3, 1e4, 1e5):
        cp = P.couplings(tb)
        c, s = cp['c'], cp['s']
        assert abs(c * tb - 1.0) < 1e-12
        assert abs(c * c + s * s - 1.0) < 1e-14
        assert abs(cp['xi_H_V'] - c) < 1e-15
        assert abs(cp['xi_h_V'] - s) < 1e-15
        # xi_H^f = c - s/tb must equal the cancellation-free 1/(tb^3(1+s))
        rel = abs(cp['xi_H_f'] - cp['xi_H_f_stable']) / cp['xi_H_f_stable']
        assert rel < 1e-5, (tb, rel)
        assert abs(cp['xi_H_f'] - (c - s / tb)) < 1e-18
        assert abs(cp['xi_h_f'] - (s + c / tb)) < 1e-15
        assert 0.0 <= c <= 1.0 and 0.0 <= s <= 1.0


def test_ckm():
    V = P.V_CKM
    assert np.max(np.abs(V.conj().T @ V - np.eye(3))) < 1e-14
    assert np.max(np.abs(V @ V.conj().T - np.eye(3))) < 1e-14
    # GIM: equal-mass-weighted sums must vanish for i != j
    for i in range(3):
        for j in range(3):
            if i != j:
                assert abs(np.sum(np.conj(V[:, i]) * V[:, j])) < 1e-15
    mag = np.abs(V)
    # DESIGN nominal magnitudes (the only sizeable difference is |Vtd|, documented)
    nom = {(0, 0): 0.97373, (0, 1): 0.2243, (0, 2): 0.00382, (1, 0): 0.221,
           (1, 1): 0.975, (1, 2): 0.0408, (2, 0): 0.0086, (2, 1): 0.0415, (2, 2): 0.999}
    for (i, j), m in nom.items():
        assert abs(mag[i, j] - m) / m < 0.04, (i, j, mag[i, j], m)
    alt = P.ckm_matrix_from_magnitudes()
    for (i, j), m in nom.items():
        assert abs(abs(alt[i, j]) - m) < 1e-12


def test_g_functions_vs_wolfram():
    G = P.g_aux_functions(100.0)
    for q in ('u', 'c', 't'):
        for fn in ('g0', 'g1', 'g2', 'X1', 'X2'):
            ref = WOLFRAM_G[q][fn]
            assert abs(G[q][fn] - ref) / max(abs(ref), 1e-30) < 1e-10, (q, fn)


def test_xi_limits():
    # xi_H^{ij} -> 0 like 1/tanb as tanb -> inf
    a = abs(P.xi_H(1.0, 1e4, 's', 'b'))
    b = abs(P.xi_H(1.0, 1e5, 's', 'b'))
    assert abs(b / a - 0.1) < 1e-3, b / a
    for tb in (5.0, 1e5):
        for q in ('u', 'c', 't'):
            g = P.g_aux_functions(tb)[q]
            for fn in ('g0', 'g1', 'g2'):
                assert math.isfinite(g[fn])
    # magnitude of xi at large tanb is < 1e-8
    assert abs(P.xi_H(1.0, 1e5, 's', 'b')) < 1e-8
    # all three down-type pairs finite/positive magnitude
    for (i, j) in (('d', 's'), ('d', 'b'), ('s', 'b')):
        assert abs(P.xi_H(1.0, 100.0, i, j)) > 0.0


def test_lambda_HH():
    for tb, mH in ((100.0, 1.0), (5.0, 2.0)):
        c = 1.0 / tb
        assert abs(P.lambda_HHpmHpm(mH, tb) - (-(2 * 600.0 ** 2 - mH ** 2) * c / 246.0)) < 1e-10


def test_br_B_nominal():
    for tb in (5.0, 100.0, 1e4):
        for mH in (0.021, 0.5, 1.0, 3.0):
            d = P.br_B_to_XsH(mH, tb)
            assert d['total'] > 0.0
            s = sum(d['channels'].values())
            assert abs(s - d['total']) <= 1e-12 * max(d['total'], 1.0), (mH, tb, s, d['total'])
            assert all(v >= 0.0 for v in d['channels'].values())
    # B+ and B0: identical within the meson-mass difference of the phase-space factor
    # (4e-6 at mH=1); identical normalised channel fractions in every mode
    for mode in P.B_MODES:
        a = P.br_B_to_XsH(1.0, 100.0, 'B+', mode=mode)
        b = P.br_B_to_XsH(1.0, 100.0, 'B0', mode=mode)
        # Gamma(B->XsH) meson independent -> BR(B0) = (Gamma_B+/Gamma_B0) BR(B+) in
        # every mode
        # (residuals: 8e-6 for foresee/paper_mb from the per-meson m_B in the phase-space
        # factor; <=2e-4 for exclusive from the K+- vs K0 recoil mass in the channel sum)
        assert abs(b['total'] / a['total'] - P.B0_LIFETIME_SCALE) < 5e-4, (mode, b['total'] / a['total'])
        assert abs(P.B0_LIFETIME_SCALE - 0.92735) < 1e-4
        fa = {k: v / a['total'] for k, v in a['channels'].items()}
        fb = {k: v / b['total'] for k, v in b['channels'].items()}
        # the kinematic mixture is deliberately NOT identical for B+- and B0: the lightest
        # strangeness-carrying recoil is K+- (0.493677) for B+- but K0 (0.497611) for B0
        # (review F2), which shifts the K channel phase space by ~1e-4 and renormalises all
        # fractions.  The B-/B0 recoil difference must stay far below the 1% level.
        for k in fa:
            assert abs(fa[k] - fb[k]) < 3e-4, (k, fa[k], fb[k])
        assert abs(sum(b['channels'].values()) - b['total']) <= 1e-12 * max(b['total'], 1.0)
    # mode 'paper_mb' closes at m_b, 'foresee' at the B-meson mass
    assert P.br_B_to_XsH(4.18, 100.0, mode='paper_mb')['total'] == 0.0
    assert P.br_B_to_XsH(4.5, 100.0, mode='paper_mb')['total'] == 0.0
    # mesonic endpoint in the foresee mode: mH >= mB - mK -> 0 and flagged
    endp = P.M_BP - P.M_KP
    assert abs(endp - 4.785663) < 1e-5
    assert P.br_B_to_XsH(endp - 1e-4, 100.0, mode='foresee')['total'] > 0.0
    assert P.br_B_to_XsH(endp, 100.0, mode='foresee')['total'] == 0.0
    # note: mH=4.5 is still BELOW the mesonic endpoint (4.7857), so B->K H is open
    assert P.br_B_to_XsH(4.5, 100.0, mode='foresee')['total'] > 0.0
    assert P.br_B_to_XsH(4.8, 100.0, mode='foresee')['total'] == 0.0
    assert P.br_B_to_XsH(5.28, 100.0, mode='foresee')['total'] == 0.0
    assert P.br_B_to_XsH(endp, 100.0, mode='foresee')['beyond_mesonic_endpoint'] is True
    assert P.br_B_to_XsH(4.0, 100.0, mode='foresee')['beyond_mesonic_endpoint'] is False
    assert P.br_B_to_XsH(4.5, 100.0, mode='exclusive')['total'] > 0.0
    # every mode keeps sum(channels) == total
    for mode in P.B_MODES:
        for mH in (0.1, 1.0, 4.0, 4.5, 5.0):
            d = P.br_B_to_XsH(mH, 100.0, mode=mode)
            assert abs(sum(d['channels'].values()) - d['total']) <= 1e-12 * max(d['total'], 1.0)
    # foresee and paper_mb coincide at mH -> 0 and differ at high mass
    a = P.br_B_to_XsH(0.021, 100.0, mode='foresee')['total']
    b = P.br_B_to_XsH(0.021, 100.0, mode='paper_mb')['total']
    assert abs(a / b - 1.0) < 2e-3
    a = P.br_B_to_XsH(3.0, 100.0, mode='foresee')['total']
    b = P.br_B_to_XsH(3.0, 100.0, mode='paper_mb')['total']
    assert a / b > 1.5
    # monotonic decreasing in tanb
    brs = [P.br_B_to_XsH(0.1, t)['total'] for t in (5.0, 10.0, 30.0, 100.0, 1e3, 1e4, 1e5)]
    assert all(brs[i] > brs[i + 1] for i in range(len(brs) - 1))
    # ~1/tanb^2 at large tanb
    r = P.br_B_to_XsH(0.1, 1e5, mode='foresee')['total'] / P.br_B_to_XsH(0.1, 1e4, mode='foresee')['total']
    assert abs(r - 0.01) < 1e-3 * 0.01 / 0.01


def test_br_B_exclusive_variant():
    for tb in (10.0, 100.0, 1e4):
        for mH in (0.021, 1.0, 4.5, 4.78):
            d = P.br_B_to_XsH_exclusive(mH, tb)
            assert abs(sum(d['channels'].values()) - d['total']) <= 1e-12 * max(d['total'], 1.0)
    # B0 lifetime factor 0.9274 (Gamma_B+/Gamma_B0 = 0.93 in Boiarska's caption)
    a = P.br_B_to_XsH_exclusive(1.0, 100.0, 'B+')['total']
    b = P.br_B_to_XsH_exclusive(1.0, 100.0, 'B0')['total']
    assert abs(b / a - 0.9274) < 1e-3, b / a
    # extends beyond the free-quark threshold
    assert P.br_B_to_XsH_exclusive(4.5, 100.0)['total'] > 0.0
    assert P.br_B_to_XsH_exclusive(4.79, 100.0)['total'] == 0.0
    # inclusive/exclusive ~ 2 at mH -> 0 (Boiarska eq:exclusive-to-inclusive)
    r = P.br_B_to_XsH(0.1, 100.0)['total'] / P.br_B_to_XsH_exclusive(0.1, 100.0)['total']
    assert 1.7 < r < 2.2, r


def test_K_pi():
    assert P.br_K_to_piH(0.36, 100.0, 'K+') == 0.0
    assert P.br_K_to_piH(0.37, 100.0, 'K0L') == 0.0
    assert P.br_K_to_piH(0.0, 100.0, 'K+') > 0.0
    assert P.br_K_to_piH(0.0, 100.0, 'K0S') == 0.0
    # K0L > K+ (longer lifetime): ratio ~4.3
    r = P.br_K_to_piH(0.0, 100.0, 'K0L') / P.br_K_to_piH(0.0, 100.0, 'K+')
    assert 4.0 < r < 4.6, r
    # coherent amplitude: interference term is present and negative for the paper's
    # sign conventions
    a = P.k_to_pi_amplitude(0.0, 100.0, 'K+')
    assert abs(a['interference']) > 0.0
    assert a['interference'] < 0.0
    assert abs(a['M'] ** 2 if isinstance(a['M'], float) else abs(a['M']) ** 2
               - (abs(a['four_quark']) ** 2 + abs(a['fcnc']) ** 2 + a['interference'])) < 1e-30
    # thresholds exactly at mK - mpi
    assert P.br_K_to_piH(0.35410661, 100.0, 'K+') < 1e-12
    # ~1/tanb^2 scaling at large tanb
    r2 = P.br_K_to_piH(0.0, 1e5, 'K+') / P.br_K_to_piH(0.0, 1e4, 'K+')
    assert abs(r2 - 0.01) < 1e-4


def test_semileptonic():
    assert P.br_pi_to_Henu(0.1396, 100.0) == 0.0
    assert P.br_K_to_Henu(0.5, 100.0) == 0.0
    assert P.br_K_to_Henu(0.1, 100.0, 'K0L') == 0.0
    p = P.br_pi_to_Henu(0.0, 100.0)
    k = P.br_K_to_Henu(0.0, 100.0, 'K+')
    assert p > 0 and k > 0
    # 1/tanb^2 scaling
    assert abs(P.br_K_to_Henu(0.0, 1e5) / P.br_K_to_Henu(0.0, 1e4) - 0.01) < 1e-4
    # K/p ratio at c=1 (Boiarska table normalisation sanity): ~3.9e-8 / 2e-9
    r = P.br_K_to_Henu(0.0, 1.0) / P.br_pi_to_Henu(0.0, 1.0)
    assert 15.0 < r < 25.0, r


def test_higgs_sector():
    for tb in (5.0, 10.0, 100.0, 1e5):
        hs = P.higgs_sector(1.0, tb)
        assert abs(hs['g3'] + 2.0 * hs['g_C']) < 1e-18
        assert 0.0 <= hs['BR_hHH'] < 1.0
        assert hs['Gamma_h_other'] > 0.0
        # kappa_g == kappa_f in Type-I (all fermion couplings equal)
        assert abs(hs['kappa_g'] - hs['kappa_f']) < 1e-12
        assert len(hs['sigma_process']) == 9
        assert abs(sum(hs['sigma_process'].values()) - hs['sigma_total']) < 1e-12
    # Gamma(h->HH) formula
    mH, tb = 1.0, 5.0
    hs = P.higgs_sector(mH, tb)
    beta = math.sqrt(1.0 - 4.0 * mH ** 2 / 125.0 ** 2)
    assert abs(hs['Gamma_hHH'] - hs['g3'] ** 2 / (32 * math.pi * 125.0) * beta) < 1e-15
    # BR formula: Gamma_HH/(Gamma_other+Gamma_HH)
    assert abs(hs['BR_hHH'] - hs['Gamma_hHH'] / (hs['Gamma_h_other'] + hs['Gamma_hHH'])) < 1e-15
    # closed above mh/2
    assert P.higgs_sector(62.6, 100.0)['BR_hHH'] == 0.0
    # kappa -> 1 as tanb -> infinity and sigma_total -> 55.4527 pb
    hs = P.higgs_sector(1.0, 1e5)
    for k in ('kappa_f', 'kappa_V', 'kappa_g', 'kappa_gamma'):
        assert abs(hs[k] - 1.0) < 1e-8, (k, hs[k])
    assert abs(hs['sigma_total'] - P.SIGMA_H_TOTAL_SM) < 1e-6
    assert abs(P.SIGMA_H_TOTAL_SM - 55.4527) < 1e-3
    # tHW bracket
    lo, hi = hs['sigma_tHW_range']
    assert lo <= hs['sigma_process']['tHW'] <= hi
    # 4700 (g/v)^2 relation: BR = 4700 (g/v)^2 * (Gam_SM/Gam_other) / (1+Gam_HH/Gam_other)
    mH = 1.0
    hs = P.higgs_sector(mH, 10.0)
    pred = 4700.0 * (hs['g_C'] / P.V_VEV) ** 2
    pred *= (P.GAM_H_SM / hs['Gamma_h_other']) / (1.0 + hs['Gamma_hHH'] / hs['Gamma_h_other'])
    assert abs(hs['BR_hHH'] / pred - 1.0) < 0.02, hs['BR_hHH'] / pred


def test_grid_sanity():
    pts = [(mH, tb) for mH in (0.021, 0.1, 0.5, 1.0, 2.0, 3.5, 4.0, 4.17)
           for tb in (5.0, 10.0, 100.0, 1e3, 1e5)]
    for e in P.evaluate_grid(pts):
        for k, v in e.items():
            if isinstance(v, float):
                assert math.isfinite(v), (k, v)
        for k in ('br_B_to_XsH', 'br_K_to_piH_K+', 'br_K_to_piH_K0L', 'br_K_to_Henu_K+',
                  'br_pi_to_Henu', 'sigma_total', 'BR_hHH', 'Gamma_h_other'):
            assert e[k] >= 0.0, k
        assert e['BR_hHH'] <= 1.0
        assert e['sigma_total'] > 0.0
    # evaluate equals the individual functions
    e = P.evaluate(0.5, 100.0)
    assert abs(e['br_B_to_XsH'] - P.br_B_to_XsH(0.5, 100.0)['total']) < 1e-18
    assert abs(e['sigma_total'] - P.higgs_sector(0.5, 100.0)['sigma_total']) < 1e-15
    # validity flags
    assert P.evaluate(0.1, 5.0)['br_B_to_XsH_validity'] != 'OK'
    assert P.evaluate(0.1, 1e5)['br_B_to_XsH_validity'] == 'OK'


def test_b_channels():
    w = P.b_channel_weights(1.0)
    assert abs(sum(w.values()) - 1.0) < 1e-14
    assert all(v >= 0.0 for v in w.values())
    assert set(w) == set(P.channel_masses())
    assert len(w) == 9
    assert 'pi' not in w          # b->d channel dropped
    # tanb independence is by construction (mH-only function signature)
    assert b_channels.channel_weights(1.0) == w
    # closed-weight behaviour
    assert all(v == 0.0 for v in P.b_channel_weights(4.79).values())


def main():
    tests = [(k, v) for k, v in sorted(globals().items()) if k.startswith('test_')]
    nfail = 0
    for name, fn in tests:
        try:
            fn()
            print('PASS', name)
        except AssertionError as e:
            nfail += 1
            print('FAIL', name, '->', e)
        except Exception as e:  # noqa: BLE001
            nfail += 1
            print('ERROR', name, '->', type(e).__name__, e)
    print('\n%d/%d tests passed' % (len(tests) - nfail, len(tests)))
    return 1 if nfail else 0


if __name__ == '__main__':
    sys.exit(main())
