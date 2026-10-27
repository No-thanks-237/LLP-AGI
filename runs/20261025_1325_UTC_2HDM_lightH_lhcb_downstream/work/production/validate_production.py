# -*- coding: utf-8 -*-
"""Validation driver for the 2HDM production module (DESIGN.md section 3, items a-f).

Run:  python work/production/validate_production.py
Writes: work/production/validation.json  (machine readable)
        prints every number quoted in logs/production_validation.md

Independent cross-checks used here (deliberately NOT the module's own code):
  * Wolfram Language (mcp__wolfram__WolframLanguageEvaluator) re-evaluation of App. D
    g0,g1,g2 at tanb=100 for k=u,c,t -> WOLFRAM_G (hard-coded below, 12 digits).
  * Boiarska arXiv:1904.10447 master formulas re-implemented here with Boiarska's own
    quark masses / lifetimes (0.0924, 0.00467, Gamma_K = hbar/1.238e-8) and their
    theta^2-normalised xi_ds=3.3e-6, xi_sb=3.6e-4.
  * Lead's independent bracket calculation: -0.3609 (tanb=5), -0.0180 (tanb=100).
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
import twohdm_production as P          # noqa: E402
import b_channels                     # noqa: E402

OUT = {}

# Wolfram Language re-evaluation of App. D (g0,g1,g2,X1,X2), 30-digit precision,
# tanb=100, m_Hpm=600, m_W=80.385; copied from the MCP tool output.
WOLFRAM_G = {
    'u': {'g0': -3.36549386575e-07, 'g1': -0.750000000000, 'g2': 0.0351515305077,
          'X1': 3.51515305077, 'X2': -3.19741896996e-10},
    'c': {'g0': -3.36526819507e-07, 'g1': -0.750000005029, 'g2': 0.0350966706407,
          'X1': 3.50966706910, 'X2': -5.02910117637e-05},
    't': {'g0': -2.73828701384e-07, 'g1': -0.750013977299, 'g2': 0.0105130174230,
          'X1': 1.05131571960, 'X2': -0.139772989813},
}
LEAD_BRACKET = {5.0: -0.3609, 100.0: -0.0180}


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def hdr(t):
    print('\n' + '=' * 100)
    print(t)
    print('=' * 100)


# ======================================================================================
# independent re-implementations (Boiarska master formulas, own constants)
# ======================================================================================
V_BOI = 246.0
HBAR = 6.582119569e-25
GAM_KP_BOI = HBAR / 1.2380e-8
GAM_KL_BOI = HBAR / 5.116e-8
XI_DS_BOI, XI_SB_BOI = 3.3e-6, 3.6e-4
MS_BOI, MD_BOI, MB_BOI = 0.0924, 0.00467, 4.18
B_CHANNELS_BOI = {  # Boiarska theta^2-normalised B -> X_s S channel BRs at mS->0
    'K': 4.3e-1, 'K*892': 4.7e-1, 'K*1410': 2.1e-1, 'K*1680': 1.3e-1,
    'K0*700': 7.6e-1, 'K0*1430': 8.1e-2, 'K1_1270': 9.1e-1, 'K1_1400': 1.6e-2,
    'K2*1430': 3.0e-1,
}


def pstar(mX, mP, mS):
    if mS + mP >= mX:
        return 0.0
    return math.sqrt(max((mX ** 2 - (mS + mP) ** 2) * (mX ** 2 - (mS - mP) ** 2), 0.0)) / (2.0 * mX)


def boiarska_br_K(ms, which='K+'):
    """Boiarska master formula for BR(K -> pi S), theta^2-normalised, own constants.

    uses their chiral matrix element |M_Kpi| = (mK^2-mpi^2)/(m_s-m_d) * f0(0.96).
    """
    if which == 'K+':
        mX, mP, GX = 0.493677, 0.13957039, GAM_KP_BOI
    else:
        mX, mP, GX = 0.497611, 0.1349768, GAM_KL_BOI
    p = pstar(mX, mP, ms)
    if p <= 0:
        return 0.0
    M = abs((mX ** 2 - mP ** 2) / (MS_BOI - MD_BOI)) * 0.96
    return XI_DS_BOI ** 2 * MS_BOI ** 2 * M ** 2 / (32.0 * math.pi * V_BOI ** 2 * GX) * p / mX ** 2


def foresee_f(c):
    """The paper's own FORESEE 2HDM_H bracket (references/tools/FORESEE/Models/2HDM_H/2HDM_H.ipynb).

    F(c) = -3/(4c) - 0.1415/c^3 - (1.0386/c - 0.1415/c^3) sqrt(1-1/c^2) - 0.608446664/c^3
    (c = tanb); it corresponds to our (g1 c - g2 s - g0 k) bracket with fitted numbers.
    """
    return (-3.0 / (4 * c) - 0.1415 / c ** 3
            - (1.0386 / c - 0.1415 / c ** 3) * math.sqrt(max(1.0 - 1.0 / c ** 2, 0.0))
            - 0.608446664 / c ** 3)


def foresee_br_B(mH, tanb, meson='B+'):
    """The paper's FORESEE 2HDM_H B -> X_s H branching ratio (their notebook formula).

    br = 9.68544e-13 * mB/Gamma_B * (1-mH^2/mB^2)^2 * F(tanb)^2
    with the scikit-particle PDG mB and width that FORESEE uses.
    """
    mB = 5.27934 if meson == 'B+' else 5.27965          # particle.Particle 521 / 511
    GB = HBAR / 1.638e-12 if meson == 'B+' else HBAR / 1.519e-12
    if mH >= mB:
        return 0.0
    return 9.68544e-13 * mB / GB * (1.0 - (mH / mB) ** 2) ** 2 * foresee_f(tanb) ** 2


def boiarska_exclusive_B_sum(ms, meson='B+'):
    """Sum of Boiarska exclusive b->s channels at mS: table(0) x p(mS)/p(0), theta^2-norm.

    For the shape we use the table values at mS=0 and scale by the two-body phase
    space p*/p*(0) (form factors at q^2=mS^2 are a smaller correction that we keep
    exactly in b_channels; here the point is an independent normalisation check).
    """
    mX = 5.27934 if meson == 'B+' else 5.27965
    GX = HBAR / 1.638e-12 if meson == 'B+' else HBAR / 1.519e-12
    tot = 0.0
    for name, br0 in B_CHANNELS_BOI.items():
        mP = b_channels.B_CHANNEL_MASS[name]
        p = pstar(mX, mP, ms)
        if p <= 0:
            continue
        p0 = pstar(mX, mP, 0.0)
        tot += br0 * p / p0
    return tot * (1.0 if meson == 'B+' else 0.927)


# ======================================================================================
# a.  xi_H^{bs} hand calculation at (mH, tanb) = (1, 100)
# ======================================================================================
def validation_a():
    hdr('a.  xi_H^{bs} at mH = 1 GeV, tanb = 100: hand calculation, per-term')
    mH, tb = 1.0, 100.0
    d = P.xi_H_details(mH, tb, 's', 'b')
    print('prefactor -(4 sqrt2 G_F)/(16 pi^2) = %.12e' % d['prefactor'])
    print('lambda_HH+H- = -(2 m_Hpm^2 - mH^2) c/v = %.10e GeV' % d['lambda_HH'])
    print('2v/mW^2 * lambda                  = %.10e' % d['k'])
    print('c = %.10f   s = %.12f' % (d['c'], d['s']))
    print('\npaper eq. xi_H_2HDM lower (H) row per loop quark:')
    print('  bracket_k = g1_k*c - g2_k*s - g0_k*(2v/mW^2)*lambda_HH+H-')
    print('%-3s %12s %12s %12s %14s %14s %14s' %
          ('k', 'g1_k*c', '-g2_k*s', '-g0_k*k', 'bracket', 'V*_ks V_kb', 'm_k^2*ckm*br'))
    for q in ('u', 'c', 't'):
        r = d['rows'][q]
        print('%-3s %12.6e %12.6e %12.6e %14.10e  (%9.3e%+9.3ei) %14.6e%+14.6ei' %
              (q, r['g1_c'], r['minus_g2_s'], r['minus_g0_k'], r['bracket'],
               r['ckm_product'].real, r['ckm_product'].imag,
               r['m2_ckm_bracket'].real, r['m2_ckm_bracket'].imag))
    print('sum_k = %.10e %+.10ei' % (d['sum'].real, d['sum'].imag))
    print('xi_H^{bs} = prefactor * sum = %.10e %+.10ei' % (d['xi_H'].real, d['xi_H'].imag))
    print('|xi_H^{bs}| = %.10e' % abs(d['xi_H']))
    print('\ncomparison with the LEAD independent bracket calculation:')
    b = {}
    for tb2 in (5.0, 100.0):
        t = P.xi_H_terms(1.0, tb2)
        bt = t['per_quark']['t']['bracket']
        b[tb2] = bt
        print('  tanb=%-6g bracket(t) = %.8f   lead = %.4f   delta = %+.2e' %
              (tb2, bt, LEAD_BRACKET[tb2], bt - LEAD_BRACKET[tb2]))
    # independent re-implementation in numpy from the paper text
    V = P.V_CKM
    m2 = [P.M_QUARK['u'] ** 2, P.M_QUARK['c'] ** 2, P.M_QUARK['t'] ** 2]
    xH = P.M_HPM ** 2 / P.M_W ** 2
    cotb = 1.0 / tb
    c, s = 1.0 / tb, math.sqrt(1.0 - 1.0 / tb ** 2)
    lam = -(2.0 * P.M_HPM ** 2 - mH ** 2) * c / P.V_VEV
    tot = 0.0 + 0.0j
    for iq, q in enumerate(('u', 'c', 't')):
        x = P.M_QUARK[q] ** 2 / P.M_W ** 2
        dd = xH - x
        L = math.log(xH / x)
        g0 = -cotb ** 2 * (3 * xH ** 2 - 4 * xH * x + x ** 2 - 2 * x * (2 * xH - x) * L) / (16 * dd ** 3)
        g1 = -0.75 + cotb ** 2 * x * (5 * xH ** 2 - 8 * xH * x + 3 * x ** 2
                                      - 2 * xH * (2 * xH - x) * L) / (4 * dd ** 3)
        A = xH * (3 * xH - 2 * x) / dd ** 2
        X1 = -0.25 * ((xH / dd - 6 / (x - 1) ** 2 + 3) - A * math.log(xH)
                      + (A + 3 * (x + 1) / (x - 1) ** 3) * math.log(x))
        X2 = x * (5 * xH - 3 * x) / (4 * dd ** 2) - xH * x * (2 * xH - x) / (2 * dd ** 3) * L
        g2 = cotb * X1 + cotb ** 3 * X2
        br = g1 * c - g2 * s - g0 * (2 * V_BOI / P.M_W ** 2) * lam
        tot += np.conj(V[iq, 1]) * m2[iq] * br * V[iq, 2]
    xi_ind = -4.0 * math.sqrt(2.0) * P.GF / (16 * math.pi ** 2) * tot
    print('\nindependent numpy re-implementation of the whole sum: '
          'xi_H^{bs} = %.10e %+.10ei  (delta = %.2e)'
          % (xi_ind.real, xi_ind.imag, abs(xi_ind - d['xi_H'])))
    # Wolfram check of App. D
    print('\nWolfram-Language check of App. D g0,g1,g2 (tanb=100):')
    G = P.g_aux_functions(tb)
    wmax = 0.0
    for q in ('u', 'c', 't'):
        for fn in ('g0', 'g1', 'g2', 'X1', 'X2'):
            ref = WOLFRAM_G[q][fn]
            got = G[q][fn]
            rel = abs(got - ref) / max(abs(ref), 1e-30)
            wmax = max(wmax, rel)
            print('  %-2s %-3s python %.12e  wolfram %.12e  rel.diff %.2e' % (q, fn, got, ref, rel))
    OUT['a'] = {'mH': mH, 'tanb': tb,
                'xi_H_bs': [d['xi_H'].real, d['xi_H'].imag], 'xi_H_bs_abs': abs(d['xi_H']),
                'lambda_HH': d['lambda_HH'],
                'rows': {q: {k: (r[k] if not isinstance(r[k], complex) else [r[k].real, r[k].imag])
                             for k in ('g0', 'g1', 'g2', 'g1_c', 'minus_g2_s', 'minus_g0_k',
                                       'bracket', 'ckm_product', 'm2_ckm_bracket')}
                         for q, r in d['rows'].items()},
                'bracket_t': {str(k): v for k, v in b.items()},
                'lead_bracket': {str(k): v for k, v in LEAD_BRACKET.items()},
                'independent_xi': [xi_ind.real, xi_ind.imag],
                'wolfram_g_max_rel_diff': wmax}


# ======================================================================================
# b.  Br(B -> X_s H) at mH -> 0
# ======================================================================================
def validation_b():
    hdr('b.  Br(B->XsH) at mH -> 0, tanb = 100: inclusive (nominal) vs Boiarska exclusive')
    tb = 100.0
    rows = []
    for mH in (0.021, 0.1, 1.0, 2.0, 3.0, 4.0, 4.17):
        nom = P.br_B_to_XsH(mH, tb)
        ex = P.br_B_to_XsH_exclusive(mH, tb)
        boi_ex = boiarska_exclusive_B_sum(mH) * (abs(P.xi_H(mH, tb, 's', 'b')) / XI_SB_BOI) ** 2
        rows.append({'mH': mH, 'nominal': nom['total'], 'exclusive_variant': ex['total'],
                     'boiarska_exclusive_rescaled': boi_ex,
                     'ratio_nominal_over_excl': (nom['total'] / ex['total']) if ex['total'] else None,
                     'ratio_nominal_over_boiarska': (nom['total'] / boi_ex) if boi_ex else None})
        print('mH=%-6.3f nominal=%.6e  exclusive(module)=%.6e  Boiarska-excl-rescaled=%.6e  '
              'nom/excl=%.3f  nom/boi=%.3f'
              % (mH, nom['total'], ex['total'], boi_ex,
                 rows[-1]['ratio_nominal_over_excl'] or 0, rows[-1]['ratio_nominal_over_boiarska'] or 0))
    # explicit hand calculation at (mH=0.1, tanb=100)
    mH = 0.1
    xib = abs(P.xi_H(mH, tb, 's', 'b'))
    f = P._f_phase_space((P.M_C_QUARK / P.M_B_QUARK) ** 2)
    t1 = 12.0 * math.pi ** 2
    t2 = P.V_VEV ** 2 / P.M_B_QUARK ** 2
    t3 = (1.0 - mH ** 2 / P.M_B_QUARK ** 2) ** 2
    t4 = 1.0 / f
    t5 = (xib / abs(P.V_CKM[1, 2])) ** 2
    print('\nHAND CHECK Br(B->XsH) at mH=0.1, tanb=100 (all terms printed):')
    print('  12 pi^2            = %.8f' % t1)
    print('  v^2/m_b^2          = %.8f  (v=%.1f, m_b=%.2f)' % (t2, P.V_VEV, P.M_B_QUARK))
    print('  (1-mH^2/m_b^2)^2   = %.8f' % t3)
    print('  1/f(mc^2/mb^2)     = %.8f  (f=%.8f, mc=%.2f)' % (t4, f, P.M_C_QUARK))
    print('  |xi_H^bs/Vcb|^2    = %.8e  (|xi|=%.6e, |Vcb|=%.6f)' % (t5, xib, abs(P.V_CKM[1, 2])))
    print('  product            = %.8e' % (t1 * t2 * t3 * t4 * t5))
    print('  x Br(B->Xc e nu)   = %.8e  (0.104)' % (t1 * t2 * t3 * t4 * t5 * P.BR_B_XC_ENU))
    print('  module total       = %.8e' % P.br_B_to_XsH(mH, tb)['total'])
    # ---- Lead request: the three rate modes at mH = 0.1, 1, 3, 4, 4.5 ----------
    print('\nRATE MODES (DEFAULT_B_MODE=%s), tanb=100:' % P.DEFAULT_B_MODE)
    modes = {}
    for mH in (0.1, 1.0, 3.0, 4.0, 4.5):
        d = {md: P.br_B_to_XsH(mH, tb, mode=md) for md in P.B_MODES}
        d = {k: v['total'] for k, v in d.items()}
        r1 = d['foresee'] / d['paper_mb'] if d['paper_mb'] else float('inf')
        r2 = d['foresee'] / d['exclusive'] if d['exclusive'] else float('inf')
        print('  mH=%-4.1f foresee=%.5e paper_mb=%.5e exclusive=%.5e | foresee/paper_mb=%-9s '
              'foresee/excl=%s' % (mH, d['foresee'], d['paper_mb'], d['exclusive'],
                                   ('%.3f' % r1) if np.isfinite(r1) else 'inf',
                                   ('%.3f' % r2) if np.isfinite(r2) else 'inf'))
        modes[str(mH)] = {'foresee': d['foresee'], 'paper_mb': d['paper_mb'],
                          'exclusive': d['exclusive'],
                          'ratio_foresee_paper_mb': r1, 'ratio_foresee_exclusive': r2}
    # --- Lead request: our bracket vs FORESEE's f(tanb), tanb grid ---
    print('\n(i) our top-quark bracket vs FORESEE f(c):')
    brac = {}
    for tb2 in (5.0, 10.0, 20.0, 50.0, 100.0, 300.0, 1000.0, 1e4, 1e5):
        b = P.xi_H_terms(1.0, tb2)['per_quark']['t']['bracket']
        fF = foresee_f(tb2)
        brac[str(tb2)] = {'bracket': b, 'foresee': fF, 'ratio': b / fF}
        print('  tanb=%-8g bracket=%+.8f  FORESEE f=%+.8f  ratio=%.5f' % (tb2, b, fF, b / fF))
    # --- Lead request: our BR vs FORESEE BR at (mH,tanb) nodes, mH < 4 ---
    print('\n(ii) our br_B_to_XsH(mode="foresee") vs FORESEE formula, B+ (mH<4):')
    fs = {}
    non = []
    for mH in (0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0):
        for tb2 in (5.0, 10.0, 30.0, 100.0, 300.0, 1000.0):
            ours = P.br_B_to_XsH(mH, tb2, 'B+', mode='foresee')['total']
            ref = foresee_br_B(mH, tb2, 'B+')
            non.append(ours / ref)
            fs['mH=%g_tanb=%g' % (mH, tb2)] = {'ours': ours, 'foresee': ref, 'ratio': ours / ref}
    print('  %d nodes: ratio ours/FORESEE min=%.4f max=%.4f mean=%.4f' %
          (len(non), min(non), max(non), sum(non) / len(non)))
    for k in ('mH=1_tanb=100', 'mH=3_tanb=100', 'mH=4_tanb=10'):
        print('   %-14s ours=%.5e FORESEE=%.5e ratio=%.4f' % (k, fs[k]['ours'], fs[k]['foresee'], fs[k]['ratio']))
    # B0 (lifetime factor differs in the FORESEE normalisation)
    o = P.br_B_to_XsH(1.0, 100.0, 'B0', mode='foresee')['total']
    rr = foresee_br_B(1.0, 100.0, 'B0')
    print('  B0 at (1,100): ours=%.5e FORESEE=%.5e ratio=%.4f  (both now carry '
          'BR(B0)=Gamma_B+/Gamma_B0*BR(B+), Gamma_B+/Gamma_B0=%.5f)' %
          (o, rr, o / rr, P.B0_LIFETIME_SCALE))
    print('  B0/B+ lifetime scaling (all modes): ' + ', '.join(
        '%s=%.5f' % (md, P.br_B_to_XsH(1.0, 100.0, 'B0', mode=md)['total'] /
                     P.br_B_to_XsH(1.0, 100.0, 'B+', mode=md)['total']) for md in P.B_MODES))
    # high mass reach and the mesonic endpoint (mB - mK = 4.785663 for B+)
    print('  mesonic endpoint: mB - mK = %.6f (B+), %.6f (B0); foresee mode returns 0 '
          'at/above it' % (P.M_BP - P.M_KP, P.M_B0 - P.M_KP))
    print('  foresee mode near the endpoint: ' + ', '.join(
        'mH=%.3f:total=%.4e flag=%s' % (
            mH, P.br_B_to_XsH(mH, 100.0, mode='foresee')['total'],
            P.br_B_to_XsH(mH, 100.0, mode='foresee')['beyond_mesonic_endpoint'])
        for mH in (4.5, 4.7, 4.78, 4.7857, 4.79, 5.0)))
    print('  other modes at mH=4.8: paper_mb=%.3e exclusive=%.3e (natural thresholds); '
          'beyond_mesonic_endpoint=%s'
          % (P.br_B_to_XsH(4.8, 100.0, mode='paper_mb')['total'],
             P.br_B_to_XsH(4.8, 100.0, mode='exclusive')['total'],
             P.br_B_to_XsH(4.8, 100.0)['beyond_mesonic_endpoint']))
    print('\nBoiarska inclusive/exclusive relation (their eq:exclusive-to-inclusive): 0.5 '
          '-> inclusive/exclusive = 2.  We find %.3f at mH=0.1 (1.957 expected): agreement %.1f%%'
          % (rows[1]['ratio_nominal_over_excl'], 100 * (rows[1]['ratio_nominal_over_excl'] / 2.0 - 1)))
    print('Boiarska total exclusive B->XsS at mS=0: 3.3 theta^2 (their eq:total-exclusive);')
    print('our Boiarska-normalised exclusive sum at mS=0.021: %.3f theta^2' % boiarska_exclusive_B_sum(0.021))
    print('\nHigh-mass behaviour: foresee closes at the B-meson mass 5.279, paper_mb at '
          'm_b=4.18, exclusive at each channel threshold (~4.79)')
    for mH in (4.17, 4.2, 4.5, 4.78, 4.79, 5.0):
        d = {md: P.br_B_to_XsH(mH, tb, mode=md)['total'] for md in P.B_MODES}
        print('  mH=%.2f  foresee=%.4e  paper_mb=%.4e  exclusive=%.4e  channels open: %s'
              % (mH, d['foresee'], d['paper_mb'], d['exclusive'],
                 [k for k, v in P.br_B_to_XsH_exclusive(mH, tb)['channels'].items() if v > 0]))
    print('\nvalidity boundaries (Br(B->XsH) = 0.1 and 1):')
    vb = {}
    for mH in (0.021, 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 4.17):
        v = P.validity_boundaries(mH)
        vb[str(mH)] = v
        print('  mH=%-6.3f  tanb(BR=0.1)=%s  tanb(BR=1)=%s'
              % (mH, '%.2f' % v['tanb_BR_gt_0.1'] if v['tanb_BR_gt_0.1'] else 'none',
                 '%.2f' % v['tanb_BR_gt_1'] if v['tanb_BR_gt_1'] else 'none'))
    # monotonicity of BR in tanb (bisection assumption)
    mono = True
    for mH in (0.021, 0.1, 1.0, 3.0):
        brs = [P.br_B_to_XsH(mH, t)['total'] for t in np.logspace(math.log10(5), 5, 60)]
        if not all(brs[i] >= brs[i + 1] * (1 - 1e-12) for i in range(len(brs) - 1)):
            mono = False
    print('\nBr(B->XsH) monotonically decreasing in tanb on the 60-point check: %s' % mono)
    print('\nFASER-paper B picture: the paper uses exactly this inclusive ratio with '
          'Br(B->Xc e nu)=0.104 (their Eq. after eq:Lfceff_phi) and takes B production only')
    print('for m_H above the kaon window; at tanb < %.1f (mH=1) Br(B->XsH)>0.1 -> the '
          'effective/inclusive description is saturated (flag).' % P.tanb_valid(1.0, 0.1))
    OUT['b'] = {'tanb': tb, 'rows': rows, 'validity': vb, 'monotonic': mono,
                'boiarska_ratio_expected': 2.0,
                'boiarska_exclusive_theta2_at_0p021': boiarska_exclusive_B_sum(0.021),
                'modes': modes, 'bracket_vs_foresee': brac,
                'foresee_nodes': fs,
                'foresee_ratio_min': min(non), 'foresee_ratio_max': max(non),
                'foresee_ratio_mean': sum(non) / len(non),
                'default_mode': P.DEFAULT_B_MODE}


# ======================================================================================
# c.  K -> pi H at mH -> 0, tanb = 100
# ======================================================================================
def validation_c():
    hdr('c.  K->pi H at mH -> 0, tanb = 100: FCNC vs four-quark vs interference')
    tb = 100.0
    xi_ds = P.xi_H(mH=0.0, tanb=tb, i='d', j='s')
    print('xi_H^{ds} = %.10e %+.10ei  |xi| = %.6e   (Boiarska xi_ds = %.2e)'
          % (xi_ds.real, xi_ds.imag, abs(xi_ds), XI_DS_BOI))
    out = {}
    for mH in (0.0, 0.05, 0.1396, 0.2, 0.3, 0.3541):
        a = P.k_to_pi_amplitude(mH, tb, 'K+')
        M2 = abs(a['M']) ** 2
        br = P.br_K_to_piH(mH, tb, 'K+')
        fq2 = abs(a['four_quark']) ** 2
        fc2 = abs(a['fcnc']) ** 2
        intf = a['interference']
        print('\nmH=%.4f  p*=%.6f  |M|=%.6e' % (mH, a['p_star'], abs(a['M'])))
        print('   |FCNC|^2          = %.6e   (BR_FCNC_only = %.6e)'
              % (fc2, (1 / P.GAM_KP) * (2 * a['p_star'] / a['mK']) * fc2 / (16 * math.pi * a['mK'])))
        print('   |4-quark|^2       = %.6e   (BR_4q_only   = %.6e)'
              % (fq2, (1 / P.GAM_KP) * (2 * a['p_star'] / a['mK']) * fq2 / (16 * math.pi * a['mK'])))
        print('   interference 2Re  = %.6e   (%.2f %% of |M|^2)' % (intf, 100 * intf / M2))
        print('   |M|^2 total       = %.6e   BR(K+->pi+H) = %.6e' % (M2, br))
        # Boiarska master formula rescaled by (xi_H^ds/xi_B^ds)^2
        boi = boiarska_br_K(mH) * (abs(xi_ds) / XI_DS_BOI) ** 2
        print('   Boiarska K->piS master formula rescaled: %.6e   FCNC-only/boiarska = %.4f'
              % (boi, ((1 / P.GAM_KP) * (2 * a['p_star'] / a['mK']) * fc2
                       / (16 * math.pi * a['mK'])) / boi if boi else 0))
        out[str(mH)] = {'br': br, 'fcnc2': fc2, 'fourq2': fq2, 'interference': intf,
                        'M2': M2, 'boiarska_rescaled': boi,
                        'br_fcnc_only': (1 / P.GAM_KP) * (2 * a['p_star'] / a['mK']) * fc2 / (16 * math.pi * a['mK'])}
    # K0L and K0S
    print('\nK0L->pi0 H at mH=0, tanb=100: %.6e  (per K0L; response applies 1/2 PDG-311)'
          % P.br_K_to_piH(0.0, tb, 'K0L'))
    print('K0S->pi0 H: %.1f  (CP suppressed, as in BC5)' % P.br_K_to_piH(0.0, tb, 'K0S'))
    print('K0L / K+ ratio: %.4f   (Boiarska table: 7e-3/1.7e-3 = %.2f)'
          % (P.br_K_to_piH(0.0, tb, 'K0L') / P.br_K_to_piH(0.0, tb, 'K+'), 7.0 / 1.7))
    OUT['c'] = {'tanb': tb, 'xi_H_ds': [xi_ds.real, xi_ds.imag],
                'xi_ds_abs': abs(xi_ds), 'points': out,
                'K0L_over_Kplus': P.br_K_to_piH(0.0, tb, 'K0L') / P.br_K_to_piH(0.0, tb, 'K+')}


# ======================================================================================
# d.  semileptonic pi/K -> H e nu
# ======================================================================================
def validation_d():
    hdr('d.  Br(pi->H e nu), Br(K->H e nu) at mH -> 0')
    tb = 100.0
    out = {}
    for name, val, mX, br0 in (('pi', P.br_pi_to_Henu(0.0, tb), P.M_PI_P, P.BR_PI_MUNU),
                               ('K', P.br_K_to_Henu(0.0, tb, 'K+'), P.M_KP, P.BR_K_MUNU)):
        c = 1.0 / tb
        # independent: Boiarska table normalisation = same formula with theta=1 (xi_W=1)
        x = (P.M_MU / mX) ** 2
        unit = (math.sqrt(2) * P.GF * mX ** 4 / (96 * math.pi ** 2 * P.M_MU ** 2 * (1 - x) ** 2)
                * br0 * (7.0 / 9.0) ** 2)
        print('%-3s: Br = %.6e at tanb=100 (c^2=%.2e);  at tanb=1 (c=1): %.6e' % (name, val, c * c, unit))
        out[name] = {'br_tanb100': val, 'theta1_normalisation': unit}
    print('\nBoiarska table tab:BR-3-body: BR(K->S e nu)/f(x)/theta^2 = 4.1e-8; we get %.3e '
          '(normalisation, theta=1)' % out['K']['theta1_normalisation'])
    print('Boiarska table tab:BR-3-body (their text) and our pi normalisation: %.3e' % out['pi']['theta1_normalisation'])
    print('K0L -> H e nu: %.1f (not a physical two-body leptonic channel)' % P.br_K_to_Henu(0.0, tb, 'K0L'))
    th = {}
    for mH in (0.021, 0.1, 0.13, 0.1395, 0.1396):
        th[str(mH)] = {'pi': P.br_pi_to_Henu(mH, tb), 'K': P.br_K_to_Henu(mH, tb, 'K+')}
        print('mH=%-7.4f Br(pi->He nu)=%.4e  Br(K->He nu)=%.4e' % (mH, th[str(mH)]['pi'], th[str(mH)]['K']))
    OUT['d'] = {'tanb': tb, 'items': out, 'thresholds': th}


# ======================================================================================
# e.  g_C, h -> HH
# ======================================================================================
def validation_e():
    hdr('e.  g_C, g3, Gamma(h->HH), BR(h->HH)')
    approx = lambda mH, tb: -P.M_H_SM ** 2 * (1.0 / tb) ** 2 / (4.0 * P.V_VEV)   # noqa: E731
    out = {}
    print('paper approximation g_hHH ~ -mh^2 c^2/(4v) (their text after eq:BM_H).')
    for mH in (0.021, 0.1, 1.0, 2.0, 4.0):
        for tb in (5.0, 10.0, 100.0):
            g = P.g_hHH(mH, tb)
            hs = P.higgs_sector(mH, tb)
            appr = approx(mH, tb)
            br4700 = 4700.0 * (g / P.V_VEV) ** 2
            print('mH=%-5.2f tanb=%-6g g_C=%+.6e  |g_C|/|appr|=%.4f  g3=%+.6e  '
                  'Gamma_hHH=%.4e  Gamma_h_other=%.6e  BR_hHH=%.6e  BR(4700 g^2/v^2)=%.6e  ratio=%.4f'
                  % (mH, tb, g, abs(g / appr) if appr else 0, hs['g3'], hs['Gamma_hHH'],
                     hs['Gamma_h_other'], hs['BR_hHH'], br4700,
                     hs['BR_hHH'] / br4700 if br4700 else 0))
            out['mH=%g_tanb=%g' % (mH, tb)] = {'g_C': g, 'g3': hs['g3'], 'appr': appr,
                                               'Gamma_hHH': hs['Gamma_hHH'],
                                               'Gamma_h_other': hs['Gamma_h_other'],
                                               'BR_hHH': hs['BR_hHH'], 'BR_4700': br4700,
                                               'kappa_g': hs['kappa_g'], 'kappa_f': hs['kappa_f'],
                                               'kappa_gamma': hs['kappa_gamma']}
    # mH -> 0 convergence to the paper approximation
    print('\nmH -> 0 convergence of g_C to -mh^2 c^2/(4v) at tanb=100:')
    conv = []
    for mH in (1.0, 0.1, 0.021):
        g = P.g_hHH(mH, 100.0)
        r = g / approx(mH, 100.0)
        conv.append({'mH': mH, 'ratio_to_paper_approx': r})
        print('  mH=%-6.4f g_C=%+.6e  paper approx=%+.6e  ratio=%.5f' % (mH, g, approx(mH, 100.0), r))
    print('  (the deviation at mH=1 is the exact -(mH^2 c^2/2) term of the bracket; it '
          'vanishes as mH->0)')
    print('\n4700*(g/v)^2 vs our Gamma(h->HH)/(Gamma_h^other): '
          'the paper coefficient 4700 = v^2/(8 pi mh Gamma_h^SM) = %.1f'
          % (P.V_VEV ** 2 / (8 * math.pi * P.M_H_SM * P.GAM_H_SM)))
    print('paper text also claims Br(h->HH) ~ 75 c^4; our exact formula gives ~%.2f c^4 at mH->0 '
          '(flagged as a paper inconsistency, see report)' % (P.higgs_sector(0.021, 10.0)['BR_hHH'] / (0.1 ** 4)))
    OUT['e'] = {'table': out, 'convergence': conv,
                'paper_coeff_4700': P.V_VEV ** 2 / (8 * math.pi * P.M_H_SM * P.GAM_H_SM)}


# ======================================================================================
# f.  kappa -> 1, sigma_total -> 55.45 pb; per-process at tanb=5,1000
# ======================================================================================
def validation_f():
    hdr('f.  kappa -> 1 and sigma_total -> 55.4527 pb as tanb -> infinity')
    out = {}
    for tb in (10.0, 100.0, 1000.0, 1e4, 1e5):
        hs = P.higgs_sector(1.0, tb)
        print('tanb=%-8g kappa_f-1=%+.3e kappa_V-1=%+.3e kappa_g-1=%+.3e kappa_gamma-1=%+.3e '
              'sigma_total=%.6f' % (tb, hs['kappa_f'] - 1, hs['kappa_V'] - 1, hs['kappa_g'] - 1,
                                    hs['kappa_gamma'] - 1, hs['sigma_total']))
        out['tanb=%g' % tb] = {'kappa_f': hs['kappa_f'], 'kappa_V': hs['kappa_V'],
                                'kappa_g': hs['kappa_g'], 'kappa_gamma': hs['kappa_gamma'],
                                'sigma_total': hs['sigma_total']}
    print('\nSIGMA_H sum = %.6f pb ; at tanb=1e5 sigma_total = %.6f pb (delta = %.2e)'
          % (P.SIGMA_H_TOTAL_SM, P.higgs_sector(1.0, 1e5)['sigma_total'],
             abs(P.higgs_sector(1.0, 1e5)['sigma_total'] - P.SIGMA_H_TOTAL_SM)))
    print('\nper-process at tanb=5 and tanb=1000 vs sigma_SM*kappa^2 (Lead request):')
    per = {}
    for tb in (5.0, 1000.0):
        hs = P.higgs_sector(1.0, tb)
        print('  tanb=%g   kappa_f=%.6f kappa_V=%.6f kappa_g=%.6f' % (tb, hs['kappa_f'], hs['kappa_V'], hs['kappa_g']))
        print('  %-5s %10s %10s %10s %12s' % ('proc', 'sigma_SM', 'sigma_2HDM', 'SM*kappa^2', 'ratio'))
        for p in ('ggF', 'VBF', 'WH', 'ZH', 'bbH', 'ttH', 'tHj', 'tHW', 'tHb'):
            exp = P.SIGMA_H[p] * (hs['kappa_g'] ** 2 if p == 'ggF' else
                                  hs['kappa_V'] ** 2 if p in ('VBF', 'WH', 'ZH') else
                                  hs['kappa_f'] ** 2)
            print('  %-5s %10.4f %10.4f %10.4f %12.6f' % (p, P.SIGMA_H[p], hs['sigma_process'][p], exp,
                                                          hs['sigma_process'][p] / exp if exp else 0))
            per.setdefault(p, {})[str(tb)] = {'sigma_SM': P.SIGMA_H[p], 'sigma': hs['sigma_process'][p], 'expect': exp}
        print('  total = %.4f pb' % hs['sigma_total'])
    # tHW bracket
    print('\n tHW interference bracket: sigma = sigma_SM*[kt^2+kV^2+2 kt kV cos(phi)]/4, baseline phi=0')
    thw = {}
    for tb in (5.0, 100.0, 1000.0):
        hs = P.higgs_sector(1.0, tb)
        lo, hi = hs['sigma_tHW_range']
        frac = hs['sigma_process']['tHW'] / hs['sigma_total']
        print('  tanb=%-6g kt=%.6f kV=%.6f -> sigma_tHW base=%.6e pb  range=[%.6e, %.6e]  '
              'weight in total=%.4f%%' % (tb, hs['kappa_f'], hs['kappa_V'],
                                          hs['sigma_process']['tHW'], lo, hi, 100 * frac))
        thw[str(tb)] = {'base': hs['sigma_process']['tHW'], 'lo': lo, 'hi': hi, 'frac': frac}
    OUT['f'] = {'scan': out, 'per_process': per, 'tHW': thw,
                'SIGMA_H_TOTAL_SM': P.SIGMA_H_TOTAL_SM}


# ======================================================================================
# extra: CKM choice, couplings identities, thresholds, provenance
# ======================================================================================
def validation_extra():
    hdr('extra.  CKM choice, coupling identities, sanity/provenance')
    V = P.V_CKM
    Valt = P.ckm_matrix_from_magnitudes()
    print('unitary CKM |V|:')
    for i, row in enumerate(np.abs(V)):
        print('   ' + ' '.join('%.6f' % x for x in row))
    print('max |V^dag V - 1| = %.2e' % np.max(np.abs(V.conj().T @ V - np.eye(3))))
    print('|Vtd Vts| unitary = %.6e ; nominal-list variant = %.6e ; ratio = %.4f'
          % (abs(V[2, 0] * V[2, 1]), abs(Valt[2, 0] * Valt[2, 1]),
             abs(Valt[2, 0] * Valt[2, 1]) / abs(V[2, 0] * V[2, 1])))
    br_u = P.br_K_to_piH(0.0, 100.0, 'K+')
    # recompute with the nominal-magnitude matrix by temporarily swapping V_CKM
    V_save = P.V_CKM
    P.V_CKM = Valt
    br_alt = P.br_K_to_piH(0.0, 100.0, 'K+')
    P.V_CKM = V_save
    print('BR(K+->pi+H) with unitary CKM = %.6e ; with nominal list = %.6e ; ratio = %.4f'
          % (br_u, br_alt, br_alt / br_u))
    print('\nGIM check: sum_k V*_ki V_kj for i!=j = max %.2e (must vanish)'
          % max(abs(np.sum(np.conj(V[:, i]) * V[:, j])) for i in range(3) for j in range(3) if i != j))
    print('\ncoupling identities at tanb = 100, 1e3, 1e5:')
    for tb in (100.0, 1e3, 1e5):
        cp = P.couplings(tb)
        print('  tanb=%-8g xi_H^f = c-s/tb = %.12e ; 1/(tb^3(1+s)) = %.12e ; rel %.1e'
              % (tb, cp['xi_H_f'], cp['xi_H_f_stable'],
                 abs(cp['xi_H_f'] - cp['xi_H_f_stable']) / cp['xi_H_f_stable']))
    print('\nlambda_HH+H- at lambda v^2 = 0: -(2 m_Hpm^2 - mH^2) c/v')
    for mH in (1.0, 10.0):
        print('  mH=%-5.1f tanb=100: %.8f  (c=0.01 -> -%.8f)'
              % (mH, P.lambda_HHpmHpm(mH, 100.0), (2 * 600.0 ** 2 - mH ** 2) * 0.01 / 246.0))
    print('\nBR(h->HH)=0 above mh/2: %.1f' % P.higgs_sector(62.6, 100.0)['BR_hHH'])
    print('channel weights sum: %.12f  (mH=1)' % sum(P.b_channel_weights(1.0).values()))
    print('channel weights independent of tanb: yes by construction (mH-only function)')
    print('\nprovenance sha256:')
    sha = {}
    for rel in ('inputs/references/FASER_light.tex', 'inputs/bc5_reuse/bc5_production.py',
                'work/analysis/grids/mass_grid.txt', 'work/production/twohdm_production.py',
                'work/production/b_channels.py', 'DESIGN.md'):
        p = os.path.join(RUN, rel)
        if os.path.exists(p):
            sha[rel] = sha256(p)
            print('  %-45s %s' % (rel, sha[rel]))
    OUT['extra'] = {'ckm_unitary': np.abs(V).tolist(),
                    'unitarity_dev': float(np.max(np.abs(V.conj().T @ V - np.eye(3)))),
                    'VtdVts_unitary': float(abs(V[2, 0] * V[2, 1])),
                    'VtdVts_nominal': float(abs(Valt[2, 0] * Valt[2, 1])),
                    'br_Kpi_unitary': br_u, 'br_Kpi_nominal': br_alt,
                    'br_Kpi_ratio': br_alt / br_u,
                    'sha256': sha}


def main():
    OUT['meta'] = {'module': 'work/production/twohdm_production.py',
                   'model_tag': P.MODEL_TAG, 'version': P.__version__,
                   'design': 'DESIGN.md section 3 (items a-f)',
                   'python': sys.version,
                   'numpy': np.__version__,
                   'mW': P.M_W, 'mt': P.M_QUARK['t'], 'mHpm': P.M_HPM, 'v': P.V_VEV,
                   'mb': P.M_B_QUARK, 'ms_K': P.M_S_QUARK, 'md_K': P.M_D_QUARK}
    validation_a()
    validation_b()
    validation_c()
    validation_d()
    validation_e()
    validation_f()
    validation_extra()
    with open(os.path.join(HERE, 'validation.json'), 'w', encoding='utf-8') as f:
        json.dump(OUT, f, indent=1, default=float)
    print('\nwrote', os.path.join(HERE, 'validation.json'))


if __name__ == '__main__':
    main()
