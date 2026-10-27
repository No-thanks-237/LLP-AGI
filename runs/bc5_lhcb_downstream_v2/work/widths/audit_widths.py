# -*- coding: utf-8 -*-
"""
Audit of the scalar width model against every reference data set available locally.

    python audit_widths.py

Compares, channel by channel,
  A. the absolute widths against the digitised Winkler Fig. 2 (work/widths/reference_figures)
  B. the branching ratios against the SensCalc tables (phenomenology/darkscalar/prepared)
  C. the branching ratios against the digitised branching-ratio figure of Boiarska et al.
  D. internal consistency: BR sum, 2 GeV matching, thresholds, table sanity

Writes outputs/data/widths_audit.json and prints a table.  Exit code 0 always: this is a
diagnostic, not a pass/fail test.
"""
from __future__ import annotations
import csv, json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, '..', '..'))
WORK = os.path.abspath(os.path.join(RUN, '..', '..'))
sys.path.insert(0, HERE)
import scalar_widths as sw

TAB = os.path.join(HERE, 'hadronic_table.csv')
OUT = os.path.join(RUN, 'outputs', 'data', 'widths_audit.json')

# my channel -> the corresponding curve id in the Winkler partial-width figure
WINKLER = {
    'mu': 'muons', 'tau': 'taus', 'pipi': 'pions', 'KK': 'kaons',
    'multi': 'multimeson_4pi_etaeta_rhorho',
    'ss': 'strange_quarks_ss', 'cc': 'charm_quarks_cc', 'gg': 'gluons_gg',
}
# my channel -> the columns of the SensCalc table that must be SUMMED to match it
# (the SensCalc table lists exclusive channels; mine are isospin totals)
SENSCALC = {
    'e': ['e_plus_e_minus'], 'mu': ['mu_plus_mu_minus'], 'tau': ['tau_plus_tau_minus'],
    'pipi': ['pi_plus_pi_minus', 'pi0_pi0'], 'KK': ['K_plus_K_minus', 'K_L_K_S'],
    'multi': ['2pi_plus_2pi_minus', 'pi_plus_pi_minus_2pi0'],
    'ss': ['Jets_ss'], 'cc': ['Jets_cc'], 'bb': ['Jets_bb'], 'gg': ['Jets_GG'],
}
BOIARSKA = {'e': 'BR_ee', 'mu': 'BR_mumu', 'tau': 'BR_tautau', 'pipi': 'BR_pions',
            'KK': 'BR_kaons', 'multi': 'BR_4pi', 'ss': 'BR_ss', 'cc': 'BR_cc', 'gg': 'BR_gg'}


def load_curves(path):
    d = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            d.setdefault(r['curve_id'], []).append((float(r['mass_GeV']), float(r['value'])))
    return {k: np.array(sorted(v)) for k, v in d.items()}


def widths_and_br(ms):
    w = sw.widths(ms)
    tot = sum(w.values())
    return w, tot, {k: (v / tot if tot > 0 else 0.0) for k, v in w.items()}


def compare(x, y, masses):
    ok = (np.isfinite(x)) & (np.isfinite(y)) & (x > 0) & (y > 0)
    if ok.sum() == 0:
        return None
    r = x[ok] / y[ok]
    i = int(np.argmax(np.abs(r - 1.0)))
    return {'n': int(ok.sum()), 'median': float(np.median(r)),
            'min': float(r.min()), 'max': float(r.max()),
            'worst_mass': float(masses[ok][i]), 'worst_ratio': float(r[i])}


def main():
    sw.load_hadronic_table(TAB)
    tab = np.genfromtxt(TAB, delimiter=',', names=True)
    rep = {'sources': {}, 'A_absolute_widths_vs_winkler': {}, 'B_br_vs_senscalc': {},
           'C_br_vs_boiarska': {}, 'D_consistency': {}}

    wink = load_curves(os.path.join(HERE, 'reference_figures', 'decayrates_curves.csv'))
    rep['sources']['winkler_fig2'] = os.path.join('work', 'widths', 'reference_figures', 'decayrates_curves.csv')
    sens_path = os.path.join(WORK, 'phenomenology', 'darkscalar', 'prepared', 'branching_ratios.csv')
    sens = np.genfromtxt(sens_path, delimiter=',', names=True, dtype=None, encoding=None)
    rep['sources']['senscalc'] = sens_path
    boi_path = os.path.join(HERE, 'reference_figures', 'scalardecaybr_curves.csv')
    boi = load_curves(boi_path) if os.path.exists(boi_path) else {}
    rep['sources']['boiarska_br_figure'] = boi_path

    # ---- A -------------------------------------------------------------------------------
    for key, cid in WINKLER.items():
        if cid not in wink:
            continue
        a = wink[cid]
        mine = np.array([sw.widths(m).get(key, 0.0) for m in a[:, 0]])
        rep['A_absolute_widths_vs_winkler'][key] = compare(mine, a[:, 1], a[:, 0])
        if rep['A_absolute_widths_vs_winkler'][key]:
            rep['A_absolute_widths_vs_winkler'][key]['mass_range'] = [float(a[0, 0]), float(a[-1, 0])]
    if 'hadrons' in wink:
        a = wink['hadrons']
        mine = np.array([sum(v for k, v in sw.widths(m).items() if k not in ('e', 'mu', 'tau', 'gamma'))
                         for m in a[:, 0]])
        rep['A_absolute_widths_vs_winkler']['hadronic_total'] = compare(mine, a[:, 1], a[:, 0])

    # ---- A2: absolute widths against SensCalc (its BRs times its own total width) --------
    ms = np.exp(np.linspace(math.log(max(sens['mass_GeV'].min(), 1e-3)),
                            math.log(sens['mass_GeV'].max()), 120))
    gpath = os.path.join(WORK, 'phenomenology', 'darkscalar', 'GammaS_reference.csv')
    rep['A2_absolute_widths_vs_senscalc'] = {}
    if os.path.exists(gpath) and os.path.exists(sens_path):
        g = np.genfromtxt(gpath, delimiter=',', names=True, dtype=None, encoding=None)
        rep['sources']['senscalc_width'] = gpath
        for key, cols in SENSCALC.items():
            if key == 'multi' or key == 'bb':
                continue
            mm, mine, ref = [], [], []
            for m in ms:
                w = sw.widths(m)
                gtot = float(np.interp(m, g['mass_GeV'], g['Gamma_reference_GeV']))
                br_ref = sum(float(np.interp(m, sens['mass_GeV'], sens[c])) for c in cols)
                mm.append(m); mine.append(w.get(key, 0.0)); ref.append(br_ref * gtot)
            s = compare(np.array(mine), np.array(ref), np.array(mm))
            if s:
                s['senscalc_total_width_over_mine_at_2GeV'] = float(
                    np.interp(2.0, g['mass_GeV'], g['Gamma_reference_GeV']) / sw.gamma_total(2.0))
            rep['A2_absolute_widths_vs_senscalc'][key] = s

    # ---- B -------------------------------------------------------------------------------
    ms = np.exp(np.linspace(math.log(max(sens['mass_GeV'].min(), 1e-3)),
                            math.log(sens['mass_GeV'].max()), 120))
    for key, cols in SENSCALC.items():
        mine, ref = [], []
        for m in ms:
            w, tot, br = widths_and_br(m)
            if key == 'multi':
                continue                      # my multi is Winkler's toy, not the two 4pi channels
            mine.append(br.get(key, 0.0))
            ref.append(sum(float(np.interp(m, sens['mass_GeV'], sens[c])) for c in cols))
        rep['B_br_vs_senscalc'][key] = compare(np.array(mine), np.array(ref), ms)

    # ---- C -------------------------------------------------------------------------------
    for key, cid in BOIARSKA.items():
        if cid not in boi:
            continue
        a = boi[cid]
        mine = np.array([widths_and_br(m)[2].get(key, 0.0) for m in a[:, 0]])
        rep['C_br_vs_boiarska'][key] = compare(mine, a[:, 1], a[:, 0])

    # ---- D -------------------------------------------------------------------------------
    d = rep['D_consistency']
    mss = np.exp(np.linspace(math.log(0.021), math.log(62.49), 200))
    sumdev = []
    for m in mss:
        w = sw.widths(m)
        t = sum(w.values())
        sumdev.append(abs(sum(w.values()) / t - 1.0) if t > 0 else 0.0)
    d['br_sum_max_deviation'] = float(max(sumdev))

    w1, w2 = sw.widths(2.0 - 1e-9), sw.widths(2.0 + 1e-9)
    h1 = sum(v for k, v in w1.items() if k not in ('e', 'mu', 'tau', 'gamma'))
    h2 = sum(v for k, v in w2.items() if k not in ('e', 'mu', 'tau', 'gamma'))
    d['hadronic_width_below_2GeV'] = float(h1)
    d['hadronic_width_above_2GeV'] = float(h2)
    d['hadronic_matching_ratio'] = float(h2 / h1)

    d['table'] = {
        'rows': int(len(tab)),
        'mass_range_GeV': [float(tab['mass_GeV'][0]), float(tab['mass_GeV'][-1])],
        'spacing_min_GeV': float(np.min(np.diff(tab['mass_GeV']))),
        'spacing_max_GeV': float(np.max(np.diff(tab['mass_GeV']))),
        'non_finite_rows': int((~np.isfinite(tab['gamma_pipi_GeV'])
                                | ~np.isfinite(tab['gamma_KK_GeV'])).sum()),
        'negative_rows': int((tab['gamma_pipi_GeV'] < 0).sum() + (tab['gamma_KK_GeV'] < 0).sum()),
    }
    d['thresholds'] = {}
    for key, thr in (('e', 0.0010219979), ('mu', 2 * 0.1056583755), ('tau', 2 * 1.77686),
                     ('pipi', 2 * 0.13957039), ('KK', 2 * 0.493677), ('multi', 4 * 0.13957039),
                     ('cc', 2 * 1.8697), ('bb', 2 * 5.2797)):
        d['thresholds'][key] = {
            'nominal_GeV': thr,
            'below': float(sw.widths(thr * 0.995).get(key, 0.0)),
            'above': float(sw.widths(thr * 1.02).get(key, 0.0)),
            'ok': bool(sw.widths(thr * 0.995).get(key, 0.0) == 0.0
                       and sw.widths(thr * 1.02).get(key, 0.0) > 0.0),
        }

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, 'w') as f:
        json.dump(rep, f, indent=2)

    # ---- print ---------------------------------------------------------------------------
    def line(tag, s):
        if s is None:
            print('  %-9s   no comparable points' % tag)
        else:
            print('  %-9s n=%4d  median %6.3f  range [%5.3f, %6.3f]  worst m=%.3f'
                  % (tag, s['n'], s['median'], s['min'], s['max'], s['worst_mass']))

    print('A. absolute widths / Winkler Fig. 2')
    for k, s in rep['A_absolute_widths_vs_winkler'].items():
        line(k, s)
    print()
    print('A2. absolute widths / SensCalc   (its BRs x its own total width)')
    for k, s in rep['A2_absolute_widths_vs_senscalc'].items():
        line(k, s)
    print()
    print('B. branching ratios / SensCalc   (mine are isospin totals, SensCalc exclusive summed)')
    for k, s in rep['B_br_vs_senscalc'].items():
        line(k, s)
    print()
    print('C. branching ratios / Boiarska digitised figure')
    for k, s in rep['C_br_vs_boiarska'].items():
        line(k, s)
    print()
    print('D. consistency')
    print('  BR sum max deviation from 1      : %.2e' % d['br_sum_max_deviation'])
    print('  hadronic width below/above 2 GeV : %.5e / %.5e  ratio %.4f'
          % (h1, h2, d['hadronic_matching_ratio']))
    print('  table: %d rows, m in [%.4f, %.4f], spacing %.1e - %.1e, %d NaN, %d negative'
          % (d['table']['rows'], d['table']['mass_range_GeV'][0], d['table']['mass_range_GeV'][1],
             d['table']['spacing_min_GeV'], d['table']['spacing_max_GeV'],
             d['table']['non_finite_rows'], d['table']['negative_rows']))
    bad = [k for k, v in d['thresholds'].items() if not v['ok']]
    print('  thresholds: %s' % ('all OK' if not bad else 'check ' + ', '.join(bad)
                                + ' (ss is only used above 2 GeV by construction)'))
    print()
    print('written: %s' % OUT)


if __name__ == '__main__':
    main()
