# -*- coding: utf-8 -*-
"""Validation-only width model built directly from the published reference branching ratios.

This module is NOT part of the BC5 model implementation.  `scalar_widths.py` implements Winkler
(arXiv:1809.01876) only, channel by channel.  The present module exists so that the same
sensitivity pipeline can be run with the decay description of the *other* reference
(Boiarska et al. arXiv:1904.10447 / SensCalc), which separates two questions:

    "is our implementation of Winkler right?"      -> scalar_widths.py vs the printed formulas
    "does the pipeline reproduce the reference?"   -> this module fed into run_scan.py

Construction
------------
* the table is `outputs/data/bc5_branching_ratios_senscalc.csv` (SensCalc BrRatioScalar),
  columns e, mu, pi+pi-, K+K-, Jets_ss, 2pi+2pi-, Jets_GG, tau, Jets_cc, 8115 mass rows
  from 0.021 to 62.49 GeV;
* Gamma_tot(m_S) = Gamma_l(m_S) / BR_l(m_S) with l = mu (or e below 2 m_mu).  The leptonic width
  is model independent (same exact G_F formula in every model), so it is a thermometer for the
  total width.  This is the only anchor used - no other number from the reference is imported;
* Gamma_a = BR_a * Gamma_tot for every listed channel;
* whatever is left over (their table lists only ONE charge state for pi pi / K K, and it omits
  gamma-gamma, the pi+pi-2pi0 channel and b bbar) is assigned to `multi` below 2 GeV (where the
  omitted channels are the extra hadronic ones) and to `bb` above 2 GeV (where the omitted
  piece is the b bbar channel: at 40 GeV the listed columns only sum to 0.144).

Everything is theta^2 = 1 and in GeV, exactly like `scalar_widths.widths`.
"""
from __future__ import annotations
import csv, math, os

import numpy as np

import scalar_widths as sw

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, '..', '..'))
CSV = os.path.join(RUN, 'outputs', 'data', 'bc5_branching_ratios_senscalc.csv')

# reference column -> our channel name
_COLS = {'e_plus_e_minus': 'e', 'mu_plus_mu_minus': 'mu', 'tau_plus_tau_minus': 'tau',
         'pi_plus_pi_minus': 'pipi', 'K_plus_K_minus': 'KK', '2pi_plus_2pi_minus': 'multi',
         'Jets_ss': 'ss', 'Jets_cc': 'cc', 'Jets_GG': 'gg'}
_TABLE = None


def _load():
    global _TABLE
    if _TABLE is not None:
        return _TABLE
    with open(CSV) as fh:
        rows = [r for r in csv.DictReader(l for l in fh if not l.startswith('#'))]
    mass = np.array([float(r['mass_GeV']) for r in rows])
    lm = np.log(mass)
    data = {}
    for col, name in _COLS.items():
        data[name] = np.array([float(r[col]) for r in rows])
    _TABLE = (lm, data)
    return _TABLE


def branching(ms):
    """Interpolated reference branching ratios (log-log in the mass)."""
    lm, data = _load()
    lm0 = math.log(ms)
    out = {}
    for name, vals in data.items():
        if vals.max() <= 0.0:
            out[name] = 0.0
        else:
            lo = max(int(np.searchsorted(lm, lm0)) - 1, 0)
            hi = min(lo + 1, len(lm) - 1)
            w = 0.0 if hi == lo else (lm0 - lm[lo]) / (lm[hi] - lm[lo])
            w = min(max(w, 0.0), 1.0)
            out[name] = float((1.0 - w) * vals[lo] + w * vals[hi])
    return out


def widths(ms: float):
    """Drop-in replacement for scalar_widths.widths() (validation model)."""
    br = branching(ms)
    out = {ch: 0.0 for ch in ('e', 'mu', 'tau', 'pipi', 'KK', 'multi', 'ss', 'cc', 'bb', 'gg',
                              'gamma')}
    # thermometer: the leptonic width is the same in every model
    if br['mu'] > 1e-12:
        g_tot = sw.gamma_ll(ms, sw.M_MU) / br['mu']
    elif br['e'] > 1e-12:
        g_tot = sw.gamma_ll(ms, sw.M_E) / br['e']
    else:
        return out
    for ch, b in br.items():
        out[ch] = b * g_tot
    # F13: the table lists only ONE charge state for pi pi and K K (pi+pi-, K+K-), while the
    # program's 'pipi' and 'KK' are the ISOSCALAR rates (the pi0pi0 and K0Kbar0 halves are
    # included, exactly as in the digitised Winkler curve).  An isoscalar pair is 2/3 pi+pi- +
    # 1/3 pi0pi0 and 1/2 K+K- + 1/2 K0Kbar0, so the columns must be multiplied by 1.5 and 2.0;
    # otherwise the visible validation rate is wrong by up to 16 % (and the neutral halves end up
    # double-counted inside 'multi').  At 0.5 GeV the table's own leftover is exactly
    # BR(pi+pi-)/2, which is the pi0pi0 half and confirms the identification.
    out['pipi'] *= 1.5
    out['KK'] *= 2.0
    # the remainder: extra hadronic channels below 2 GeV, b bbar above
    rest = g_tot - sum(out.values())
    if rest > 0.0:
        out['multi' if ms <= sw.PERT_THRESHOLD else 'bb'] += rest
    return out


def gamma_total(ms: float) -> float:
    return sum(widths(ms).values())


if __name__ == '__main__':
    sw.load_hadronic_table(os.path.join(HERE, 'hadronic_table.csv'))
    for m in (0.05, 0.3, 0.5, 1.0, 1.5, 1.99, 2.0, 3.0, 5.0, 10.0, 40.0, 62.0):
        w = widths(m); tot = sum(w.values())
        print('m=%7.3f  Gamma=%.5e  ' % (m, tot)
              + ' '.join('%s=%.3f' % (k, v / tot) for k, v in sorted(w.items(), key=lambda kv: -kv[1]) if v > 0))
    print()
    print('%8s %13s %13s %8s' % ('m_S', 'ref-BR model', 'Winkler', 'ratio'))
    for m in (0.3, 0.5, 1.0, 1.5, 1.99, 2.1, 3.0, 5.0, 10.0, 40.0):
        print('%8.3f %13.4e %13.4e %8.3f' % (m, gamma_total(m), sum(sw.widths(m).values()),
                                             gamma_total(m) / sum(sw.widths(m).values())))
