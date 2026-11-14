#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Standalone driver for the BC5 dark-scalar decay widths.

    python widths_cli.py --mass 1.0
    python widths_cli.py --mass 1.0 --theta2 1e-9
    python widths_cli.py --mass 2.5 --json
    python widths_cli.py --scan 0.021 62.49 52
    python widths_cli.py --check
    python widths_cli.py --mass 1.0 --provenance

See README.md for the model, the two mass regimes, the provenance of every channel, and the
validation results.  Requires Python 3 and numpy only.
"""
from __future__ import annotations
import argparse, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import scalar_widths as sw

# channel -> (label, how it is obtained, is it computed here?)
PROVENANCE = {
    'e':     ('e+e-',        'closed form, G_F m_S m_l^2 beta^3 / (4 sqrt2 pi)', True),
    'mu':    ('mu+mu-',      'closed form, G_F m_S m_l^2 beta^3 / (4 sqrt2 pi)', True),
    'tau':   ('tau+tau-',    'closed form, G_F m_S m_l^2 beta^3 / (4 sqrt2 pi)', True),
    'gamma': ('gamma gamma', 'SM h -> gamma gamma loop functions, m_h -> m_S', True),
    'pipi':  ('pi pi',       'digitised Winkler Fig. 2 (dispersive)', False),
    'KK':    ('K K',         'digitised Winkler Fig. 2 (dispersive)', False),
    'multi': ('multi-meson', 'C m_S^3 beta_2pi, C from Winkler', None),
    'ss':    ('s sbar',      'N_c m_S mbar_q^2(m_S) beta^3 / (8 pi v^2), LO spectator', True),
    'cc':    ('c cbar',      'N_c m_S mbar_q^2(m_S) beta^3 / (8 pi v^2), LO spectator', True),
    'bb':    ('b bbar',      'N_c m_S mbar_q^2(m_S) beta^3 / (8 pi v^2), LO spectator', True),
    'gg':    ('g g',         'alpha_s^2 m^3 |sum F_1/2|^2 / (128 pi^3 v^2)', True),
}
ORDER = ['e', 'mu', 'tau', 'gamma', 'pipi', 'KK', 'multi', 'ss', 'cc', 'bb', 'gg']


def select_hadronic(mode):
    if mode == 'table':
        path = os.path.join(HERE, 'hadronic_table.csv')
        if not os.path.exists(path):
            raise SystemExit('hadronic_table.csv is missing next to widths_cli.py')
        sw.load_hadronic_table(path)
        return 'digitised Winkler Fig. 2  (work/widths/hadronic_table.csv)'
    if mode == 'chpt':
        sw.HADRONIC_TABLE = None
        return 'lowest-order ChPT form factors (valid below ~0.5 GeV only)'
    raise SystemExit('unknown --hadronic mode %r' % mode)


def report(ms, theta2=1.0, gg_k=1.0):
    w = sw.widths(ms, gg_k)
    tot = sum(w.values())
    return {
        'm_S_GeV': ms,
        'theta2': theta2,
        'regime': 'dispersive hadrons' if ms < sw.PERT_THRESHOLD else 'perturbative quarks/gluons',
        'Gamma_GeV_at_theta2_1': tot,
        'ctau_m': sw.HBAR_C / (tot * theta2) if tot > 0 else float('inf'),
        'widths_GeV': {k: w.get(k, 0.0) for k in ORDER},
        'branching_ratios': {k: (w.get(k, 0.0) / tot if tot > 0 else 0.0) for k in ORDER},
        'quark_nlo': sw.QUARK_NLO,
    }


def print_report(r, provenance=False):
    print('m_S = %.4f GeV,  theta^2 = %.3g' % (r['m_S_GeV'], r['theta2']))
    print('  regime   : %s' % r['regime'])
    print('  Gamma_S  : %.6e GeV' % r['Gamma_GeV_at_theta2_1'])
    print('  c*tau    : %.6e m' % r['ctau_m'])
    print()
    print('  %-14s %14s %10s   %s' % ('channel', 'Gamma [GeV]', 'BR',
                                      'source' if provenance else ''))
    for k in ORDER:
        g = r['widths_GeV'][k]
        if g <= 0:
            continue
        lab, src, computed = PROVENANCE[k]
        if provenance:
            mark = 'computed' if computed else ('external' if computed is False else 'formula')
            src = '%-8s %s' % (mark, src)
        print('  %-14s %14.6e %10.4f   %s' % (lab, g, r['branching_ratios'][k], src))
    print()
    print('  widths are at theta^2 = 1; Gamma scales as theta^2')


def main():
    ap = argparse.ArgumentParser(description='BC5 dark scalar decay widths (see README.md)')
    ap.add_argument('--mass', type=float, help='scalar mass in GeV')
    ap.add_argument('--scan', nargs=3, metavar=('MIN', 'MAX', 'N'), help='log-spaced mass scan')
    ap.add_argument('--theta2', type=float, default=1.0, help='mixing angle squared')
    ap.add_argument('--gg-k', type=float, default=1.0, help='NLO K factor for S -> g g')
    ap.add_argument('--hadronic', default='table', choices=('table', 'chpt'),
                    help='source of the pi pi / K K widths')
    ap.add_argument('--json', action='store_true', help='machine-readable output')
    ap.add_argument('--provenance', action='store_true', help='show where each channel comes from')
    ap.add_argument('--check', action='store_true', help='validate against the reference curves')
    args = ap.parse_args()

    if args.check:
        import validate_widths
        validate_widths.main()
        return

    if not args.json:
        print('# pi pi / K K from: %s' % select_hadronic(args.hadronic))
        print()
        if args.theta2 == 1.0:
            print('# theta^2 not set: widths below are at theta^2 = 1 (Gamma scales as theta^2)')
            print()
    else:
        select_hadronic(args.hadronic)

    if args.scan:
        lo, hi, n = float(args.scan[0]), float(args.scan[1]), int(args.scan[2])
        masses = [lo * (hi / lo) ** (i / max(n - 1, 1)) for i in range(n)]
    elif args.mass is not None:
        masses = [args.mass]
    else:
        ap.error('one of --mass, --scan or --check is required')

    rows = [report(m, args.theta2, args.gg_k) for m in masses]
    if args.json:
        print(json.dumps(rows if len(rows) > 1 else rows[0], indent=2))
    else:
        for r in rows:
            print_report(r, provenance=args.provenance)


if __name__ == '__main__':
    main()
