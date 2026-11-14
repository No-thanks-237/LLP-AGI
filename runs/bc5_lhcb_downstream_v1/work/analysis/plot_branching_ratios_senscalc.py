# -*- coding: utf-8 -*-
"""
BC5 dark scalar branching ratios from the SensCalc tables in this working directory.

DATA SOURCE (stated on the figure):
    phenomenology/darkscalar/prepared/branching_ratios.csv
    converted from SensCalc  phenomenology/Scalar/decay widths/BrRatioScalar.mx
    interpretation: codes/llp-decays.nb
    the tables are built on the Winkler 1809.01876 scalar model.
Nothing on this figure is computed by the present analysis: every curve is SensCalc's own
number, including the pi0 pi0, K_L K_S and pi+pi-2pi0 channels that the model figure derives
by isospin instead of listing.

Output: outputs/figures/bc5_branching_ratios_senscalc.{png,pdf}
        outputs/data/bc5_branching_ratios_senscalc.csv
"""
from __future__ import annotations
import csv, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, '..', '..'))
WORK = os.path.abspath(os.path.join(RUN, '..', '..'))
SRC = os.path.join(WORK, 'phenomenology', 'darkscalar', 'prepared', 'branching_ratios.csv')
FIG = os.path.join(RUN, 'outputs', 'figures')
DATA = os.path.join(RUN, 'outputs', 'data')

# SensCalc column, label, colour -- same legend as the model figure
CHANNELS = [
    ('e_plus_e_minus',        r'$e^+e^-$',               'black',   'e+e-'),
    ('mu_plus_mu_minus',      r'$\mu^+\mu^-$',           'red',     'mu+mu-'),
    ('pi_plus_pi_minus',      r'$\pi^+\pi^-$',           '#2ca02c', 'pi+pi-'),
    ('K_plus_K_minus',        r'$K^+K^-$',               'blue',    'K+K-'),
    ('Jets_ss',               r'$s\bar s$',              '#7B68EE', 'ssbar (jets)'),
    ('2pi_plus_2pi_minus',    r'$\pi^+\pi^-\pi^+\pi^-$', '#FFA500', '2pi+2pi-'),
    ('Jets_GG',               r'$gg$',                   'magenta', 'gg (jets)'),
    ('tau_plus_tau_minus',    r'$\tau^+\tau^-$',         'cyan',    'tau+tau-'),
    ('Jets_cc',               r'$c\bar c$',              '#9ACD32', 'ccbar (jets)'),
]
EXTRA = [('Jets_bb', r'$b\bar b$', '#8B4513', 'bbbar (jets)'),
         ('pi0_pi0', r'$\pi^0\pi^0$', '#808080', 'pi0pi0'),
         ('K_L_K_S', r'$K_LK_S$', '#00BFFF', 'KLKS'),
         ('pi_plus_pi_minus_2pi0', r'$\pi^+\pi^-\pi^0\pi^0$', '#FFD700', 'pi+pi-2pi0')]

XMIN, XMAX = 0.0, 4.5
YMIN, YMAX = 1e-3, 1.1
PERT = 2.0

SOURCE_NOTE = ('data source: SensCalc\n'
               'phenomenology/Scalar/decay widths/BrRatioScalar.mx\n'
               '(Winkler 1809.01876 scalar model),\n'
               'as converted in phenomenology/darkscalar/')


def main():
    df = pd.read_csv(SRC)
    m = df['mass_GeV'].values

    def make(name, chans, source_note=None):
        fig, ax = plt.subplots(figsize=(8.0, 6.0))
        for col, lab, c, _ in chans:
            y = df[col].values.astype(float)
            ax.plot(m, np.where(y > 0, y, np.nan), color=c, lw=1.8, label=lab)
        ax.axvline(PERT, color='0.55', lw=0.9, ls=':')
        ax.text(PERT * 1.02, 1.02, r'$m_S=2$ GeV: scheme switch', fontsize=7, color='0.35')
        ax.set_xscale('linear'); ax.set_yscale('log')
        ax.set_xlim(XMIN, XMAX); ax.set_ylim(YMIN, YMAX)
        ax.set_xlabel(r'$m_S$ [GeV]')
        ax.set_ylabel(r'$\mathrm{BR}(S\to X)$')
        ax.set_title('BC5 dark scalar branching ratios  —  data source: SensCalc')
        ax.grid(alpha=0.3, which='both')
        ax.legend(loc='lower right', fontsize=8.5, ncol=1, framealpha=0.95)
        if source_note:
            ax.text(0.015, 0.02, source_note, transform=ax.transAxes, fontsize=7.0, va='bottom',
                    bbox=dict(boxstyle='round', fc='#FFF8DC', ec='#B8860B', alpha=0.95))
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, name + '.png'), dpi=150)
        fig.savefig(os.path.join(FIG, name + '.pdf'))
        plt.close(fig)

    make('bc5_branching_ratios_senscalc', CHANNELS, SOURCE_NOTE)
    make('bc5_branching_ratios_senscalc_all', CHANNELS + EXTRA, SOURCE_NOTE)

    # dump the plotted numbers
    with open(os.path.join(DATA, 'bc5_branching_ratios_senscalc.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['# data source: SensCalc BrRatioScalar.mx (Winkler 1809.01876 model)'])
        w.writerow(['mass_GeV'] + [c for c, _, _, _ in CHANNELS])
        for i, mm in enumerate(m):
            w.writerow(['%.6g' % mm] + ['%.8e' % df[c].values[i] for c, _, _, _ in CHANNELS])
    print('wrote SensCalc figures and CSV from %s' % SRC)
    print('mass range in file: %.4f - %.4f GeV, %d rows' % (m.min(), m.max(), len(m)))
    for mm in (0.5, 1.0, 2.0, 3.0, 4.5):
        i = int(np.argmin(np.abs(m - mm)))
        print('  m=%6.3f  ' % m[i] + '  '.join('%s=%.3g' % (c, df[c].values[i]) for c, _, _, _ in CHANNELS))
    print()
    print('SensCalc hadronic total across the 2 GeV switch:')
    for mm in (1.99, 1.999, 2.001, 2.01):
        i = int(np.argmin(np.abs(m - mm)))
        dis = sum(df[c].values[i] for c in ('pi_plus_pi_minus', 'pi0_pi0', 'K_plus_K_minus',
                                            'K_L_K_S', '2pi_plus_2pi_minus', 'pi_plus_pi_minus_2pi0'))
        per = sum(df[c].values[i] for c in ('Jets_ss', 'Jets_cc', 'Jets_bb', 'Jets_GG'))
        print('  m=%8.4f  dispersive=%.4f  perturbative=%.4f  ratio=%.3f'
              % (m[i], dis, per, per / dis if dis else float('nan')))


if __name__ == '__main__':
    main()
