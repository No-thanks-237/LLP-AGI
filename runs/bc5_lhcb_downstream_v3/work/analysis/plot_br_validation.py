# -*- coding: utf-8 -*-
"""Validation figure: our Winkler branching ratios vs the published reference table.

Two things are kept separate in this project:

  * `work/widths/scalar_widths.py`  -- the model IMPLEMENTATION, Winkler (arXiv:1809.01876) only;
  * `work/widths/reference_br_model.py` -- a VALIDATION model built from the published reference
    branching-ratio table, used only to cross-check the pipeline.

This script is the validation half: it plots the two sets of branching ratios side by side, per
channel, without feeding one into the other.

Caveats drawn on the figure:
  * the reference table lists only ONE charge state for pi pi and K K (pi+pi-, K+K-), while ours
    are isospin sums (pi+pi- + pi0pi0, K+K- + K0K0bar), so those two panels carry a constant
    factor and are marked as such;
  * the reference table has no gamma-gamma and no b bbar column, so those two channels are drawn
    from our model only.

    python work/analysis/plot_br_validation.py
"""
from __future__ import annotations
import csv, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(RUN, 'work', 'widths'))
import scalar_widths as sw

FIG = os.path.join(RUN, 'outputs', 'figures')
REF_CSV = os.path.join(RUN, 'outputs', 'data', 'bc5_branching_ratios_senscalc.csv')

# our channel -> reference column (None = not in the reference table)
PANELS = [
    ('e',         'e_plus_e_minus',       r'$e^+e^-$'),
    ('mu',        'mu_plus_mu_minus',     r'$\mu^+\mu^-$'),
    ('tau',       'tau_plus_tau_minus',   r'$\tau^+\tau^-$'),
    ('pipi',      'pi_plus_pi_minus',     r'$\pi\pi$  ([W]: isoscalar, i.e. $\pi^+\pi^-+\pi^0\pi^0$)'),
    ('KK',        'K_plus_K_minus',       r'$K\bar K$  ([W]: isoscalar, i.e. $K^+K^-+K^0\bar K^0$)'),
    ('multi',     '2pi_plus_2pi_minus',   r'$4\pi,\ \eta\eta,\ \rho\rho,\dots$'),
    ('gg',        'Jets_GG',              r'$gg$'),
    ('ss',        'Jets_ss',              r'$s\bar s$'),
    ('cc',        'Jets_cc',              r'$c\bar c$'),
    ('bb',        None,                   r'$b\bar b$  (not in the reference table)'),
    ('gamma',     None,                   r'$\gamma\gamma$  (not in the reference table)'),
]


def main():
    sw.load_hadronic_table(os.path.join(RUN, 'work', 'widths', 'hadronic_table.csv'))
    with open(REF_CSV) as fh:
        ref = [r for r in csv.DictReader(l for l in fh if not l.startswith('#'))]
    ref_m = np.array([float(r['mass_GeV']) for r in ref])
    ref_br = {c: np.array([float(r[c]) for r in ref]) for c in ref[0] if c != 'mass_GeV'}

    grid = np.exp(np.linspace(np.log(0.021), np.log(62.49), 400))
    ours = {}
    for i, m in enumerate(grid):
        w = sw.widths(float(m))
        tot = sum(w.values())
        for ch, _col, _lab in PANELS:
            ours.setdefault(ch, np.zeros(len(grid)))[i] = (w.get(ch, 0.0) / tot) if tot > 0 else 0.0

    fig, axes = plt.subplots(3, 4, figsize=(16.5, 9.6), sharex=True)
    axes = axes.ravel()
    for k, (ch, col, label) in enumerate(PANELS):
        ax = axes[k]
        ax.plot(grid, np.maximum(ours[ch], 1e-12), '-', color='#1f77b4', lw=1.8,
                label='this work (Winkler)')
        if col is not None:
            ax.plot(ref_m, np.maximum(ref_br[col], 1e-12), '--', color='0.35', lw=1.4,
                    label='reference table')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(1e-8, 2.0)
        ax.set_title(label, fontsize=9)
        ax.grid(alpha=0.3, which='both')
        if k % 4 == 0:
            ax.set_ylabel('branching ratio', fontsize=9)
        if k >= 7:
            ax.set_xlabel(r'$m_S$ [GeV]', fontsize=9)
        ax.tick_params(labelsize=8)
    ax = axes[len(PANELS)]
    ax.axis('off')
    handles, labels = axes[0].get_legend_handles_labels()
    ax.legend(handles, labels, loc='center', fontsize=10, frameon=False)
    ax.text(0.0, 0.32,
            'Pi pi / K K: this work uses the digitised Winkler curve,\n'
            'which is the isoscalar rate (the CG coefficient of the\n'
            'isoscalar projection is absorbed into his form factors),\n'
            'i.e. pi+pi- + pi0pi0 and K+K- + K0K0bar.  The reference\n'
            'columns are the single charge states pi+pi- and K+K-,\n'
            'so those two panels differ by 1.5 and 2.0 by construction.\n'
            'The reference table has neither gamma-gamma nor b bbar.',
            fontsize=7.5, va='top', transform=ax.transAxes)
    fig.suptitle('Branching ratios: our Winkler implementation vs the published reference table',
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = os.path.join(FIG, 'br_validation')
    fig.savefig(out + '.png', dpi=140)
    fig.savefig(out + '.pdf')
    plt.close(fig)
    print('wrote', out + '.png / .pdf')


if __name__ == '__main__':
    main()
