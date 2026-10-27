# -*- coding: utf-8 -*-
"""
Overlay our BC5 downstream exclusion contour on the reference curve of
references/model/BC5.pdf (the BC5 benchmark figure, Gorkavenko et al.).

    python plot_comparison_reference.py                  # s = 2.3 (90% CL)
    python plot_comparison_reference.py --s 3 --suffix _s3

Inputs
  outputs/data/exclusion_boundary<suffix>.csv               our contour
  work/widths/reference_figures/bc5_reference_curves.csv   digitised reference,
                                                           columns curve_id,mass_GeV,theta2
Output
  outputs/figures/bc5_vs_reference<suffix>.{png,pdf}

The reference curves were digitised from the PDF's vector path geometry, so the solid 25 fb-1
and dashed 300 fb-1 curves are separated exactly by their dash patterns.  Each of them is the
outline of a closed excluded band, which inside the axes appears as two branches joined by a
near-vertical segment at the right-hand edge; the branches are split here and drawn separately.
"""
from __future__ import annotations
import argparse, csv, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, '..', '..'))
DATA = os.path.join(RUN, 'outputs', 'data')
FIG = os.path.join(RUN, 'outputs', 'figures')
REF = os.path.join(RUN, 'work', 'widths', 'reference_figures', 'bc5_reference_curves.csv')
REF_REPORT = os.path.join(RUN, 'work', 'widths', 'reference_figures', 'bc5_digitisation_report.json')

MINE = {25.0: ('#1f77b4', r'ours, 25 fb$^{-1}$'), 300.0: ('#d62728', r'ours, 300 fb$^{-1}$')}
REF_COLOUR = {25.0: '0.0', 300.0: '0.45'}


def load_reference(path, report_path):
    """curve_id -> list of (mass, theta2) branches, using the row ranges recorded by the
    digitisation.  Each blue curve in BC5.pdf is the outline of a closed excluded band, so the
    CSV stores its two edges as contiguous blocks; slicing by those ranges avoids drawing a
    fake connector between them."""
    import json
    with open(report_path) as f:
        rep = json.load(f)
    rows = []
    with open(path) as f:
        for row in csv.DictReader(f):
            rows.append((row['curve_id'], float(row['mass_GeV']), float(row['theta2'])))
    # row_start / row_end are 1-based CSV line numbers with the header on line 1, so the data
    # rows of a segment are rows[row_start-2 : row_end-1].  Getting this off by one splices the
    # first point of the next branch onto the end of the previous one and draws a fake connector.
    out = {}
    for cid, info in rep['curves'].items():
        out[cid] = []
        for seg in info['segments']:
            a = rows[seg['row_start'] - 2: seg['row_end'] - 1]
            out[cid].append((np.array([r[1] for r in a]), np.array([r[2] for r in a])))
            print('  %-18s %-10s %4d points  m %.3f - %.3f'
                  % (cid, seg.get('branch', '?'), len(a), a[0][1], a[-1][1]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--s', type=float, default=2.3)
    ap.add_argument('--suffix', type=str, default='')
    ap.add_argument('--s3', action='store_true', help='also write the s = 3 comparison')
    args = ap.parse_args()

    jobs = [(args.s, args.suffix)]
    if args.s3:
        jobs.append((3.0, args.suffix + '_s3'))
    for s_cut, sfx in jobs:
        make(s_cut, sfx)


def make(s_cut, sfx):
    bnd = np.genfromtxt(os.path.join(DATA, 'exclusion_boundary%s.csv' % sfx),
                        delimiter=',', names=True)
    fig, ax = plt.subplots(figsize=(7.8, 5.8))

    for L in (25.0, 300.0):
        c = bnd[bnd['lumi_fb'] == L]
        colour, label = MINE[L]
        ax.plot(c['mass_GeV'], c['theta2_long_lived_edge'], color=colour, lw=2.0, ls='--',
                label='%s, long-lived edge' % label)
        ax.plot(c['mass_GeV'], c['theta2_short_lived_edge'], color=colour, lw=2.0, ls='-',
                label='%s, short-lived edge' % label)

    if os.path.exists(REF) and os.path.exists(REF_REPORT):
        for cid, branches in load_reference(REF, REF_REPORT).items():
            L = 300.0 if '300' in cid else 25.0
            ls = '--' if L == 300.0 else '-'
            for j, (mm, tt) in enumerate(branches):
                ax.plot(mm, tt, color=REF_COLOUR[L], lw=1.5, ls=ls, zorder=1,
                        label=('reference, %g fb$^{-1}$ (digitised from BC5.pdf)' % L) if j == 0 else None)
    else:
        print('WARNING: %s not found; plotting our contour only' % REF)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$m_S$ [GeV]')
    ax.set_ylabel(r'$\theta^2$')
    ax.set_title(r'BC5 dark scalar, LHCb downstream: this analysis vs the reference  '
                 r'($s=%.4g$)' % s_cut)
    ax.set_xlim(0.25, 70)
    ax.set_ylim(1e-20, 3e-5)
    ax.grid(alpha=0.3, which='both')
    ax.legend(loc='lower left', fontsize=7.2, ncol=1, framealpha=0.95)
    fig.tight_layout()
    out = os.path.join(FIG, 'bc5_vs_reference%s' % sfx)
    fig.savefig(out + '.png', dpi=150)
    fig.savefig(out + '.pdf')
    plt.close(fig)
    print('wrote %s.png / .pdf' % out)


if __name__ == '__main__':
    main()
