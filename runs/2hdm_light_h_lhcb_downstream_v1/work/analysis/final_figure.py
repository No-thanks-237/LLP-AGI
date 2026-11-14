# -*- coding: utf-8 -*-
"""Final figures for the 2HDM light-H LHCb-downstream scan.

Main exclusion figure: ONE closed outline per luminosity connecting the two limits of the
excluded band (lower = short-lived crossing / scan floor, upper = long-lived edge), with
L = 300 fb^-1 as a SOLID line and L = 25 fb^-1 as a DASHED line.  No vertical threshold
lines, no shading, no data points, no MC band on the main figure.  The EFT-validity flags
and the short-crossing MC caveat live in outputs/data/exclusion_boundary.csv,
logs/run_report.md and exclusion_crossings_mc.png.
"""
from __future__ import annotations
import json, math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import final_analysis as FA

FIG = os.path.join(FA.RUN, 'outputs', 'figures')
DATA = os.path.join(FA.RUN, 'outputs', 'data')


def _save_fig(fig, path_png, path_pdf, dpi=180):
    """Save PNG/PDF with retries; fall back to a .new file if an external viewer locks it."""
    import time, os
    for attempt in range(4):
        try:
            fig.savefig(path_png, dpi=dpi)
            break
        except (PermissionError, OSError) as exc:
            print('  save retry %d/4 %s: %s' % (attempt + 1, path_png, exc))
            time.sleep(8)
    else:
        alt = path_png.replace('.png', '.new.png')
        fig.savefig(alt, dpi=dpi)
        print('  WARNING: %s locked; wrote %s instead' % (path_png, alt))
    for attempt in range(4):
        try:
            fig.savefig(path_pdf)
            break
        except (PermissionError, OSError) as exc:
            print('  save retry %d/4 %s: %s' % (attempt + 1, path_pdf, exc))
            time.sleep(8)
    else:
        alt = path_pdf.replace('.pdf', '.new.pdf')
        fig.savefig(alt)
        print('  WARNING: %s locked; wrote %s instead' % (path_pdf, alt))


def contour_xy(sub, floor=5.0, endpoint=4.78566):
    """Closed outline of the excluded band at one luminosity.

    sub : boundary rows (dicts) of ONE luminosity with a long-lived edge.
    The outline runs along the long-lived edge from low to high mass, closes at the
    kinematic endpoint (mB - mK, where the band collapses to the scan floor), returns
    along the short-lived crossing (floor where no short crossing exists) and closes on
    the left edge with a vertical segment up to the first long-lived point.
    """
    up = sorted((float(r['mass_GeV']), float(r['tanb_long'])) for r in sub
                if r['tanb_long'] != '')
    if not up:
        return [], []
    low = {float(r['mass_GeV']): (float(r['tanb_short']) if r['tanb_short'] != '' else floor)
           for r in sub if r['tanb_long'] != ''}
    xs = [m for m, _ in up]
    ys = [y for _, y in up]
    if xs[-1] < endpoint:
        xs.append(endpoint)
        ys.append(floor)
    for m, _y in reversed(up):
        xs.append(m)
        ys.append(max(low[m], floor))
    xs.append(xs[0])
    ys.append(up[0][1])
    return xs, ys


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    scan = FA.load_scan(os.path.join(DATA, 'signal_scan.json'))
    rows = FA.boundary_rows(scan)
    colors = {25.0: '#1f77b4', 300.0: '#d62728'}
    lumis = [float(L) for L in scan['lumis']]

    # ---------------- main exclusion figure: one closed outline per luminosity
    style = {300.0: ('-', '#b2182b', r'$\mathcal{L}=300$ fb$^{-1}$'),
             25.0: ('--', '#2166ac', r'$\mathcal{L}=25$ fb$^{-1}$')}
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    for L in (300.0, 25.0):
        sub = [r for r in rows if float(r['lumi_fb']) == L]
        x, y = contour_xy(sub)
        if not x:
            continue
        ls, col, lab = style[L]
        ax.plot(x, y, ls, color=col, lw=2.0, label=lab, solid_capstyle='round')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlim(0.021, 10.0); ax.set_ylim(5.0, 1e5)
    ax.set_xlabel(r'$m_H$ [GeV]')
    ax.set_ylabel(r'$\tan\beta$')
    ax.set_title('Type-I 2HDM light $H$, LHCb downstream (zero background, $s=2.30$, 90% CL)\n'
                 r'$B\to X_s H$: FORESEE ($m_B$) prescription, $m_A=m_{H^\pm}=600$ GeV, $\lambda v^2=0$',
                 fontsize=10.5)
    ax.text(1.05, 210.0, 'excluded (90 % CL)', fontsize=11, color='0.35', style='italic')
    ax.legend(fontsize=9, loc='upper left', framealpha=0.95)
    ax.grid(True, which='both', alpha=0.22)
    fig.tight_layout()
    _save_fig(fig, os.path.join(FIG, 'exclusion_2hdm_lhcb_downstream.png'),
              os.path.join(FIG, 'exclusion_2hdm_lhcb_downstream.pdf'), dpi=180)
    plt.close(fig)

    # ---------------- diagnostic figure: crossings, their tanb error and s-level MC error
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4))
    for L in lumis:
        sub = [r for r in rows if float(r['lumi_fb']) == L]
        for tag, mk, ls, name in (('long', 'o', '--', 'long-lived edge'),
                                  ('short', 's', ':', 'short-lived crossing')):
            pts = [(r['mass_GeV'], float(r['tanb_' + tag]),
                    float(r['tanb_' + tag + '_err']) if r['tanb_' + tag + '_err'] not in ('', 'nan') else np.nan,
                    float(r['rel_ds_at_' + tag]) if r['rel_ds_at_' + tag] not in ('', 'nan') else np.nan)
                   for r in sub if r['tanb_' + tag] != '']
            if not pts:
                continue
            x = np.array([p[0] for p in pts]); y = np.array([p[1] for p in pts])
            e = np.array([p[2] for p in pts]); rel = np.array([p[3] for p in pts])
            axes[0].errorbar(x, y, yerr=e, fmt=mk + ls, ms=2.6, lw=1.0, color=colors[L],
                             alpha=0.85, label='%s, %g fb$^{-1}$' % (name, L))
            axes[1].plot(x, rel, mk + ls, ms=2.6, lw=1.0, color=colors[L], alpha=0.85,
                         label='%s, %g fb$^{-1}$' % (name, L))
    axes[0].set_xscale('log'); axes[0].set_yscale('log')
    axes[0].set_xlim(0.021, 10); axes[0].set_ylim(5, 1e5)
    axes[0].set_xlabel('m_H [GeV]'); axes[0].set_ylabel('tan(beta) at crossing')
    axes[0].set_title('crossings and their tan(beta) MC error')
    axes[0].legend(fontsize=6.5, ncol=2, loc='upper left'); axes[0].grid(True, which='both', alpha=0.2)
    axes[1].axhline(0.05, color='k', lw=1.0, ls='--')
    axes[1].text(0.023, 0.055, '5 % MC target', fontsize=7)
    axes[1].set_xscale('log'); axes[1].set_yscale('log')
    axes[1].set_xlim(0.021, 10); axes[1].set_ylim(0.02, 2.0)
    axes[1].set_xlabel('m_H [GeV]'); axes[1].set_ylabel('delta_s / s at crossing')
    axes[1].set_title('event-level MC error of s at each crossing')
    axes[1].legend(fontsize=6.5, ncol=2, loc='upper left'); axes[1].grid(True, which='both', alpha=0.2)
    fig.tight_layout()
    _save_fig(fig, os.path.join(FIG, 'exclusion_crossings_mc.png'),
              os.path.join(FIG, 'exclusion_crossings_mc.pdf'), dpi=170)
    plt.close(fig)
    print('figures written (main figure: curves only, no vertical lines)')


if __name__ == '__main__':
    main()
