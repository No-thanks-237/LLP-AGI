# -*- coding: utf-8 -*-
"""Signal significance map log10(s/s90) in the (mH, tanb) plane (25 fb^-1) with the contour."""
from __future__ import annotations
import math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import final_analysis as FA

FIG = os.path.join(FA.RUN, 'outputs', 'figures')
DATA = os.path.join(FA.RUN, 'outputs', 'data')


def _save(fig, png, pdf, dpi=170):
    import time
    for attempt in range(4):
        try:
            fig.savefig(png, dpi=dpi)
            break
        except (PermissionError, OSError) as exc:
            print('  save retry %d/4 %s: %s' % (attempt + 1, png, exc))
            time.sleep(8)
    else:
        alt = png.replace('.png', '.new.png')
        fig.savefig(alt, dpi=dpi)
        print('  WARNING: locked, wrote', alt)
    for attempt in range(4):
        try:
            fig.savefig(pdf)
            break
        except (PermissionError, OSError):
            time.sleep(8)


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    scan = FA.load_scan(os.path.join(DATA, 'signal_scan.json'))
    by = FA.group_by_mass(scan)
    masses = sorted(by)
    tbs = sorted(set(float(p['tanb']) for p in scan['points']))
    z = np.full((len(tbs), len(masses)), np.nan)
    for j, m in enumerate(masses):
        for p in by[m]:
            i = tbs.index(float(p['tanb']))
            z[i, j] = math.log10(max(FA.svalue(p, 25.0), 1e-30) / FA.S90)
    fig, ax = plt.subplots(figsize=(8.0, 5.4))
    vmax = 6.0
    im = ax.pcolormesh(masses, tbs, np.clip(z, -8, vmax), shading='auto', cmap='viridis')
    cb = fig.colorbar(im, ax=ax, pad=0.015)
    cb.set_label('log10( s(25 fb-1) / s90 )   [clipped at 6]')
    cs = ax.contour(masses, tbs, z, levels=[0.0], colors='w', linewidths=1.6)
    tbg = np.array(tbs)
    for L, ls in ((25.0, '-'), (300.0, '--')):
        rows = FA.boundary_rows(scan)
        pts = [(r['mass_GeV'], float(r['tanb_long'])) for r in rows
               if float(r['lumi_fb']) == L and r['tanb_long'] != '']
        if pts:
            x, y = zip(*pts)
            ax.plot(x, y, ls, color='w', lw=2.4 if L == 25.0 else 1.4,
                    label='%g fb-1 long-lived edge' % L)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlim(0.021, 10); ax.set_ylim(5, 1e5)
    ax.set_xlabel('m_H [GeV]'); ax.set_ylabel('tan(beta)')
    ax.set_title('Expected signal significance in the (m_H, tan beta) plane')
    ax.legend(fontsize=7.5, loc='lower left', framealpha=0.85)
    fig.tight_layout()
    _save(fig, os.path.join(FIG, 'signal_significance_map.png'), os.path.join(FIG, 'signal_significance_map.pdf'))
    plt.close(fig)
    print('heatmap written')


if __name__ == '__main__':
    main()
