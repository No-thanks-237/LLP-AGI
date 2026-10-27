# -*- coding: utf-8 -*-
"""Final contour extraction, tables and figures for the 2HDM light-H LHCb-downstream scan.

Reads outputs/data/signal_scan.json (schema of work/analysis/aggregate_scan.py) which carries,
per mass, the refined crossings from the engine in mass_files[mass]['crossings'].
"""
from __future__ import annotations
import csv, json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, '..', '..'))
DATA = os.path.join(RUN, 'outputs', 'data')
FIGS = os.path.join(RUN, 'outputs', 'figures')
S90 = -math.log(0.1)


def get_float(d, key):
    if not isinstance(d, dict):
        return None
    if key in d:
        return float(d[key])
    for k, v in d.items():
        try:
            if abs(float(k) - float(key)) < 1e-9:
                return float(v)
        except (TypeError, ValueError):
            continue
    return None


def svalue(pt, lumi):
    v = get_float(pt.get('s', {}), lumi)
    return float(v) if v is not None else float('nan')


def dsvalue(pt, lumi):
    v = get_float(pt.get('ds', {}), lumi)
    return float(v) if v is not None else 0.0


def load_scan(path):
    with open(path) as f:
        return json.load(f)


def group_by_mass(scan):
    out = {}
    for p in scan['points']:
        out.setdefault(float(p['mH']), []).append(p)
    for m in out:
        out[m].sort(key=lambda p: p['tanb'])
    return out


def relerr_at(pts, i, lumi):
    s = svalue(pts[i], lumi)
    d = dsvalue(pts[i], lumi)
    return d / s if s > 0 and d > 0 else float('nan')


def crossings_grid(pts, lumi):
    tb = np.array([p['tanb'] for p in pts], float)
    sv = np.array([svalue(p, lumi) for p in pts], float)
    out = []
    for i in range(len(tb) - 1):
        if not (np.isfinite(sv[i]) and np.isfinite(sv[i + 1])) or sv[i] <= 0 or sv[i + 1] <= 0:
            continue
        a, b = sv[i] - S90, sv[i + 1] - S90
        if a * b < 0:
            t = (math.log(S90) - math.log(sv[i])) / (math.log(sv[i + 1]) - math.log(sv[i]))
            x = math.exp(math.log(tb[i]) + t * (math.log(tb[i + 1]) - math.log(tb[i])))
            out.append({'param_refined': x, 'edge': ('short_lived' if sv[i + 1] > sv[i] else 'long_lived'),
                        'outside_validity_B': None, 's_at_refined': S90, '_grid_i': i})
    return out


def boundary_rows(scan):
    by = group_by_mass(scan)
    mf = scan.get('mass_files', {})
    rows = []
    for m in sorted(by):
        pts = by[m]
        meta = mf.get('%.6g' % m, {})
        cross = meta.get('crossings', {})
        for L in scan['lumis']:
            cr = []
            for k in cross:
                if abs(float(k) - float(L)) < 1e-9:
                    cr = cross[k]
            if not cr:
                cr = crossings_grid(pts, L)
            tbg = np.array([p['tanb'] for p in pts], float)
            rec = {'mass_GeV': m, 'lumi_fb': L, 'n_crossings': len(cr), 'b_mode': scan.get('b_mode'),

                   'tanb_short': '', 'tanb_short_err': '', 'short_outside_validity': '',
                   'rel_ds_at_short': '', 'tanb_long': '', 'tanb_long_err': '',
                   'long_outside_validity': '', 'rel_ds_at_long': '',
                   's_at_tanb_5': svalue(pts[0], L),
                   's_at_tanb_10': svalue(pts[int(np.argmin(np.abs(np.log(tbg) - math.log(10.0))))], L)}
            all_cr = []
            for c in cr:
                tag = 'short' if c.get('edge') == 'short_lived' else 'long'
                tb_b = float(c.get('param_refined', c.get('param')))
                i = int(np.argmin(np.abs(np.log(tbg) - math.log(tb_b))))
                i0, i1 = max(i - 1, 0), min(i + 1, len(pts) - 1)
                sv0, sv1 = svalue(pts[i0], L), svalue(pts[i1], L)
                d0, d1 = dsvalue(pts[i0], L), dsvalue(pts[i1], L)
                # local significance of the separation from the threshold on both sides:
                # a crossing is only "resolved" (usable as a plot boundary) when at least one
                # neighbouring grid point is more than 2 sigma away from s90.  Near the B
                # endpoint (mH=4.784) the signal curve grazes the threshold and the extra
                # crossings are 0.1-0.2 sigma excursions, i.e. not resolved (review F4).
                sig_lo = (sv0 - S90) / d0 if d0 > 0 else float('nan')
                sig_hi = (sv1 - S90) / d1 if d1 > 0 else float('nan')
                resolved = (max(abs(sig_lo) if np.isfinite(sig_lo) else 0.0,
                                abs(sig_hi) if np.isfinite(sig_hi) else 0.0) > 2.0)
                tb0, tb1 = pts[i0]['tanb'], pts[i1]['tanb']
                slope = ((math.log(sv1) - math.log(sv0)) / (math.log(tb1) - math.log(tb0))
                         if sv0 > 0 and sv1 > 0 and tb1 != tb0 else float('nan'))
                re = relerr_at(pts, i, L)
                err = (tb_b * re / abs(slope)) if (np.isfinite(slope) and abs(slope) > 1e-12
                                                   and np.isfinite(re)) else float('nan')
                if tag == 'short' and not resolved:
                    # keep the raw value for the CSV but do not use it as a plot boundary
                    rec['tanb_short_unresolved'] = tb_b
                    all_cr.append({'mass_GeV': m, 'lumi_fb': L, 'b_mode': scan.get('b_mode'),
                                   'index': len(all_cr), 'edge': tag, 'tanb': tb_b,
                                   'tanb_err': err, 'rel_ds_s': re,
                                   'outside_validity': c.get('outside_validity_B'),
                                   's_at_crossing': S90, 'sig_lo': sig_lo, 'sig_hi': sig_hi,
                                   'resolved': resolved})
                    continue
                rec['tanb_' + tag] = tb_b
                rec['tanb_' + tag + '_err'] = err
                rec['rel_ds_at_' + tag] = re
                rec[tag + '_outside_validity'] = c.get('outside_validity_B')
                all_cr.append({'mass_GeV': m, 'lumi_fb': L, 'b_mode': scan.get('b_mode'),
                               'index': len(all_cr), 'edge': tag, 'tanb': tb_b,
                               'tanb_err': err, 'rel_ds_s': re,
                               'outside_validity': c.get('outside_validity_B'),
                               's_at_crossing': S90, 'sig_lo': sig_lo, 'sig_hi': sig_hi,
                               'resolved': resolved})
            # complete crossing list: the two-slot columns above keep only the OUTERMOST long
            # edge and short crossing (for plotting); nothing is lost here (review F4)
            rec['all_crossings'] = ';'.join('%s@%.7g' % (c['edge'], c['tanb']) for c in all_cr)
            rec['_all'] = all_cr
            rows.append(rec)
    return rows


def write_crossings(rows, path):
    """One row per crossing (complete list, unlike the two-slot boundary columns)."""
    cols = ['mass_GeV', 'lumi_fb', 'b_mode', 'index', 'edge', 'tanb', 'tanb_err',
            'rel_ds_s', 'outside_validity', 's_at_crossing', 'sig_lo', 'sig_hi', 'resolved']
    n = 0
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            for c in r.get('_all', []):
                w.writerow({k: c[k] for k in cols})
                n += 1
    print('wrote %s (%d crossings)' % (path, n))
    return n


def write_boundary(rows, path):
    cols = ['mass_GeV', 'lumi_fb', 'n_crossings', 'b_mode',
            'tanb_short', 'tanb_short_err', 'short_outside_validity', 'rel_ds_at_short',
            'tanb_short_unresolved',
            'tanb_long', 'tanb_long_err', 'long_outside_validity', 'rel_ds_at_long',
            'all_crossings', 's_at_tanb_5', 's_at_tanb_10']
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            out = []
            for c in cols:
                v = r.get(c, '')
                if isinstance(v, (float, np.floating)):
                    v = 'nan' if not np.isfinite(v) else '%.6g' % float(v)
                out.append(v)
            w.writerow(out)


def write_s_curves(scan, path):
    lumis = [float(L) for L in scan['lumis']]
    with open(path, 'w') as f:
        f.write('mH_GeV,tanb,' + ','.join('s_%gfb' % L for L in lumis) + ','
                + ','.join('ds_%gfb' % L for L in lumis) + '\n')
        for p in sorted(scan['points'], key=lambda q: (q['mH'], q['tanb'])):
            f.write('%.6g,%.6g' % (p['mH'], p['tanb']))
            for L in lumis:
                f.write(',%.6g' % svalue(p, L))
            for L in lumis:
                f.write(',%.6g' % dsvalue(p, L))
            f.write('\n')


def excluded_intervals(scan, lumi):
    by = group_by_mass(scan)
    out = []
    for m in sorted(by):
        pts = by[m]
        tb = np.array([p['tanb'] for p in pts], float)
        sv = np.array([svalue(p, lumi) for p in pts], float)
        above = np.isfinite(sv) & (sv > S90)
        i = 0
        while i < len(tb):
            if not above[i]:
                i += 1
                continue
            j = i
            while j + 1 < len(tb) and above[j + 1]:
                j += 1
            lo = tb[i] if i == 0 else _cross(tb, sv, i - 1)
            hi = tb[j] if j == len(tb) - 1 else _cross(tb, sv, j)
            out.append({'mass': m, 'tb_lo': lo, 'tb_hi': hi})
            i = j + 1
    return out


def _cross(tb, sv, i):
    a, b = sv[i], sv[i + 1]
    if a <= 0 or b <= 0:
        return float(np.exp(0.5 * (math.log(tb[i]) + math.log(tb[i + 1]))))
    t = (math.log(S90) - math.log(a)) / (math.log(b) - math.log(a))
    return float(math.exp(math.log(tb[i]) + t * (math.log(tb[i + 1]) - math.log(tb[i]))))


def validity_curves(scan):
    by = group_by_mass(scan)
    out = {}
    for label, key in (('0.1', 'tanb_BR_gt_0.1'), ('1', 'tanb_BR_gt_1')):
        xs, ys = [], []
        for m in sorted(by):
            for p in by[m]:
                v = p.get('validity')
                if isinstance(v, dict) and v.get(key):
                    xs.append(m); ys.append(float(v[key])); break
        out[label] = (np.array(xs), np.array(ys))
    return out


def _save_fig(fig, path_png, path_pdf):
    """Save PNG/PDF; tolerate a PDF locked by an external viewer (WPS etc.)."""
    import time
    fig.savefig(path_png, dpi=180)
    for attempt in range(3):
        try:
            fig.savefig(path_pdf)
            return
        except PermissionError:
            print('  WARNING: %s is locked by another program; retry %d/3' % (path_pdf, attempt + 1))
            time.sleep(10)
    print('  WARNING: kept the previous %s (locked)' % path_pdf)


def plot_exclusion(scan, rows, path_png, path_pdf, variants=None):
    """LEGACY BC5-style diagnostic (shaded bands + threshold lines).  Not called by main()."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    fig, ax = plt.subplots(figsize=(7.6, 6.0))
    colors = {25.0: '#1f77b4', 300.0: '#d62728'}
    lumis = [float(L) for L in scan['lumis']]
    for L in lumis:
        for seg in excluded_intervals(scan, L):
            ax.fill_between([seg['mass'] * 0.99, seg['mass'] * 1.01], seg['tb_lo'], seg['tb_hi'],
                            color=colors.get(L, 'grey'), alpha=0.10, lw=0)
    for L in lumis:
        sub = [r for r in rows if float(r['lumi_fb']) == L]
        up = [(r['mass_GeV'], float(r['tanb_long']), float(r['tanb_long_err']))
              for r in sub if r['tanb_long'] != '']
        lo = [(r['mass_GeV'], float(r['tanb_short']), float(r['tanb_short_err']))
              for r in sub if r['tanb_short'] != '']
        if up:
            x, y, e = zip(*up)
            ax.plot(x, y, '-', color=colors.get(L), lw=2.1, label='%g fb$^{-1}$: long-lived edge' % L)
            es, xs, ys = np.array(e, float), np.array(x), np.array(y)
            good = np.isfinite(es) & (es > 0)
            if good.sum() >= 2:
                ax.fill_between(xs[good], ys[good] - es[good], ys[good] + es[good],
                                color=colors.get(L), alpha=0.28, lw=0)
        if lo:
            x, y, e = zip(*lo)
            ax.plot(x, y, '--', color=colors.get(L), lw=1.6, label='%g fb$^{-1}$: short-lived edge' % L)
    vc = validity_curves(scan)
    if len(vc['0.1'][0]):
        ax.plot(vc['0.1'][0], vc['0.1'][1], ':', color='0.35', lw=1.5)
        ax.plot(vc['1'][0], vc['1'][1], ':', color='k', lw=1.8)
    for x, lab in ((0.354, 'K->pi H closes'), (0.996, '2 m_K'), (3.553, '2 m_tau'),
                   (3.729, '2 m_D'), (4.786, 'm_B - m_K')):
        ax.axvline(x, color='0.75', lw=0.8, zorder=0)
        ax.text(x, 5.3e4, lab, rotation=90, fontsize=6.5, color='0.4', ha='right', va='top')
    if variants:
        style = {'paper_mb': ('--', '#2ca02c'), 'exclusive': ('-.', '#9467bd')}
        for name, sc in variants.items():
            ls, col = style.get(name, ('-', 'grey'))
            rr = boundary_rows(sc)
            pts = [(r['mass_GeV'], float(r['tanb_long'])) for r in rr
                   if r['lumi_fb'] == 25.0 and r['tanb_long'] != '']
            if pts:
                x, y = zip(*pts)
                ax.plot(x, y, ls, color=col, lw=1.0, alpha=0.9, label='25 fb$^{-1}$: %s' % name)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlim(0.021, 10.0); ax.set_ylim(5.0, 1e5)
    ax.set_xlabel('m_H [GeV]'); ax.set_ylabel('tan(beta)')
    ax.set_title('Type-I 2HDM light H: LHCb downstream 90% CL (zero background)')
    handles, labels = ax.get_legend_handles_labels()
    extra = [Line2D([], [], color='0.35', ls=':', lw=1.5, label='Br(B->X_s H)=0.1 (EFT validity)'),
             Line2D([], [], color='k', ls=':', lw=1.8, label='Br(B->X_s H)=1 (unphysical)')]
    ax.legend(handles=handles + extra, fontsize=7.5, loc='upper left', ncol=1, framealpha=0.92)
    ax.grid(True, which='both', alpha=0.22)
    fig.tight_layout(); fig.savefig(path_png, dpi=170); fig.savefig(path_pdf)
    plt.close(fig)


def plot_systematics(scan, variants, path_png, path_pdf):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    styles = {'foresee': ('-', '#1f77b4'), 'paper_mb': ('--', '#2ca02c'), 'exclusive': ('-.', '#9467bd')}
    allscans = {'foresee': scan}
    allscans.update(variants)
    for name, sc in allscans.items():
        for L, lw in ((25.0, 1.8), (300.0, 1.0)):
            rr = boundary_rows(sc)
            pts = [(r['mass_GeV'], float(r['tanb_long'])) for r in rr
                   if float(r['lumi_fb']) == L and r['tanb_long'] != '']
            if not pts:
                continue
            x, y = zip(*pts)
            ls, col = styles.get(name, ('-', 'grey'))
            ax.plot(x, y, ls, color=col, lw=lw, label='%s, %g fb$^{-1}$' % (name, L))
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlim(0.021, 10.0); ax.set_ylim(5.0, 1e5)
    ax.set_xlabel('m_H [GeV]'); ax.set_ylabel('tan(beta)')
    ax.set_title('B -> X_s H prescription systematics: long-lived edge')
    ax.legend(fontsize=7, loc='upper left', ncol=2)
    ax.grid(True, which='both', alpha=0.22)
    fig.tight_layout(); fig.savefig(path_png, dpi=170); fig.savefig(path_pdf)
    plt.close(fig)


def write_widths_table():
    src = os.path.join(DATA, 'decay_widths_grid.csv')
    if not os.path.isfile(src):
        return
    lines = [ln for ln in open(src) if not ln.startswith('#')]
    rows = list(csv.DictReader(lines))
    tbs = sorted(set(float(r['tanb']) for r in rows))
    picked = set()
    for t in (10.0, 100.0, 1000.0):
        picked.add(min(tbs, key=lambda x: abs(math.log(x) - math.log(t))))
    keep = [r for r in rows if float(r['tanb']) in picked]
    with open(os.path.join(DATA, 'twohdm_widths_br.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in keep:
            w.writerow(r)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--scan', default=os.path.join(DATA, 'signal_scan.json'))
    ap.add_argument('--variant-dir', default=os.path.join(RUN, 'work', 'analysis'))
    args = ap.parse_args()
    scan = load_scan(args.scan)
    rows = boundary_rows(scan)
    write_boundary(rows, os.path.join(DATA, 'exclusion_boundary.csv'))
    write_s_curves(scan, os.path.join(DATA, 'exclusion_s_curves.csv'))
    write_widths_table()
    variants = {}
    for name in ('paper_mb', 'exclusive'):
        p = os.path.join(args.variant_dir, 'scan_2hdm_%s' % name, 'signal_scan.json')
        if os.path.isfile(p):
            variants[name] = load_scan(p)
    # complete crossing list: nominal + every prescription variant (review F4: the variants
    # carry three-crossing points such as paper_mb 3.992 GeV and exclusive 4.784 GeV)
    rows_all = list(rows)
    for name, sc in variants.items():
        try:
            rows_all += boundary_rows(sc)
        except Exception as exc:      # a broken variant must not kill the nominal products
            print('warning: could not derive crossings for %s (%s)' % (name, exc))
    write_crossings(rows_all, os.path.join(DATA, 'exclusion_crossings.csv'))
    # The main exclusion figure is produced by final_figure.py (ONE closed outline per
    # luminosity: 300 fb^-1 solid, 25 fb^-1 dashed, upper and lower limits connected).
    # plot_exclusion() below is kept only as the legacy BC5-style diagnostic and is not
    # called any more: it would overwrite the main figure.
    if variants:
        plot_systematics(scan, variants, os.path.join(FIGS, 'exclusion_systematics_bmode.png'),
                         os.path.join(FIGS, 'exclusion_systematics_bmode.pdf'))
    n2 = sum(1 for r in rows if r['n_crossings'] == 2)
    print('boundary rows %d (%d with two crossings); b_mode=%s' % (len(rows), n2, scan.get('b_mode')))
    for L in (25.0, 300.0):
        sub = [r for r in rows if float(r['lumi_fb']) == L and r['tanb_long'] != '']
        if sub:
            rel = [float(r['rel_ds_at_long']) for r in sub if np.isfinite(float(r['rel_ds_at_long']))]
            print(' L=%g: long-lived edge %.1f..%.1f (median %.1f), MC rel median %.3f max %.3f'
                  % (L, min(float(r['tanb_long']) for r in sub), max(float(r['tanb_long']) for r in sub),
                     float(np.median([float(r['tanb_long']) for r in sub])),
                     np.median(rel), np.max(rel)))


if __name__ == '__main__':
    main()
