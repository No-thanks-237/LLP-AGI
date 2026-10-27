# -*- coding: utf-8 -*-
"""
Extract the s = s_cut contour (zero background, 95% CL for s_cut = 3) in the (m_S, theta^2)
plane and plot it.  The threshold is a CLI parameter / environment variable:
    python make_curves.py --s 2.3 [--suffix _s23]
Outputs: outputs/figures/bc5_lhcb_downstream_exclusion<suffix>.{png,pdf}
         outputs/data/exclusion_boundary<suffix>.csv, exclusion_s_curves<suffix>.csv,
         exclusion_report<suffix>.json

With zero background, P(n = 0 | s) = exp(-s), so
    s_cut = 3.0   -> 95% CL   (default of prompt.md)
    s_cut = 2.3   -> 90% CL   (exp(-2.3) = 0.1003)
"""
from __future__ import annotations
import csv, json, math, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, '..', '..'))
DATA = os.path.join(RUN, 'outputs', 'data')
FIG = os.path.join(RUN, 'outputs', 'figures')
os.makedirs(FIG, exist_ok=True)
S95 = -math.log(0.05)          # 2.9957; overridden by --s / BC5_S_THRESHOLD


def load(scan_path):
    with open(scan_path) as f:
        d = json.load(f)
    return d


def s_curve(row):
    """row: list over theta^2 of per-channel expected events (B, K, H)."""
    return np.array([r.get('B', 0.0) + r.get('K', 0.0) + r.get('H', 0.0) for r in row])


def s_error(row):
    """MC standard error of the total, combining the channels in quadrature."""
    out = []
    for r in row:
        v = 0.0
        for c in ('B', 'K', 'H'):
            if r.get(c, 0.0) > 0:
                v += (r[c] * r.get('relerr_' + c, 0.0)) ** 2
        out.append(math.sqrt(v))
    return np.array(out)


def crossings(x, y, level):
    """All x where y = level, with linear interpolation in log10(x)-log10(y)."""
    out = []
    lx = np.log10(x); ly = np.log10(np.maximum(y, 1e-300))
    for i in range(len(x) - 1):
        f0, f1 = ly[i] - math.log10(level), ly[i + 1] - math.log10(level)
        if f0 == 0:
            out.append(x[i])
        elif f0 * f1 < 0:
            t = f0 / (f0 - f1)
            out.append(10 ** (lx[i] + t * (lx[i + 1] - lx[i])))
    return out


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--s', type=float,
                    default=float(os.environ.get('BC5_S_THRESHOLD', '2.3')),
                    help='expected-event threshold s_cut (3.0 = 95%% CL, 2.3 = 90%% CL)')
    ap.add_argument('--suffix', type=str, default='', help='suffix for the output file names')
    args = ap.parse_args()
    global S95
    S95 = args.s
    sfx = args.suffix
    print('s_cut = %.4f  (P(n=0|s) = %.4f)  output suffix %r' % (S95, math.exp(-S95), sfx))

    scan = load(os.path.join(DATA, os.environ.get('BC5_SCAN_IN', 'signal_scan.json')))
    masses = np.array(scan['masses'])
    thetas = np.array(scan['theta2'])
    sig = scan['signal_per_fb_per_channel']
    LUMI = [25.0, 300.0]

    rows = []
    curves = {}
    for L in LUMI:
        curves[L] = {'m': [], 'hi': [], 'lo': [], 'hi_err': [], 'lo_err': [],
                     'max_rel_s_err': []}
    for m in masses:
        row = sig[str(m)]
        s1 = s_curve(row)                # per fb^-1
        e1 = s_error(row)                # per fb^-1
        for L in LUMI:
            s = s1 * L
            ds = e1 * L
            cr = crossings(thetas, s, S95)
            if len(cr) >= 2:
                lo, hi = min(cr), max(cr)
                open_band = False
            elif len(cr) == 1:
                # only the long-lifetime edge lies inside the scanned theta^2 window:
                # the excluded band continues beyond THETA2_MAX, so it is left open
                lo, hi = cr[0], float('nan')
                open_band = True
            else:
                lo = hi = float('nan')
                open_band = False
            # boundary uncertainty: crossing of s -/+ ds
            cr_p = crossings(thetas, s + ds, S95)
            cr_m = crossings(thetas, np.maximum(s - ds, 0.0), S95)
            # lo is the long-lifetime edge (small theta^2) and hi the short-lifetime edge
            # (large theta^2), so the crossing that belongs to lo is the SMALLEST one.
            # Taking max for lo and min for hi assigned each error bar to the other branch.
            lo_e = abs(min(cr_m) - min(cr_p)) / 2.0 if (cr_p and cr_m) else float('nan')
            hi_e = abs(max(cr_m) - max(cr_p)) / 2.0 if (cr_p and cr_m) else float('nan')
            rel = np.where(s > 1e-30, ds / np.maximum(s, 1e-30), 0.0)
            rel_at_edge = float(np.interp(math.log10(hi), np.log10(thetas), rel)) if np.isfinite(hi) else float('nan')
            rel_at_lo = float(np.interp(math.log10(lo), np.log10(thetas), rel)) if np.isfinite(lo) else float('nan')
            curves[L]['m'].append(m)
            curves[L]['lo'].append(lo)   # long-lifetime boundary (small theta^2)
            curves[L]['hi'].append(hi)   # short-lifetime boundary (large theta^2)
            curves[L]['lo_err'].append(lo_e)
            curves[L]['hi_err'].append(hi_e)
            curves[L]['max_rel_s_err'].append(float(np.nanmax(rel)))
            rows.append({'mass_GeV': m, 'lumi_fb': L,
                         's_at_theta2_1e-8': float(np.interp(math.log10(1e-8), np.log10(thetas), s1) * L),
                         'rel_s_mc_err_at_1e-8': float(np.interp(math.log10(1e-8), np.log10(thetas), np.where(s1 > 0, e1 / np.maximum(s1, 1e-30), 0))),
                         'max_rel_s_mc_err': float(np.nanmax(np.where(s1 > 0, e1 / np.maximum(s1, 1e-30), 0))),
                         'rel_s_mc_err_at_long_edge': rel_at_lo,
                         'rel_s_mc_err_at_short_edge': rel_at_edge,
                         'theta2_long_lived_edge': lo, 'theta2_long_lived_edge_err': lo_e,
                         'theta2_short_lived_edge': hi, 'theta2_short_lived_edge_err': hi_e,
                         'n_crossings': len(cr), 'band_open_above': open_band})

    with open(os.path.join(DATA, 'exclusion_boundary%s.csv' % sfx), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    with open(os.path.join(DATA, 'exclusion_s_curves%s.csv' % sfx), 'w', newline='') as f:
        f.write('mass_GeV,theta2,s_per_fb,total_25fb,total_300fb\n')
        for m in masses:
            s1 = s_curve(sig[str(m)])
            for t, v in zip(thetas, s1):
                f.write('%.6g,%.6g,%.6g,%.6g,%.6g\n' % (m, t, v, v * 25, v * 300))

    # ---- figure ------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    colors = {25.0: '#1f77b4', 300.0: '#d62728'}
    for L in LUMI:
        c = curves[L]
        m = np.array(c['m']); hi = np.array(c['hi']); lo = np.array(c['lo'])
        good = np.isfinite(hi)
        ax.plot(m[good], hi[good], '-', color=colors[L], lw=2,
                label=r'$\mathcal{L}=%d\ \mathrm{fb}^{-1}$' % L)
        he = np.array(c['hi_err'])
        if np.any(np.isfinite(he)):
            ax.fill_between(m[good], hi[good] * np.exp(-np.nan_to_num(he[good]) / np.maximum(hi[good], 1e-300)),
                            hi[good] * np.exp(+np.nan_to_num(he[good]) / np.maximum(hi[good], 1e-300)),
                            color=colors[L], alpha=0.20, lw=0)
        good2 = np.isfinite(lo)
        if np.any(good2):
            ax.plot(m[good2], lo[good2], '--', color=colors[L], lw=1.6)
        if bool(np.any(~good)):
            ax.plot(m[~good], np.full((~good).sum(), 1e-3), 'v', color=colors[L], ms=5,
                    label=None)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlim(0.021, 62.49); ax.set_ylim(1e-20, 1e-3)
    ax.set_xlabel(r'$m_S$ [GeV]'); ax.set_ylabel(r'$\theta^2$')
    ax.set_title(r'BC5 dark scalar, LHCb downstream (zero background, $s=%.4g$, %.0f%% CL)'
                 % (S95, 100.0 * (1.0 - math.exp(-S95))))
    ax.grid(alpha=0.3, which='both')
    ax.legend(loc='upper right', fontsize=9)
    ax.text(0.03, 0.04, 'solid: short-lifetime edge\ndashed: long-lifetime edge\n(excluded region between them)',
            transform=ax.transAxes, fontsize=7.5, va='bottom')
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'bc5_lhcb_downstream_exclusion%s.png' % sfx), dpi=150)
    fig.savefig(os.path.join(FIG, 'bc5_lhcb_downstream_exclusion%s.pdf' % sfx))
    print('wrote figures')

    rep = {'s_cut': S95, 'cl': 1.0 - math.exp(-S95), 'luminosities_fb': LUMI,
           'n_mass_points': len(masses),
           'points_with_two_edges_25fb': int(np.isfinite(curves[25.0]['hi']).sum()),
           'points_with_two_edges_300fb': int(np.isfinite(curves[300.0]['hi']).sum())}
    with open(os.path.join(DATA, 'exclusion_report%s.json' % sfx), 'w') as f:
        json.dump(rep, f, indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == '__main__':
    main()
