# -*- coding: utf-8 -*-
"""BC5 v3.6 end-to-end regression: reproduce the published 25 fb^-1 long-lived boundary.

This is DESIGN.md section 4, cross-check (iv) -- the strongest check of the whole chain
(parent samples -> production -> per-channel response -> lifetime integral -> event OR ->
normalisation -> threshold crossing).  The SAME engine (work/analysis/event_engine.py) is
run with the BC5 model adapter (BC5 v3.6 scalar_widths + BC5 production + theta^2 mixing),
the BC5 v3.6 parent samples and the BC5 hadronisation library.  The engine is NOT allowed any
BC5-specific special case beyond the model adapter; the comparison target is

    <bc5 run>/outputs/data/exclusion_boundary.csv : theta2_long_lived_edge (L = 25 fb^-1)
    <bc5 run>/outputs/data/signal_scan.json       : the full s(mS, theta^2) curve

Usage:
    python bc5_regression.py [--nsamp 3000,1500,3000] [--masses m1,m2,...] [--out ...]
"""
from __future__ import annotations
import argparse, json, math, os, sys, time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, '..', '..'))
BC5_RUN = os.path.abspath(os.path.join(RUN, '..', 'bc5_lhcb_downstream_v3'))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(RUN, 'work', 'response'))

import twohdm_lib as L
import models as M
import response_table as RT
import event_engine as E

BC5_SAMPLES = os.path.join(BC5_RUN, 'work', 'magnus', 'samples')
BC5_LIB = os.path.join(BC5_RUN, 'work', 'magnus', 'hadron_lib')
BC5_DATA = os.path.join(BC5_RUN, 'outputs', 'data')
SAMPLE_FILES = {'B': 'bbar_heavy.txt.gz', 'inel': 'minbias_light.txt.gz', 'h': 'higgs.txt.gz'}

# 16 masses on the BC5 grid spanning the region where the long-lived edge is well defined
DEFAULT_MASSES = [0.483517, 0.661654, 0.905421, 0.9884, 1.0, 1.2, 1.44937, 1.98335,
                  2.2, 2.6834, 3.17488, 3.7579, 4.34457, 5.08225, 6.95466, 10.5834]


def load_published():
    import csv
    edges = {}
    with open(os.path.join(BC5_DATA, 'exclusion_boundary.csv')) as fh:
        for row in csv.DictReader(fh):
            if float(row['lumi_fb']) != 25.0:
                continue
            m = round(float(row['mass_GeV']), 6)
            edges[m] = {'theta2_long': float(row['theta2_long_lived_edge'])
                        if row['theta2_long_lived_edge'] not in ('', 'nan') else None,
                        'theta2_short': float(row['theta2_short_lived_edge'])
                        if row['theta2_short_lived_edge'] not in ('', 'nan') else None,
                        'relerr_long': float(row['rel_s_mc_err_at_long_edge'])
                        if row['rel_s_mc_err_at_long_edge'] not in ('', 'nan') else None}
    with open(os.path.join(BC5_DATA, 'signal_scan.json')) as fh:
        scan = json.load(fh)
    return edges, scan


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--nsamp', type=str, default='3000,1500,3000',
                    help='B,inel,h parents sampled (default 3000,1500,3000 = 2x BC5)')
    ap.add_argument('--ndecay', type=str, default=None)
    ap.add_argument('--masses', type=str, default=None)
    ap.add_argument('--seed', type=int, default=RT.DEFAULT_SEED)
    ap.add_argument('--out', type=str, default=os.path.join(HERE, 'bc5_regression.json'))
    ap.add_argument('--no-cache', action='store_true')
    ap.add_argument('--force', action='store_true')
    ap.add_argument('--stride', type=int, default=2,
                    help='use every Nth published theta^2 grid point (crossings are refined)')
    args = ap.parse_args(argv)
    nb, nk, nh = (int(x) for x in args.nsamp.split(','))
    nsamp = {'B': nb, 'inel': nk, 'h': nh}
    ndecay = json.loads(args.ndecay) if args.ndecay else None
    masses = ([float(x) for x in args.masses.split(',')] if args.masses
              else list(DEFAULT_MASSES))
    # snap every requested mass onto the exact published grid value (the CSV carries 6-7
    # significant digits; a 1e-6 mismatch would also break the edges look-up)
    _pub = sorted(load_published()[0].keys())
    masses = sorted(set(min(_pub, key=lambda x: abs(math.log(x / m))) for m in masses))
    for p in (BC5_SAMPLES, BC5_LIB, BC5_DATA):
        if not os.path.isdir(p):
            raise SystemExit('missing BC5 reference input: ' + p)
    if not L.load_hadron_lib(BC5_LIB):
        raise SystemExit('could not load the BC5 hadron library from ' + BC5_LIB)
    model = M.BC5Model(verbose=True)
    samples = RT.Samples(directory=BC5_SAMPLES, files=SAMPLE_FILES, verbose=True,
                         inel_ev_sub=None)   # BC5 used the full 200k-event minbias sample
    edges, scan = load_published()
    grid_full = np.array(scan['theta2'])
    grid = grid_full[::max(1, int(args.stride))]
    result = {'model': 'bc5', 'regression_of': BC5_RUN, 'nsamp': nsamp, 'seed': args.seed,
              'theta2_grid': len(grid), 'theta2_stride': int(args.stride), 'lumi_fb': 25.0,
              's90': E.S90, 'code_sha256': RT._code_sha_now(),
              'dll_sha256': RT.sha256_file(os.path.join(RUN, 'work', 'response',
                                                        'downstream_kernel.dll')),
              'masses': {}, 'summary': {}}
    t_all = time.time()
    ratios, diffs = [], []
    for mS in masses:
        t0 = time.time()
        tab = RT.get_table(model, mS, samples, seed=args.seed, nsamp=nsamp, ndecay=ndecay,
                           cache=not args.no_cache, force=args.force, verbose=True)
        rows = []
        t1 = time.time()
        for t2 in grid:
            r = E.evaluate(tab, model, mS, float(t2), lumi_fb=25.0)
            rows.append(r)
        t_eval = time.time() - t1
        s = np.array([r['s_total'] for r in rows])
        se = np.array([r['s_total_se'] for r in rows])
        cross = E.find_crossings(grid, s, lifetime_sign=model.param_lifetime_sign)
        # refine every crossing by a local bisection so the comparison is grid-independent
        for c in cross:
            j = int(np.argmin(np.abs(np.log(grid / c['param']))))
            got = E.refine_crossing(tab, model, mS, 25.0,
                                    grid[max(0, j - 1)], grid[min(len(grid) - 1, j + 1)])
            if got:
                c.update(got)
        long_c = next((c.get('param_refined', c['param']) for c in cross
                       if c['edge'] == 'long_lived'), None)
        short_c = next((c.get('param_refined', c['param']) for c in cross
                        if c['edge'] == 'short_lived'), None)
        pub = edges.get(round(float(mS), 6), {})
        pub_long = pub.get('theta2_long')
        rec = {'theta2_long_mine': long_c, 'theta2_long_published': pub_long,
               'theta2_short_mine': short_c, 'theta2_short_published': pub.get('theta2_short'),
               't_table_s': tab['meta'].get('t_build_s'), 't_eval_s': t_eval,
               'n_scalar': int(len(tab['scal_hist'])),
               'n_scalar_nonzero_response': int(tab['meta'].get('n_scalar_nonzero_response', 0)),
               'channels': list(tab['channels']),
               's_at_long_edge_GeV': float(2.302585093)}
        rec['crossings'] = [{k: v for k, v in c.items() if k != 'bracket'} for c in cross]
        if long_c and pub_long:
            rec['boundary_ratio'] = long_c / pub_long
            ratios.append(rec['boundary_ratio'])
        if short_c and pub.get('theta2_short'):
            rec['boundary_ratio_short'] = short_c / pub['theta2_short']
        # MC error of s at the long-lived crossing and the implied boundary error
        if long_c:
            i = int(np.argmin(np.abs(np.log(grid / long_c))))
            rec['s_at_edge'] = float(s[i])
            rec['rel_s_mc_err_at_edge'] = float(se[i] / s[i]) if s[i] > 0 else None
            # local slope dlns/dln(theta2) around the crossing -> boundary MC error
            j = max(0, i - 2); k = min(len(grid) - 1, i + 2)
            if k > j and s[k] > 0 and s[j] > 0:
                slope = (math.log(s[k]) - math.log(s[j])) / (math.log(grid[k]) - math.log(grid[j]))
                rec['dlns_dlntheta2'] = float(slope)
                if slope != 0:
                    rec['rel_boundary_mc_err'] = float(rec['rel_s_mc_err_at_edge'] / abs(slope))
        # curve comparison on the published grid (where the published value is non-negligible)
        pub_row = scan['signal_per_fb_per_channel'].get(str(round(float(mS), 6)))
        if pub_row is None:
            pub_row = scan['signal_per_fb_per_channel'].get(str(float(mS)))
        if pub_row is not None:
            ref_full = 25.0 * np.array([r['B'] + r['K'] + r['H'] for r in pub_row])
            ref = ref_full[::max(1, int(args.stride))][:len(s)]
            mm = (ref > 1e-3) & (s > 1e-12) & np.isfinite(ref)
            if mm.sum() > 5:
                rr = s[mm] / ref[mm]
                rec['curve_ratio_median'] = float(np.median(rr))
                rec['curve_max_abs_ln_ratio'] = float(np.max(np.abs(np.log(rr))))
                diffs.append(rec['curve_ratio_median'])
        result['masses'][str(mS)] = rec
        print('  mS=%-9.6g table %5.1fs eval %4.1fs  edge mine/published = %s  ratio %s  '
              'delta_s/s@edge %.3f' % (
                  mS, rec['t_table_s'] or 0.0, t_eval,
                  ('%.4g' % long_c) if long_c else 'none',
                  ('%.4f' % rec['boundary_ratio']) if 'boundary_ratio' in rec else '-',
                  rec.get('rel_s_mc_err_at_edge') or float('nan')), flush=True)
        del tab
    if ratios:
        ratios = np.array(ratios)
        result['summary'] = {
            'n_masses': int(len(ratios)),
            'boundary_ratio_median': float(np.median(ratios)),
            'boundary_ratio_min': float(ratios.min()),
            'boundary_ratio_max': float(ratios.max()),
            'max_abs_relative_deviation': float(np.max(np.abs(ratios - 1.0))),
            'median_abs_relative_deviation': float(np.median(np.abs(ratios - 1.0))),
            'fraction_within_5pct': float(np.mean(np.abs(ratios - 1.0) <= 0.05)),
            'curve_ratio_median': float(np.median(diffs)) if diffs else None,
        }
    result['wall_s'] = time.time() - t_all
    with open(args.out, 'w') as fh:
        json.dump(result, fh, indent=1)
    print('\nBC5 REGRESSION SUMMARY')
    print(json.dumps(result['summary'], indent=1))
    print('wrote ' + args.out)
    return result


if __name__ == '__main__':
    main()
