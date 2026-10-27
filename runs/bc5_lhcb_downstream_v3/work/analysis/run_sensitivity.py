#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
One-command runner for the BC5 downstream sensitivity.

    python run_sensitivity.py                 # full scan (52 masses) then the s = 2.3 curve
    python run_sensitivity.py --quick          # 8 masses, 200 samples: smoke test, ~2 min
    python run_sensitivity.py --curves-only    # re-plot from the existing signal_scan.json
    python run_sensitivity.py --s 3            # 95% CL curve, written with the _s3 suffix
    python run_sensitivity.py --s 2.3 --s3     # both thresholds in one go

The scan stage is the expensive part: it Monte-Carlos the scalar decays through the LHCb
downstream detector kernel for every mass and theta^2 point (about 20 minutes at the default
settings).  The curve stage is seconds and only needs outputs/data/signal_scan.json, so
--curves-only re-plots without re-running anything.

Everything it needs lives under the run directory:

    work/magnus/samples/            bbar_heavy.txt, minbias_light.txt, higgs.txt + normalisation.json
    work/magnus/hadron_lib/         the Pythia hadronisation libraries
    work/response/downstream_kernel.dll
    work/widths/hadronic_table.csv
    work/analysis/                  this pipeline
"""
from __future__ import annotations
import argparse, json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(RUN, 'work', 'widths'))


def preflight():
    """Check that every input exists; return True if the scan stage can run."""
    import json
    need = [
        ('samples/bbar_heavy.txt', os.path.join(RUN, 'work', 'magnus', 'samples', 'bbar_heavy.txt')),
        ('samples/minbias_light.txt', os.path.join(RUN, 'work', 'magnus', 'samples', 'minbias_light.txt')),
        ('samples/higgs.txt', os.path.join(RUN, 'work', 'magnus', 'samples', 'higgs.txt')),
        ('samples/normalisation.json', os.path.join(RUN, 'work', 'magnus', 'samples', 'normalisation.json')),
        ('response/downstream_kernel.dll', os.path.join(RUN, 'work', 'response', 'downstream_kernel.dll')),
        ('widths/hadronic_table.csv', os.path.join(RUN, 'work', 'widths', 'hadronic_table.csv')),
    ]
    libdir = os.path.join(RUN, 'work', 'magnus', 'hadron_lib')
    for ch in ('gg', 'ss', 'cc', 'bb', 'tautau'):
        need.append(('hadron_lib/hadron_lib_%s.txt.gz' % ch,
                     os.path.join(libdir, 'hadron_lib_%s.txt.gz' % ch)))
    ok = True
    print('preflight:')
    for tag, path in need:
        # load_sample() also accepts the gzipped form, so the preflight must too
        exists = os.path.exists(path) or os.path.exists(path + '.gz')
        ok = ok and exists
        print('  [%s] %s' % ('ok' if exists else 'MISSING', tag))
    npath = os.path.join(RUN, 'work', 'magnus', 'samples', 'normalisation.json')
    if os.path.exists(npath):
        n = json.load(open(npath))
        print('  normalisation: n_gen B=%s  K=%s  H=%s' % (
            n.get('B', {}).get('n_gen'), n.get('K', {}).get('n_gen'), n.get('H', {}).get('n_gen')))
    else:
        print('  WARNING: normalisation.json missing, the MC normalisation will use fallbacks')
    return ok


def _crossings(x, y, level):
    """All x where y = level, linear in log10(x)-log10(y).  Same rule as make_curves.py."""
    out = []
    lx = np.log10(x)
    ly = np.log10(np.maximum(y, 1e-300))
    for i in range(len(x) - 1):
        f0 = ly[i] - math.log10(level)
        f1 = ly[i + 1] - math.log10(level)
        if f0 * f1 < 0:
            t = f0 / (f0 - f1)
            out.append(10 ** (lx[i] + t * (lx[i + 1] - lx[i])))
    return out


def crossing_errors(scan, level, lumi=1.0):
    """For every mass: (theta2_lo, theta2_hi, rel MC error of s at lo, rel MC error at hi).

    The relative MC error of s at a FIXED theta^2 is luminosity independent (both s and its
    standard error scale with L), but the CROSSING is not: the stored signal is per fb^-1, so the
    exclusion at luminosity L sits at s = level / L.  F11: the old code searched for crossings at
    level = 2.3 on the per-fb curve, i.e. it validated the 1 fb^-1 contour instead of the 25/300
    fb^-1 contours that are drawn, and the adaptive rounds therefore optimised the wrong point.
    """
    th = np.array(scan['theta2'], dtype=float)
    out = {}
    for mstr, row in scan['signal_per_fb_per_channel'].items():
        s = np.array([(r['B'] + r['K'] + r['H']) for r in row], dtype=float)
        var = np.array([(r['B'] * r.get('relerr_B', 0.0)) ** 2
                        + (r['K'] * r.get('relerr_K', 0.0)) ** 2
                        + (r['H'] * r.get('relerr_H', 0.0)) ** 2 for r in row], dtype=float)
        rel = np.where(s > 0, np.sqrt(var) / np.maximum(s, 1e-300), 0.0)
        cr = _crossings(th, s, level / lumi)
        lo = min(cr) if cr else float('nan')
        hi = max(cr) if len(cr) >= 2 else float('nan')

        def err_at(x):
            if not np.isfinite(x):
                return float('nan')
            return float(np.interp(math.log10(x), np.log10(th), rel))

        out[float(mstr)] = (lo, hi, err_at(lo), err_at(hi))
    return out


LUMIS = (25.0, 300.0)


def adaptive_plan(scan, level, target, nsamp_base, nsamp_max, growth_cap, lumis=LUMIS):
    """{mass: nsamp} for the masses whose crossing is noisier than target.

    rel ~ nsamp^(-1/2) in the long-lifetime regime, so the new sample size is
    nsamp * (rel/target)^2, capped by growth_cap per round and by nsamp_max.
    """
    cur = {float(k): float(v) for k, v in scan.get('nsamp_by_mass', {}).items()}
    plan = {}
    # F11: optimise the LONG-LIVED edge of every luminosity that is drawn, i.e. take the worst
    # of the 25 and 300 fb^-1 contours (300 fb^-1 is the more demanding one).
    worst = {}
    for lum in lumis:
        for m, (lo, hi, rel_lo, rel_hi) in crossing_errors(scan, level, lum).items():
            rel = rel_lo if np.isfinite(rel_lo) else rel_hi
            if np.isfinite(rel):
                worst[m] = max(worst.get(m, 0.0), rel)
    for m, rel in worst.items():
        if not np.isfinite(rel) or rel <= target:
            continue
        n_now = cur.get(m, float(nsamp_base))
        factor = min((rel / target) ** 2, growth_cap)
        n_new = int(math.ceil(n_now * factor / 100.0) * 100)
        if n_new > n_now:
            plan[m] = int(min(n_new, nsamp_max))
    return plan


def _run_scan(scan_name, nsamp, ndecay, masses, nmass, mmin, mmax, only=None, nsamp_map=None):
    import run_scan
    os.environ['BC5_SCAN_OUT'] = scan_name
    argv = ['run_scan.py', '--nsamp', str(nsamp), '--ndecay', str(ndecay)]
    # NOTE: '%.6f', never '%g'.  run_scan keys every mass by round(m, 6), so a 6-significant-
    # digit '%g' (which is what '%g' gives for masses above 1 GeV, e.g. 1.98335 instead of
    # 1.983347) silently creates a NEW mass point instead of replacing the old one.  That bug
    # turned a 67-point grid into 88 points in the first v3 production run.
    argv += ['--masses', ','.join('%.6f' % m for m in masses)] if masses else []
    if not masses:
        argv += ['--nmass', str(nmass), '--mmin', str(mmin), '--mmax', str(mmax)]
    if only:
        argv += ['--only-masses', ','.join('%.6f' % m for m in only)]
    if nsamp_map:
        path = os.path.join(RUN, 'outputs', 'data', '_nsamp_map.json')
        with open(path, 'w') as fh:
            json.dump({('%.6f' % m): int(n) for m, n in sorted(nsamp_map.items())}, fh, indent=2)
        argv += ['--nsamp-map', path]
    sys.argv = argv
    run_scan.main()


def main():
    ap = argparse.ArgumentParser(description='BC5 downstream sensitivity, full pipeline')
    ap.add_argument('--quick', action='store_true',
                    help='small settings for a smoke test (8 masses, 200 samples, ndecay 5)')
    ap.add_argument('--curves-only', action='store_true',
                    help='skip the scan and only re-plot from the existing signal_scan.json')
    ap.add_argument('--nmass', type=int, help='number of mass points')
    ap.add_argument('--nsamp', type=int, help='scalars sampled per mass point per channel')
    ap.add_argument('--ndecay', type=int, help='decay samplings per scalar for A(z)')
    ap.add_argument('--mmin', type=float, default=0.021)
    ap.add_argument('--mmax', type=float, default=62.49)
    ap.add_argument('--masses', type=str, default=None,
                    help='explicit comma-separated mass list in GeV, e.g. 0.5,1,2,5,10')
    ap.add_argument('--s', type=float, default=2.3, help='threshold for the main curve')
    ap.add_argument('--suffix', type=str, default='', help='suffix for the output file names')
    ap.add_argument('--s3', action='store_true', help='also write the s = 3 (95%% CL) curves')
    ap.add_argument('--adaptive', action='store_true',
                    help='per-mass adaptive parent sample size: run every mass at --nsamp, then '
                         're-run only the masses whose crossing is noisier than --target-rel-err '
                         'with nsamp * (rel/target)^2, up to --adaptive-rounds times')
    ap.add_argument('--target-rel-err', type=float, default=0.05,
                    help='target relative MC error of s at the exclusion crossing (default 5%%)')
    ap.add_argument('--adaptive-rounds', type=int, default=3)
    ap.add_argument('--nsamp-max', type=int, default=25000,
                    help='cap on the per-mass parent sample size in the adaptive rounds')
    ap.add_argument('--force', action='store_true',
                    help='with --adaptive, re-run the base pass even if the scan already exists')
    args = ap.parse_args()

    if args.masses:
        masses = [float(x) for x in args.masses.split(',') if x.strip()]
        nmass = len(masses)
    else:
        masses = None
        nmass = args.nmass or (8 if args.quick else 52)
    nsamp = args.nsamp or (200 if args.quick else 1200)
    ndecay = args.ndecay or (5 if args.quick else 8)

    if not preflight():
        raise SystemExit('preflight failed: fix the missing inputs listed above')

    # anything other than a plain production run writes to its own files, so a test can never
    # overwrite the production scan or the production figures
    sfx = args.suffix + ('_quick' if args.quick else '')
    scan_name = ('signal_scan%s.json' % sfx) if sfx else 'signal_scan.json'

    if not args.curves_only:
        print()
        scan_path = os.path.join(RUN, 'outputs', 'data', scan_name)
        base_needed = True
        if args.adaptive and os.path.exists(scan_path) and not args.force:
            base_needed = False
            print('base pass skipped: %s already exists (use --force to redo it)' % scan_name)
        if base_needed:
            print('scan stage: nsamp=%d  ndecay=%d  s_cut=%.4g  -> %s'
                  % (nsamp, ndecay, args.s, scan_name))
            _run_scan(scan_name, nsamp, ndecay, masses, nmass,
                      args.mmin, args.mmax)
        if args.adaptive:
            for rnd in range(1, args.adaptive_rounds + 1):
                with open(scan_path) as fh:
                    scan = json.load(fh)
                plan = adaptive_plan(scan, args.s, args.target_rel_err, nsamp,
                                     args.nsamp_max, growth_cap=8.0)
                if not plan:
                    print('adaptive: every crossing is within %.1f%% after round %d'
                          % (100.0 * args.target_rel_err, rnd - 1))
                    break
                print()
                print('adaptive round %d: re-running %d masses, nsamp %d..%d  (total %.0f'
                      % (rnd, len(plan), min(plan.values()), max(plan.values()),
                         sum(plan.values())))
                _run_scan(scan_name, nsamp, ndecay, masses, nmass, args.mmin, args.mmax,
                          only=sorted(plan.keys()), nsamp_map=plan)
    else:
        scan = os.path.join(RUN, 'outputs', 'data', 'signal_scan.json')
        if not os.path.exists(scan):
            raise SystemExit('--curves-only needs %s, which does not exist' % scan)
        print('scan stage skipped, re-plotting from %s' % scan)

    print()
    os.environ['BC5_SCAN_IN'] = scan_name
    import make_curves
    jobs = [(args.s, sfx)]
    if args.s3:
        jobs.append((3.0, sfx + '_s3'))
    for s, s_sfx in jobs:
        print('curve stage: s = %.4g, suffix %r' % (s, s_sfx))
        sys.argv = ['make_curves.py', '--s', str(s)] + (['--suffix', s_sfx] if s_sfx else [])
        make_curves.main()
    # The BC5-versus-reference comparison is part of the standard output set: every run must
    # produce it, so that the agreement with the published curve is never checked by hand.
    print()
    print('comparison stage: our contour vs the digitised BC5 reference')
    import plot_comparison_reference
    for s, s_sfx in jobs:
        print('  comparison for s = %.4g, suffix %r' % (s, s_sfx))
        sys.argv = ['plot_comparison_reference.py', '--s', str(s)] + (
            ['--suffix', s_sfx] if s_sfx else [])
        plot_comparison_reference.main()
    # ---- achieved MC error at the crossings (prompt.md 4.3 requires <= 5%) ------------------
    print()
    with open(os.path.join(RUN, 'outputs', 'data', scan_name)) as fh:
        scan = json.load(fh)
    nsamp_by = {float(k): float(v) for k, v in scan.get('nsamp_by_mass', {}).items()}
    def _stat(v):
        if not v:
            return 'none'
        return 'median %.1f%%  max %.1f%%  above 5%%: %d/%d' % (
            100 * float(np.median(v)), 100 * float(np.max(v)),
            sum(1 for x in v if x > args.target_rel_err), len(v))
    # F11: every luminosity has its own contour (s = level / L), so report them separately.
    for lum in LUMIS:
        info = crossing_errors(scan, args.s, lum)
        lo_rels = [v[2] for v in info.values() if np.isfinite(v[2])]
        hi_rels = [v[3] for v in info.values() if np.isfinite(v[3])]
        print('[L = %.0f fb^-1]' % lum)
        print('  MC error of s  delta_s/s   at the long-lived edge: ' + _stat(lo_rels))
        print('  MC error of s  delta_s/s   at the short-lived edge: ' + _stat(hi_rels))
    # The boundary POSITION error is a different quantity: it is obtained by crossing the
    # s -/+ ds bands, and on the steep short-lifetime wall it is much smaller than delta_s/s.
    bnd_path = os.path.join(RUN, 'outputs', 'data', 'exclusion_boundary%s.csv' % sfx)
    if os.path.exists(bnd_path):
        import csv as _csv
        lo_b, hi_b = [], []
        with open(bnd_path) as fh:
            for row in _csv.DictReader(fh):
                for col, val, bucket in (('theta2_long_lived_edge', 'theta2_long_lived_edge_err', lo_b),
                                         ('theta2_short_lived_edge', 'theta2_short_lived_edge_err', hi_b)):
                    try:
                        v = float(row[col]); e = float(row[val])
                    except (KeyError, ValueError):
                        continue
                    if v and abs(v) > 0 and abs(e) == abs(e):
                        bucket.append(abs(e) / abs(v))
        print('boundary-position error  dtheta2/theta2  long-lived edge: ' + _stat(lo_b))
        print('boundary-position error  dtheta2/theta2  short-lived edge: ' + _stat(hi_b))
    if nsamp_by:
        print('parent sample size per mass: min %d  median %d  max %d'
              % (min(nsamp_by.values()), int(np.median(list(nsamp_by.values()))),
                 max(nsamp_by.values())))
    print()
    print('done.  figures in outputs/figures, data in outputs/data')


if __name__ == '__main__':
    main()
