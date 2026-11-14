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
import argparse, os, sys

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
        exists = os.path.exists(path)
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
        print('scan stage: nMASS=%d  nsamp=%d  ndecay=%d  m_S in [%.4f, %.4f]  -> %s'
              % (nmass, nsamp, ndecay, args.mmin, args.mmax, scan_name))
        import run_scan
        os.environ['BC5_SCAN_OUT'] = scan_name
        sys.argv = ['run_scan.py', '--nsamp', str(nsamp), '--ndecay', str(ndecay)]
        if masses:
            sys.argv += ['--masses', ','.join('%g' % m for m in masses)]
        else:
            sys.argv += ['--nmass', str(nmass), '--mmin', str(args.mmin), '--mmax', str(args.mmax)]
        run_scan.main()
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
    print()
    print('done.  figures in outputs/figures, data in outputs/data')


if __name__ == '__main__':
    main()
