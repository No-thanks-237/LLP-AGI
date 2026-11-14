# -*- coding: utf-8 -*-
"""Aggregate per-mass event_engine JSONs into outputs/data/signal_scan.json (Lead).

Input : work/analysis/scan_2hdm/m_<mass>.json produced by
        python work/analysis/event_engine.py --model 2hdm --mass <m> --nparam 40 ...
Output: outputs/data/signal_scan.json in the schema consumed by work/analysis/final_analysis.py
Every input file's sha256 is recorded; the merge refuses to mix different b_modes / model tags
unless --allow-mixed is given.
"""
from __future__ import annotations
import argparse, glob, hashlib, json, math, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, '..', '..'))
S90 = -math.log(0.1)


def sha256(path, n=1 << 20):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while True:
            b = f.read(n)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float('nan')


def load_mass_file(path):
    with open(path) as f:
        return json.load(f)


def rows_to_points(res, key_s, key_ds):
    """One mass file -> {param: point dict} (each point already carries s for all lumis)."""
    by_param = {}
    for lname, blk in res.get('curves', {}).items():
        for row in blk.get('rows', []):
            par = _f(row.get('param', row.get('tanb')))
            p = by_param.setdefault(par, {'tanb': par})
            p['Gamma_H_GeV'] = _f(row.get('Gamma_H_GeV'))
            p['ctau_m'] = _f(row.get('ctau_m'))
            p.setdefault('s', {})['%g' % _f(row.get('lumi_fb', lname))] = _f(row.get(key_s))
            p.setdefault('ds', {})['%g' % _f(row.get('lumi_fb', lname))] = _f(row.get(key_ds))
            comp = row.get('s', {}) if isinstance(row.get('s'), dict) else {}
            p.setdefault('s_comp', {})['%g' % _f(row.get('lumi_fb', lname))] = {
                k: _f(v) for k, v in comp.items()}
            if 'eps' in row:
                p['eps'] = {k: _f(v) for k, v in row['eps'].items()}
            if 'sigma_group' in row:
                p['sigma_group_pb'] = {k: _f(v) for k, v in row['sigma_group'].items()}
            for k in ('sigma_h_pb', 'BR_hHH', 'Gamma_h_other_GeV', 'g_C', 'g3', 'b_mode',
                      'validity', 'br'):
                if k in row:
                    p[k] = row[k]
    return by_param


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--indir', default=os.path.join(RUN, 'work', 'analysis', 'scan_2hdm'))
    ap.add_argument('--out', default=os.path.join(RUN, 'outputs', 'data', 'signal_scan.json'))
    ap.add_argument('--b-mode', default=None, help='keep only this b_mode')
    ap.add_argument('--mass-grid', default=os.path.join(RUN, 'work', 'analysis', 'grids',
                                                       'scan_mass_grid.txt'),
                    help='expected mass grid (empty string to skip the check)')
    ap.add_argument('--allow-mixed', action='store_true')
    args = ap.parse_args()
    files = sorted(glob.glob(os.path.join(args.indir, 'm_*.json')))
    if not files:
        raise SystemExit('no per-mass JSONs in ' + args.indir)
    points = []
    masses = set()
    tanb_union = set()
    bmodes = set()
    mass_meta = {}
    fingerprints = {}          # dependency fingerprint -> [mass files]  (review F5)
    for path in files:
        res = load_mass_file(path)
        m = _f(res.get('mass_GeV'))
        masses.add(m)
        bm = None
        by = rows_to_points(res, 's_total', 's_total_se')
        met = res.get('table_meta', {})
        mass_meta['%.6g' % m] = {
            'file': os.path.basename(path), 'sha256': sha256(path),
            'nsamp': met.get('nsamp'), 'ndecay': met.get('ndecay'),
            't_response_s': met.get('t_response_s'), 't_build_s': met.get('t_build_s'),
            'n_scalar': met.get('n_scalar'),
            'n_scalar_nonzero_response': met.get('n_scalar_nonzero_response'),
            'samples': met.get('samples'), 'samples_n_gen': met.get('samples_n_gen'),
            'dll_sha256': met.get('dll_sha256'), 'code_sha256': met.get('code_sha256'),
            'crossings': {str(k): v.get('crossings', []) for k, v in res.get('curves', {}).items()},
        }
        # dependency fingerprint: code / DLL / parent samples / nsamp / ndecay must agree
        # across every per-mass file, otherwise the aggregate mixes incompatible responses
        fp = (met.get('code_sha256'), met.get('dll_sha256'),
              tuple(sorted((met.get('samples') or {}).items())),
              tuple(sorted((met.get('nsamp') or {}).items())),
              tuple(sorted((met.get('ndecay') or {}).items())))
        fingerprints.setdefault(fp, []).append(os.path.basename(path))
        for par, p in sorted(by.items()):
            bm_ = p.get('b_mode')
            if bm_ is not None:
                bmodes.add(bm_)
            if args.b_mode is not None and bm_ is not None and bm_ != args.b_mode:
                continue
            pt = {'mH': m, 'tanb': par}
            pt.update(p)
            points.append(pt)
            tanb_union.add(par)
    if len(bmodes) > 1 and not args.allow_mixed:
        raise SystemExit('mixed b_modes in input: %s (use --b-mode or --allow-mixed)' % sorted(bmodes))
    if len(fingerprints) > 1 and not args.allow_mixed:
        detail = '; '.join('%d files: %s' % (len(v), ', '.join(v[:3]) + ('...' if len(v) > 3 else ''))
                           for v in fingerprints.values())
        raise SystemExit('input files were produced by DIFFERENT code/detector/sample '
                         'configurations (use --allow-mixed to override): ' + detail)
    # expected mass grid (review F5: the aggregate must cover exactly the requested grid)
    if args.mass_grid and os.path.isfile(args.mass_grid):
        with open(args.mass_grid) as fh:
            want = sorted(float(l.split()[0]) for l in fh
                          if l.strip() and not l.lstrip().startswith('#'))
        got = sorted(masses)
        if len(want) != len(got) or any(abs(a - b) > 1e-6 for a, b in zip(want, got)):
            missing = sorted(set(round(x, 6) for x in want) - set(round(x, 6) for x in got))
            extra = sorted(set(round(x, 6) for x in got) - set(round(x, 6) for x in want))
            raise SystemExit('mass grid mismatch: %d expected, %d found; missing %s extra %s'
                             % (len(want), len(got), missing[:5], extra[:5]))
    points.sort(key=lambda q: (q['mH'], q['tanb']))
    out = {
        'model': '2HDM-TypeI-lightH-LHCb-downstream',
        'version': '2HDM-lightH-scan-1.0.0',
        'b_mode': (sorted(bmodes)[0] if len(bmodes) == 1 else sorted(bmodes)),
        'mass_grid': sorted(masses),
        'tanb_grid': sorted(tanb_union),
        'lumis': [25.0, 300.0],
        's_threshold': S90,
        'n_points': len(points),
        'mass_files': mass_meta,
        'points': points,
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, 'w') as f:
        json.dump(out, f, indent=1)
    print('wrote %s : %d points, %d masses, %d tanb, b_mode=%s'
          % (args.out, len(points), len(masses), len(tanb_union), out['b_mode']))


if __name__ == '__main__':
    main()
