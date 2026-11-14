# -*- coding: utf-8 -*-
"""Parallel launcher for the 2HDM 67->66 mass scan (Lead).

Runs one event_engine.py subprocess per mass, N workers in parallel, writing
work/analysis/scan_2hdm/m_<mass>.json.  Skips masses whose JSON exists and matches
the requested nsamp/b_mode.  Logs each mass to work/analysis/scan_2hdm/logs/.
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, '..', '..'))
ENGINE = os.path.join(RUN, 'work', 'analysis', 'event_engine.py')


def masses_from(path):
    out = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            out.append(float(line.split()[0]))
    return out


def _current_code_sha():
    """Response-code fingerprint of the present source tree (same definition as the engine)."""
    sys.path.insert(0, HERE)
    import response_table as RT
    return RT._code_sha_now()


def done_ok(path, nsamp, b_mode, nparam=None, ndecay=None, code_sha=None):
    """True only if the existing per-mass JSON was produced by the SAME configuration.

    Checks the b_mode, nsamp, nparam, ndecay and the response-code fingerprint; a change of
    any of them (e.g. a new estimator) forces a rebuild instead of a silent skip.
    """
    if not os.path.isfile(path):
        return False
    try:
        d = json.load(open(path))
    except Exception:
        return False
    if b_mode and d.get('b_mode') != b_mode:
        return False
    if nsamp and d.get('table_nsamp') != nsamp and d.get('table_meta', {}).get('nsamp') != nsamp:
        return False
    met = d.get('table_meta', {})
    if nparam is not None:
        try:
            if any(int(blk.get('nparam', -1)) != int(nparam) for blk in d.get('curves', {}).values()):
                return False
        except (TypeError, ValueError):
            return False
    if ndecay:
        if dict(met.get('ndecay', {})) != dict(ndecay):
            return False
    if code_sha is not None and met.get('code_sha256') != code_sha:
        return False
    return 'curves' in d and len(d['curves']) == 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--masses', default=os.path.join(RUN, 'work', 'analysis', 'grids', 'scan_mass_grid.txt'))
    ap.add_argument('--outdir', default=os.path.join(RUN, 'work', 'analysis', 'scan_2hdm'))
    ap.add_argument('--nsamp', default='{"B":17293,"inel":2500,"h":4620}')
    ap.add_argument('--ndecay', default=None)
    ap.add_argument('--b-mode', default='foresee')
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--nparam', type=int, default=40)
    ap.add_argument('--force', action='store_true')
    args = ap.parse_args()
    ms = masses_from(args.masses)
    od = args.outdir
    ld = os.path.join(od, 'logs')
    os.makedirs(od, exist_ok=True)
    os.makedirs(ld, exist_ok=True)
    nsamp = json.loads(args.nsamp)
    ndecay = json.loads(args.ndecay) if args.ndecay else None
    code_sha = _current_code_sha()

    def one(m):
        outp = os.path.join(od, 'm_%g.json' % m)
        if not args.force and done_ok(outp, nsamp, args.b_mode, nparam=args.nparam,
                                      ndecay=ndecay, code_sha=code_sha):
            return (m, 'skip', 0.0, '')
        cmd = [sys.executable, ENGINE, '--model', '2hdm', '--mass', repr(m),
               '--nparam', str(args.nparam), '--param-lo', '5', '--param-hi', '1e5',
               '--lumi', 'both', '--b-mode', args.b_mode, '--nsamp', json.dumps(nsamp),
               '--out', outp]
        if args.ndecay:
            cmd += ['--ndecay', args.ndecay]
        t0 = time.time()
        with open(os.path.join(ld, 'm_%g.log' % m), 'w') as lf:
            p = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
        return (m, 'ok' if p.returncode == 0 else 'FAIL', time.time() - t0, 'rc=%d' % p.returncode)

    t0 = time.time()
    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(one, m): m for m in ms}
        for fut in as_completed(futs):
            m, st, dt, note = fut.result()
            results.append((m, st, dt, note))
            print('%-4s m=%-10.6g %7.1fs %s' % (st, m, dt, note), flush=True)
    nfail = sum(1 for r in results if r[1] == 'FAIL')
    print('DONE %d masses, %d failed, wall %.1f s' % (len(results), nfail, time.time() - t0), flush=True)
    sys.exit(1 if nfail else 0)


if __name__ == '__main__':
    main()
