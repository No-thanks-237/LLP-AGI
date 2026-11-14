# -*- coding: utf-8 -*-
"""Regenerate outputs/sha256_manifest.txt (delivered outputs + code) and inputs/source.sha256.

Both manifests are sha256sum-compatible and are EXCLUDING themselves and any .pyc/__pycache__
entries, so that
    cd <run> && sha256sum -c outputs/sha256_manifest.txt && sha256sum -c inputs/source.sha256
verifies completely (the earlier versions recorded their own empty-file hash and two deleted
pyc files, which made the check fail by construction).

Run:  python work/analysis/make_manifest.py
"""
import glob, hashlib, os

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, '..', '..'))


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def _write(out, rel_files):
    rel_files = sorted(set(rel_files))
    with open(out, 'w', newline='\n') as fh:
        for rel in rel_files:
            fh.write('%s *%s\n' % (sha256(os.path.join(RUN, rel)), rel.replace(os.sep, '/')))
    print('wrote %s with %d entries' % (out, len(rel_files)))


def main():
    # ---- delivered outputs + the code behind them
    files = []
    for root, dirs, names in os.walk(os.path.join(RUN, 'outputs')):
        dirs[:] = sorted(d for d in dirs if d != '__pycache__')
        for n in sorted(names):
            if n.startswith('.') or n.endswith('.new.png') or n.endswith('.new.pdf'):
                continue
            files.append(os.path.relpath(os.path.join(root, n), RUN))
    files += [os.path.relpath(p, RUN) for p in sorted(glob.glob(os.path.join(HERE, '*.py')))]
    files += [os.path.relpath(p, RUN) for p in sorted(glob.glob(os.path.join(HERE, 'grids', '*')))]
    for rel in ('work/production/twohdm_production.py', 'work/production/b_channels.py',
                'work/decay/decay_widths.py', 'work/decay/make_decay_figures.py',
                'work/response/twohdm_lib.py', 'work/response/downstream_ctypes.py'):
        p = os.path.join(RUN, rel)
        if os.path.isfile(p):
            files.append(rel)
    out = os.path.join(RUN, 'outputs', 'sha256_manifest.txt')
    self_rel = os.path.relpath(out, RUN)
    _write(out, [x for x in files if x != self_rel])

    # ---- frozen inputs tree (no .pyc / __pycache__ / self-reference)
    inp = []
    for root, dirs, names in os.walk(os.path.join(RUN, 'inputs')):
        dirs[:] = sorted(d for d in dirs if d != '__pycache__')
        for n in sorted(names):
            if n.endswith('.pyc'):
                continue
            inp.append(os.path.relpath(os.path.join(root, n), RUN))
    out_in = os.path.join(RUN, 'inputs', 'source.sha256')
    self_in = os.path.relpath(out_in, RUN)
    _write(out_in, [x for x in inp if x != self_in])


if __name__ == '__main__':
    main()
