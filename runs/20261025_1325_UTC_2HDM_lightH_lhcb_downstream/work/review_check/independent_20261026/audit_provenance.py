"""Check physical input fingerprints and the frozen input manifest, without modifying them."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
RUN = HERE.parents[2]


def sha(path):
    with open(path, 'rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def main():
    checks = {}
    norm = json.loads((RUN / 'work/magnus/samples/normalisation.json').read_text())
    for name in ('B', 'inel', 'H'):
        d = norm[name]
        path = RUN / 'work/magnus/samples' / d['file']
        checks[name] = {'matches': sha(path) == d['sha256'], 'bytes': path.stat().st_size}
    man = json.loads((RUN / 'work/magnus/hadron_lib/manifest.json').read_text())
    for name, d in man['channels'].items():
        checks[name] = {'matches': sha(RUN / 'work/magnus/hadron_lib' / d['file']) == d['sha256']}
    result = json.loads((RUN / 'work/analysis/scan_2hdm/m_1.json').read_text())
    checks['DLL'] = {'matches': sha(RUN / 'work/response/downstream_kernel.dll') == result['table_meta']['dll_sha256']}
    first = result['curves']['25.0']['rows'][0]
    for key, file in (('production_sha256', 'work/production/twohdm_production.py'),
                      ('decay_sha256', 'work/decay/decay_widths.py'),
                      ('decay_binary_sha256', 'work/decay/build/light_scalar_decay.exe')):
        checks[key] = {'present': key in first, 'matches': sha(RUN / file) == first.get(key)}
    bad = []
    for line in (RUN / 'inputs/source.sha256').read_text().splitlines():
        expected, rel = line.split(' *', 1)
        path = RUN / rel
        if not path.is_file() or sha(path) != expected:
            bad.append(rel)
    checks['inputs_manifest_mismatches'] = bad
    (HERE / 'provenance_results.json').write_text(json.dumps(checks, indent=2), encoding='utf-8')
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
