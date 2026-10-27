# -*- coding: utf-8 -*-
"""Collect the T4 validation numbers into logs/response_validation.json.

Reads the unit-test output (re-run), the BC5 regression JSON and the response-table metadata,
and writes a machine-readable summary next to the human-readable logs/response_validation.md.
"""
from __future__ import annotations
import json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(RUN, 'work', 'response'))

import response_table as RT


def main():
    out = {'stage': 'T4 detector response and event engine', 'owner': 'response',
           'design': 'DESIGN.md section 4', 'generated_at': None,
           'provenance': {
               'code_sha256': RT._code_sha_now(),
               'dll_sha256': RT.sha256_file(os.path.join(RUN, 'work', 'response',
                                                         'downstream_kernel.dll')),
               'files': {}},
           'kernel_test': {}, 'unit_tests': {}, 'bc5_regression': {}, 'timing': {},
           'mc': {}, 'limitations': []}
    import datetime
    out['generated_at'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    for rel in ('work/response/twohdm_lib.py', 'work/analysis/models.py',
                'work/analysis/response_table.py', 'work/analysis/event_engine.py',
                'work/analysis/test_engine.py', 'work/analysis/bc5_regression.py'):
        p = os.path.join(RUN, rel)
        if os.path.isfile(p):
            out['provenance']['files'][rel] = RT.sha256_file(p)
    # unit tests
    p = subprocess.run([sys.executable, os.path.join(HERE, 'test_engine.py')],
                       capture_output=True, text=True, cwd=HERE)
    lines = p.stdout.splitlines()
    summ = [l for l in lines if l.startswith('SUMMARY:')]
    out['unit_tests'] = {'summary': summ[-1] if summ else None,
                         'exit_code': p.returncode,
                         'checks': [l[6:] for l in lines if l.startswith('PASS: ')]}
    # kernel suite
    pk = subprocess.run([sys.executable, os.path.join(RUN, 'work', 'response', 'test_kernel.py')],
                        capture_output=True, text=True,
                        cwd=os.path.join(RUN, 'work', 'response'))
    ks = [l for l in pk.stdout.splitlines() if 'SUMMARY:' in l]
    out['kernel_test'] = {'summary': ks[-1] if ks else None,
                          'failed': [l.split('FAILED: ')[1] for l in pk.stdout.splitlines()
                                     if 'FAILED: ' in l]}
    # BC5 regression
    reg = os.path.join(HERE, 'bc5_regression.json')
    if os.path.isfile(reg):
        with open(reg) as fh:
            r = json.load(fh)
        out['bc5_regression'] = {'file': 'work/analysis/bc5_regression.json',
                                 'summary': r.get('summary'),
                                 'nsamp': r.get('nsamp'), 'stride': r.get('theta2_stride'),
                                 'n_masses': len(r.get('masses', {})),
                                 'per_mass': {k: {'ratio': v.get('boundary_ratio'),
                                                  'ratio_short': v.get('boundary_ratio_short'),
                                                  'rel_s_mc_err_at_edge': v.get('rel_s_mc_err_at_edge'),
                                                  'curve_ratio_median': v.get('curve_ratio_median')}
                                              for k, v in r.get('masses', {}).items()}}
    # timing: from the 2HDM test masses (measured; see logs/response_validation.md)
    out['timing'] = {
        'response_table_2hdm_nsamp1500': {
            'mH_0.3': {'table_s': 9.4, 'response_s': 5.1, 'n_scalar': 13458},
            'mH_1.0': {'table_s': 14.0, 'response_s': 9.8, 'n_scalar': 4500},
            'mH_3.07067_nsamp500': {'table_s': 6.4, 'response_s': 2.6, 'n_scalar': 1500},
        },
        'cached_40tanb_2lumi_s': 5.8,
        'per_process_memory_MB': 250,
    }
    out['mc'] = {
        'mH_1.0_boundary_25fb': {'tanb': 357.86, 'rel_s_mc_err': 0.065,
                                 'dlns_dlntanb': -8.0, 'rel_boundary_mc_err': 0.008},
        'recommendation': 'nsamp B=2500-3000, h=2500 gives delta_s/s ~4.5% at the boundary',
    }
    out['limitations'] = [
        'minbias_pik read as a 100k/500k event subsample (N_gen_sub=100000)',
        'T4 review fixes applied (sorted histories/scalars, unique parent ids, sorted sampled events, model-aware edge labels, n_tot/n_smp factor)',
        'A_ch sampled with ndecay=10 (analytic) / 6 (hadron library) decays per channel',
    ]
    path = os.path.join(RUN, 'logs', 'response_validation.json')
    with open(path, 'w') as fh:
        json.dump(out, fh, indent=1)
    print('wrote', path)
    print('unit tests :', out['unit_tests']['summary'])
    print('kernel     :', out['kernel_test']['summary'], out['kernel_test']['failed'])
    print('regression :', json.dumps(out['bc5_regression'].get('summary'), indent=1))


if __name__ == '__main__':
    main()
