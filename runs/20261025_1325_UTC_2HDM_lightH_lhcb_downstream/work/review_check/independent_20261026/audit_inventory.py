"""Read-only audit of the published light-H scan; writes only into this audit folder."""
from pathlib import Path
import csv
import hashlib
import json
import math
import sys
from collections import Counter

import numpy as np

HERE = Path(__file__).resolve().parent
RUN = HERE.parents[2]
sys.path.insert(0, str(RUN / 'work/analysis'))
import response_table as RT
import event_engine as E
import models as M
import final_analysis as FA
import launch_scan as LS
import twohdm_lib as L


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def check_scan(mode):
    folder = RUN / 'work/analysis' / ('scan_2hdm' if mode == 'foresee' else 'scan_2hdm_' + mode)
    aggregated = RUN / 'outputs/data/signal_scan.json' if mode == 'foresee' else folder / 'signal_scan.json'
    scan = read(aggregated)
    files = sorted(folder.glob('m_*.json'))
    ds = [read(p) for p in files]
    version = Counter()
    nsamp = Counter()
    code = Counter()
    mismatch = []
    residual = []
    errors_exact, errors_grid = [], []
    crossings = Counter()
    mc_by_m = []
    lumi_dev = []
    brsum_dev = []
    recomposition = []
    missing_mass_files = []
    grid = np.loadtxt(RUN / 'work/analysis/grids/scan_mass_grid.txt')
    model = M.TwoHDMModel(verbose=False, b_mode=mode)
    point_map = {(p['mH'], p['tanb']): p for p in scan['points']}
    all_relative_s = []
    channel_components = {k: [] for k in ('B', 'inel', 'h')}
    for path, d in zip(files, ds):
        m = d['mass_GeV']
        met = d['table_meta']
        code[met.get('code_sha256')] += 1
        nsamp[json.dumps(d['table_nsamp'], sort_keys=True)] += 1
        a = scan['mass_files'].get('%.6g' % m)
        if not a or a['sha256'] != sha(path):
            mismatch.append(path.name)
        tabpath = Path(RT.cache_path(model, m))
        if not tabpath.is_file():
            missing_mass_files.append(str(tabpath))
            continue
        tab = RT.load_table(tabpath)
        version[tab['version']] += 1
        for lname, blk in d['curves'].items():
            Lfb = float(lname)
            crossings[len(blk['crossings'])] += 1
            for row in blk['rows']:
                par = row['param']
                p = point_map[(m, par)]
                old = p['s']['%g' % Lfb]
                if old != row['s_total']:
                    mismatch.append([m, par, lname, 'aggregated point differs'])
                if Lfb == 25:
                    all_relative_s.append(row['s_total_se']/row['s_total'] if row['s_total'] > 0 else 0)
                    for k in channel_components:
                        channel_components[k].append(row['s'][k]/row['s_total'] if row['s_total'] > 0 else 0)
                den = row['s_total']
                rec = 1e3*Lfb*sum(row['sigma_group'][k]*row['eps'][k] for k in M.GROUP_NAMES)
                recomposition.append(abs(rec/den-1) if den else abs(rec))
                brsum_dev.append(abs(sum(row['br'].values())-1))
                if Lfb == 300:
                    lo = point_map[(m, par)]['s']['25']
                    lumi_dev.append(abs(den/(12*lo)-1) if lo else abs(den))
            for c in blk['crossings']:
                val = E.evaluate(tab, model, m, c.get('param_refined', c['param']), Lfb)
                residual.append(abs(val['s_total']/E.S90-1))
                errors_exact.append(val['s_total_se']/val['s_total'])
                pp = np.array(d['param_grid'])
                i = np.argmin(abs(np.log(pp/c.get('param_refined', c['param']))))
                near = blk['rows'][i]
                errors_grid.append(near['s_total_se']/near['s_total'])
                mc_by_m.append([m, Lfb, errors_exact[-1], errors_grid[-1]])
    boundary = FA.boundary_rows(scan)
    nonempty = [r for r in boundary if r['tanb_long'] != '']
    return {'n_mass_files': len(ds), 'n_points': len(scan['points']),
            'mass_grid': scan['mass_grid'], 'expected_grid_equal': bool(np.array_equal(sorted(grid),scan['mass_grid'])),
            'table_versions': dict(version), 'nsamp': dict(nsamp), 'table_code_shas': dict(code),
            'current_response_sha': RT._code_sha_now(), 'aggregate_mismatches': mismatch,
            'missing_tables': missing_mass_files, 'crossings_count': dict(crossings),
            'max_crossing_relative_residual': max(residual, default=None),
            'exact_crossing_mc_median_max': [float(np.median(errors_exact)), max(errors_exact)] if errors_exact else [],
            'nearest_grid_mc_median_max': [float(np.median(errors_grid)), max(errors_grid)] if errors_grid else [],
            'mc_by_mass': mc_by_m,
            'max_signal_recomposition_error': max(recomposition), 'max_lumi_scaling_error': max(lumi_dev),
            'max_br_sum_error': max(brsum_dev),
            'grid_signal_relative_mc_median_max': [float(np.median(all_relative_s)),max(all_relative_s)],
            'largest_signal_fraction': {k:max(v) for k,v in channel_components.items()},
            'long_edge_range': {str(l):[min(r['tanb_long'] for r in nonempty if r['lumi_fb']==l),
                                      max(r['tanb_long'] for r in nonempty if r['lumi_fb']==l)] for l in (25,300)}}


def main():
    out = {}
    man = RUN / 'outputs/sha256_manifest.txt'
    bad, missing = [], []
    for line in man.read_text().splitlines():
        expected, rel = line.split(' *', 1)
        path = RUN / rel
        if not path.is_file():
            missing.append(rel)
        elif sha(path) != expected:
            bad.append({'file':rel, 'expected':expected, 'actual':sha(path)})
    out['manifest']={'entries':len(man.read_text().splitlines()),'mismatches':bad,'missing':missing}
    for mode in ('foresee','paper_mb','exclusive'):
        print('checking',mode,flush=True)
        out[mode]=check_scan(mode)
    nominal = read(RUN / 'work/analysis/scan_2hdm/m_1.json')
    fpath = RUN / 'work/analysis/scan_2hdm/m_1.json'
    out['resume_ignores_nparam_and_ndecay'] = {
        'done_ok':LS.done_ok(str(fpath),nominal['table_nsamp'],'foresee'),
        'existing_nparam':len(nominal['param_grid']), 'existing_ndecay':nominal['table_ndecay']}
    out['cache_validation'] = {}
    model = M.TwoHDMModel(verbose=False)
    tab = RT.load_table(RT.cache_path(model,1.0))
    L.load_hadron_lib(str(RUN/'work/magnus/hadron_lib'))
    for key, val in (('dll_sha256','intentionally altered'),('samples',{'B':'intentionally altered'}),
                     ('detector_config',{'by':99}),('birth_cuts',{'minE':999})):
        import copy
        copytab = copy.deepcopy(tab)
        copytab['meta'][key]=val
        out['cache_validation'][key] = RT._cache_matches(copytab,tab['nsamp'],tab['meta']['ndecay'],tab['meta']['seed'])
    outpath = HERE / 'inventory_results.json'
    outpath.write_text(json.dumps(out,indent=2,allow_nan=True),encoding='utf-8')
    compact={k:v for k,v in out.items() if k not in ('foresee','paper_mb','exclusive')}
    compact['scans']={k:{kk:vv for kk,vv in out[k].items() if kk not in ('mass_grid','mc_by_mass')}
                      for k in ('foresee','paper_mb','exclusive')}
    print(json.dumps(compact,indent=2),flush=True)


if __name__ == '__main__':
    main()
