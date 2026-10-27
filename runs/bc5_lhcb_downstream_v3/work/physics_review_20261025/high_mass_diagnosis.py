"""Read-only diagnosis of the high-mass reference-BR residual; audit outputs only."""
from pathlib import Path
import csv, json, math, sys, hashlib
import numpy as np

HERE = Path(__file__).resolve().parent
RUN = HERE.parent.parent
ROOT = RUN.parent.parent
sys.path.insert(0, str(RUN / 'work/analysis'))
sys.path.insert(0, str(RUN / 'work/widths'))
import reference_br_model as ref
import scalar_widths as sw
import run_sensitivity as rs

paths = [RUN/'outputs/data/signal_scan_refbr.json', RUN/'outputs/data/exclusion_boundary_refbr.csv',
         RUN/'work/widths/reference_br_model.py']
out = {'sha256': {str(p.relative_to(RUN)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
scan = json.loads(paths[0].read_text())
theta = np.array(scan['theta2'])
bnd = list(csv.DictReader(paths[1].open()))
repro = RUN/'work/widths/reference_figures'
rep = json.loads((repro/'bc5_digitisation_report.json').read_text())
rr = list(csv.DictReader((repro/'bc5_reference_curves.csv').open()))
branches = {}
for cid, info in rep['curves'].items():
    lum = 300 if '300' in cid else 25
    for seg in info['segments']:
        r = rr[seg['row_start']-2:seg['row_end']-1]
        a = np.array([(float(x['mass_GeV']), float(x['theta2'])) for x in r])
        branches[(lum, seg['branch'])] = a

life = np.genfromtxt(ROOT/'phenomenology/darkscalar/prepared/lifetime.csv', delimiter=',', names=True)
rows = []
for b in bnd:
    m, lum = float(b['mass_GeV']), int(float(b['lumi_fb']))
    if m < 5 or m > 61: continue  # exclude the reference's artificially closed endpoint
    row = scan['signal_per_fb_per_channel'][str(m)]
    signals = np.array([sum(x[c] for c in ['B','K','H']) for x in row])
    h = np.array([x['H'] for x in row])
    y = float(b['theta2_long_lived_edge'])
    yr = branches[(lum,'lower')]
    yref = float(np.exp(np.interp(np.log(m), np.log(yr[:,0]), np.log(yr[:,1]))))
    syref = float(np.exp(np.interp(np.log(yref), np.log(theta), np.log(np.maximum(signals, 1e-300))))) * lum
    s0 = np.interp(np.log(y), np.log(theta), signals)
    hh = np.interp(np.log(y), np.log(theta), h)
    sl = np.gradient(np.log(np.maximum(signals, 1e-300)), np.log(theta))
    slope = float(np.interp(np.log(y), np.log(theta), sl))
    corrected = rs._crossings(theta, signals * .87, 2.3/lum)[0]
    phi = branches[(lum, 'upper')]
    yyref = float(np.exp(np.interp(np.log(m), np.log(phi[:,0]), np.log(phi[:,1]))))
    gtab = float(np.exp(np.interp(np.log(m), np.log(life['mass_GeV']), np.log(life['Gamma_reference_GeV']))))
    rows.append({'mass':m, 'lumi':lum, 'theta2_ours':y,'theta2_ref':yref,
                 'long_ratio':y/yref, 'short_ratio':float(b['theta2_short_lived_edge'])/yyref,
                 'H_fraction':float(hh/s0), 'local_log_slope':slope,
                 'signal_at_paper_edge':syref, 'signal_at_paper_edge_times087':syref*.87,
                 'long_ratio_after087':corrected/yref,
                 'total_width_thermometer_over_supplied_lifetime':ref.gamma_total(m)/gtab,
                 'reference_BRs':{k:v/sum(ref.widths(m).values()) for k,v in ref.widths(m).items()}})
out['residuals'] = rows
for b in rows:
    if b['lumi']==25:
        print('m=%8.4f long=%6.3f after087=%6.3f short=%6.3f H=%7.4f slope=%5.3f G/table=%6.4f required_signal=%5.2f' %
              (b['mass'], b['long_ratio'],b['long_ratio_after087'],b['short_ratio'],b['H_fraction'],b['local_log_slope'],b['total_width_thermometer_over_supplied_lifetime'],b['signal_at_paper_edge']))
out['high_mass_median_long_ratio']={str(l):float(np.median([r['long_ratio'] for r in rows if r['lumi']==l])) for l in [25,300]}
out['high_mass_median_after087']={str(l):float(np.median([r['long_ratio_after087'] for r in rows if r['lumi']==l])) for l in [25,300]}
(HERE/'high_mass_diagnosis.json').write_text(json.dumps(out, indent=2), encoding='utf-8')
print('medians',out['high_mass_median_long_ratio'],out['high_mass_median_after087'])
