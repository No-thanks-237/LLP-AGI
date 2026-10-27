"""Read-only physics diagnostics; writes only in this audit directory.

Run from any directory with python -B audit_physics.py. No production scan is run.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import random
import sys
import tarfile
import numpy as np

HERE = Path(__file__).resolve().parent
RUN = HERE.parent.parent
sys.path.insert(0, str(RUN / 'work' / 'analysis'))
import bc5_lib as L
import bc5_production as P
import run_scan as R
import scalar_widths as W

report = {'date': '2026-10-24', 'model_tag': R.MODEL_TAG, 'production_files_modified': False}
report['source_sha256'] = {str(p.relative_to(RUN)): hashlib.sha256(p.read_bytes()).hexdigest()
    for p in [RUN/'work/analysis/bc5_lib.py', RUN/'work/analysis/bc5_production.py',
              RUN/'work/analysis/run_scan.py', RUN/'work/widths/scalar_widths.py']}

# Retain the exact archived generator source so findings can be linked to line numbers.
for archive, member, target in [('hadronlib2.tar.gz', 'shadlib.cc', 'shadlib_archived.cc'),
                                ('higgs_only2.tar.gz', 'parents.cc', 'higgs_parents_archived.cc')]:
    with tarfile.open(RUN/'work/magnus/bundles'/archive) as t:
        (HERE/target).write_bytes(t.extractfile(member).read())

# All four massless daughters must have the same mean E=M/4 under uniform phase space.
rng = random.Random(24001)
n = 30000
arr = np.array([L.nbody(1.0, [0.0]*4, rng) for _ in range(n)])
report['massless_four_body'] = {'n': n, 'mean_E_over_M': arr[:,:,1].mean(axis=0).tolist(),
    'standard_error': (arr[:,:,1].std(axis=0)/math.sqrt(n)).tolist(),
    'expected_each': 0.25,
    'max_energy_residual': float(np.max(np.abs(arr[:,:,1].sum(axis=1)-1))),
    'max_momentum_residual': float(np.max(np.abs(arr[:,:,2:5].sum(axis=1))))}
print('four-body', report['massless_four_body'], flush=True)

# Uniform massive 3-body phase space dPhi3/ds23 is proportional to both two-body factors.
M, masses = 1.0, [0.13957039]*3
rng = random.Random(24002)
arr = np.array([L.nbody(M, masses, rng) for _ in range(n)])
s23 = (M-arr[:,0,1])**2 - np.sum(arr[:,0,2:5]**2, axis=1)
s = np.linspace((masses[1]+masses[2])**2, (M-masses[0])**2, 20001)
lam1 = (M*M-(masses[0]+np.sqrt(s))**2)*(M*M-(masses[0]-np.sqrt(s))**2)
lam2 = (s-(masses[1]+masses[2])**2)*(s-(masses[1]-masses[2])**2)
weight = np.sqrt(np.maximum(lam1*lam2, 0.0))/s
expected = np.trapezoid(s*weight, s)/np.trapezoid(weight, s)
report['massive_three_body'] = {'n': n, 'mean_s23': float(s23.mean()),
    'se_mean_s23': float(s23.std()/math.sqrt(n)), 'expected_s23': float(expected),
    'mean_E': arr[:,:,1].mean(axis=0).tolist(), 'expected_each_E': M/3}
print('three-body', report['massive_three_body'], flush=True)

report['widths'] = {}
for m in [0.021, 0.1, 0.3, 0.5, 0.99, 1, 1.3, 1.5, 1.99, 2, 2.000001, 3, 5, 10, 30, 62.49]:
    w = W.widths(m)
    report['widths'][str(m)] = {'partial_GeV': w, 'total_GeV': sum(w.values()),
        'BR': {k:v/sum(w.values()) for k,v in w.items()},
        'missing_probability_in_channel_branching': 1-sum(L.channel_branching(m)[0].values())}
report['width_anchors'] = {'alpha_s_2': W.alpha_s(2), 'alpha_s_5': W.alpha_s(5),
    'Gamma_gg_125': W.gamma_gg(125.25), 'Gamma_gammagamma_125': W.gamma_gammagamma(125.25),
    'junction_ratio': W.gamma_total(2.000001)/W.gamma_total(2)}

# Compare the production code with Winkler's actual inclusive prescription, not with itself.
def w_inclusive(m, parent_mass=P.M_BP, gamma=P.GAM_BP):
    return (P.M_B_QUARK/P.V)**2*P.XI_SB**2*(parent_mass**2-m*m)**2/(32*math.pi*parent_mass**3*gamma)
report['mixing_B_production'] = {str(m): {'exclusive_code':sum(P.br_B_to_S(m).values()),
    'winkler_inclusive': w_inclusive(m),
    'code_over_inclusive':sum(P.br_B_to_S(m).values())/w_inclusive(m)}
    for m in [0, 0.3, 1, 2, 3, 4, 4.346, 4.348, 4.69]}
report['Bs_cutoff'] = {str(m):P.br_Bs_inclusive(m) for m in [4.347, 4.348]}

# Exact lifetime integral for A(z)=1 can be evaluated analytically.
report['lifetime_integral'] = []
for z0 in [0, 1.01, 1.10, 1.249]:
    for decay_length in [0.1, 0.3, 1, 10]:
        ms, pz = 1.0, np.array([10.0])
        theta2 = W.HBAR_C/(W.gamma_total(ms)*decay_length*ms/pz[0])
        got = R.per_scalar_p(ms, {}, theta2, pz, np.array([z0]), np.ones((1,len(L.ZGRID))))[0]
        exact = math.exp(-(max(1,z0)-z0)/decay_length)-math.exp(-(2.5-z0)/decay_length)
        report['lifetime_integral'].append({'z0_m':z0, 'decay_length_z_m':decay_length,
            'code':float(got), 'exact':exact, 'ratio':float(got/exact)})

L.load_hadron_lib(str(RUN/'work/magnus/hadron_lib'))
report['hadron_library'] = {}
for ch, c in L.HADRON_LIB.chan.items():
    # Read PDG IDs separately: HadronLib does not retain them.
    path = RUN/'work/magnus/hadron_lib'/('hadron_lib_'+ch+'.txt.gz')
    pdgs = np.loadtxt(path, comments='#', usecols=(2,), dtype=int)
    ids, cnt = np.unique(pdgs, return_counts=True)
    stats = {'nodes':c['mass_nodes'].tolist(), 'n_attempt':int(c['n_attempt']),
        'pid_counts':{str(k):int(v) for k,v in zip(ids,cnt)}, 'node_checks':{}}
    for wanted in [2,5,10,30,62.49]:
        node = float(c['mass_nodes'][np.argmin(abs(c['mass_nodes']-wanted))])
        sel = c['mass']==node
        p2 = c['px'][sel]**2+c['py'][sel]**2+c['pz'][sel]**2
        cs2 = c['pz'][sel]**2/np.maximum(p2,1e-30)
        # Any charged subset of a physical decay must obey E_ch+|p_ch|<=M.
        # The unseen neutral system has nonnegative invariant mass and positive energy.
        starts = c['start'][c['mass'][c['start']]==node]
        energy = np.add.reduceat(c['e'],c['start'])[c['mass'][c['start']]==node]
        qx = np.add.reduceat(c['px'],c['start'])[c['mass'][c['start']]==node]
        qy = np.add.reduceat(c['py'],c['start'])[c['mass'][c['start']]==node]
        qz = np.add.reduceat(c['pz'],c['start'])[c['mass'][c['start']]==node]
        lower = energy+np.sqrt(qx*qx+qy*qy+qz*qz)
        stats['node_checks'][str(wanted)] = {'actual_node':node,
            'tracks':int(sel.sum()), 'mean_cos2_z':float(cs2.mean()),
            'fraction_abs_cos_gt_0p8':float(np.mean(cs2>0.64)),
            'mean_Echarged_over_requested_mass':float(energy.sum()/c['n_attempt']/wanted),
            'nonempty_attempts':len(starts),
            'fraction_charged_subset_impossible_at_requested_mass':float(np.mean(lower>wanted+1e-4)),
            'fraction_charged_subset_impossible_at_node':float(np.mean(lower>node+1e-4)),
            'max_Echarged_plus_p_over_node':float(lower.max()/node)}
    report['hadron_library'][ch]=stats
    print('library',ch,stats['node_checks'],flush=True)

# Pairwise comparison with a random SO(3) rotation of the same library event.
def rotate(d, rng):
    q=np.array([rng.gauss(0,1) for _ in range(4)])
    q/=np.linalg.norm(q)
    w,x,y,z=q
    mat=np.array([[1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)],
                  [2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)],
                  [2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)]])
    return [(a[0],a[1],*list(mat@np.array(a[2:5])),a[5]) for a in d]
report['fixed_axis_acceptance'] = []
for ch,m in [('gg',10.0),('bb',30.0),('tau',5.0)]:
    pz=100.0
    px=8.0
    e=math.sqrt(m*m+px*px+pz*pz)
    rng=random.Random(24003)
    rot_rng=random.Random(24004)
    raw=[];rotated=[]
    for i in range(5000):
        d=L.s_decay_products(m,ch,rng)
        dr=rotate(d,rot_rng)
        vals=[]
        for a in [d,dr]:
            lab=L.boost([v[:5] for v in a],e,(px,0,pz))
            lab=[(*b,a[k][5]) for k,b in enumerate(lab)]
            vals.append(L._passes_pair(lab,2*px/pz,0,2.0))
        raw.append(vals[0]);rotated.append(vals[1])
    delta=np.array(rotated)-np.array(raw)
    report['fixed_axis_acceptance'].append({'channel':ch,'mass':m,'p4':[e,px,0,pz],
        'n':len(raw),'raw':float(np.mean(raw)),'rotated':float(np.mean(rotated)),
        'paired_delta':float(delta.mean()),'paired_delta_se':float(delta.std()/math.sqrt(len(delta)))})
print('axis acceptance',report['fixed_axis_acceptance'],flush=True)

# Read the real Higgs sample and compare its first and last copy.
d=np.loadtxt(RUN/'work/magnus/samples/higgs.txt',comments='#')
ev=d[:,0].astype(int)
unique,first,count=np.unique(ev,return_index=True,return_counts=True)
last=first+count-1
pt=np.hypot(d[:,2],d[:,3])
delta=np.linalg.norm(d[first,2:6]-d[last,2:6],axis=1)
minv2=d[:,5]**2-np.sum(d[:,2:5]**2,axis=1)
hm=np.sqrt(np.maximum(minv2,0))
report['higgs_copies']={'rows':len(d),'unique_events':len(unique),
    'fraction_first_last_differ':float(np.mean(delta>1e-4)),
    'first_zero_pT_fraction':float(np.mean(pt[first]<1e-5)),
    'last_zero_pT_fraction':float(np.mean(pt[last]<1e-5)),
    'mean_first_pT':float(pt[first].mean()),'mean_last_pT':float(pt[last].mean()),
    'mass_first_quantiles':np.quantile(hm[first],[0,.01,.5,.99,1]).tolist(),
    'mass_last_quantiles':np.quantile(hm[last],[0,.01,.5,.99,1]).tolist()}
print('Higgs',report['higgs_copies'],flush=True)

# B and K sample species/yields/eta windows.
report['parents']={}
for channel,name in [('B','bbar_heavy.txt.gz'),('K','minbias_light.txt')]:
    d=np.loadtxt(RUN/'work/magnus/samples'/name,comments='#')
    ids,counts=np.unique(d[:,1].astype(int),return_counts=True)
    pt=np.hypot(d[:,2],d[:,3]); eta=np.arcsinh(d[:,4]/np.maximum(pt,1e-20))
    k={'rows':len(d),'max_event_id':int(d[:,0].max()),'unique_events':len(np.unique(d[:,0])),
       'pid_counts':{str(k):int(v) for k,v in zip(ids,counts)},
       'eta_min':float(eta.min()),'eta_max':float(eta.max()),
       'fraction_eta_2_5':float(np.mean((eta>2)&(eta<5))),
       'fraction_eta_1p5_5p5':float(np.mean((eta>1.5)&(eta<5.5)))}
    report['parents'][channel]=k
    print('parents',channel,k,flush=True)

# Existing scan identity: review current artifacts separately from code.
scan=json.loads((RUN/'outputs/data/signal_scan.json').read_text())
report['scan_metadata']={k:v for k,v in scan.items() if k not in ['signal_per_fb_per_channel','theta2']}
out=HERE/'audit_results.json'
out.write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print('Wrote',out,flush=True)
