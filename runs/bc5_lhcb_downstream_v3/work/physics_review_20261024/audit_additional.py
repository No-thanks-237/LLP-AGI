"""Further read-only checks; temporary substitutions are in-memory diagnostics only."""
from pathlib import Path
import json, math, random, sys
import numpy as np
HERE=Path(__file__).resolve().parent
RUN=HERE.parent.parent
sys.path.insert(0,str(RUN/'work/analysis'))
import run_scan as R
import bc5_lib as L
import bc5_production as P
import scalar_widths as W
import run_sensitivity as RS

out={}
scan=json.loads((RUN/'outputs/data/signal_scan.json').read_text())
out['mc_at_actual_crossings']={}
for lumi in [1,25,300]:
    errors=RS.crossing_errors(scan,2.3/lumi)
    lo=[v[2] for v in errors.values() if np.isfinite(v[2])]
    hi=[v[3] for v in errors.values() if np.isfinite(v[3])]
    out['mc_at_actual_crossings'][str(lumi)]={
        'long_median':float(np.median(lo)), 'long_max':float(np.max(lo)),
        'long_over_5pct':sum(x>.05 for x in lo),'long_count':len(lo),
        'short_median':float(np.median(hi)), 'short_max':float(np.max(hi)),
        'short_over_5pct':sum(x>.05 for x in hi),'short_count':len(hi)}
print('mc errors',out['mc_at_actual_crossings'],flush=True)

# Gamma is in total widths but absent from CHANNELS: its leftover draw is sent to gg.
out['photon_leftover']={}
for m in [.021,.1,.3,1,1.5,2.2,5,30]:
    w=W.widths(m)
    out['photon_leftover'][str(m)]={'BR_gamma':w['gamma']/sum(w.values()),
        'sum_sampled_BR':sum(L.channel_branching(m)[0].values())}

# The strict-Winkler fixed-mass s/c prescription does not reach gg loop masses.
orig=W._Q_LOOP
out['loop_mass_mismatch']={}
for m in [2.00001,3,5,10]:
    old=W.gamma_gg(m)
    W._Q_LOOP=(W.M_U_MSBAR,W.M_D_MSBAR,W.W_M_S,W.W_M_C,W.W_M_B,W.M_T_POLE)
    new=W.gamma_gg(m)
    W._Q_LOOP=orig
    out['loop_mass_mismatch'][str(m)]={'actual':old,'fixed_Winkler_masses':new,'ratio':old/new}

# K_SS and K_SS_pi are mutually exclusive but have the same hid in run_mass.
tab={'chan':np.array(['K']*4),'ev':np.array([0]*4),'par':np.array([0]*4),
     'hid':np.array([1]*4),'kq':np.array([0.0]*4),
     'w':np.array([4.4e-13]*2+[6.6e-15]*2),'counts':{'K':[1,1]},
     'samp_ev':{'K':np.array([0])}}
p=np.array([.1]*4)
actual=R.channel_efficiency(tab,p,'K',1,1)[0]
tab['hid']=np.array([1,1,2,2])
correct=R.channel_efficiency(tab,p,'K',1,1)[0]
out['K_quartic_history_collision']={'actual':actual,'separate_histories':correct,'ratio':actual/correct}

# Quantify two-body-channel normalization from independent Winkler expressions.
def wk(m,parent=P.M_BP,mfin=P.M_KP):
    lam=(parent**2-(mfin+m)**2)*(parent**2-(mfin-m)**2)/parent**4
    f=.33/(1-m*m/37.5)
    mat=.25*(parent**2-mfin**2)**2/(P.M_B_QUARK-P.M_S_QUARK)**2*f*f
    g2=(P.M_B_QUARK/P.V)**2*P.XI_SB**2
    return g2*mat*math.sqrt(max(lam,0))/(16*math.pi*parent*P.GAM_BP)
out['exclusive_K_vs_Winkler']={str(m):{'code':P.br_B_to_S(m).get('K',0),'Winkler':wk(m)}
    for m in [0,1,3,4.7,4.78]}

# Conditional H acceptance: current first record vs last record; literal (unrotated)
# and isotropically rotated library. Same parent selection seed in all configurations.
L.load_hadron_lib(str(RUN/'work/magnus/hadron_lib'))
d=np.loadtxt(RUN/'work/magnus/samples/higgs.txt',comments='#')
_,first,count=np.unique(d[:,0],return_index=True,return_counts=True)
last=first+count-1
def sample(rows):
    a=d[rows]
    return {k:a[:,i].astype(int) if i in [0,1] else a[:,i]
        for i,k in enumerate(['ev','pdg','px','py','pz','e','x','y','z','t'])}
R._SAMPLE_CACHE[R.SAMPLE_B]=None
R._SAMPLE_CACHE[R.SAMPLE_K]=None
first_smp=sample(first); last_smp=sample(last)
original_decay=L.s_decay_products
rot_rng=random.Random(2)
def rotated_decay(ms,ch,rng):
    rest=original_decay(ms,ch,rng)
    # Uniform SO(3), rotating the whole correlated final state together.
    q=np.array([rot_rng.gauss(0,1) for _ in range(4)])
    q/=np.linalg.norm(q)
    w,x,y,z=q
    mat=np.array([[1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)],
                  [2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)],
                  [2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)]])
    return [(a[0],a[1],*(mat@np.array(a[2:5])),a[5]) for a in rest]
out['higgs_acceptance_diagnostic']=[]
for ms in [2.0,5.082255,28.527411,62.49]:
    values={}
    for mode in ['first','last','first_rotated']:
        R._SAMPLE_CACHE[R.SAMPLE_H]=last_smp if mode=='last' else first_smp
        L.s_decay_products=rotated_decay if mode=='first_rotated' else original_decay
        R._np_rng=np.random.default_rng(24006)
        rot_rng=random.Random(24007)
        tab=R.run_mass(ms,10000,random.Random(24005),ndecay=5,verbose=False)
        results={}
        for ctau in [.001,.01,.1,1,10]:
            t2=W.HBAR_C/(W.gamma_total(ms)*ctau)
            p=R.per_scalar_p(ms,tab,t2,tab['pz'],tab['z'],tab['A'])
            eps,se=R.channel_efficiency(tab,p,'H',R.N_GEN_HIGGS,t2)
            results[str(ctau)]={'epsilon':eps,'se':se}
        values[mode]=results
    out['higgs_acceptance_diagnostic'].append({'mass':ms,'nparents':10000,'ndecay':5,'results':values})
    print('Higgs diag',ms,values,flush=True)
L.s_decay_products=original_decay
(HERE/'additional_results.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print('wrote additional_results.json',flush=True)
