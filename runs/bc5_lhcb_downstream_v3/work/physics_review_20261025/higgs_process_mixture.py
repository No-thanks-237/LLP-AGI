"""Diagnose the Pythia Higgs process mixture using a common reference response.

This only writes audit outputs. Reweighting is diagnostic, not a production patch.
The acceptance interpolator was independently compared to the current transport MC.
"""
from pathlib import Path
import json,math,sys,hashlib
import numpy as np
from scipy.interpolate import RegularGridInterpolator
HERE=Path(__file__).resolve().parent;RUN=HERE.parent.parent
sys.path.insert(0,str(RUN/'work/analysis'))
import bc5_production as P

groups={'ggF':[902], 'VBF':[906,907], 'WH':[905], 'ZH':[904],
        'bb_proxy':[901], 'ttH':[908,909], 'gamma_gamma':[903]}
targets={'ggF':P.SIGMA_H['ggF'],'VBF':P.SIGMA_H['VBF'],'WH':P.SIGMA_H['WH'],
         'ZH':P.SIGMA_H['ZH'],'bb_proxy':P.SIGMA_H['bbH'],'ttH':P.SIGMA_H['ttH'],
         'gamma_gamma':0.0}
h=np.loadtxt(RUN/'work/magnus/samples/higgs.txt.gz');h=h[abs(h[:,-1])==62]
codes=h[:,1].astype(int);E=h[:,5];hp=h[:,2:5]
hm=np.sqrt(np.maximum(E*E-(hp*hp).sum(1),1e-12));be=hp/E[:,None];g=E/hm;b2=(be*be).sum(1)
count={k:int(np.isin(codes,ids).sum()) for k,ids in groups.items()}
weights={k:(targets[k]/P.SIGMA_H_TOTAL)/(count[k]/len(h)) for k in groups}
w=np.zeros(len(h))
out={'method':'current Higgs sample, common tabulated conditional response; linear long-lifetime coefficient',
     'normalization_note':'relative group frequencies estimated among retained status-62 Higgs; residual process-dependent generation/filter losses not separately measured',
     'bb_note':'901 is fermion annihilation, mainly bbar; its shape is not the requested 4FS bbH shape',
     'missing_target_processes':{k:v for k,v in P.SIGMA_H.items() if k.startswith('tH')},
     'groups':{},'masses':{}}
for key,ids in groups.items():
    sel=np.isin(codes,ids);w[sel]=weights[key]
    out['groups'][key]={'codes':ids,'count':count[key], 'sample_fraction':count[key]/len(h),
         'target_fraction':targets[key]/P.SIGMA_H_TOTAL,'weight':weights[key],
         'mean_E':float(E[sel].mean()),'mean_pT':float(np.hypot(hp[sel,0],hp[sel,1]).mean())}
a=np.loadtxt(HERE/'reference_acceptance_high.csv',delimiter=',')
a=a[np.lexsort(tuple(a[:,i] for i in [3,2,1,0]))]
axes=[np.unique(a[:,i]) for i in range(4)];eps=(a[:,4]*a[:,5]).reshape(tuple(len(x) for x in axes))
ip=RegularGridInterpolator(tuple(np.log(x) for x in axes),np.log(np.maximum(eps,1e-90)),bounds_error=False,fill_value=math.log(1e-90))
zs=np.linspace(1,2.5,31)
rng=np.random.default_rng(25301);ct=rng.uniform(-1,1,len(h));ph=rng.uniform(-math.pi,math.pi,len(h));st=np.sqrt(1-ct*ct)
for ms in [5.082255,13.023095,33.371213,53.419632,62.49]:
    p=np.sqrt(np.maximum(hm*hm/4-ms*ms,0));v=np.c_[p*st*np.cos(ph),p*st*np.sin(ph),p*ct];bp=(be*v).sum(1)
    # Coefficient per produced Higgs, summing both scalars. Event OR cross terms vanish
    # in this long-lifetime limit. This keeps the independent unit at the Higgs level.
    coeff=np.zeros(len(h))
    for sign in [1,-1]:
        kk=(g-1)*sign*bp/np.maximum(b2,1e-30)+g*hm/2;mom=sign*v+kk[:,None]*be;en=g*(hm/2+sign*bp)
        theta=np.arctan2(np.hypot(mom[:,0],mom[:,1]),mom[:,2])
        sel=(theta>=axes[1][0])&(theta<=axes[1][-1]);t=theta[sel];ee=en[sel]
        profiles=[]
        for z in zs:
            q=np.c_[np.full_like(t,ms),t,ee,np.full_like(t,z)]
            profiles.append(np.exp(ip(np.log(q))))
        prof=np.array(profiles).T
        coeff[sel]+=ms/mom[sel,2]*np.trapezoid(prof,zs,axis=1)
    current=coeff.sum();weighted=(coeff*w).sum();ratio=weighted/current
    # Paired MC uncertainty for the ratio; finite response-table uncertainty is excluded.
    residual=coeff*w-ratio*coeff
    err=math.sqrt(len(h)*np.var(residual,ddof=1))/current
    row={'weighted_over_current_signal':float(ratio),'paired_MC_standard_error':err,
         'long_edge_weighted_over_current':float(1/ratio),
         'per_scalar_K_current':float(current/(2*200000)),
         'per_scalar_K_weighted':float(weighted/(2*200000)),
         'current_signal_fractions':{k:float(coeff[np.isin(codes,ids)].sum()/current) for k,ids in groups.items()}}
    out['masses'][str(ms)]=row
    (HERE/'higgs_process_mixture.json').write_text(json.dumps(out,indent=2),encoding='utf8')
    print(ms,row,flush=True)
print('DONE',flush=True)
