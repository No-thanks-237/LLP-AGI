"""Check precision and two-step vertex approximation without changing scan inputs."""
from pathlib import Path
import json
import math
import random
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
RUN = HERE.parents[2]
sys.path.insert(0,str(RUN/'work/analysis'))
import models as M
import response_table as RT
import twohdm_lib as L
import event_engine as E


def stats(v):
    a=np.asarray(v,dtype=float)
    return dict(zip(('min','p01','median','p99','max'), map(float,np.quantile(a,[0,.01,.5,.99,1]))))


def stable_boost(parts, e_b, p_b, mass=None):
    p=np.asarray(p_b,dtype=float)
    if mass is None:
        mass=math.sqrt(max(float(e_b)**2-float(p@p),1e-30))
    energy=math.sqrt(float(p@p)+mass*mass)
    gamma=energy/mass
    out=[]
    for m,e,px,py,pz in parts:
        q=np.array([px,py,pz],dtype=float)
        pq=float(p@q)
        qq=q+p*(e/mass+pq/(mass*(energy+mass)))
        out.append((m,gamma*e+pq/mass,*qq))
    return out


def main():
    out={}
    samples=RT.Samples(verbose=False)
    b=samples.get('B')
    sel=np.isin(b['pdg'],[521,-521,511,-511])
    inds=np.flatnonzero(sel)
    p=np.stack([b[k][inds].astype(float) for k in ('px','py','pz')],axis=1)
    energy=b['e'][inds].astype(float)
    m2=energy**2-np.sum(p*p,axis=1)
    b2_float32=sum((b[k][inds]/b['e'][inds])**2 for k in ('px','py','pz'))
    out['B_precision']={'n':len(inds),'non_timelike_float64_count':int((m2<=0).sum()),
        'non_timelike_beta_float32_count':int((b2_float32>=1).sum()),
        'mass_GeV':stats(np.sqrt(np.maximum(m2,0))),
        'mass_GeV_pz_gt_1500':stats(np.sqrt(np.maximum(m2[np.abs(p[:,2])>1500],0)))}
    # Replay EXACT candidate random seeds.  Full candidate builders are used, so the
    # first B histories coincide with the cache; inelastic and h load lazily later.
    model=M.TwoHDMModel(verbose=False)
    L.load_hadron_lib(str(RUN/'work/magnus/hadron_lib'))
    allp={}
    # Only B sample is required for the precision check; preserve candidate B random stream.
    class BOnly:
        n_gen={'B':samples.n_gen['B']}
        def get(self,key): return b if key=='B' else None
    for mh in (.021,.10137,1.,4.784):
        seed=(RT.DEFAULT_SEED+round(mh*1e6))%1000003
        rng=random.Random(seed)
        nr=np.random.default_rng(seed)
        cand=model.build_candidates(BOnly(),mh,rng,nr,{'B':16000,'inel':2500,'h':10000})
        p4=cand['scal_p4']
        mm2=p4[:,0]**2-np.sum(p4[:,1:]**2,axis=1)
        pt=np.hypot(p4[:,1],p4[:,2])
        eta=np.arcsinh(p4[:,3]/np.maximum(pt,1e-100))
        acc=(eta>2)&(eta<5)
        ratios=[]
        en_ratios=[]
        dzxy=[]
        # old boost of a rest daughter pair vs stable boost using the KNOWN scalar mass
        for v in p4[acc]:
            parts=L.two_body(mh,L.M_E,L.M_E,random.Random(42))
            old=L.boost(parts,v[0],v[1:])
            new=stable_boost(parts,v[0],v[1:],mh)
            en_ratios.append(sum(a[1] for a in old)/v[0])
            ratios.append(old[0][1]/new[0][1])
        allp[str(mh)]={'n_candidate':len(p4),'n_eta_accepted':int(acc.sum()),
            'scalar_m2_ratio_to_nominal':stats(mm2[acc]/mh**2),
            'non_timelike_eta_count':int((mm2[acc]<=0).sum()),
            'ee_energy_sum_over_H_energy':stats(en_ratios),
            'ee_one_daughter_old_over_known_mass_boost':stats(ratios),
            'ee_energy_nonconservation_gt_1pct_count':int((abs(np.array(en_ratios)-1)>.01).sum())}
        # Conditioned response difference at selected points: same candidates and RNG,
        # replace only the H -> daughters boost with a mass-shell-conserving boost.
        chosen=np.flatnonzero(acc)
        chosen=chosen[:min(1500,len(chosen))]
        ao=[];an=[]
        oldboost=L.boost
        for j in chosen:
            pp=tuple(p4[j]); vv=tuple(cand['scal_v0'][j])
            arng=random.Random(20261026+int(j))
            o=L.acceptance_profiles(mh,pp,vv,arng,channels=['ee'],ndecay=30)
            try:
                L.boost=lambda parts,e,p:stable_boost(parts,e,p,mh)
                n=L.acceptance_profiles(mh,pp,vv,random.Random(20261026+int(j)),channels=['ee'],ndecay=30)
            finally:
                L.boost=oldboost
            ao.append(o[0]);an.append(n[0])
        pz=p4[chosen,3]
        lam=cand['hist_lam'][chosen]
        scale=cand['hist_scale'][chosen]
        for tb in (5.,25.,100.,350.):
            ct=model.ctau_m(mh,tb)
            po=L.p_window_vertex(mh,pz,np.array(ao),ct,lam)
            pn=L.p_window_vertex(mh,pz,np.array(an),ct,lam)
            allp[str(mh)]['ee_acceptance_ratio_new_over_old_tb_%g'%tb]=float(np.dot(pn,scale)/np.dot(po,scale)) if np.dot(po,scale)>0 else None
    out['scalar_boost']=allp
    # Inspect cached minbias subsample rather than rereading the 520 MB compressed source.
    for f in (RUN/'work/analysis/cache').glob('sample_minbias*'):
        d=np.load(f)
        q=np.stack([d[k].astype(float) for k in ('px','py','pz')],axis=1)
        eb=d['e'].astype(float)
        mb2=eb*eb-(q*q).sum(axis=1)
        nonpositive=mb2<=0
        out['inel_precision']={'file':f.name,'rows':len(eb),'n_non_timelike':int(nonpositive.sum()),
                              'n_pz_gt100_non_timelike':int((nonpositive&(abs(q[:,2])>100)).sum()),
                              'max_momentum_GeV':float(np.sqrt((q*q).sum(axis=1)).max())}
        break
    (HERE/'kinematics_results.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2),flush=True)


if __name__=='__main__':main()
