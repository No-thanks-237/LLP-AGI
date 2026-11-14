"""Paired low-mass response experiment; all modifications are confined to memory."""
from pathlib import Path
import json, random, sys
import numpy as np
HERE=Path(__file__).resolve().parent
RUN=HERE.parents[2]
sys.path.insert(0,str(HERE))
from audit_kinematics import stable_boost
sys.path.insert(0,str(RUN/'work/analysis'))
import models as M
import response_table as RT
import event_engine as E
import twohdm_lib as L


def main():
    samples=RT.Samples(verbose=False)
    b=samples.get('B')
    class BOnly:
        n_gen={'B':samples.n_gen['B']}
        def get(self,key):return b if key=='B' else None
    model=M.TwoHDMModel(verbose=False)
    L.load_hadron_lib(str(RUN/'work/magnus/hadron_lib'))
    out={}
    for mh in (.021,.10137):
        seed=(RT.DEFAULT_SEED+round(mh*1e6))%1000003
        cand=model.build_candidates(BOnly(),mh,random.Random(seed),np.random.default_rng(seed),
                                    {'B':16000,'inel':2500,'h':10000})
        tab=RT.load_table(RT.cache_path(model,mh))
        hist=tab['scal_hist']
        inds=np.flatnonzero(tab['hist_group'][hist]==0)
        parents=tab['hist_parent'][hist[inds]]
        p4=cand['scal_p4'][parents]
        assert np.array_equal(p4[:,3],tab['pz'][inds]),'replay does not match published cache'
        oldA=tab['A'].copy(); newA=tab['A'].copy()
        ic=tab['channels'].index('ee')
        oldboost=L.boost
        used=0
        for k,(idx,par,pp) in enumerate(zip(inds,parents,p4)):
            if pp[3]<=0:continue
            eta=np.arcsinh(pp[3]/max(np.hypot(pp[1],pp[2]),1e-30))
            if not 2<eta<5:continue
            vv=tuple(cand['scal_v0'][par])
            seedA=20261026+int(par)
            o=L.acceptance_profiles(mh,tuple(pp),vv,random.Random(seedA),channels=['ee'],ndecay=100)
            try:
                L.boost=lambda parts,e,p:stable_boost(parts,e,p,mh)
                n=L.acceptance_profiles(mh,tuple(pp),vv,random.Random(seedA),channels=['ee'],ndecay=100)
            finally:L.boost=oldboost
            oldA[idx,ic,:]=o[0];newA[idx,ic,:]=n[0];used+=1
        r={'replay_pz_bit_identical':True,'n_B':len(inds),'n_eta_accepted':used,
           'ndecay':100,'scope':'Only H-to-ee boost changed; parent boost/momenta unchanged.'}
        # Branching fractions in this range contain only visible ee, with gamma gamma
        # and closed-channel 1e-40 floors; all other table rows are held fixed.
        for Lfb,tb in ((25,22.011742571604433 if mh==.021 else 32.674),
                       (300,30.025649091092642 if mh==.021 else 44.6)):
            vals={}
            for label,A in (('published',tab['A']),('paired_old',oldA),('paired_known_mass',newA)):
                copytab=dict(tab,A=A)
                v=E.evaluate(copytab,model,mh,tb,Lfb)
                vals[label]={'s':v['s_total'],'ds':v['s_total_se'],'s_B':v['s']['B']}
            vals['signal_ratio_new_over_old']=vals['paired_known_mass']['s']/vals['paired_old']['s']
            vals['boundary_ratio_eighth_root_estimate']=vals['signal_ratio_new_over_old']**(1/8)
            r[str(Lfb)]=vals
        out[str(mh)]=r
        print(mh,json.dumps(r),flush=True)
    (HERE/'low_mass_response_results.json').write_text(json.dumps(out,indent=2),encoding='utf-8')


if __name__=='__main__':main()
