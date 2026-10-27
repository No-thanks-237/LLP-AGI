# -*- coding: utf-8 -*-
"""F13: verify the reference-BR validation model's isospin factors and the running-vs-fixed
mass check in validate_widths.py.  Writes f13_results.json."""
import sys, os, math, random, json
import numpy as np
RUN=r'E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3'
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(RUN,'work','analysis'))
sys.path.insert(0, os.path.join(RUN,'work','widths'))
import bc5_lib as L, scalar_widths as sw, reference_br_model as RB
L.load_hadron_lib(os.path.join(RUN,'work','magnus','hadron_lib'))

def make_p4(ms, eta=3.0, pmod=60.0):
    pt=pmod/math.cosh(eta); pz=pmod*math.tanh(eta)
    return (math.sqrt(pmod*pmod+ms*ms), pt*0.94, pt*0.34, pz)

def accept_chan(ms,p4,chan,N=4000,seed=3):
    rng=random.Random(seed); e,px,py,pz=p4; z0=1.5; tot=0
    for _ in range(N):
        rest=L.s_decay_products(ms,chan,rng)
        if not rest: continue
        lab5=L.boost([(a[0],a[1],a[2],a[3],a[4]) for a in rest],e,(px,py,pz))
        chg=[a[5] for a in rest]
        lab=[(b[0],b[1],b[2],b[3],b[4],c) for b,c in zip(lab5,chg)]
        tot+=L._passes_pair(lab,0.0,0.0,z0)
    return tot/N, math.sqrt(max(tot/N*(1-tot/N),0)/N)

def main():
    res={}
    # ---- 1. reference table vs production isoscalar widths -----------------------------
    lm,data=RB._load()
    rows=[]
    for ms in (0.5,0.8,1.0,1.3,1.6,1.99,2.5,3.0,4.0,5.0):
        br=RB.branching(ms)
        if br['mu']>1e-12: gtot=sw.gamma_ll(ms,sw.M_MU)/br['mu']
        elif br['e']>1e-12: gtot=sw.gamma_ll(ms,sw.M_E)/br['e']
        else: continue
        w=sw.widths(ms)
        r_pipi=(w.get('pipi',0.0)/gtot/br['pipi']) if br['pipi']>0 else float('nan')
        r_KK=(w.get('KK',0.0)/gtot/br['KK']) if br['KK']>0 else float('nan')
        rows.append({'m':ms,'BR_pipi':br['pipi'],'BR_KK':br['KK'],
                     'G_pipi_ours_over_Gref_pipi':r_pipi,'G_KK_ours_over_Gref_KK':r_KK})
    res['isospin_ratios']=rows
    # ---- 2. what reference_br_model actually assigns -----------------------------------
    res['ref_model']=[]
    for ms in (1.0,1.3,1.6,1.99):
        w=RB.widths(ms); tot=sum(w.values())
        br=RB.branching(ms)
        res['ref_model'].append({'m':ms,'Gamma_tot':tot,
            'frac':{k:(w[k]/tot if tot>0 else 0) for k in ('pipi','KK','multi','e','mu','tau')},
            'BRref_pipi':br['pipi'],'BRref_KK':br['KK'],
            'BRref_4pi':br['multi']})
    # ---- 3. visible-rate bias of the reference model at fixed parent kinematics --------
    res['visible_bias']=[]
    for ms in (1.0,1.3,1.6,1.99):
        p4=make_p4(ms)
        A={}; Ae={}
        for c in ('e','mu','tau','pipi','KK','multi'):
            A[c],Ae[c]=accept_chan(ms,p4,c)
        w=RB.widths(ms); tot=sum(w.values())
        # current validation model: pipi/KK widths are the single-charge-state ref values
        # visible = sum_c Gamma_c * A_c  (Gamma from the validation model)
        vis_cur=sum(w[c]*A[c] for c in A)
        # corrected: enforce the isospin relation between the reference single-charge-state
        # columns and our isoscalar widths, then take the leftover
        gtot=tot
        br=RB.branching(ms)
        gp=1.5*br['pipi']*gtot
        gk=2.0*br['KK']*gtot
        # the validation model's leftover (assigned to multi below 2 GeV)
        rest_raw=1.0-sum(br.values())   # genuinely missing below 2 GeV (pi0pi0, other hadrons)
        multi_raw=gtot*br['multi']      # the reference 2pi+2pi- column alone
        rest_new=rest_raw-0.5*br['pipi']-1.0*br['KK']
        w2=dict(w); w2['pipi']=gp; w2['KK']=gk
        w2['multi']=max(multi_raw+rest_new*gtot,0.0)
        res.setdefault('rest_check',[]).append({'m':ms,'rest_raw':rest_raw,'rest_new':rest_new,'multi_old':w['multi'],'multi_new':w2['multi']})
        vis_new=sum(w2[c]*A[c] for c in A)
        res['visible_bias'].append({'m':ms,'visible_current':vis_cur,'visible_corrected':vis_new,
            'ratio_cur_over_new':vis_cur/vis_new if vis_new>0 else float('nan'),
            'A':A,'A_se':Ae,'W_cur':{k:w[k]/tot for k in A},'W_new':{k:w2[k]/tot for k in A}})
    # ---- 4. validate_widths running vs fixed mass --------------------------------------
    res['running_vs_fixed']=[]
    for ms in (2.2,2.5,3.0,4.0,5.0,7.0):
        run_s=sw.gamma_qq(ms,sw.MSB_S_2GEV,2.0,sw.M_K)
        run_c=sw.gamma_qq(ms,sw.MSB_C_2GEV,2.0,sw.M_D)
        f=sw.widths(ms)
        res['running_vs_fixed'].append({'m':ms,'ss_running':run_s,'ss_fixed':f.get('ss',0),
            'ratio_ss':run_s/f['ss'] if f.get('ss',0)>0 else float('nan'),
            'cc_running':run_c,'cc_fixed':f.get('cc',0),
            'ratio_cc':run_c/f['cc'] if f.get('cc',0)>0 else float('nan')})
    with open(os.path.join(HERE,'f13_results.json'),'w') as fh: json.dump(res,fh,indent=1)
    print(json.dumps(res['isospin_ratios'],indent=1))
    for r in res['visible_bias']:
        print('m=%.2f  visible_cur/visible_corrected = %.4f  A_pipi=%.3f A_KK=%.3f A_multi=%.3f'%(
            r['m'],r['ratio_cur_over_new'],r['A']['pipi'],r['A']['KK'],r['A']['multi']))
    for r in res['running_vs_fixed']:
        print('m=%.1f  ss run/fixed=%.3f  cc run/fixed=%.3f'%(r['m'],r['ratio_ss'],r['ratio_cc']))
if __name__=='__main__': main()
