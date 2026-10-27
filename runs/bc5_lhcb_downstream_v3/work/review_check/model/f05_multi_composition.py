# -*- coding: utf-8 -*-
"""F05: visible-rate effect of the arbitrary 'multi' composition split at fixed parent
kinematics, with the real transport kernel.  Writes f05_results.json."""
import sys, os, time, math, random, json
import numpy as np
RUN=r'E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3'
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(RUN,'work','analysis'))
sys.path.insert(0, os.path.join(RUN,'work','widths'))
import bc5_lib as L
MPI=L.M_PI; MPI0=L.M_PI0; META=L.M_ETA

def comp_current(ms, rng):
    x=rng.random()
    if x<1/3.: masses,charges=[MPI]*4,[1,1,-1,-1]
    elif x<2/3.: masses,charges=[MPI,MPI,MPI0,MPI0],[1,-1,0,0]
    else: masses,charges=[META,META,MPI,MPI],[0,0,1,-1]
    if ms<sum(masses):
        masses,charges=[MPI]*4,[1,1,-1,-1]
        if ms<sum(masses): return []
    d=L.nbody(ms,masses,rng)
    if d is None: return []
    return [(d[i][0],d[i][1],d[i][2],d[i][3],d[i][4],charges[i]) for i in range(len(d))]

def comp_all4pi(ms, rng):
    if ms < 4*MPI: return []
    d=L.nbody(ms,[MPI]*4,rng)
    if d is None: return []
    ch=[1,1,-1,-1]
    return [(d[i][0],d[i][1],d[i][2],d[i][3],d[i][4],ch[i]) for i in range(4)]

def comp_rho(ms, rng):
    if ms < 2*0.775: return []
    d=L.two_body(ms,0.775,0.775,rng)
    if d is None: return []
    out=[]
    for i,sg in ((0,1),(1,-1)):
        dd=L.two_body(0.775,MPI,MPI,rng)
        if dd is None: return []
        lab=L.boost([(a[0],a[1],a[2],a[3],a[4]) for a in dd],d[i][1],d[i][2:5])
        out.append((lab[0][0],lab[0][1],lab[0][2],lab[0][3],lab[0][4],sg))
        out.append((lab[1][0],lab[1][1],lab[1][2],lab[1][3],lab[1][4],-sg))
    return out

def _eta_daughters(rng):
    r=rng.random()
    if r<0.394: return []          # gamma gamma
    if r<0.721: return []          # 3 pi0
    if r<0.950:                    # pi+ pi- pi0 : two charged
        d=L.nbody(META,[MPI,MPI,MPI0],rng)
        if d is None: return []
        return [(d[0][0],d[0][1],d[0][2],d[0][3],d[0][4],1),
                (d[1][0],d[1][1],d[1][2],d[1][3],d[1][4],-1)]
    return []

def comp_eta_decays(ms, rng):
    x=rng.random()
    if x<1/3.: return comp_all4pi(ms,rng)
    if x<2/3.:
        if ms < 2*MPI+2*MPI0: return comp_all4pi(ms,rng)
        d=L.nbody(ms,[MPI,MPI,MPI0,MPI0],rng)
        if d is None: return []
        ch=[1,-1,0,0]
        return [(d[i][0],d[i][1],d[i][2],d[i][3],d[i][4],ch[i]) for i in range(4)]
    if ms < 2*META+2*MPI: return comp_all4pi(ms,rng)
    d=L.nbody(ms,[META,META,MPI,MPI],rng)
    if d is None: return []
    out=[]
    for i in range(2):
        sub=_eta_daughters(rng)
        if sub:
            lab=L.boost([(a[0],a[1],a[2],a[3],a[4]) for a in sub],d[i][1],(d[i][2],d[i][3],d[i][4]))
            out.append((lab[0][0],lab[0][1],lab[0][2],lab[0][3],lab[0][4],1))
            out.append((lab[1][0],lab[1][1],lab[1][2],lab[1][3],lab[1][4],-1))
    out.append((d[2][0],d[2][1],d[2][2],d[2][3],d[2][4],1))
    out.append((d[3][0],d[3][1],d[3][2],d[3][3],d[3][4],-1))
    return out

COMPS=[('current 3-way (eta stable)',comp_current),
       ('all 4 charged pions',comp_all4pi),
       ('rho rho -> 4 charged pions',comp_rho),
       ('current 3-way + eta decays',comp_eta_decays)]

def make_p4(ms, eta, pmod=60.0):
    pt=pmod/math.cosh(eta); pz=pmod*math.tanh(eta)
    e=math.sqrt(pmod*pmod+ms*ms)
    return (e, pt*0.94, pt*0.34, pz)

def accept_profile(ms,p4,v0_mm,decay_fn,N,rng,zgrid=L.ZGRID):
    e,px,py,pz=p4
    tot=np.zeros(len(zgrid)); sq=np.zeros(len(zgrid)); z0=v0_mm[2]*1e-3
    dirx,diry=px/pz,py/pz
    for _ in range(N):
        rest=decay_fn(ms,rng)
        if not rest: continue
        lab5=L.boost([(a[0],a[1],a[2],a[3],a[4]) for a in rest],e,(px,py,pz))
        chg=[a[5] for a in rest]
        lab=[(b[0],b[1],b[2],b[3],b[4],c) for b,c in zip(lab5,chg)]
        v=np.array([L._passes_pair(lab,v0_mm[0]*1e-3+(zd-z0)*dirx,
                                   v0_mm[1]*1e-3+(zd-z0)*diry,zd) for zd in zgrid])
        tot+=v; sq+=v*v
    mean=tot/N
    var=np.maximum(sq/N-mean*mean,0.0)*N/max(N-1,1)
    return mean, np.sqrt(var/N)

def main():
    res={}; t0=time.time()
    for eta in (3.0,4.0):
        for ms in (1.0,1.3,1.6,1.99):
            key='m%.3f_eta%.1f'%(ms,eta)
            p4=make_p4(ms,eta); v0=(0.0,0.0,1500.0)
            br,tot=L.channel_branching(ms)
            def _base_decay(m,r,br=br):
                names=list(br.keys()); cum=np.cumsum([br[n] for n in names])
                u=r.random(); idx=min(int(np.searchsorted(cum,u)),len(names)-1)
                return L.s_decay_products(m,names[idx],r)
            base,base_se=accept_profile(ms,p4,v0,_base_decay,20000,random.Random(7001))
            A_base=float(np.mean(base)); A_base_se=float(np.mean(base_se))
            out={'A_base':A_base,'A_base_se':A_base_se,'BR_multi':float(br['multi']),
                 'BR':{c:float(v) for c,v in br.items()}}
            cur,cur_se=accept_profile(ms,p4,v0,comp_current,6000,random.Random(7002))
            out['A_multi_current']=float(np.mean(cur)); out['A_multi_current_se']=float(np.mean(cur_se))
            for name,fn in COMPS:
                m,se=accept_profile(ms,p4,v0,fn,6000,random.Random(7003))
                Am=float(np.mean(m)); Amse=float(np.mean(se))
                dA=Am-out['A_multi_current']
                dse=math.hypot(Amse,out['A_multi_current_se'])
                out[name]={'A_multi':Am,'se':Amse,'dA':dA,'dA_se':dse,
                           'delta_R_over_R':float(br['multi']*dA/A_base),
                           'delta_R_over_R_se':float(br['multi']*dse/A_base)}
            res[key]=out
            print('m=%.2f eta=%.1f BRmulti=%.3f  '%(ms,eta,br['multi'])
                  +' | '.join('%s A=%.4f+-%.4f dR/R=%+.1f%%'%(n[:14],out[n]['A_multi'],out[n]['se'],100*out[n]['delta_R_over_R']) for n,_ in COMPS)
                  ,flush=True)
    with open(os.path.join(HERE,'f05_results.json'),'w') as f: json.dump(res,f,indent=1)
    print('done %.1fs'%(time.time()-t0))
if __name__=='__main__': main()
