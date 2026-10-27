# -*- coding: utf-8 -*-
"""F05 addendum: the two missing compositions - pure pi+pi-2pi0 and the isospin-weighted 4pi
(2/3 2pi+2pi- + 1/3 pi+pi-2pi0) - at the same fixed kinematics as f05_multi_composition.py."""
import sys, os, math, random, json
import numpy as np
RUN=r'E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3'
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(RUN,'work','analysis'))
sys.path.insert(0, os.path.join(RUN,'work','widths'))
import bc5_lib as L
MPI=L.M_PI; MPI0=L.M_PI0

def comp_pipi2pi0(ms,rng):
    if ms<2*MPI+2*MPI0: return []
    d=L.nbody(ms,[MPI,MPI,MPI0,MPI0],rng)
    if d is None: return []
    ch=[1,-1,0,0]
    return [(d[i][0],d[i][1],d[i][2],d[i][3],d[i][4],ch[i]) for i in range(4)]

def comp_all4pi(ms,rng):
    if ms<4*MPI: return []
    d=L.nbody(ms,[MPI]*4,rng)
    if d is None: return []
    ch=[1,1,-1,-1]
    return [(d[i][0],d[i][1],d[i][2],d[i][3],d[i][4],ch[i]) for i in range(4)]

def comp_iso4pi(ms,rng):
    # isoscalar 4pi: 2/3 (2pi+2pi-) + 1/3 (pi+pi-2pi0)
    if rng.random()<2/3.:
        return comp_all4pi(ms,rng)
    return comp_pipi2pi0(ms,rng)

def make_p4(ms,eta,pmod=60.0):
    pt=pmod/math.cosh(eta); pz=pmod*math.tanh(eta)
    return (math.sqrt(pmod*pmod+ms*ms),pt*0.94,pt*0.34,pz)

def accept_profile(ms,p4,v0,fn,N,rng,zgrid=L.ZGRID):
    e,px,py,pz=p4; tot=np.zeros(len(zgrid)); sq=np.zeros(len(zgrid)); z0=v0[2]*1e-3
    for _ in range(N):
        rest=fn(ms,rng)
        if not rest: continue
        lab5=L.boost([(a[0],a[1],a[2],a[3],a[4]) for a in rest],e,(px,py,pz))
        chg=[a[5] for a in rest]
        lab=[(b[0],b[1],b[2],b[3],b[4],c) for b,c in zip(lab5,chg)]
        v=np.array([L._passes_pair(lab,0.0,0.0,zd) for zd in zgrid])
        tot+=v; sq+=v*v
    mean=tot/N; var=np.maximum(sq/N-mean*mean,0)*N/max(N-1,1)
    return mean,np.sqrt(var/N)

def main():
    base=json.load(open(os.path.join(HERE,'f05_results.json')))
    out={}
    for eta in (3.0,4.0):
        for ms in (1.0,1.3,1.6,1.99):
            key='m%.3f_eta%.1f'%(ms,eta)
            p4=make_p4(ms,eta); v0=(0.0,0.0,1500.0)
            br=base[key]['BR']['multi']; Ab=base[key]['A_base']
            cur=base[key]['A_multi_current']; cur_se=base[key]['A_multi_current_se']
            r={}
            for name,fn in (('pi pi 2pi0',comp_pipi2pi0),('isospin 4pi',comp_iso4pi)):
                m,se=accept_profile(ms,p4,v0,fn,6000,random.Random(7005))
                Am=float(np.mean(m)); Ase=float(np.mean(se))
                r[name]={'A':Am,'se':Ase,'delta_R_over_R':float(br*(Am-cur)/Ab),
                         'delta_R_over_R_se':float(br*math.hypot(Ase,cur_se)/Ab)}
            out[key]=r
            print('m=%.2f eta=%.1f Acur=%.4f  '%(ms,eta,cur)
                  +' '.join('%s A=%.4f dR/R=%+.1f%%'%(k,v['A'],100*v['delta_R_over_R']) for k,v in r.items()),flush=True)
    json.dump(out,open(os.path.join(HERE,'f05_addendum.json'),'w'),indent=1)
    print('done')
if __name__=='__main__': main()
