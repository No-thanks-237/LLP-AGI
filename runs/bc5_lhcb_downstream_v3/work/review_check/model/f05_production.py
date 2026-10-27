# -*- coding: utf-8 -*-
"""F05 production cross-check: run the real B/K/H parent pipeline at m_S=1.6 and 1.99 with
three 'multi' compositions and compare the expected signal at the 25 fb crossing theta^2."""
import sys, os, math, random, json, time
import numpy as np
RUN=r'E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3'
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(RUN,'work','analysis'))
sys.path.insert(0, os.path.join(RUN,'work','widths'))
import bc5_lib as L, bc5_production as P, run_scan as R
import scalar_widths as sw
L.load_hadron_lib(os.path.join(RUN,'work','magnus','hadron_lib'))
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
    if ms<4*MPI: return []
    d=L.nbody(ms,[MPI]*4,rng)
    if d is None: return []
    ch=[1,1,-1,-1]
    return [(d[i][0],d[i][1],d[i][2],d[i][3],d[i][4],ch[i]) for i in range(4)]

def comp_pipi2pi0(ms,rng):
    if ms<2*MPI+2*MPI0: return []
    d=L.nbody(ms,[MPI,MPI,MPI0,MPI0],rng)
    if d is None: return []
    ch=[1,-1,0,0]
    return [(d[i][0],d[i][1],d[i][2],d[i][3],d[i][4],ch[i]) for i in range(4)]

def comp_iso4pi(ms,rng):
    if rng.random()<2/3.: return comp_all4pi(ms,rng)
    return comp_pipi2pi0(ms,rng)

COMPS={'current':comp_current,'all4pi':comp_all4pi,'iso4pi':comp_iso4pi}
_ORIG=L.s_decay_products
def patch(fn):
    def wrapped(ms,ch,rng):
        if ch=='multi': return fn(ms,rng)
        return _ORIG(ms,ch,rng)
    L.s_decay_products=wrapped

def main():
    sc=json.load(open(os.path.join(RUN,'outputs','data','signal_scan.json')))
    masses=np.array(sc['masses']); th=np.array(sc['theta2']); sig=sc['signal_per_fb_per_channel']
    def crossings(x,y,level):
        out=[]; lx=np.log10(x); ly=np.log10(np.maximum(y,1e-300))
        for i in range(len(x)-1):
            f0=ly[i]-math.log10(level); f1=ly[i+1]-math.log10(level)
            if f0*f1<0: out.append(10**(lx[i]+f0/(f0-f1)*(lx[i+1]-lx[i])))
        return out
    res={'settings':{'nsamp':2000,'ndecay':8},'rows':{}}
    for ms in (1.6,1.99):
        i=np.argmin(abs(masses-ms)); s=np.array([r['B']+r['K']+r['H'] for r in sig[str(masses[i])]])
        cr=crossings(th,s,2.3/25.0); t2=min(cr)
        out={'m':ms,'scan_mass':float(masses[i]),'theta2_25fb':t2,'comps':{}}
        for name,fn in COMPS.items():
            patch(fn)
            R._np_rng=np.random.default_rng(31337)
            t0=time.time()
            tab=R.run_mass(ms,2000,random.Random(777),ndecay=8,verbose=False)
            p=R.per_scalar_p(ms,tab,t2,tab['pz'],tab['z'],tab['A'])
            d={}; se={}
            for chan,ng in (('B',R.N_GEN_BBAR),('K',R.N_GEN_MINBIAS),('H',R.N_GEN_HIGGS)):
                e,ee=R.channel_efficiency(tab,p,chan,ng,t2); d[chan]=e; se[chan]=ee
            s25=1e3*25.0*(R.SIGMA_BB_PB*d['B']+R.SIGMA_INEL_PB*d['K']+P.SIGMA_H_TOTAL*P.B_H*d['H'])
            s25_se=1e3*25.0*math.sqrt((R.SIGMA_BB_PB*se['B'])**2+(R.SIGMA_INEL_PB*se['K'])**2+(P.SIGMA_H_TOTAL*P.B_H*se['H'])**2)
            out['comps'][name]={'eps':d,'s25':s25,'s25_se':s25_se,'time':time.time()-t0}
            print('m=%.2f %-8s s25=%.4f+-%.4f (B=%.3f K=%.3f H=%.3f) %.0fs'%(
                ms,name,s25,s25_se,d['B'],d['K'],d['H'],time.time()-t0),flush=True)
        ref=out['comps']['current']['s25']
        for name in COMPS:
            out['comps'][name]['delta_R_over_R']=out['comps'][name]['s25']/ref-1.0
        res['rows']['m%.2f'%ms]=out
        L.s_decay_products=_ORIG
    json.dump(res,open(os.path.join(HERE,'f05_production.json'),'w'),indent=1)
    print('done')
if __name__=='__main__': main()
