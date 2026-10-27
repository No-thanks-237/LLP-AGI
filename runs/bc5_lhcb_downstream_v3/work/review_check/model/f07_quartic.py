# -*- coding: utf-8 -*-
"""F07: is the quartic (X->SS / X->X'SS) rate integral consistent with the generated
kinematics, and how big is the quartic contribution to the total signal at m_S=1,2,4 GeV?
Writes f07_results.json."""
import sys, os, math, random, json, time
import numpy as np
RUN=r'E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3'
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(RUN,'work','analysis'))
sys.path.insert(0, os.path.join(RUN,'work','widths'))
import bc5_lib as L, bc5_production as P, run_scan as R
import scalar_widths as sw
L.load_hadron_lib(os.path.join(RUN,'work','magnus','hadron_lib'))

def subset(tab,chan,hid):
    sel=(tab['chan']==chan)&(tab['hid']==hid)
    out={k:(v[sel] if isinstance(v,np.ndarray) and len(v)==len(tab['chan']) else v) for k,v in tab.items()}
    out['chan']=tab['chan'][sel]
    return out

def main():
    sc=json.load(open(os.path.join(RUN,'outputs','data','signal_scan.json')))
    masses=np.array(sc['masses']); th=np.array(sc['theta2']); sig=sc['signal_per_fb_per_channel']
    def crossings(x,y,level):
        out=[]; lx=np.log10(x); ly=np.log10(np.maximum(y,1e-300))
        for i in range(len(x)-1):
            f0=ly[i]-math.log10(level); f1=ly[i+1]-math.log10(level)
            if f0*f1<0: out.append(10**(lx[i]+f0/(f0-f1)*(lx[i+1]-lx[i])))
        return out
    # B_YIELD_CORR from the sample (same rule as run_scan.main)
    smpB=R.load_sample(R.SAMPLE_B)
    sel=np.isin(smpB['pdg'],[521,-521,511,-511,531,-531])
    per_ev=float(sel.sum())/float(R.N_GEN_BBAR)
    R.B_YIELD_CORR=R.REF_B_MESON_YIELD/per_ev
    print('B_YIELD_CORR %.5f (b mesons/event %.4f)'%(R.B_YIELD_CORR,per_ev),flush=True)
    # parent composition (rows) used by the estimator
    ids,cts=np.unique(smpB['pdg'][sel],return_counts=True)
    comp={int(k):int(v) for k,v in zip(ids,cts)}
    print('B parent rows',comp,flush=True)
    res={'settings':{'nsamp':1200,'ndecay':8,'B_YIELD_CORR':R.B_YIELD_CORR,'parent_rows':comp},'rows':{}}
    for ms in (1.0,2.0,4.0):
        t0=time.time()
        R._np_rng=np.random.default_rng(999)
        tab=R.run_mass(ms,1200,random.Random(555),ndecay=8,verbose=True)
        i=np.argmin(abs(masses-ms)); mrow=sig[str(masses[i])]
        s=np.array([r['B']+r['K']+r['H'] for r in mrow])
        cr=crossings(th,s,2.3/25.0); t2_lo=min(cr) if cr else 1e-9
        cr3=crossings(th,s,2.3/300.0); t2_lo3=min(cr3) if cr3 else 1e-9
        out={'m':ms,'scan_mass':float(masses[i]),'theta2_lo_25fb':t2_lo,'theta2_lo_300fb':t2_lo3,'t2':{},'brs':{}}
        # BR decomposition at this mass
        brB={bm:sum(P.br_B_to_S(ms,bm).values()) for bm in ('B+','B0','Bs')}
        qB={bm:sum(q[3] for q in P.quartic_channels(ms) if q[0]==bm) for bm in ('B+','B0','Bs')}
        brK=P.br_K_to_S(ms,'K+')+P.br_K_to_S(ms,'K0L')
        qK=sum(q[2] for q in P.quartic_kaon_channels(ms))
        out['brs']={'mixB':brB,'quarticB':qB,'mixK':brK,'quarticK':qK}
        # |M(q^2)|^2 variation for the 3-body channels: ratio of M^2 at the q^2 endpoints
        mm={}
        for (mes,kind,mfin,br0,mclose,chan) in P.QUARTIC:
            if kind!='3b': continue
            mX={'B+':P.M_BP,'B0':P.M_B0,'Bs':P.M_BS}[mes]
            lo=4*ms*ms; hi=(mX-mfin)**2
            if hi<=lo: continue
            v=[]
            for q2 in (lo+1e-9*(hi-lo), hi-1e-9*(hi-lo)):
                m=P.matrix_element(chan,math.sqrt(q2),mX,mfin); v.append(m*m)
            mm[mes+'_'+chan]={'M2_lo':v[0],'M2_hi':v[1],'ratio_hi_over_lo':(v[1]/v[0] if v[0]>0 else float('nan'))}
        out['M2_variation']=mm
        # contributions
        for t2name,t2 in (('lo_25fb',t2_lo),('lo_300fb',t2_lo3),('1e-10',1e-10),('1e-8',1e-8),('1e-6',1e-6)):
            p=R.per_scalar_p(ms,tab,t2,tab['pz'],tab['z'],tab['A'])
            d={}
            for chan in ('B','K','H'):
                e,se=R.channel_efficiency(tab,p,chan,R.N_GEN_BBAR if chan=='B' else (R.N_GEN_MINBIAS if chan=='K' else R.N_GEN_HIGGS),t2)
                d[chan]=e
            emix={}; eq={}
            for chan,ng in (('B',R.N_GEN_BBAR),('K',R.N_GEN_MINBIAS)):
                for hid,store in ((0,emix),(1,eq)):
                    sub=subset(tab,chan,hid)
                    if len(sub['chan'])==0:
                        store[chan]=0.0; continue
                    psub=R.per_scalar_p(ms,sub,t2,sub['pz'],sub['z'],sub['A'])
                    ev,_=R.channel_efficiency(sub,psub,chan,ng,t2)
                    store[chan]=ev
            sb=1e3*25.0*R.SIGMA_BB_PB*d['B']*R.B_YIELD_CORR
            sk=1e3*25.0*R.SIGMA_INEL_PB*d['K']
            sh=1e3*25.0*P.SIGMA_H_TOTAL*P.B_H*d['H']
            qb=1e3*25.0*R.SIGMA_BB_PB*eq['B']*R.B_YIELD_CORR
            qk=1e3*25.0*R.SIGMA_INEL_PB*eq['K']
            tot=sb+sk+sh
            out['t2'][t2name]={'theta2':t2,'eps':d,'eps_mix':emix,'eps_quartic':eq,
                               's25_B':sb,'s25_K':sk,'s25_H':sh,'s25_total':tot,
                               's25_quartic':qb+qk,'quartic_share':(qb+qk)/tot if tot>0 else 0.0,
                               'quartic_share_B_channel':qb/sb if sb>0 else 0.0}
            print('m=%.2f t2=%.2e  s25=%.4g  quartic=%.3g share=%.3f  Bshare=%.4f (%.0fs)'%(
                ms,t2,tot,qb+qk,(qb+qk)/tot,qb/sb if sb>0 else 0,time.time()-t0),flush=True)
        res['rows']['m%.1f'%ms]=out
    with open(os.path.join(HERE,'f07_results.json'),'w') as f: json.dump(res,f,indent=1)
    print('done')
if __name__=='__main__': main()
