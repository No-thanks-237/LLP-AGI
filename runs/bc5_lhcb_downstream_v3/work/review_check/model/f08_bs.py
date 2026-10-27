# -*- coding: utf-8 -*-
"""F08: B_s-channel share of the signal vs mass, the visibility of the 4.347 GeV cutoff in
the crossing curve, and the continuous B_s -> phi S (K* form factor) estimate.
Writes f08_results.json."""
import sys, os, math, json
import numpy as np
RUN=r'E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3'
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(RUN,'work','analysis'))
sys.path.insert(0, os.path.join(RUN,'work','widths'))
import bc5_production as P
import run_scan as R

M_PHI=1.019461
M_ETA=0.547862
M_ETAP=0.95778

def bs_phi_S(ms):
    """BR(B_s -> phi S)/theta^2 with Winkler Eq.(BKphi) and the K* form factor, m_B->m_Bs."""
    mX=P.M_BS; mXp=M_PHI
    p=P.pS(mX,mXp,ms)
    if p<=0: return 0.0,0.0
    q2=ms*ms
    M=abs(-2.0*mX*p/(P.M_B_QUARK+P.M_S_QUARK))*P.ff_vec892(q2)
    br=P.XI_SB**2*P.M_B_QUARK**2*M**2/(32.0*math.pi*P.V**2*P.GAM_BS)*p/mX**2
    return br,M

def main():
    res={}
    # parent composition of the B estimator (rows, as used by run_scan)
    smp=R.load_sample(R.SAMPLE_B)
    sel=np.isin(smp['pdg'],[521,-521,511,-511,531,-531])
    ids,cts=np.unique(smp['pdg'][sel],return_counts=True)
    comp={int(k):int(v) for k,v in zip(ids,cts)}
    nBp=comp.get(521,0)+comp.get(-521,0); nB0=comp.get(511,0)+comp.get(-511,0); nBs=comp.get(531,0)+comp.get(-531,0)
    res['parent_rows']={'B+':nBp,'B0':nB0,'Bs':nBs,'Bs_fraction_of_B_rows':nBs/(nBp+nB0+nBs)}
    # scan signals for the channel fractions
    sc=json.load(open(os.path.join(RUN,'outputs','data','signal_scan.json')))
    masses=np.array(sc['masses']); th=np.array(sc['theta2']); sig=sc['signal_per_fb_per_channel']
    def crossing_lo(s,L):
        out=[]; lx=np.log10(th); ly=np.log10(np.maximum(s,1e-300)); lev=2.3/L
        for i in range(len(th)-1):
            f0=ly[i]-math.log10(lev); f1=ly[i+1]-math.log10(lev)
            if f0*f1<0: out.append(10**(lx[i]+f0/(f0-f1)*(lx[i+1]-lx[i])))
        return min(out) if out else float('nan')
    def row_at(m):
        i=np.argmin(abs(masses-m)); return i
    res['mass_table']=[]
    for ms in (3.5,4.0,4.3,4.4,5.0,6.0,8.0):
        i=row_at(ms); m=masses[i]; row=sig[str(m)]
        s=np.array([r['B']+r['K']+r['H'] for r in row])
        thlo=crossing_lo(s,25.0)
        j=np.argmin(abs(np.log(th)-np.log(thlo)))
        r=row[j]; tot=r['B']+r['K']+r['H']
        Bfrac=r['B']/tot if tot>0 else 0.0
        # B_s share of the B channel at this mass
        brB={bm:sum(P.br_B_to_S(ms,bm).values()) for bm in ('B+','B0','Bs')}
        n={'B+':nBp,'B0':nB0,'Bs':nBs}
        num=n['Bs']*brB['Bs']; den=sum(n[k]*brB[k] for k in brB)
        xBs=num/den if den>0 else 0.0
        res['mass_table'].append({'m':ms,'scan_mass':float(m),'theta2_lo_25fb':thlo,
            'B_frac_at_crossing':Bfrac,'brB':brB,'Bs_share_of_B':xBs,
            'Bs_share_of_total':xBs*Bfrac,
            'bs_phi_S':bs_phi_S(ms)[0],
            'bs_inclusive':list(P.br_Bs_inclusive(ms).values())[0] if P.br_Bs_inclusive(ms) else 0.0})
    # cliff in the crossing curve: at the scan point nearest 4.3474, remove B_s from B
    row=res['mass_table'][2]   # m=4.3
    i=row_at(4.3446); m=masses[i]; rr=sig[str(m)]
    s=np.array([r['B']+r['K']+r['H'] for r in rr])
    brB={bm:sum(P.br_B_to_S(m,bm).values()) for bm in ('B+','B0','Bs')}
    n={'B+':nBp,'B0':nB0,'Bs':nBs}; den=sum(n[k]*brB[k] for k in brB)
    xBs=n['Bs']*brB['Bs']/den if den>0 else 0.0
    base25=crossing_lo(s,25.0); base300=crossing_lo(s,300.0)
    s_nobs=s*(1-xBs*(np.array([r['B'] for r in rr])/np.maximum(s,1e-300)))
    res['cliff']={'scan_mass':float(m),'theta2_lo_25fb':base25,'theta2_lo_300fb':base300,
                  'xBs':xBs,'B_frac_at_crossing_25fb':float(np.interp(math.log(base25),np.log(th),np.array([r['B'] for r in rr])/np.maximum(s,1e-300))),
                  'theta2_lo_25fb_without_Bs':crossing_lo(s_nobs,25.0),
                  'theta2_lo_300fb_without_Bs':crossing_lo(s_nobs,300.0)}
    res['cliff']['rel_shift_25fb']=res['cliff']['theta2_lo_25fb_without_Bs']/base25-1
    res['cliff']['rel_shift_300fb']=res['cliff']['theta2_lo_300fb_without_Bs']/base300-1
    # the continuous phi-S curve and the endpoint behaviour
    res['phi_curve']=[]
    for ms in (0.5,1.0,2.0,3.0,4.0,4.3,4.34,4.347,4.3474,4.3475,4.4,4.5):
        b_phi,_=bs_phi_S(ms)
        binc=list(P.br_Bs_inclusive(ms).values())[0] if P.br_Bs_inclusive(ms) else 0.0
        res['phi_curve'].append({'m':ms,'BR_bs_phi_S':b_phi,'BR_bs_inclusive':binc,
                                 'phi_over_inclusive':(b_phi/binc if binc>0 else float('nan'))})
    # thresholds of the remaining channels
    res['thresholds']={'Bs_to_phiS':P.M_BS-M_PHI,'Bs_to_etaS':P.M_BS-M_ETA,
                       'Bs_to_etapS':P.M_BS-M_ETAP,
                       'Bs_to_f0S_approx':P.M_BS-0.990}
    with open(os.path.join(HERE,'f08_results.json'),'w') as f: json.dump(res,f,indent=1)
    print(json.dumps(res['parent_rows'],indent=1))
    for r in res['mass_table']:
        print('m=%.2f  Bfrac=%.4f  Bs/B=%.3f  Bs/total=%.4f  BR_phiS=%.3e  BR_incl=%.3e'%(
            r['m'],r['B_frac_at_crossing'],r['Bs_share_of_B'],r['Bs_share_of_total'],
            r['bs_phi_S'],r['bs_inclusive']))
    print('cliff',json.dumps(res['cliff'],indent=1))
    for r in res['phi_curve']:
        print('phi: m=%.4f BR_phi=%.4e BR_incl=%.4e ratio=%.4f'%(r['m'],r['BR_bs_phi_S'],r['BR_bs_inclusive'],r['phi_over_inclusive']))
if __name__=='__main__': main()
