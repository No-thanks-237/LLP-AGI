# -*- coding: utf-8 -*-
"""F06: effect of using the first vs the last Pythia event-record copy of the Higgs on the
Higgs-channel efficiency, and of using the sample invariant mass instead of P.M_H.
3000 parents, ndecay=8, real transport kernel.  Writes f06_results.json."""
import sys, os, math, random, json, time
import numpy as np
RUN=r'E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3'
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(RUN,'work','analysis'))
sys.path.insert(0, os.path.join(RUN,'work','widths'))
import bc5_lib as L, bc5_production as P, run_scan as R
import scalar_widths as sw
L.load_hadron_lib(os.path.join(RUN,'work','magnus','hadron_lib'))

NP=3000; NDECAY=8

def load_higgs():
    d=np.loadtxt(os.path.join(RUN,'work','magnus','samples','higgs.txt'),comments='#')
    ev=d[:,0].astype(int)
    u,first,counts=np.unique(ev,return_index=True,return_counts=True)
    last=first+counts-1
    def smp(rows):
        a=d[rows]
        return {'ev':a[:,0].astype(np.int64),'pdg':a[:,1].astype(np.int32),
                'px':a[:,2],'py':a[:,3],'pz':a[:,4],'e':a[:,5],
                'x':a[:,6],'y':a[:,7],'z':a[:,8],'t':a[:,9]}
    return smp(first), smp(last)

def build_tab(ms, sample, idxs, m_eff_mode, seed=12345):
    """Replicates run_scan.run_mass's Higgs block for a fixed parent index list."""
    rng=random.Random(seed)
    cols={k:[] for k in ('ev','e','px','py','pz','x','y','z','A','w','par','hid','kq')}
    n=len(idxs)
    for ip,i in enumerate(idxs):
        e=sample['e'][i]; px=sample['px'][i]; py=sample['py'][i]; pz=sample['pz'][i]
        if m_eff_mode=='inv':
            m_eff=math.sqrt(max(e*e-px*px-py*py-pz*pz,0.0))
        else:
            m_eff=P.M_H
        d=L.two_body(m_eff,ms,ms,rng)
        if d is None: continue
        pair=L.boost(d,e,(px,py,pz))
        for s4 in pair:
            r=R.decay_vertex_and_selection(ms,(s4[1],s4[2],s4[3],s4[4],sample['x'][i],sample['y'][i],sample['z'][i]),rng,ndecay=NDECAY)
            if r is None: continue
            p4,v0,a=r
            cols['ev'].append(i); cols['e'].append(p4[0]); cols['px'].append(p4[1]); cols['py'].append(p4[2]); cols['pz'].append(p4[3])
            cols['x'].append(v0[0]*1e-3); cols['y'].append(v0[1]*1e-3); cols['z'].append(v0[2]*1e-3)
            cols['A'].append(a); cols['w'].append(1.0); cols['par'].append(ip); cols['hid'].append(0); cols['kq'].append(0.0)
    tab={k:np.asarray(v) for k,v in cols.items()}
    tab['chan']=np.array(['H']*len(tab['ev']))
    tab['counts']={'H':[199105,n]}
    tab['samp_ev']={'H':np.asarray(idxs,dtype=np.int64)}
    return tab

def channel_efficiency_events(tab,p,chan,n_gen=200000):
    sel=tab['chan']==chan
    if not np.any(sel) or n_gen<=0: return 0.0,0.0,np.zeros(0),np.zeros(0,dtype=np.int64)
    c=tab['counts'][chan]; n_tot,n_smp=int(c[0]),int(c[1])
    if n_smp<=0: return 0.0,0.0,np.zeros(0),np.zeros(0,dtype=np.int64)
    ev=tab['ev'][sel].astype(np.int64); par=tab['par'][sel].astype(np.int64); hid=tab['hid'][sel].astype(np.int64)
    kq=tab['kq'][sel]; br=tab['w'][sel]
    pe=p[sel]*np.where(kq>0.5,1.0,1.0)
    log1mp=np.log1p(-np.clip(pe,0.0,1.0-1e-15))
    order=np.lexsort((hid,par,ev))
    ev_s,par_s,hid_s,log_s=ev[order],par[order],hid[order],log1mp[order]
    br_s=br[order]
    new_h=np.ones(len(ev_s),dtype=bool)
    if len(ev_s)>1:
        new_h[1:]=(ev_s[1:]!=ev_s[:-1])|(par_s[1:]!=par_s[:-1])|(hid_s[1:]!=hid_s[:-1])
    start_h=np.flatnonzero(new_h)
    o_inner=-np.expm1(np.add.reduceat(log_s,start_h))
    o_hist=br_s[start_h]*o_inner
    ev_h,par_h=ev_s[start_h],par_s[start_h]
    new_p=np.ones(len(o_hist),dtype=bool)
    if len(o_hist)>1:
        new_p[1:]=(ev_h[1:]!=ev_h[:-1])|(par_h[1:]!=par_h[:-1])
    start_p=np.flatnonzero(new_p)
    p_par=np.add.reduceat(o_hist,start_p); ev_p=ev_h[start_p]
    new_e=np.ones(len(p_par),dtype=bool)
    if len(p_par)>1: new_e[1:]=ev_p[1:]!=ev_p[:-1]
    start_e=np.flatnonzero(new_e)
    o_found=-np.expm1(np.add.reduceat(np.log1p(-np.clip(p_par,0.0,1.0-1e-15)),start_e))
    samp=tab['samp_ev'][chan]
    all_ev=np.unique(samp)
    o_all=np.zeros(len(all_ev)); have=ev_p[start_e]
    if len(have):
        pos=np.searchsorted(all_ev,have); o_all[pos]=o_found
    norm=(n_tot/(n_smp*n_gen))*L.RECO_EFF
    eps=float(np.sum(o_all))*norm
    se=norm*math.sqrt(max(float(np.var(o_all)),0.0)*len(o_all))
    return eps,se,o_all,all_ev

def main():
    first_smp,last_smp=load_higgs()
    sc=json.load(open(os.path.join(RUN,'outputs','data','signal_scan.json')))
    masses=np.array(sc['masses']); th=np.array(sc['theta2']); sig=sc['signal_per_fb_per_channel']
    def crossings(x,y,level):
        out=[]; lx=np.log10(x); ly=np.log10(np.maximum(y,1e-300))
        for i in range(len(x)-1):
            f0=ly[i]-math.log10(level); f1=ly[i+1]-math.log10(level)
            if f0*f1<0: out.append(10**(lx[i]+f0/(f0-f1)*(lx[i+1]-lx[i])))
        return out
    res={'settings':{'nparents':NP,'ndecay':NDECAY,'n_gen':200000,'n_tot':199105},
         'rows':{}, 'copy_stats':{}}
    d=np.loadtxt(os.path.join(RUN,'work','magnus','samples','higgs.txt'),comments='#')
    ev=d[:,0].astype(int); u,first,counts=np.unique(ev,return_index=True,return_counts=True)
    last=first+counts-1
    pt=np.hypot(d[:,2],d[:,3])
    res['copy_stats']={'rows':int(len(d)),'events':int(len(u)),
        'copies_min':int(counts.min()),'copies_max':int(counts.max()),'copies_mean':float(counts.mean()),
        'frac_first_ne_last':float(np.mean(np.linalg.norm(d[first,2:5]-d[last,2:5],axis=1)>1e-6)),
        'mean_pT_first':float(pt[first].mean()),'mean_pT_last':float(pt[last].mean()),
        'rows_with_pT_below_1e-4':int((pt<1e-4).sum()),
        'first_copy_pT_median':float(np.median(pt[first])),
        'mass_mean':float(np.mean(np.sqrt(np.maximum(d[:,5]**2-np.sum(d[:,2:5]**2,axis=1),0)))),
        'mass_std':float(np.std(np.sqrt(np.maximum(d[:,5]**2-np.sum(d[:,2:5]**2,axis=1),0))))}
    for ms in (5.0,30.0,62.0):
        i=np.argmin(abs(masses-ms))
        row=sig[str(masses[i])]
        s=np.array([r['B']+r['K']+r['H'] for r in row])
        cr=crossings(th,s,2.3/25.0)
        t2_cross=min(cr) if cr else None
        ctau_cross=(sw.HBAR_C/(sw.gamma_total(ms)*t2_cross)) if t2_cross else None
        # parent subset identical for all variants
        R._np_rng=np.random.default_rng(424242)
        nrows=len(u)            # first/last arrays are indexed by row position 0..nrows-1
        sub=R._np_rng.choice(nrows,size=min(NP,nrows),replace=False)
        sub=np.sort(sub)
        out={'m':ms,'theta2_cross_25fb':t2_cross,'ctau_at_cross_m':ctau_cross,'variants':{}}
        tabs={}
        for name,smp,mode in (('first',first_smp,'MH'),('last',last_smp,'MH'),
                              ('first_invmass',first_smp,'inv'),('last_invmass',last_smp,'inv')):
            t0=time.time()
            tab=build_tab(ms,smp,sub,mode)
            tabs[name]=tab
            rowres={}
            for ctau in (0.1,1.0,10.0,100.0):
                t2=sw.HBAR_C/(sw.gamma_total(ms)*ctau)
                p=R.per_scalar_p(ms,tab,t2,tab['pz'],tab['z'],tab['A'])
                eps,se,o,aev=channel_efficiency_events(tab,p,'H')
                rowres['ctau_%.3g'%ctau]={'theta2':t2,'eps':eps,'se':se}
            out['variants'][name]={'result':rowres,'nrows':int(len(tab['ev'])),'time':time.time()-t0}
        # paired comparisons (same parent list, same seeds -> same sampled parents)
        for other in ('last','first_invmass'):
            for ctau in (0.1,1.0,10.0,100.0):
                k='ctau_%.3g'%ctau
                t2=sw.HBAR_C/(sw.gamma_total(ms)*ctau)
                p1=R.per_scalar_p(ms,tabs['first'],t2,tabs['first']['pz'],tabs['first']['z'],tabs['first']['A'])
                p2=R.per_scalar_p(ms,tabs[other],t2,tabs[other]['pz'],tabs[other]['z'],tabs[other]['A'])
                e1,s1,o1,ae1=channel_efficiency_events(tabs['first'],p1,'H')
                e2,s2,o2,ae2=channel_efficiency_events(tabs[other],p2,'H')
                ratio=o2.sum()/o1.sum() if o1.sum()>0 else float('nan')
                # paired bootstrap over the common event list
                n=len(o1); dr=np.zeros(n)
                if n>0 and len(o2)==n:
                    dr=o2-o1
                rngb=np.random.default_rng(7)
                bs=[]
                for _ in range(400):
                    idx=rngb.integers(0,n,n)
                    bs.append(o2[idx].sum()/o1[idx].sum() if o1[idx].sum()>0 else np.nan)
                bs=np.array(bs); bs=bs[np.isfinite(bs)]
                out.setdefault('ratios',{})[other+'_'+k]={'ratio':float(ratio),'boot_se':float(np.std(bs)),
                    'eps1':e1,'eps2':e2,'se1':s1,'se2':s2,'n':n}
        res['rows']['m%.1f'%ms]=out
        print('m=%.1f th2cross=%.3e  '%(ms,t2_cross if t2_cross else float('nan'))
              +' '.join('%s: eps(ctau10)=%.4e ratio(last/first)=%.4f'%(
                  v,out['variants'][v]['result']['ctau_10']['eps'],
                  out['ratios'].get(v+'_ctau_10',{}).get('ratio',float('nan')))
                  for v in ('first','last','first_invmass','last_invmass')),flush=True)
    with open(os.path.join(HERE,'f06_results.json'),'w') as f: json.dump(res,f,indent=1)
    print('done')
if __name__=='__main__': main()
