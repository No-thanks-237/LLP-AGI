# -*- coding: utf-8 -*-
import sys, math, json
import numpy as np
sys.path.insert(0, 'runs/bc5_lhcb_downstream_v3/work/analysis')
import run_scan as R
import bc5_production as P
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read_foresee(energy):
    c=[]
    for line in open('temp/foresee/files/hadrons/%sTeV/Pythia8/Pythia8_%sTeV_25.txt'%(energy,energy)):
        if line.startswith('#'): continue
        c.append([float(x) for x in line.split()])
    a=np.array(c); return a[:,0],a[:,1],a[:,2]

def cells(energy, nsub=4):
    lth,lp,w = read_foresee(energy); nz=w>1e-6
    lth,lp,w = lth[nz],lp[nz],w[nz]
    g=np.linspace(-0.025,0.025,nsub); a,b=np.meshgrid(g,g,indexing='ij')
    th = 10**(lth[:,None]+a.ravel()[None,:]).ravel(); pp = 10**(lp[:,None]+b.ravel()[None,:]).ravel()
    ww = np.repeat(w/(nsub*nsub), nsub*nsub)
    return th,pp,ww

smp=R.load_sample('higgs.txt.gz'); cand=np.where(np.abs(smp['status'])==62)[0]
E,px,py,pz=smp['e'][cand],smp['px'][cand],smp['py'][cand],smp['pz'][cand]
code=smp['pdg'][cand]
eta=np.arcsinh(pz/np.maximum(np.hypot(px,py),1e-300)); pt=np.hypot(px,py); p=np.sqrt(px*px+py*py+pz*pz)
grp=np.array([P.HIGGS_PROC_GROUP.get(int(c),'other') for c in code])
wH=np.ones(len(cand))
for g_ in sorted(set(grp.tolist())):
    if g_=='other' or g_ not in P.SIGMA_H: continue
    m_=grp==g_; wH[m_]=P.SIGMA_H[g_]/(P.SIGMA_H_TOTAL*float(m_.sum())/len(cand))
ggf=grp=='ggF'; sel=eta>0
th,pp,ww=cells('14'); eta_f=-np.log(np.tan(th/2.0)); pt_f=pp*np.sin(th)

def cdf_ratio(xs, v, wv):
    out=[]
    for x in xs:
        f=float(np.sum(wv[v<x])/np.sum(wv)); o=float(np.sum(v[xs is None]))
        out.append(f)
    return np.array(out)
def frac(v,x):
    return float(np.mean(v<x))
def wfrac(v,w,x):
    return float(np.sum(w[v<x])/np.sum(w))

fig, ax = plt.subplots(2,3, figsize=(15,8))
# (a) eta distribution
eb=np.arange(0.,5.51,0.5); cen=0.5*(eb[1:]+eb[:-1]); wd=eb[1:]-eb[:-1]
def hn(v,w,b):
    h,_=np.histogram(v,bins=b,weights=w); return h/h.sum()
ax[0,0].step(cen, hn(eta_f,ww,eb)/wd, where='mid', color='C3', lw=1.8, label='FORESEE 14 TeV h(25)')
ax[0,0].step(cen, hn(eta[sel],None,eb)/wd, where='mid', color='C0', lw=1.4, label='ours all (unweighted)')
ax[0,0].step(cen, hn(eta[ggf&sel],None,eb)/wd, where='mid', color='C2', lw=1.2, ls='--', label='ours ggF only')
ax[0,0].set_xlabel(r'$\eta_h$'); ax[0,0].set_ylabel(r'$1/N\,dN/d\eta$'); ax[0,0].legend(fontsize=7)
ax[0,0].set_title(r'(a) Higgs $\eta$ (forward hemisphere)')
# (b) eta CDF ratio
xs=np.array([1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.])
fo=np.array([wfrac(eta[sel],np.ones(sel.sum()),x) for x in xs])
fg=np.array([wfrac(eta[ggf&sel],np.ones((ggf&sel).sum()),x) for x in xs])
ff=np.array([wfrac(eta_f,ww,x) for x in xs])
ax[0,1].plot(xs, ff/fo, 'o-', color='C0', label='FORESEE / ours all')
ax[0,1].plot(xs, ff/fg, 's--', color='C2', label='FORESEE / ours ggF')
ax[0,1].axhline(1,color='k',lw=0.6); ax[0,1].set_xlabel(r'$\eta$'); ax[0,1].set_ylabel('ratio of CDF')
ax[0,1].legend(fontsize=7); ax[0,1].set_title(r'(b) cumulative $dN/d\eta$ ratio')
# (c) pT distribution
pb=np.array([0,2,5,8,12,20,30,40,60,80,120,200,400]); cp=0.5*(pb[1:]+pb[:-1]); wp=pb[1:]-pb[:-1]
ax[0,2].step(cp, hn(pt_f,ww,pb)/wp, where='mid', color='C3', lw=1.8, label='FORESEE 14 TeV')
ax[0,2].step(cp, hn(pt[sel],None,pb)/wp, where='mid', color='C0', lw=1.4, label='ours all')
ax[0,2].step(cp, hn(pt[ggf&sel],None,pb)/wp, where='mid', color='C2', lw=1.2, ls='--', label='ours ggF')
ax[0,2].set_xscale('log'); ax[0,2].set_xlabel(r'$p_{T,h}$ [GeV]'); ax[0,2].set_ylabel(r'$1/N\,dN/dp_T$ [GeV$^{-1}$]')
ax[0,2].legend(fontsize=7); ax[0,2].set_title(r'(c) Higgs $p_T$')
# (d) pT CDF ratio
xs2=np.array([5.,10.,15.,20.,30.,50.,80.,120.,200.])
fo2=np.array([wfrac(pt[sel],np.ones(sel.sum()),x) for x in xs2])
fg2=np.array([wfrac(pt[ggf&sel],np.ones((ggf&sel).sum()),x) for x in xs2])
ff2=np.array([wfrac(pt_f,ww,x) for x in xs2])
ax[1,0].plot(xs2, ff2/fo2, 'o-', color='C0', label='FORESEE / ours all')
ax[1,0].plot(xs2, ff2/fg2, 's--', color='C2', label='FORESEE / ours ggF')
ax[1,0].axhline(1,color='k',lw=0.6); ax[1,0].set_xlabel(r'$p_T$ [GeV]'); ax[1,0].set_ylabel('ratio of CDF')
ax[1,0].legend(fontsize=7); ax[1,0].set_title(r'(d) cumulative $dN/dp_T$ ratio')
# (e) K ratio
tab=json.load(open('runs/bc5_lhcb_downstream_v3/work/review_check/foresee/results_final_table.json'))
ms=np.array([5.08,13.02,33.37,53.42])
for key,lab,col in (('ratio','all / FORESEE','C0'),('ratio_mix','SIGMA_H mix / FORESEE','C1'),('ratio_ggf','ggF / FORESEE','C2')):
    y=np.array([tab[str(m)][key] for m in ms]); ey=np.array([tab[str(m)]['SE_ratio'] for m in ms])
    ax[1,1].errorbar(ms,y,yerr=ey,fmt='o-',ms=4,color=col,label=lab,lw=1.2)
ax[1,1].axhline(1,color='k',lw=0.6); ax[1,1].set_xscale('log'); ax[1,1].set_xlabel(r'$m_S$ [GeV]')
ax[1,1].set_ylabel(r'$K_{\rm ours}/K_{\rm FORESEE}$'); ax[1,1].legend(fontsize=7)
ax[1,1].set_title(r'(e) $K=\langle 1_{2<\eta_S<5}\,m_S/p_z\rangle$'); ax[1,1].set_ylim(0.94,1.03)
# (f) K_FORESEE vs energy
tr=json.load(open('runs/bc5_lhcb_downstream_v3/work/review_check/foresee/results_K.json'))['trend']
ens=[14,27,100]
for m,col in zip(ms,['C0','C1','C2','C3']):
    y=[tr[str(e)][str(m)] for e in ens]
    ax[1,2].plot(ens,y,'o-',color=col,label='$m_S$=%.2f'%m)
ax[1,2].set_xscale('log'); ax[1,2].set_xticks([14,27,100]); ax[1,2].set_xticklabels(['14','27','100'])
ax[1,2].minorticks_off(); ax[1,2].set_xlabel(r'$\sqrt{s}$ [TeV]'); ax[1,2].set_ylabel(r'$K_{\rm FORESEE}$')
ax[1,2].legend(fontsize=7); ax[1,2].set_title(r'(f) energy trend of $K$ (FORESEE h spectra)')
for a in ax.ravel(): a.grid(alpha=0.25)
plt.tight_layout()
plt.savefig('runs/bc5_lhcb_downstream_v3/work/review_check/foresee/fig_foresee_compare.png', dpi=120)
print('wrote figure')
