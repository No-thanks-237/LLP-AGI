
import sys, math
import numpy as np
sys.path.insert(0, 'runs/bc5_lhcb_downstream_v3/work/analysis')
import run_scan as R
import bc5_production as P

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
    th = 10**(lth[:,None]+a.ravel()[None,:]).ravel()
    pp = 10**(lp[:,None]+b.ravel()[None,:]).ravel()
    ww = np.repeat(w/(nsub*nsub), nsub*nsub)
    return th,pp,ww

# ours
smp=R.load_sample('higgs.txt.gz'); cand=np.where(np.abs(smp['status'])==62)[0]
E,px,py,pz=smp['e'][cand],smp['px'][cand],smp['py'][cand],smp['pz'][cand]
code=smp['pdg'][cand]
p=np.sqrt(px*px+py*py+pz*pz)
eta=np.arcsinh(pz/np.maximum(np.hypot(px,py),1e-300)); pt=np.hypot(px,py)
grp=np.array([P.HIGGS_PROC_GROUP.get(int(c),'other') for c in code])
wH=np.ones(len(cand))
for g_ in sorted(set(grp.tolist())):
    if g_=='other' or g_ not in P.SIGMA_H: continue
    m_=grp==g_; wH[m_]=P.SIGMA_H[g_]/(P.SIGMA_H_TOTAL*float(m_.sum())/len(cand))
ggf=grp=='ggF'; sel=eta>0

th,pp,ww = cells('14')
eta_f=-np.log(np.tan(th/2.0)); pt_f=pp*np.sin(th)
# full-solid-angle FORESEE: add mirrored hemisphere (only needed for eta, not for |eta|-CDF)
eta_f_full=np.concatenate([eta_f,-eta_f]); pt_f_full=np.concatenate([pt_f,pt_f])
ww_full=np.concatenate([ww,ww])

def frac_lt(vals, wts, x):
    return float(np.sum(wts[vals< x])/np.sum(wts))

print('=== CDF comparison, forward hemisphere (eta>0), FORESEE vs ours ===')
print('eta cut  < eta  :   ours all   ours ggF   ours mix   FORESEE   F/all   F/ggF')
for x in (1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0):
    fo = frac_lt(eta[sel], np.ones(sel.sum()), x)
    fg = frac_lt(eta[ggf&sel], np.ones((ggf&sel).sum()), x)
    fm = frac_lt(eta[sel], wH[sel], x)
    ff = frac_lt(eta_f, ww, x)
    print('eta<%4.1f        : %9.4f %9.4f %9.4f %9.4f %7.3f %7.3f' % (x,fo,fg,fm,ff,ff/fo,ff/fg))
print()
print('pT cut  < pT  :   ours all   ours ggF   ours mix   FORESEE   F/all   F/ggF')
for x in (5.,10.,15.,20.,30.,50.,80.,120.,200.):
    fo = frac_lt(pt[sel], np.ones(sel.sum()), x)
    fg = frac_lt(pt[ggf&sel], np.ones((ggf&sel).sum()), x)
    fm = frac_lt(pt[sel], wH[sel], x)
    ff = frac_lt(pt_f, ww, x)
    print('pT <%6.1f    : %9.4f %9.4f %9.4f %9.4f %7.3f %7.3f' % (x,fo,fg,fm,ff,ff/fo,ff/fg))
print()
# median etc from CDF
for name, v, w in (('ours all', eta[sel], np.ones(sel.sum())), ('ours ggF', eta[ggf&sel], np.ones((ggf&sel).sum())),
                   ('FORESEE', eta_f, ww)):
    c = np.cumsum(w)/np.sum(w); med = v[np.searchsorted(c,0.5)]
    print('%-10s eta median %.3f' % (name, med))
for name, v, w in (('ours all', pt[sel], np.ones(sel.sum())), ('ours ggF', pt[ggf&sel], np.ones((ggf&sel).sum())),
                   ('FORESEE', pt_f, ww)):
    c = np.cumsum(w)/np.sum(w); med = v[np.searchsorted(c,0.5)]
    print('%-10s pT  median %.2f' % (name, med))
