# -*- coding: utf-8 -*-
"""Final table + figure. FORESEE's histogram-stat error on K is the SPREAD of K over
Poisson-resampled histograms (the shipped file is one realisation), not the SE of the mean."""
import sys, math, json
import numpy as np
sys.path.insert(0, 'runs/bc5_lhcb_downstream_v3/work/analysis')
import run_scan as R
import bc5_production as P

M_S_LIST = [5.08, 13.02, 33.37, 53.42]
FS = 'temp/foresee/files/hadrons/%sTeV/Pythia8/Pythia8_%sTeV_25.txt'

def read_foresee(energy):
    c = []
    for line in open(FS % (energy, energy)):
        if line.startswith('#'): continue
        c.append([float(x) for x in line.split()])
    a = np.array(c); return a[:,0], a[:,1], a[:,2]

def cell_geometry(energy, nsub=4, seed=1):
    lth, lp, w = read_foresee(energy)
    nz = w > 1e-6
    lth, lp, w = lth[nz], lp[nz], w[nz]
    g = np.linspace(-0.025, 0.025, nsub)
    a, b = np.meshgrid(g, g, indexing='ij')
    th = 10**(lth[:, None] + a.ravel()[None, :]).ravel()
    pp = 10**(lp[:, None] + b.ravel()[None, :]).ravel()
    rng = np.random.default_rng(seed)
    phi = rng.uniform(-math.pi, math.pi, len(th))
    th2 = np.concatenate([th, math.pi-th]); p2 = np.concatenate([pp, pp])
    phi2 = np.concatenate([phi, phi+math.pi])
    wsub = np.repeat(w/(nsub*nsub), nsub*nsub)
    wsub = np.concatenate([wsub, wsub])
    return th2, p2, phi2, wsub, w, float(w.min())

def p4_from(th, pp, phi, mh):
    return (np.sqrt(pp*pp+mh*mh), pp*np.sin(th)*np.cos(phi), pp*np.sin(th)*np.sin(phi), pp*np.cos(th))

def scalar_R(E, px, py, pz, mH, ms, rng, ndecay=4):
    N = len(E)
    pstar = np.sqrt(np.maximum(mH*mH/4.0-ms*ms, 0.0))
    u = rng.normal(size=(N*ndecay, 3)); u /= np.maximum(np.linalg.norm(u, axis=1, keepdims=True), 1e-300)
    d = np.repeat(pstar, ndecay)[:, None]*u
    idx = np.repeat(np.arange(N), ndecay)
    Eb, pxb, pyb, pzb = E[idx], px[idx], py[idx], pz[idx]
    bvec = np.stack([pxb, pyb, pzb], axis=1)/np.maximum(Eb, 1e-300)[:, None]
    b2 = np.sum(bvec*bvec, axis=1); gamma = Eb/mH[idx]
    per = np.zeros(N); Ed = np.repeat(mH/2.0, ndecay)
    for dd in (d, -d):
        bp = np.sum(dd*bvec, axis=1)
        k = (gamma-1.0)*bp/np.maximum(b2, 1e-300) + gamma*Ed
        dp = dd + k[:, None]*bvec
        pz2 = dp[:, 2]; pt2 = np.hypot(dp[:, 0], dp[:, 1])
        e2 = np.arcsinh(pz2/np.maximum(pt2, 1e-300))
        good = (e2 > 2.0) & (e2 < 5.0) & (pz2 > 0.0) & np.isfinite(e2)
        r = np.zeros(len(dp)); r[good] = ms/pz2[good]
        per += 0.5*r.reshape(N, ndecay).mean(axis=1)
    return per

def wmean(v, w):
    W = w/w.sum(); m = float(np.sum(W*v))
    return m, float(math.sqrt(np.sum(W*W*(v-m)**2)*len(v)/(len(v)-1)))

smp = R.load_sample('higgs.txt.gz')
cand = np.where(np.abs(smp['status']) == 62)[0]
cand = cand[np.argsort(smp['ev'][cand], kind='stable')]
E0, px0, py0, pz0 = smp['e'][cand], smp['px'][cand], smp['py'][cand], smp['pz'][cand]
code = smp['pdg'][cand]
p0 = np.sqrt(px0*px0+py0*py0+pz0*pz0)
mH0 = np.sqrt(np.maximum(E0*E0-p0*p0, 0.0))
etaH0 = np.arcsinh(pz0/np.maximum(np.hypot(px0, py0), 1e-300)); ptH0 = np.hypot(px0, py0)
grp = np.array([P.HIGGS_PROC_GROUP.get(int(c), 'other') for c in code])
wH = np.ones(len(cand))
for g_ in sorted(set(grp.tolist())):
    if g_ == 'other' or g_ not in P.SIGMA_H: continue
    m_ = grp == g_
    wH[m_] = P.SIGMA_H[g_]/(P.SIGMA_H_TOTAL*float(m_.sum())/len(cand))
ggf = grp == 'ggF'; sel = etaH0 > 0

th, pp, phi, wsub, wraw, q = cell_geometry('14', nsub=4, seed=1)
Ef, pxf, pyf, pzf = p4_from(th, pp, phi, P.M_H)
eta_f = np.arcsinh(pzf/np.maximum(np.hypot(pxf, pyf), 1e-300)); pt_f = np.hypot(pxf, pyf)

# ---------------- K table with proper errors
print('=== FINAL K TABLE ===')
tab = {}
for ms in M_S_LIST:
    rng = np.random.default_rng(int(round(1000*ms))+3)
    pa = scalar_R(E0, px0, py0, pz0, mH0, ms, rng)
    pg = scalar_R(E0[ggf], px0[ggf], py0[ggf], pz0[ggf], mH0[ggf], ms, rng)
    K_all, SE_all = pa.mean(), pa.std(ddof=1)/math.sqrt(len(pa))
    K_g, SE_g = pg.mean(), pg.std(ddof=1)/math.sqrt(len(pg))
    K_mix, SE_mix = wmean(pa, wH)
    vals = []
    for irep in range(24):
        rg = np.random.default_rng(9000+irep)
        npois = rg.poisson(wraw[wraw > 1e-6]/q).astype(float)
        wres = np.repeat(npois*q/16.0, 16); wres = np.concatenate([wres, wres])
        pf = scalar_R(Ef, pxf, pyf, pzf, np.full(len(Ef), P.M_H), ms, rg)
        vals.append(wmean(pf, wres)[0])
    vals = np.array(vals)
    K_f = float(vals.mean()); SE_f_hist = float(vals.std(ddof=1))
    ratio = K_all/K_f
    SE_ratio = ratio*math.sqrt((SE_all/K_all)**2 + (SE_f_hist/K_f)**2)
    print()
    print('m_S=%5.2f  K_all=%.5f+-%.5f  K_ggf=%.5f  K_mix=%.5f  K_F=%.5f+-%.5f  all/F=%.4f+-%.4f  mix/F=%.4f  ggf/F=%.4f'
          % (ms, K_all, SE_all, K_g, K_mix, K_f, SE_f_hist, ratio, SE_ratio, K_mix/K_f, K_g/K_f))
    tab[str(ms)] = dict(K_all=K_all, SE_all=SE_all, K_ggf=K_g, SE_ggf=SE_g, K_mix=K_mix,
                        SE_mix=SE_mix, K_f=K_f, SE_f_hist=SE_f_hist, ratio=ratio,
                        SE_ratio=SE_ratio, ratio_mix=K_mix/K_f, ratio_ggf=K_g/K_f)
json.dump(tab, open('runs/bc5_lhcb_downstream_v3/work/review_check/foresee/results_final_table.json','w'), indent=1)

# ---------------- figure
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 2, figsize=(11, 8.4))
eta_bins = np.array([0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5])
pt_bins = np.array([2,5,8,12,20,30,40,60,80,120,200,400])

def hnorm(vals, wts, bins):
    h,_ = np.histogram(vals, bins=bins, weights=wts); return h/h.sum()
def hcount(vals, bins):
    return np.histogram(vals, bins=bins)[0]
def hw(vals, wts, bins):
    return np.histogram(vals, bins=bins, weights=wts)[0]

h_f = hnorm(eta_f, wsub, eta_bins)
n_f = hw(eta_f, wsub, eta_bins)/q
h_all = hnorm(etaH0[sel], None, eta_bins)
n_all = hcount(etaH0[sel], eta_bins)
h_ggf = hnorm(etaH0[ggf & sel], None, eta_bins)
n_ggf = hcount(etaH0[ggf & sel], eta_bins)
cen = 0.5*(eta_bins[1:]+eta_bins[:-1]); wd = eta_bins[1:]-eta_bins[:-1]
a0 = ax[0,0]
a0.step(cen, h_f/wd, where='mid', label='FORESEE 14 TeV Pythia8 h(25)', color='C3', lw=1.8)
a0.step(cen, h_all/wd, where='mid', label='ours 13 TeV, all HiggsSM (unweighted)', color='C0', lw=1.4)
a0.step(cen, h_ggf/wd, where='mid', label='ours 13 TeV, ggF only', color='C2', lw=1.2, ls='--')
a0.set_xlabel(r'$\eta_h$'); a0.set_ylabel(r'$1/N\ dN/d\eta$'); a0.legend(fontsize=7)
a0.set_title('(a) Higgs rapidity (forward hemisphere)')

a1 = ax[0,1]
for tag, h, n, col in (('all', h_all, n_all, 'C0'), ('ggF', h_ggf, n_ggf, 'C2')):
    r = h_f/h; e = r*np.sqrt(1/np.maximum(n_f,1) + 1/np.maximum(n,1))
    a1.errorbar(cen, r, yerr=e, fmt='o-', ms=3, color=col, label='FORESEE / ours (%s)'%tag, lw=1.2)
a1.axhline(1.0, color='k', lw=0.6)
a1.set_xlabel(r'$\eta_h$'); a1.set_ylabel('ratio'); a1.legend(fontsize=7)
a1.set_title(r'(b) dN/d$\eta$ ratio')

a2 = ax[1,0]
h_fp = hnorm(pt_f, wsub, pt_bins); n_fp = hw(pt_f, wsub, pt_bins)/q
h_ap = hnorm(ptH0[sel], None, pt_bins); n_ap = hcount(ptH0[sel], pt_bins)
h_gp = hnorm(ptH0[ggf & sel], None, pt_bins); n_gp = hcount(ptH0[ggf & sel], pt_bins)
cp = 0.5*(pt_bins[1:]+pt_bins[:-1]); wp = pt_bins[1:]-pt_bins[:-1]
a2.step(cp, h_fp/wp, where='mid', color='C3', lw=1.8, label='FORESEE 14 TeV')
a2.step(cp, h_ap/wp, where='mid', color='C0', lw=1.4, label='ours all')
a2.step(cp, h_gp/wp, where='mid', color='C2', lw=1.2, ls='--', label='ours ggF')
a2.set_xscale('log'); a2.set_xlabel(r'$p_{T,h}$ [GeV]'); a2.set_ylabel(r'$1/N\ dN/dp_T$ [GeV$^{-1}$]')
a2.legend(fontsize=7); a2.set_title(r'(c) Higgs $p_T$ (forward hemisphere)')

a3 = ax[1,1]
ms_arr = np.array(M_S_LIST)
for key, lab, col in (('ratio','all / FORESEE','C0'), ('ratio_mix','SIGMA_H mix / FORESEE','C1'),
                      ('ratio_ggf','ggF / FORESEE','C2')):
    y = np.array([tab[str(m)][key] for m in ms_arr]); ey = np.array([tab[str(m)]['SE_ratio'] for m in ms_arr])
    a3.errorbar(ms_arr, y, yerr=ey, fmt='o-', ms=4, color=col, label=lab, lw=1.2)
a3.axhline(1.0, color='k', lw=0.6); a3.set_xscale('log'); a3.set_xlabel(r'$m_S$ [GeV]')
a3.set_ylabel(r'$K_{\rm ours}/K_{\rm FORESEE}$')
a3.set_title(r'(d) long-lived-boundary coefficient $K=\langle 1_{2<\eta_S<5}\,m_S/p_z\rangle$')
a3.legend(fontsize=7); a3.set_ylim(0.94, 1.03)
for a in ax.ravel():
    a.grid(alpha=0.25)
plt.tight_layout()
plt.savefig('runs/bc5_lhcb_downstream_v3/work/review_check/foresee/fig_foresee_compare.png', dpi=130)
print('\nwrote fig_foresee_compare.png')
