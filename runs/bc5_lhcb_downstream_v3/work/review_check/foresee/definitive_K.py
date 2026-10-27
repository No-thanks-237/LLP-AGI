# -*- coding: utf-8 -*-
"""Definitive K table: our 13 TeV sample (all / ggF / SIGMA_H mix) vs FORESEE 14 TeV,
with ALL statistical errors: our sample + FORESEE cell sub-sampling + decay sampling
+ FORESEE histogram-Poisson (one event = q)."""
import sys, os, math, json
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
    lth_m = (lth[:, None] + a.ravel()[None, :]).ravel()
    lp_m = (lp[:, None] + b.ravel()[None, :]).ravel()
    th = 10**lth_m; p = 10**lp_m
    rng = np.random.default_rng(seed)
    phi = rng.uniform(-math.pi, math.pi, len(th))
    # mirror (pp symmetry): forward file -> both hemispheres
    th2 = np.concatenate([th, math.pi-th])
    p2 = np.concatenate([p, p])
    phi2 = np.concatenate([phi, phi+math.pi])
    wsub = np.repeat(w/(nsub*nsub), nsub*nsub)
    return th2, p2, phi2, wsub, w, float(w.min())

def p4_from(th, p, phi, mh):
    return (np.sqrt(p*p+mh*mh), p*np.sin(th)*np.cos(phi), p*np.sin(th)*np.sin(phi), p*np.cos(th))

def scalar_R(E, px, py, pz, mH, ms, rng, ndecay=4):
    N = len(E)
    pstar = np.sqrt(np.maximum(mH*mH/4.0-ms*ms, 0.0))
    u = rng.normal(size=(N*ndecay, 3)); u /= np.maximum(np.linalg.norm(u, axis=1, keepdims=True), 1e-300)
    d = np.repeat(pstar, ndecay)[:, None]*u
    idx = np.repeat(np.arange(N), ndecay)
    Eb, pxb, pyb, pzb = E[idx], px[idx], py[idx], pz[idx]
    bvec = np.stack([pxb, pyb, pzb], axis=1)/np.maximum(Eb, 1e-300)[:, None]
    b2 = np.sum(bvec*bvec, axis=1); gamma = Eb/mH[idx]
    per = np.zeros(N)
    for dd, Ed in ((d, np.repeat(mH/2.0, ndecay)), (-d, np.repeat(mH/2.0, ndecay))):
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
    W = w/w.sum(); return float(np.sum(W*v)), float(math.sqrt(np.sum(W*W*(v-np.sum(W*v))**2)*len(v)/(len(v)-1)))

# our sample
smp = R.load_sample('higgs.txt.gz')
cand = np.where(np.abs(smp['status']) == 62)[0]
cand = cand[np.argsort(smp['ev'][cand], kind='stable')]
E0, px0, py0, pz0 = smp['e'][cand], smp['px'][cand], smp['py'][cand], smp['pz'][cand]
code = smp['pdg'][cand]
p0 = np.sqrt(px0*px0+py0*py0+pz0*pz0)
mH0 = np.sqrt(np.maximum(E0*E0-p0*p0, 0.0))
grp = np.array([P.HIGGS_PROC_GROUP.get(int(c), 'other') for c in code])
wH = np.ones(len(cand))
for g_ in sorted(set(grp.tolist())):
    if g_ == 'other' or g_ not in P.SIGMA_H: continue
    m_ = grp == g_
    wH[m_] = P.SIGMA_H[g_]/(P.SIGMA_H_TOTAL*float(m_.sum())/len(cand))
ggf = grp == 'ggF'

print('=== DEFINITIVE K TABLE:  K = <1_{2<eta_S<5} m_S/p_z> per scalar ===')
print('(the common detector factor int_1^2.5 A(z)dz cancels in every ratio)')
tab = {}
for ms in M_S_LIST:
    rng = np.random.default_rng(int(round(1000*ms))+3)
    per_all = scalar_R(E0, px0, py0, pz0, mH0, ms, rng)
    per_g = scalar_R(E0[ggf], px0[ggf], py0[ggf], pz0[ggf], mH0[ggf], ms, rng)
    K_all, SE_all = per_all.mean(), per_all.std(ddof=1)/math.sqrt(len(per_all))
    K_g, SE_g = per_g.mean(), per_g.std(ddof=1)/math.sqrt(len(per_g))
    K_mix, SE_mix = wmean(per_all, wH)
    # FORESEE: 20 realisations, each = fresh Poisson resample of the cell weights + fresh decay
    th, p, phi, wsub, wraw, q = cell_geometry('14', nsub=4, seed=1)
    Ef, pxf, pyf, pzf = p4_from(th, p, phi, P.M_H)
    vals = []
    ncell = int((wraw > 1e-6).sum())
    for irep in range(20):
        rg = np.random.default_rng(9000+irep)
        npois = rg.poisson(wraw[wraw > 1e-6]/q).astype(float)
        wres = np.repeat(npois*q/16.0, 16)
        wres = np.concatenate([wres, wres])
        per_f = scalar_R(Ef, pxf, pyf, pzf, np.full(len(Ef), P.M_H), ms, rg)
        Kf, _ = wmean(per_f, wres)
        vals.append(Kf)
    vals = np.array(vals)
    K_f = float(vals.mean()); SE_f = float(vals.std(ddof=1)/math.sqrt(len(vals)))
    ratio = K_all/K_f
    SE_ratio = ratio*math.sqrt((SE_all/K_all)**2 + (SE_f/K_f)**2)
    print()
    print('m_S = %5.2f GeV' % ms)
    print('  K_ours all           = %.5f +- %.5f' % (K_all, SE_all))
    print('  K_ours ggF           = %.5f +- %.5f' % (K_g, SE_g))
    print('  K_ours SIGMA_H mix   = %.5f +- %.5f' % (K_mix, SE_mix))
    print('  K_FORESEE 14 TeV     = %.5f +- %.5f  (20 realisations; includes Poisson)' % (K_f, SE_f))
    print('  ratio all/FORESEE    = %.4f +- %.4f' % (ratio, SE_ratio))
    print('  ratio mix/FORESEE    = %.4f' % (K_mix/K_f))
    print('  ratio ggF/FORESEE    = %.4f' % (K_g/K_f))
    tab[str(ms)] = dict(K_all=K_all, SE_all=SE_all, K_g=K_g, SE_g=SE_g, K_mix=K_mix, SE_mix=SE_mix,
                        K_f=K_f, SE_f=SE_f, ratio=ratio, SE_ratio=SE_ratio,
                        ratio_mix=K_mix/K_f, ratio_ggf=K_g/K_f)
json.dump(tab, open('runs/bc5_lhcb_downstream_v3/work/review_check/foresee/results_definitive.json','w'), indent=1)
print('\nwrote results_definitive.json')
try:
    import matplotlib
    print('matplotlib', matplotlib.__version__)
except Exception as e:
    print('no matplotlib:', e)
