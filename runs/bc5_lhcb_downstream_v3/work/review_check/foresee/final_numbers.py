# -*- coding: utf-8 -*-
"""Final numbers for the FORESEE report: spectra tables, FORESEE histogram-Poisson error,
energy trend of the mother moments, and the 13-vs-14 TeV bias estimate."""
import sys, os, math, json
import numpy as np
sys.path.insert(0, 'runs/bc5_lhcb_downstream_v3/work/analysis')
import run_scan as R
import bc5_production as P

FS = 'temp/foresee/files/hadrons/%sTeV/Pythia8/Pythia8_%sTeV_25.txt'
M_S_LIST = [5.08, 13.02, 33.37, 53.42]

def read_foresee(energy):
    c = []
    for line in open(FS % (energy, energy)):
        if line.startswith('#'): continue
        c.append([float(x) for x in line.split()])
    a = np.array(c); return a[:,0], a[:,1], a[:,2]

def cell_points(energy, nsub=4):
    lth, lp, w = read_foresee(energy)
    nz = w > 1e-6
    lth, lp, w = lth[nz], lp[nz], w[nz]
    g = np.linspace(-0.025, 0.025, nsub)
    a, b = np.meshgrid(g, g, indexing='ij')
    lth_m = (lth[:, None] + a.ravel()[None, :]).ravel()
    lp_m = (lp[:, None] + b.ravel()[None, :]).ravel()
    w_m = np.repeat(w/(nsub*nsub), nsub*nsub)
    th, pp = 10**lth_m, 10**lp_m
    eta = -np.log(np.tan(th/2.0)); pt = pp*np.sin(th)
    return eta, pt, pp, w_m, w

# ---------- moments trend
print('=== FORESEE mother moments vs collision energy (forward hemisphere, weighted) ===')
print('%-6s %8s %8s %8s %9s %12s' % ('E[TeV]', '<eta>', '<pT>', '<p>', 'f(2<eta<5)','sum w [pb]'))
mom = {}
for e in ('14','27','100'):
    eta, pt, pp, w_m, w = cell_points(e)
    W = w_m/w_m.sum()
    mom[e] = dict(eta=float(np.sum(W*eta)), pt=float(np.sum(W*pt)), p=float(np.sum(W*pp)),
                  f=float(np.sum(w_m[(eta>2)&(eta<5)])/w_m.sum()), sig=float(w_m.sum()),
                  q=float(w.min()))
    print('%-6s %8.3f %8.2f %8.1f %9.4f %12.4f' % (e, mom[e]['eta'], mom[e]['pt'], mom[e]['p'],
                                                     mom[e]['f'], mom[e]['sig']))
print()

# ---------- our sample
smp = R.load_sample('higgs.txt.gz')
cand = np.where(np.abs(smp['status']) == 62)[0]
cand = cand[np.argsort(smp['ev'][cand], kind='stable')]
E0, px0, py0, pz0 = smp['e'][cand], smp['px'][cand], smp['py'][cand], smp['pz'][cand]
code = smp['pdg'][cand]
p0 = np.sqrt(px0*px0+py0*py0+pz0*pz0)
mH0 = np.sqrt(np.maximum(E0*E0-p0*p0, 0.0))
etaH0 = np.arcsinh(pz0/np.maximum(np.hypot(px0,py0),1e-300)); ptH0 = np.hypot(px0,py0)
grp = np.array([P.HIGGS_PROC_GROUP.get(int(c),'other') for c in code])
wH = np.ones(len(cand))
for g_ in sorted(set(grp.tolist())):
    if g_ == 'other' or g_ not in P.SIGMA_H: continue
    m_ = grp == g_
    wH[m_] = P.SIGMA_H[g_]/(P.SIGMA_H_TOTAL*float(m_.sum())/len(cand))
ggf = grp == 'ggF'
sel = etaH0 > 0

# ---------- dN/deta table
eta_bins = np.array([0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5])
eta_f, pt_f, pp_f, w_f, w_raw = cell_points('14')
def hist_norm(vals, wts, bins):
    h,_ = np.histogram(vals, bins=bins, weights=wts)
    return h/np.sum(h)
h_ours = hist_norm(etaH0[sel], wH[sel], eta_bins)
h_ggf  = hist_norm(etaH0[ggf & sel], wH[ggf & sel], eta_bins)
h_all_u = hist_norm(etaH0[sel], None, eta_bins)
h_f    = hist_norm(eta_f, w_f, eta_bins)
n_ours = np.histogram(etaH0[sel], bins=eta_bins)[0]
n_all_u = np.histogram(etaH0[sel], bins=eta_bins)[0]
n_ggf = np.histogram(etaH0[ggf & sel], bins=eta_bins)[0]
q = w_raw.min()
n_f = np.histogram(eta_f, bins=eta_bins, weights=w_f)[0]/q
print('=== dN/deta (forward hemisphere, normalised over 0.5<eta<5.5) ===')
print('%-10s %9s %9s %9s %9s' % ('eta bin', 'ours all', 'ours ggF', 'ours mix', 'FORESEE'))
for i in range(len(eta_bins)-1):
    print('%-10s %9.4f %9.4f %9.4f %9.4f' % ('%.1f-%.1f'%(eta_bins[i],eta_bins[i+1]),
          h_all_u[i], h_ggf[i], h_ours[i], h_f[i]))
print('%-10s %9s %9s %9s %9s' % ('ratio F/all','','','', '%.3f'%(h_f.sum()/h_all_u.sum())))
print()
print(' per-bin ratio FORESEE / ours(all), error from Poisson of both (FORESEE n_b = w/q):')
for i in range(len(eta_bins)-1):
    a = h_f[i]; b = h_all_u[i]
    ea = a*math.sqrt(1.0/max(n_f[i],1e-9)); eb = b*math.sqrt(1.0/max(n_all_u[i],1e-9))
    print('  %-8s %.3f +- %.3f     (FORESEE n_b ~ %.0f, ours N=%d)' % (
          '%.1f-%.1f'%(eta_bins[i],eta_bins[i+1]), a/b, (a/b)*math.sqrt((ea/a)**2+(eb/b)**2),
          n_f[i], n_all_u[i]))
print()
print(' dN/deta values per unit eta (shape normalised over 0.5-5.5):')
print('%-10s %9s %9s %9s %9s' % ('eta bin','ours all','ours ggF','ours mix','FORESEE'))
for i in range(len(eta_bins)-1):
    w = eta_bins[i+1]-eta_bins[i]
    print('%-10s %9.4f %9.4f %9.4f %9.4f' % ('%.1f-%.1f'%(eta_bins[i],eta_bins[i+1]),
          h_all_u[i]/w, h_ggf[i]/w, h_ours[i]/w, h_f[i]/w))
print()

# ---------- dN/dpT table
pt_bins = np.array([2,5,8,12,20,30,40,60,80,120,200,400])
h_pt_ours = hist_norm(ptH0[sel], wH[sel], pt_bins)
h_pt_ggf  = hist_norm(ptH0[ggf & sel], wH[ggf & sel], pt_bins)
h_pt_all  = hist_norm(ptH0[sel], None, pt_bins)
h_pt_f    = hist_norm(pt_f, w_f, pt_bins)
npt_ours = np.histogram(ptH0[sel], bins=pt_bins)[0]
npt_ggf = np.histogram(ptH0[ggf & sel], bins=pt_bins)[0]
npt_f = np.histogram(pt_f, bins=pt_bins, weights=w_f)[0]/q
print('=== dN/dpT (forward hemisphere, normalised over 2<pT<400 GeV) ===')
print('%-12s %9s %9s %9s %9s %9s' % ('pT bin','ours all','ours ggF','ours mix','FORESEE','F/all'))
for i in range(len(pt_bins)-1):
    r = h_pt_f[i]/max(h_pt_all[i],1e-12)
    print('%-12s %9.4f %9.4f %9.4f %9.4f %9.3f' % ('%.0f-%.0f'%(pt_bins[i],pt_bins[i+1]),
          h_pt_all[i], h_pt_ggf[i], h_pt_ours[i], h_pt_f[i], r))
print()
print(' per-bin ratio FORESEE / ours(all):')
for i in range(len(pt_bins)-1):
    a = h_pt_f[i]; b = h_pt_all[i]
    ea = a*math.sqrt(1.0/max(npt_f[i],1e-9)); eb = b*math.sqrt(1.0/max(npt_ours[i],1e-9))
    print('  %-10s %.3f +- %.3f' % ('%.0f-%.0f'%(pt_bins[i],pt_bins[i+1]), a/b,
          (a/b)*math.sqrt((ea/a)**2+(eb/b)**2)))
print()

# ---------- FORESEE histogram-Poisson error on K
def K_foresee(energy, ms, wts, seed, nsub=4):
    eta, pt, pp, w_m, _ = cell_points(energy, nsub)
    lth, lp, _ = read_foresee(energy)
    # rebuild theta/phi from eta? simpler: recompute theta
    th = 2.0*np.arctan(np.exp(-eta))
    # p from pt/sin(th)
    p = pt/np.sin(th)
    rng = np.random.default_rng(seed)
    phi = rng.uniform(-math.pi, math.pi, len(th))
    th = np.concatenate([th, math.pi-th]); p = np.concatenate([p,p])
    phi = np.concatenate([phi, phi+math.pi]); ww = np.concatenate([wts, wts])
    px = p*np.sin(th)*np.cos(phi); py = p*np.sin(th)*np.sin(phi); pz = p*np.cos(th)
    E = np.sqrt(p*p+P.M_H**2)
    return px, py, pz, E, ww

def K_from_p4(ms, px, py, pz, E, ww, seed, ndecay=2):
    rng = np.random.default_rng(seed)
    N = len(E)
    pstar = math.sqrt(P.M_H**2/4 - ms*ms)
    u = rng.normal(size=(N*ndecay,3)); u /= np.linalg.norm(u,axis=1,keepdims=True)
    d = pstar*u; idx = np.repeat(np.arange(N), ndecay)
    bvec = np.stack([px[idx],py[idx],pz[idx]],axis=1)/np.maximum(E[idx],1e-300)[:,None]
    b2 = np.sum(bvec*bvec,axis=1); gamma = E[idx]/P.M_H
    per = np.zeros(N)
    for dd in (d,-d):
        bp = np.sum(dd*bvec,axis=1)
        k = (gamma-1.0)*bp/np.maximum(b2,1e-300) + gamma*(P.M_H/2)
        dp = dd + k[:,None]*bvec
        pz2 = dp[:,2]; pt2 = np.hypot(dp[:,0],dp[:,1])
        eta2 = np.arcsinh(pz2/np.maximum(pt2,1e-300))
        good = (eta2>2)&(eta2<5)&(pz2>0)&np.isfinite(eta2)
        r = np.zeros(len(dp)); r[good] = ms/pz2[good]
        per += 0.5*r.reshape(N,ndecay).mean(axis=1)
    W = ww/ww.sum()
    return float(np.sum(W*per))

print('=== FORESEE K: histogram-Poisson uncertainty (one event = q = %.3e pb) ===' % q)
out_pois = {}
for ms in M_S_LIST:
    px,py,pz,E,ww = K_foresee('14', ms, w_f, 1)
    base = K_from_p4(ms, px,py,pz,E,ww, 11)
    # Poisson resampling of the cell weights
    vals = []
    for r in range(40):
        rg = np.random.default_rng(300+r)
        # resample per-cell weights; the cell_points array has nsub*nsub entries per cell
        # -> draw one Poisson per cell and spread it
        ncell = len(w_raw[w_raw>1e-6])
        nb = w_raw[w_raw>1e-6]/q
        npois = rg.poisson(nb).astype(float)
        # spread each cell's resampled weight over its 16 sub-points
        wres = np.repeat(npois*q/16.0, 16)
        # mirror handled inside K_from_p4 call below: rebuild with wres
        # reuse geometry but replace ww (ww = concat(wres,wres))
        vals.append(K_from_p4(ms, px,py,pz,E,np.concatenate([wres,wres]), 400+r))
    vals = np.array(vals)
    out_pois[str(ms)] = dict(mean=float(vals.mean()), std=float(vals.std(ddof=1)), base=base)
    print('  m_S=%5.2f  K_F = %.5f   Poisson-resample mean %.5f, sd %.5f  (%.2f%% of K)' % (
          ms, base, vals.mean(), vals.std(ddof=1), 100*vals.std(ddof=1)/base))
json.dump({'mom': mom, 'poisson': out_pois, 'q': float(q)}, open(
    'runs/bc5_lhcb_downstream_v3/work/review_check/foresee/results_final.json','w'), indent=1)
print('\nwrote results_final.json')
