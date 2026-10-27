# -*- coding: utf-8 -*-
"""Long-lived-boundary coefficient K for h -> SS, our 13 TeV Pythia8 sample vs FORESEE 14 TeV.

K(m_S) = < 1_{2<eta_S<5} * m_S/p_z >   (per scalar)

The detector integral int_1^2.5 A(z) dz is a common, event-independent factor at fixed m_S,
so it cancels in the ratio K_ours/K_FORESEE.
FORESEE ships the forward hemisphere only (14.73 pb at 14 TeV); the backward hemisphere is the
exact mirror (pp symmetry), so every forward mother is also entered mirrored (theta -> pi-theta).
"""
import sys, os, math, json, random
import numpy as np

sys.path.insert(0, 'runs/bc5_lhcb_downstream_v3/work/analysis')
import run_scan as R
import bc5_lib as L
import bc5_production as P

M_S_LIST = [5.08, 13.02, 33.37, 53.42]
ETA_LO, ETA_HI = 2.0, 5.0
FS_DIR = 'temp/foresee/files/hadrons/%sTeV/Pythia8/Pythia8_%sTeV_25.txt'

def eta_of(px, py, pz):
    return np.arcsinh(pz / np.maximum(np.hypot(px, py), 1e-300))

def scalars_from_mother(mH, E, px, py, pz, mS, rng, ndecay=2):
    pstar = np.sqrt(np.maximum(mH*mH/4.0 - mS*mS, 0.0))
    N = len(mH)
    u = rng.normal(size=(N*ndecay, 3))
    u /= np.maximum(np.linalg.norm(u, axis=1, keepdims=True), 1e-300)
    d = np.repeat(pstar, ndecay)[:, None] * u
    Estar = np.repeat(mH/2.0, ndecay)
    idx = np.repeat(np.arange(N), ndecay)
    Eb, pxb, pyb, pzb = E[idx], px[idx], py[idx], pz[idx]
    bvec = np.stack([pxb, pyb, pzb], axis=1) / np.maximum(Eb, 1e-300)[:, None]
    b2 = np.sum(bvec*bvec, axis=1)
    gamma = Eb / mH[idx]
    outs = []
    for dd, Ed in ((d, Estar), (-d, Estar)):
        bp = np.sum(dd*bvec, axis=1)
        k = (gamma - 1.0)*bp/np.maximum(b2, 1e-300) + gamma*Ed
        outs.append(dd + k[:, None]*bvec)
    return outs[0], outs[1], idx

def R_of(dp, mS):
    pz = dp[:, 2]; pt = np.hypot(dp[:, 0], dp[:, 1])
    eta = np.arcsinh(pz/np.maximum(pt, 1e-300))
    good = (eta > ETA_LO) & (eta < ETA_HI) & (pz > 0.0) & np.isfinite(eta)
    r = np.zeros(len(dp)); r[good] = mS/pz[good]
    return r

def K_of_mothers(mH, E, px, py, pz, mS, rng, weights=None, ndecay=2):
    d1, d2, idx = scalars_from_mother(mH, E, px, py, pz, mS, rng, ndecay)
    per = 0.5*(R_of(d1, mS) + R_of(d2, mS))
    N = len(mH)
    M = per.reshape(N, ndecay).mean(axis=1)
    if weights is None:
        K = float(M.mean()); SE = float(M.std(ddof=1)/math.sqrt(N))
    else:
        W = weights/weights.sum()
        K = float(np.sum(W*M))
        SE = float(math.sqrt(np.sum(W*W*(M-K)**2) * N/(N-1)))
    return K, SE, M

# ---------------------------------------------------------------- validation of the boost
def validate_boost():
    rr = random.Random(7); worst = 0.0
    for k in range(10):
        mH = 125.25
        E = rr.uniform(127.0, 500.0); th = rr.uniform(0.05, 3.0); ph = rr.uniform(-3.1, 3.1)
        p = math.sqrt(E*E - mH*mH)
        px, py, pz = p*math.sin(th)*math.cos(ph), p*math.sin(th)*math.sin(ph), p*math.cos(th)
        ms = 13.02
        ps = math.sqrt(mH*mH/4 - ms*ms)
        parts = [(ms, mH/2, ps, 0.0, 0.0), (ms, mH/2, -ps, 0.0, 0.0)]
        ref = L.boost(parts, E, (px, py, pz))
        bvec = np.array([[px, py, pz]])/E; b2 = float(np.sum(bvec*bvec)); gamma = E/mH
        for dd, refp in ((np.array([[ps, 0.0, 0.0]]), ref[0]), (np.array([[-ps, 0.0, 0.0]]), ref[1])):
            bp = float(np.sum(dd*bvec))
            kk = (gamma-1.0)*bp/b2 + gamma*(mH/2)
            got = dd[0] + kk*bvec[0]
            Ep = gamma*(mH/2 + bp)
            worst = max(worst, abs(got[0]-refp[2]) + abs(got[1]-refp[3]) + abs(got[2]-refp[4])
                        + abs(Ep-refp[1]))
    return worst

print('boost validation: max |delta| = %.3e   (bc5_lib.boost vs vectorised formula)' % validate_boost())
print()

# ---------------------------------------------------------------- our sample
smp = R.load_sample('higgs.txt.gz')
cand = np.where(np.abs(smp['status']) == 62)[0]
cand = cand[np.argsort(smp['ev'][cand], kind='stable')]
E0, px0, py0, pz0 = smp['e'][cand], smp['px'][cand], smp['py'][cand], smp['pz'][cand]
code = smp['pdg'][cand]
p0 = np.sqrt(px0*px0 + py0*py0 + pz0*pz0)
mH0 = np.sqrt(np.maximum(E0*E0 - p0*p0, 0.0))
etaH0 = eta_of(px0, py0, pz0); ptH0 = np.hypot(px0, py0)
grp = np.array([P.HIGGS_PROC_GROUP.get(int(c), 'other') for c in code])
wH = np.ones(len(cand))
for g_ in sorted(set(grp.tolist())):
    if g_ == 'other' or g_ not in P.SIGMA_H: continue
    m_ = grp == g_
    wH[m_] = P.SIGMA_H[g_] / (P.SIGMA_H_TOTAL * float(m_.sum())/len(cand))
ggf = grp == 'ggF'

# ---------------------------------------------------------------- FORESEE reader / sampling
def read_foresee(path):
    cols = []
    for line in open(path):
        if line.startswith('#'): continue
        cols.append([float(x) for x in line.split()])
    a = np.array(cols)
    return a[:,0], a[:,1], a[:,2]

def foresee_mothers(energy, nsub=4, seed=1):
    p = FS_DIR % (energy, energy)
    lth, lp, w = read_foresee(p)
    nz = w > 1e-6
    lth, lp, w = lth[nz], lp[nz], w[nz]
    g = np.linspace(-0.025, 0.025, nsub)
    a, b = np.meshgrid(g, g, indexing='ij')
    dl, dp = a.ravel(), b.ravel()
    lth_m = (lth[:, None] + dl[None, :]).ravel()
    lp_m = (lp[:, None] + dp[None, :]).ravel()
    w_m = np.repeat(w/(nsub*nsub), nsub*nsub)
    th = 10**lth_m; pp = 10**lp_m
    rng = np.random.default_rng(seed)
    phi = rng.uniform(-math.pi, math.pi, len(th))
    th = np.concatenate([th, math.pi - th])
    pp = np.concatenate([pp, pp])
    phi = np.concatenate([phi, phi + math.pi])
    w_m = np.concatenate([w_m, w_m])
    return th, pp, phi, w_m

def mothers_to_p4(th, pp, phi, mh):
    px = pp*np.sin(th)*np.cos(phi); py = pp*np.sin(th)*np.sin(phi); pz = pp*np.cos(th)
    return np.sqrt(pp*pp + mh*mh), px, py, pz

def moments_foresee(energy):
    lth, lp, w = read_foresee(FS_DIR % (energy, energy))
    nz = w > 1e-6
    th = 10**lth[nz]; pp = 10**lp[nz]; ww = w[nz]
    eta = -np.log(np.tan(th/2.0)); pt = pp*np.sin(th)
    E = np.sqrt(pp*pp + P.M_H**2)
    Wt = ww/ww.sum()
    return dict(eta=float(np.sum(Wt*eta)), pT=float(np.sum(Wt*pt)), p=float(np.sum(Wt*pp)),
                E=float(np.sum(Wt*E)), frac2_5=float(np.sum(ww[(eta>2)&(eta<5)])/ww.sum()),
                sig=float(ww.sum()))

print('=== moments of the mother Higgs (forward hemisphere) ===')
print('%-30s %8s %8s %8s %8s %10s' % ('sample', '<eta>', '<pT>', '<p>', '<E>', 'f(2<eta<5)'))
def row(tag, e, pt, p, E, f):
    print('%-30s %8.3f %8.2f %8.1f %8.1f %10.4f' % (tag, e, pt, p, E, f))
for tag, sel in (('ours 13TeV all (unweighted)', etaH0 > 0), ('ours 13TeV ggF only', ggf & (etaH0 > 0))):
    row(tag, etaH0[sel].mean(), ptH0[sel].mean(), p0[sel].mean(), E0[sel].mean(),
        np.mean((etaH0[sel]>2)&(etaH0[sel]<5)))
# mixture-weighted moments
W = wH[etaH0 > 0]/wH[etaH0 > 0].sum()
row('ours 13TeV SIGMA_H mix', np.sum(W*etaH0[etaH0>0]), np.sum(W*ptH0[etaH0>0]),
    np.sum(W*p0[etaH0>0]), np.sum(W*E0[etaH0>0]),
    np.sum(W[((etaH0>2)&(etaH0<5))[etaH0>0]]))
for e in ('14',):
    m = moments_foresee(e)
    row('FORESEE 14TeV Pythia8', m['eta'], m['pT'], m['p'], m['E'], m['frac2_5'])
print()

# ---------------------------------------------------------------- K
print('=== K = <1_{2<eta_S<5} m_S/p_z>  (per scalar; common detector factor cancelled) ===')
results = {}
for ms in M_S_LIST:
    rng = np.random.default_rng(int(round(ms*1000)) + 11)
    K_ours, SE_ours, _ = K_of_mothers(mH0, E0, px0, py0, pz0, ms, rng)
    K_g, SE_g, _ = K_of_mothers(mH0[ggf], E0[ggf], px0[ggf], py0[ggf], pz0[ggf], ms, rng)
    K_w, SE_w, _ = K_of_mothers(mH0, E0, px0, py0, pz0, ms, rng, weights=wH)
    fs_vals = []
    for irep in range(10):
        th, pp, ph, ww = foresee_mothers('14', nsub=4, seed=1000+irep)
        Ef, pxf, pyf, pzf = mothers_to_p4(th, pp, ph, P.M_H)
        Kf, SEf, _ = K_of_mothers(np.full(len(Ef), P.M_H), Ef, pxf, pyf, pzf, ms,
                                  np.random.default_rng(5000+irep), weights=ww)
        fs_vals.append(Kf)
    fs_vals = np.array(fs_vals)
    K_f = fs_vals.mean(); SE_f = fs_vals.std(ddof=1)/math.sqrt(len(fs_vals))
    ratio = K_ours/K_f
    SE_ratio = ratio*math.sqrt((SE_ours/K_ours)**2 + (SE_f/K_f)**2)
    print()
    print('m_S = %5.2f GeV' % ms)
    print('  K_ours (all)         = %.5f +- %.5f' % (K_ours, SE_ours))
    print('  K_ours (ggF only)    = %.5f +- %.5f' % (K_g, SE_g))
    print('  K_ours (SIGMA_H mix) = %.5f +- %.5f' % (K_w, SE_w))
    print('  K_FORESEE            = %.5f +- %.5f   (10 realisations, rms %.5f)' % (K_f, SE_f, fs_vals.std(ddof=1)))
    print('  ratio all/FORESEE    = %.4f +- %.4f' % (ratio, SE_ratio))
    print('  ratio ggF/FORESEE    = %.4f' % (K_g/K_f))
    print('  ratio mix/FORESEE    = %.4f' % (K_w/K_f))
    results[str(ms)] = dict(K_ours=K_ours, SE_ours=SE_ours, K_ggf=K_g, SE_ggf=SE_g, K_mix=K_w,
                            SE_mix=SE_w, K_f=K_f, SE_f=SE_f, ratio=ratio, SE_ratio=SE_ratio,
                            ratio_ggf=K_g/K_f, ratio_mix=K_w/K_f, fs=fs_vals.tolist())

# ---------------------------------------------------------------- energy trend 14 -> 27 -> 100 TeV
print()
print('=== energy trend of K_FORESEE (same h(125) spectrum machinery) ===')
trend = {}
for energy in ('14', '27', '100'):
    try:
        vals = []
        for irep in range(6):
            th, pp, ph, ww = foresee_mothers(energy, nsub=4, seed=2000+irep)
            Ef, pxf, pyf, pzf = mothers_to_p4(th, pp, ph, P.M_H)
            for ms in M_S_LIST:
                Kf, _, _ = K_of_mothers(np.full(len(Ef), P.M_H), Ef, pxf, pyf, pzf, ms,
                                        np.random.default_rng(7000+irep+int(ms)), weights=ww)
                vals.append((ms, Kf))
        trend[energy] = {}
        for ms in M_S_LIST:
            v = np.array([x[1] for x in vals if x[0] == ms])
            trend[energy][str(ms)] = float(v.mean())
            print('  %3s TeV  m_S=%5.2f  K=%.5f  (rms over realisations %.5f)' % (energy, ms, v.mean(), v.std(ddof=1)))
    except Exception as ex:
        print('  %s TeV: %s' % (energy, ex))
json.dump({'K': results, 'trend': trend}, open('runs/bc5_lhcb_downstream_v3/work/review_check/foresee/results_K.json','w'), indent=1)
print('\nwrote results_K.json')
