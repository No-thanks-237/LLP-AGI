# -*- coding: utf-8 -*-
"""FORESEE h(125) spectrum vs our Pythia8 Higgs sample.

Reads:
  * our sample  : runs/bc5_lhcb_downstream_v3/work/magnus/samples/higgs.txt.gz  (|status|==62)
  * FORESEE     : temp/foresee/files/hadrons/14TeV/Pythia8/Pythia8_14TeV_25.txt
Writes JSON with all numbers to work/review_check/foresee/results.json
"""
import sys, os, math, json
import numpy as np

sys.path.insert(0, 'runs/bc5_lhcb_downstream_v3/work/analysis')
import run_scan as R
import bc5_lib as L
import bc5_production as P

M_S_LIST = [5.08, 13.02, 33.37, 53.42]
ETA_LO, ETA_HI = 2.0, 5.0

def eta_from_p(px, py, pz):
    pt = np.hypot(px, py)
    return np.arcsinh(pz / np.maximum(pt, 1e-300))

def two_body_scalars(mH, E, px, py, pz, mS, rng):
    """Isotropic h -> S S (p* = sqrt(mH^2/4 - mS^2)) boosted to the lab.
    Vectorised, algebraically identical to bc5_lib.two_body + bc5_lib.boost."""
    N = len(mH)
    pstar2 = mH * mH / 4.0 - mS * mS
    ok = pstar2 > 0.0
    pstar = np.sqrt(np.maximum(pstar2, 0.0))
    # isotropic direction in the Higgs rest frame
    u = rng.normal(size=(N, 3))
    u /= np.maximum(np.linalg.norm(u, axis=1, keepdims=True), 1e-300)
    d1 = pstar[:, None] * u
    d2 = -d1
    Estar = mH / 2.0
    # boost velocity of the Higgs
    bvec = np.stack([px, py, pz], axis=1) / np.maximum(E, 1e-300)[:, None]
    b2 = np.sum(bvec * bvec, axis=1)
    gamma = E / mH
    out = []
    for d, Ed in ((d1, Estar), (d2, Estar)):
        bp = np.sum(d * bvec, axis=1)
        k = (gamma - 1.0) * bp / np.maximum(b2, 1e-300) + gamma * Ed
        dp = d + k[:, None] * bvec
        out.append(dp)
    return out, ok

def response(dp, mS):
    """R = 1_{ETA_LO<eta<ETA_HI} * mS/pz  (0 if pz<=0 or eta outside)."""
    pt = np.hypot(dp[:, 0], dp[:, 1])
    pz = dp[:, 2]
    eta = np.arcsinh(pz / np.maximum(pt, 1e-300))
    good = (eta > ETA_LO) & (eta < ETA_HI) & (pz > 0.0) & np.isfinite(eta)
    r = np.zeros(len(dp))
    r[good] = mS / pz[good]
    return r, eta

def per_higgs_R(mH, E, px, py, pz, mS, rng):
    (d1, d2), ok = two_body_scalars(mH, E, px, py, pz, mS, rng)
    r1, e1 = response(d1, mS)
    r2, e2 = response(d2, mS)
    return 0.5 * (r1 + r2), 0.5 * (e1 + e2), (r1, r2), (e1, e2)

# ======================================================================================
# 0. validation of the vectorised decay/boost against the production code
# ======================================================================================
def validate(rng):
    import random
    rr = random.Random(12345)
    ms = 5.08
    bad = 0.0
    for k in range(20):
        mH = 125.25 + 0.001 * (k - 10)
        E = rr.uniform(126.0, 400.0); th = rr.uniform(0.05, 2.9); ph = rr.uniform(-3.1, 3.1)
        p = math.sqrt(E * E - mH * mH)
        px, py, pz = p * math.sin(th) * math.cos(ph), p * math.sin(th) * math.sin(ph), p * math.cos(th)
        parts = L.two_body(mH, ms, ms, rr)
        ref = L.boost(parts, E, (px, py, pz))
        (d1, d2), ok = two_body_scalars(np.array([mH]), np.array([E]), np.array([px]),
                                        np.array([py]), np.array([pz]), ms, rng)
        for a, b in ((d1[0], ref[0]), (d2[0], ref[1])):
            bad = max(bad, abs(a[0]-b[2]) + abs(a[1]-b[3]) + abs(a[2]-b[4]))
    return bad

# ======================================================================================
# 1. our sample
# ======================================================================================
smp = R.load_sample('higgs.txt.gz')
cand = np.where(np.abs(smp['status']) == 62)[0]
cand = cand[np.argsort(smp['ev'][cand], kind='stable')]
E0, px0, py0, pz0 = smp['e'][cand], smp['px'][cand], smp['py'][cand], smp['pz'][cand]
code = smp['pdg'][cand]
p0 = np.sqrt(px0*px0 + py0*py0 + pz0*pz0)
mH0 = np.sqrt(np.maximum(E0*E0 - p0*p0, 0.0))
etaH0 = eta_from_p(px0, py0, pz0)
ptH0 = np.hypot(px0, py0)

print('validating vectorised decay/boost against bc5_lib: max |dp| = %.3e' % validate(np.random.default_rng(1)))
print()
print('=== OUR SAMPLE (13 TeV Pythia8, status-62) ===')
print('Higgs rows            : %d   (one per event: %s)' % (len(cand), np.all(np.diff(smp['ev'][cand]) > 0)))
print('<m_H>                 : %.5f GeV   (rms %.2e)' % (mH0.mean(), mH0.std()))
print('<p>                   : %.2f GeV   <p_T> %.2f GeV   <E> %.2f GeV' % (p0.mean(), ptH0.mean(), E0.mean()))
print('<eta_H>               : %.4f    <|eta_H|> %.4f' % (etaH0.mean(), np.abs(etaH0).mean()))
print('frac 2<eta_H<5        : %.4f' % np.mean((etaH0 > 2) & (etaH0 < 5)))
print('eta_H pct [5,25,50,75,95]:', np.round(np.percentile(etaH0, [5,25,50,75,95]), 3))
print('p_T_H pct [5,25,50,75,95]:', np.round(np.percentile(ptH0, [5,25,50,75,95]), 2))
print('process code counts   :', dict(zip(*[a.tolist() for a in np.unique(code, return_counts=True)])))
grp = np.array([P.HIGGS_PROC_GROUP.get(int(c), 'other') for c in code])
ug, cg = np.unique(grp, return_counts=True)
print('process group counts  :', dict(zip(ug.tolist(), cg.tolist())))
print()

# process reweighting exactly as run_scan.py:501-509
wH = np.ones(len(cand))
for g_ in sorted(set(grp.tolist())):
    if g_ == 'other' or g_ not in P.SIGMA_H:
        continue
    m_ = grp == g_
    frac = float(m_.sum()) / float(len(cand))
    wH[m_] = P.SIGMA_H[g_] / (P.SIGMA_H_TOTAL * frac)
print('reweight wH: <w> %.4f  min %.3f max %.3f' % (wH.mean(), wH.min(), wH.max()))
print()

# --- dN/deta and dN/dpT in the FORWARD hemisphere (eta>0), for shape comparison
eta_bins = np.arange(0.0, 6.01, 0.5)
pt_bins = np.array([0, 2, 5, 8, 12, 20, 30, 40, 60, 80, 120, 200, 400, 1e4])
self = etaH0 > 0
h_eta_ours, _ = np.histogram(etaH0[self], bins=eta_bins, weights=wH[self])
h_pt_ours, _ = np.histogram(ptH0[self], bins=pt_bins, weights=wH[self])
# unweighted versions too
h_eta_ours_u, _ = np.histogram(etaH0[self], bins=eta_bins)
h_pt_ours_u, _ = np.histogram(ptH0[self], bins=pt_bins)

print('=== OUR dN/deta (forward hemisphere, per unit eta; shape normalised) ===')
area = np.sum(h_eta_ours)
for i in range(len(eta_bins)-1):
    print('  eta %.1f-%.1f  N=%9.1f  1/N dN/deta=%.4f' % (eta_bins[i], eta_bins[i+1],
          h_eta_ours[i], h_eta_ours[i]/area/(eta_bins[i+1]-eta_bins[i])))
print()
print('=== OUR dN/dpT (forward hemisphere, per GeV; shape normalised) ===')
area_p = np.sum(h_pt_ours)
for i in range(len(pt_bins)-1):
    print('  pT %6.1f-%7.1f  N=%9.1f  1/N dN/dpT=%.5f' % (pt_bins[i], pt_bins[i+1],
          h_pt_ours[i], h_pt_ours[i]/area_p/(pt_bins[i+1]-pt_bins[i])))

out = {'our': {'N': int(len(cand)), 'mH_mean': float(mH0.mean()), 'mH_std': float(mH0.std()),
               'p_mean': float(p0.mean()), 'pt_mean': float(ptH0.mean()),
               'eta_mean': float(etaH0.mean()), 'abseta_mean': float(np.abs(etaH0).mean()),
               'frac_eta_2_5': float(np.mean((etaH0>2)&(etaH0<5)))},
       'eta_bins': eta_bins.tolist(), 'pt_bins': pt_bins.tolist(),
       'h_eta_ours': h_eta_ours.tolist(), 'h_pt_ours': h_pt_ours.tolist(),
       'h_eta_ours_u': h_eta_ours_u.tolist(), 'h_pt_ours_u': h_pt_ours_u.tolist(),
       'group_counts': dict(zip(ug.tolist(), cg.tolist()))}

# ======================================================================================
# 2. FORESEE spectrum
# ======================================================================================
FS = 'temp/foresee/files/hadrons/14TeV/Pythia8/Pythia8_14TeV_25.txt'
raw = []
hdr = {}
for line in open(FS):
    if line.startswith('#'):
        if 'Normalization' in line:
            hdr['norm'] = float(line.split(':')[1].split()[0])
        continue
    raw.append([float(x) for x in line.split()])
raw = np.array(raw)
lth, lp, w = raw[:,0], raw[:,1], raw[:,2]
nz = w > 1e-6              # same cut as foresee.convert_list_to_momenta (src/foresee.py:274)
print()
print('=== FORESEE 14 TeV Pythia8 h(125) spectrum (%s) ===' % FS)
print('grid: %d x %d (log10 theta x log10 p), step 0.05 dex' % (len(np.unique(lth)), len(np.unique(lp))))
print('normalization header : %.6f pb  (forward hemisphere)' % hdr.get('norm', float('nan')))
print('sum of all bin weights: %.6f pb  (= header: %s)' % (w.sum(), abs(w.sum()-hdr['norm'])<1e-6))
print('non-zero bins (w>1e-6): %d / %d   (sum %.6f pb)' % (nz.sum(), len(w), w[nz].sum()))
print('theta range %.3e .. %.4f rad   eta range %.3f .. %.3f' % (
      10**lth[nz].min(), 10**lth[nz].max(), eta_from_p(0,0,0) if False else -math.log(math.tan(10**lth[nz].min()/2)),
      -math.log(math.tan(10**lth[nz].max()/2))))
print('p range p_min %.4f .. %.1f GeV' % (10**lp[nz].min(), 10**lp[nz].max()))
print('weight range %.3e .. %.4e pb/bin' % (w[nz].min(), w[nz].max()))

# forward-hemisphere eta / pT histograms
th_f = 10**lth[nz]; p_f = 10**lp[nz]
eta_f = -np.log(np.tan(th_f/2.0)); pt_f = p_f*np.sin(th_f)
h_eta_f, _ = np.histogram(eta_f, bins=eta_bins, weights=w[nz])
h_pt_f, _ = np.histogram(pt_f, bins=pt_bins, weights=w[nz])
area_f = np.sum(h_eta_f)
print()
print('=== FORESEE dN/deta (forward hemisphere, per unit eta; shape normalised) ===')
for i in range(len(eta_bins)-1):
    print('  eta %.1f-%.1f  w=%9.4f pb  1/N dN/deta=%.4f' % (eta_bins[i], eta_bins[i+1],
          h_eta_f[i], h_eta_f[i]/area_f/(eta_bins[i+1]-eta_bins[i])))
print()
print('=== FORESEE dN/dpT (forward hemisphere, per GeV; shape normalised) ===')
area_pf = np.sum(h_pt_f)
for i in range(len(pt_bins)-1):
    print('  pT %6.1f-%7.1f  w=%9.4f pb  1/N dN/dpT=%.5f' % (pt_bins[i], pt_bins[i+1],
          h_pt_f[i], h_pt_f[i]/area_pf/(pt_bins[i+1]-pt_bins[i])))

print()
print('=== ratio FORESEE/ours, dN/deta (both normalised over the same bins) ===')
for i in range(len(eta_bins)-1):
    a = h_eta_f[i]/area_f; b = h_eta_ours[i]/area
    print('  eta %.1f-%.1f  ratio %.3f +- %.3f' % (eta_bins[i], eta_bins[i+1], a/b,
          a/b*math.sqrt(1.0/max(h_eta_f[i],1e-12) + 1.0/max(h_eta_ours[i],1e-12))))
print()
print('=== ratio FORESEE/ours, dN/dpT ===')
for i in range(len(pt_bins)-1):
    a = h_pt_f[i]/area_pf; b = h_pt_ours[i]/area_p
    print('  pT %6.1f-%7.1f  ratio %.3f +- %.3f' % (pt_bins[i], pt_bins[i+1], a/b,
          a/b*math.sqrt(1.0/max(h_pt_f[i],1e-12) + 1.0/max(h_pt_ours[i],1e-12))))

out['foresee'] = {'norm_pb': float(hdr.get('norm', float('nan'))), 'sum_w': float(w.sum()),
                  'nz': int(nz.sum()), 'n_rows': int(len(w)),
                  'h_eta_f': h_eta_f.tolist(), 'h_pt_f': h_pt_f.tolist(),
                  'eta_min': float(eta_f.min()), 'eta_max': float(eta_f.max()),
                  'p_min': float(p_f.min()), 'p_max': float(p_f.max())}
json.dump(out, open('runs/bc5_lhcb_downstream_v3/work/review_check/foresee/results_part1.json','w'), indent=1)
print('\nwrote results_part1.json')
