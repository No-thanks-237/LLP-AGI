
# -*- coding: utf-8 -*-
"""F09-2 / F09-3 verification with the production probability algebra."""
import sys, os, math, json, random
import numpy as np
R = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, R+'/work/analysis'); sys.path.insert(0, R+'/work/response'); sys.path.insert(0, R+'/work/widths')
import bc5_lib as L, bc5_production as P, run_scan as RS

# ---------- F09-3: the production channel_efficiency ORs K->SS and K->SSpi as one history
br_SS, br_SSpi = 4.4e-13, 6.6e-15      # Boiarska BR(alpha=1) before alpha^2
p = 0.1
tab = {'chan': np.array(['K']*4), 'ev': np.array([0]*4), 'par': np.array([0]*4),
       'hid': np.array([1]*4), 'kq': np.array([0.0]*4),
       'w': np.array([br_SS, br_SS, br_SSpi, br_SSpi]),
       'pz': np.array([1.0]), 'z': np.array([0.0]), 'A': np.zeros((1,13))}
# channel_efficiency takes p per row; emulate by giving p=0.1 to every row
pvec = np.array([p]*4)
eps_wrong, _ = RS.channel_efficiency({'chan': tab['chan'], 'ev': tab['ev'], 'par': tab['par'],
    'hid': tab['hid'], 'kq': tab['kq'], 'w': tab['w'], 'counts': {'K': [1, 1]}}, pvec, 'K', 1, 1.0)
correct = 0.4*(br_SS*(1-(1-p)**2) + br_SSpi*(1-(1-p)**2))   # exclusive histories ADD
print('F09-3: eps(wrong hid merge)/0.4 = %.6e ; correct = %.6e ; ratio = %.4f' % (eps_wrong/0.4, correct, eps_wrong/(0.4*correct)))
print('       analytic max ratio for p->0 is 2*br_SS/(br_SS+br_SSpi) = %.4f' % (2*br_SS/(br_SS+br_SSpi)))

# ---------- F09-2: K+ vs K0L split of the K->pi S mixing signal (production machinery)
smp = RS.load_sample(RS.SAMPLE_K)
p3all = np.sqrt(smp['px']**2+smp['py']**2+smp['pz']**2)
etaall = 0.5*np.log((p3all+smp['pz'])/np.maximum(p3all-smp['pz'],1e-300))
ms = 0.1; t2 = 1.390e-6
ctau = L.sw.HBAR_C/(L.sw.gamma_total(ms)*t2)
rng = random.Random(5)
def meanp(pdg_ids, which, npar=400):
    sel = np.isin(smp['pdg'], pdg_ids); idxs = np.where(sel)[0]
    tot = 0.0; n = 0
    for _ in range(npar):
        i = int(idxs[rng.randrange(len(idxs))])
        mes = {'kind':'K','kaon':which,'x':smp['x'][i],'y':smp['y'][i],'z':smp['z'][i],
               'e':smp['e'][i],'px':smp['px'][i],'py':smp['py'][i],'pz':smp['pz'][i]}
        for s in RS.sample_meson_decay_S(ms, mes, rng) or []:
            e,px,py,pz,x,y,z = s
            pmod = math.sqrt(px*px+py*py+pz*pz)
            if pz<=0 or pmod<=abs(pz): continue
            eta = 0.5*math.log((pmod+pz)/(pmod-pz))
            if not (L.ETA_S_MIN<eta<L.ETA_S_MAX): continue
            A = L.acceptance_profile(ms, (e,px,py,pz), (x,y,z), rng, ndecay=4)
            z0 = z*1e-3; alpha = ms/(ctau*pz)
            w = np.zeros(len(L.ZGRID)); mm = L.ZGRID>z0
            w[mm] = alpha*np.exp(-(L.ZGRID[mm]-z0)*alpha)
            tot += float(np.trapezoid(w*A, L.ZGRID)); n += 1
    return tot/npar, npar
for tag, ids, which, br in (('K+', [321,-321], 'K+', P.br_K_to_S(ms,'K+')),
                            ('K0L', [311,-311], 'K0L', P.br_K_to_S(ms,'K0L'))):
    mp, n = meanp(ids, which)
    print('F09-2: m_S=0.1 %s  <p>=%.4e  n_tot=%d  br/theta2=%.3e  rate~%.3e' % (tag, mp, int(np.isin(smp['pdg'],ids).sum()), br, br*mp*int(np.isin(smp['pdg'],ids).sum())))
print('   (the mixing rate per parent is br*<p>; K0L carries the full K_L BR although a K0')
print('    flavour state has only |K_L|^2 = 1/2 K_L component -> factor 2 too large)')
