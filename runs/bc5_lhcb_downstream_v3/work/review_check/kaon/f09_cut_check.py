
# -*- coding: utf-8 -*-
"""F09-1: verify the 1.5<eta_K<5.5 pre-cut of the production kaon sample and estimate the
signal lost by it (K -> pi S and K -> SS quartic).  READ-ONLY."""
import sys, os, math, gzip, json
import numpy as np
R = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, R+'/work/analysis'); sys.path.insert(0, R+'/work/response'); sys.path.insert(0, R+'/work/widths')
import bc5_lib as L, bc5_production as P
PROD = R+'/work/magnus/samples/minbias_light.txt'
NOCUT = R+'/work/magnus/outputs/nocut2/out/minbias_light.txt.gz'

def load(path, cols=('ev','pdg','px','py','pz','e')):
    d = {c: [] for c in cols}
    op = gzip.open if path.endswith('.gz') else open
    with op(path,'rt') as f:
        for line in f:
            if not line or line[0]=='#': continue
            v = line.split()
            if len(v) < 10: continue
            d['ev'].append(int(v[0])); d['pdg'].append(int(v[1]))
            d['px'].append(float(v[2])); d['py'].append(float(v[3])); d['pz'].append(float(v[4])); d['e'].append(float(v[5]))
    return {k: np.asarray(v) for k,v in d.items()}

def kin(d):
    p3 = np.sqrt(d['px']**2+d['py']**2+d['pz']**2)
    eta = 0.5*np.log((p3+d['pz'])/np.maximum(p3-d['pz'],1e-300))
    pt = np.hypot(d['px'],d['py'])
    return p3, eta, pt

prod = load(PROD); nc = load(NOCUT)
for nm, d in (('PROD', prod), ('NOCUT', nc)):
    p3, eta, pt = kin(d)
    print('%s: rows=%d events=%d  eta=[%.5f,%.5f]' % (nm, len(d['pdg']), len(np.unique(d['ev'])), eta.min(), eta.max()))
    for pdg in sorted(set(d['pdg'].tolist())):
        s = d['pdg']==pdg
        print('   pdg %5d n=%7d eta=[%.3f,%.3f] med pT=%.3f' % (pdg, s.sum(), eta[s].min(), eta[s].max(), np.median(pt[s])))
# in-window selection of the no-cut sample
p3, eta, pt = kin(nc)
inw = (eta>1.5)&(eta<5.5)
print('\nNOCUT in-window rows: %d (%.1f%%); outside: %d' % (inw.sum(), 100*inw.mean(), (~inw).sum()))
for pdg in (321,-321,311,-311):
    s = nc['pdg']==pdg
    print('  pdg %5d: total %7d  in-window %7d (%.1f%%)' % (pdg, s.sum(), (s&inw).sum(), 100*(s&inw).sum()/s.sum()))
# compare in-window spectra with the production sample (same n_gen=200k)
pp3, peta, ppt = kin(prod)
for q in (10,25,50,75,90):
    print('  eta quantiles  PROD %s | NOCUT-in %s' % (np.percentile(peta,[q])[0].round(3), np.percentile(eta[inw],[q])[0].round(3)))
for q in (10,25,50,75,90):
    print('  pT  quantiles  PROD %s | NOCUT-in %s' % (np.percentile(ppt,[q])[0].round(3), np.percentile(pt[inw],[q])[0].round(3)))
# 2-body kinematics: maximum |eta_S - eta_K| for K->pi S
def max_deta(mK, ms, p):
    if p <= 0: return 0.0
    E = math.sqrt(p*p+mK*mK); beta = p/E; g = E/mK
    ps = math.sqrt(max((mK*mK-(mK if False else 0.13957039+ms)**2)*(mK*mK-(0.13957039-ms)**2),0.0))/(2*mK) if mK > 0.13957039+ms else 0.0
    Es = math.sqrt(ps*ps+ms*ms)
    best = 0.0
    for t in np.linspace(0, math.pi, 721):
        ct = math.cos(t); st = math.sin(t)
        ppar = g*(ps*ct+beta*Es); pper = ps*st
        # eta relative to the K direction (K along z here)
        pl = math.sqrt(max(ppar*ppar+pper*pper,1e-30))
        if pl-abs(ppar) <= 0: continue
        etaS = 0.5*math.log((pl+ppar)/max(pl-ppar,1e-300))
        best = max(best, abs(etaS))
    return best
print('\nmax |eta_S-eta_K| for K->pi S (GeV):')
for ms in (0.1, 0.3):
    row = []
    for p in (0.5,1.0,2.0,5.0,10.0,30.0,100.0):
        row.append('p=%.0f:%.2f' % (p, max_deta(0.493677, ms, p)))
    print('  m_S=%.1f  %s' % (ms, '  '.join(row)))
