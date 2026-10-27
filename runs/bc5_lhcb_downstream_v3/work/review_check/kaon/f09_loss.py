
# -*- coding: utf-8 -*-
"""F09-1 direct estimate: fraction of the K->pi S / K->SS signal lost by the 1.5<eta_K<5.5
pre-cut.  Uses the full (no-cut) minbias sample  work/magnus/outputs/nocut2/out/minbias_light.txt.gz
(200k events, same Pythia settings; its in-window part reproduces the production sample)."""
import sys, os, math, gzip, json
import numpy as np
R = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, R+'/work/analysis'); sys.path.insert(0, R+'/work/response'); sys.path.insert(0, R+'/work/widths')
import bc5_lib as L, bc5_production as P
NOCUT = R+'/work/magnus/outputs/nocut2/out/minbias_light.txt.gz'

def load(path):
    ev, pdg, px, py, pz, e = [], [], [], [], [], []
    with gzip.open(path,'rt') as f:
        for line in f:
            if not line or line[0]=='#': continue
            v = line.split()
            if len(v) < 10: continue
            ev.append(int(v[0])); pdg.append(int(v[1]))
            px.append(float(v[2])); py.append(float(v[3])); pz.append(float(v[4])); e.append(float(v[5]))
    a = {k: np.asarray(v, dtype=float if k not in ('ev','pdg') else np.int64)
         for k,v in dict(ev=ev,pdg=pdg,px=px,py=py,pz=pz,e=e).items()}
    return a

d = load(NOCUT)
px,py,pz,e = d['px'],d['py'],d['pz'],d['e']
p3 = np.sqrt(px*px+py*py+pz*pz); m = np.sqrt(np.maximum(e*e-p3*p3,0))
eta = 0.5*np.log((p3+pz)/np.maximum(p3-pz,1e-300))
inw = (eta>1.5)&(eta<5.5)
print('no-cut sample rows %d, in-window %d' % (len(pdg), inw.sum()) if False else '')

M_DIR = 12
def prob_eta(mK, mP, ms, sel, chunk=100000):
    """P(2<eta_S<5 and pz_S>0) for K->P S, for the parents selected by 'sel'."""
    out = np.zeros(len(px)); n = 0
    pstar = (math.sqrt(max((mK*mK-(mP+ms)**2)*(mK*mK-(mP-ms)**2),0.0))/(2*mK)) if mK > mP+ms else 0.0
    Estar = math.sqrt(pstar*pstar+ms*ms)
    for lo in range(0, len(px), chunk):
        hi = min(lo+chunk, len(px)); sl = slice(lo,hi)
        E = e[sl]; P3 = p3[sl]; M = m[sl]
        beta = np.clip(P3/np.maximum(E,1e-30), 0, 1-1e-12); gam = E/np.maximum(M,1e-30)
        ux, uy, uz = px[sl]/P3, py[sl]/P3, pz[sl]/P3
        # orthonormal basis
        ax = np.where(np.abs(uz)<0.9, 0.0, 1.0); ay = np.where(np.abs(uz)<0.9, 0.0, 0.0); az = np.where(np.abs(uz)<0.9, 1.0, 0.0)
        ex = np.cross(np.stack([ux,uy,uz],1), np.stack([ax,ay,az],1))
        exn = np.linalg.norm(ex,axis=1); ex = ex/exn[:,None]
        ey = np.cross(np.stack([ux,uy,uz],1), ex)
        ct = np.random.default_rng(7).uniform(-1,1,(hi-lo, M_DIR))
        ph = np.random.default_rng(11).uniform(0,2*math.pi,(hi-lo, M_DIR))
        st = np.sqrt(np.maximum(1-ct*ct,0))
        ppar = gam[:,None]*(pstar*ct + beta[:,None]*Estar)
        pper = pstar*st
        vx = ppar*ux[:,None] + pper*(np.cos(ph)*ex[:,0:1] + np.sin(ph)*ey[:,0:1])
        vy = ppar*uy[:,None] + pper*(np.cos(ph)*ex[:,1:2] + np.sin(ph)*ey[:,1:2])
        vz = ppar*uz[:,None] + pper*(np.cos(ph)*ex[:,2:3] + np.sin(ph)*ey[:,2:3])
        pm = np.sqrt(vx*vx+vy*vy+vz*vz)
        ok = (vz>0) & (pm-np.abs(vz) > 0)
        eS = np.where(ok, 0.5*np.log(np.maximum(pm+vz,1e-300)/np.maximum(pm-vz,1e-300)), -99)
        keep = ok & (eS>2.0) & (eS<5.0)
        out[sl] = keep.mean(axis=1)
    return out

def loss_for(ms, mK=0.493677, mP=0.13957039, tag='K+'):
    pe = prob_eta(mK, mP, ms, None)
    for pdg in (321,-321,311,-311):
        s = d['pdg']==pdg
        if s.sum()==0: continue
    return pe

res = {}
rng = np.random.default_rng(7)
for ms in (0.1, 0.2, 0.3):
    row = {}
    for tag, (pdg_ids, mK, mP) in {
            'K+': ([321,-321], 0.493677, 0.13957039),
            'K0L': ([311,-311], 0.497611, 0.1349768)}.items():
        pe = prob_eta(mK, mP, ms, None)
        sel_species = np.isin(d['pdg'], pdg_ids)
        sum_all = pe[sel_species].sum(); sum_in = pe[sel_species & inw].sum()
        row[tag] = dict(n=int(sel_species.sum()), n_in=int((sel_species&inw).sum()),
                        sumP_all=float(sum_all), sumP_in=float(sum_in),
                        loss=float(1-sum_in/sum_all) if sum_all>0 else float('nan'),
                        br=float(P.br_K_to_S(ms, tag)))
    # combine species with the production branching ratio weights
    w1 = row['K+']['br']*row['K+']['sumP_all']; w2 = row['K0L']['br']*row['K0L']['sumP_all']
    w1i = row['K+']['br']*row['K+']['sumP_in']; w2i = row['K0L']['br']*row['K0L']['sumP_in']
    row['combined_loss'] = float(1-(w1i+w2i)/(w1+w2)) if (w1+w2)>0 else float('nan')
    res['%.2f'%ms] = row
    print('\nm_S=%.2f'%ms)
    for tag in ('K+','K0L'):
        r = row[tag]
        print('  %-4s n=%7d (in %7d)  <P_eta>_all=%.5f  <P_eta>_in=%.5f  loss=%.5f  br=%.3e' % (
              tag, r['n'], r['n_in'], r['sumP_all']/r['n'], r['sumP_in']/r['n_in'], r['loss'], r['br']))
    print('  combined loss of the K->piS channel = %.5f' % row['combined_loss'])

# kinematic reach: range of eta_S for a K at eta_K, m_S=0.1
def reach(mK=0.493677, mP=0.13957039, ms=0.1, etaK=1.4):
    pstar = math.sqrt(max((mK*mK-(mP+ms)**2)*(mK*mK-(mP-ms)**2),0.0))/(2*mK)
    Estar = math.sqrt(pstar*pstar+ms*ms)
    thK = 2*math.atan(math.exp(-etaK))
    out = []
    for p in (0.5,1,2,3,5,10,30,100):
        E = math.sqrt(p*p+mK*mK); beta=p/E; gam=E/mK
        lo, hi = 1e9, -1e9
        for t in np.linspace(0,math.pi,400):
            ct, st = math.cos(t), math.sin(t)
            for ph in np.linspace(0,2*math.pi,72,endpoint=False):
                ppar = gam*(pstar*ct + beta*Estar); pper = pstar*st
                # K along (sin thK, 0, cos thK); e1=(cos thK,0,-sin thK), e2=(0,1,0)
                vx = ppar*math.sin(thK) + pper*math.cos(ph)*math.cos(thK)
                vy = pper*math.sin(ph)
                vz = ppar*math.cos(thK) - pper*math.cos(ph)*math.sin(thK)
                pm = math.sqrt(vx*vx+vy*vy+vz*vz)
                if vz<=0 or pm-abs(vz)<=0: continue
                eS = 0.5*math.log((pm+vz)/(pm-vz))
                lo=min(lo,eS); hi=max(hi,eS)
        out.append((p, lo, hi))
    return out
print('\nkinematic reach for K->pi S, m_S=0.1, eta_K=1.4 (S eta range over all emission angles):')
for p,lo,hi in reach():
    print('   p_K=%.1f GeV : eta_S in [%.2f, %.2f]  (can reach eta>2: %s)' % (p,lo,hi,hi>2.0))
