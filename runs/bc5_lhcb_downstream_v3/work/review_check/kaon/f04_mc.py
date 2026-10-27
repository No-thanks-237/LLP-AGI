
# -*- coding: utf-8 -*-
"""F04 audit: neutral-kaon flight in S -> K0 K0bar (LHCb downstream). READ-ONLY audit script.

(a) = visible KK rate as the production code computes it: all charged daughters born at the
      S decay vertex (on the S flight line, evaluated on the z grid of the acceptance table).
(b) = same S, same final state, but each neutral kaon's flight sampled from Exp(beta*gamma*c*tau)
      (c*tau_KS = 2.684 cm, c*tau_KL = 15.337 m); daughters transported from the true decay
      point with the real kernel. A kaon decaying beyond the tracking volume (z_dec > frontZ =
      7.7 m, or outside the SciFi envelope) yields no reconstructable pair and is discarded.

Outputs f04_m<m>.npz next to this file.
"""
import sys, os, math, json, time, random, csv
import numpy as np

RUN = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, RUN+'/work/analysis'); sys.path.insert(0, RUN+'/work/response'); sys.path.insert(0, RUN+'/work/widths')
import bc5_lib as L, bc5_production as P, run_scan as RS
import downstream_ctypes as dk

M_K0 = L.M_K0L
CTAU_KS = 0.026844
CTAU_KL = L.sw.HBAR_C / P.GAM_KL
ZG = np.array(L.ZGRID)
_track = dk.transport_track
_birth = dk.passes_birth_kinematics
MINE, MINPT = dk.BIRTH_CUTS['minE'], dk.BIRTH_CUTS['minPT']
FRONTZ = dk.DEFAULT_CONFIG['frontZ']


def _split_kk(ms, rng):
    if rng.random() < 0.5:
        d = L.two_body(ms, L.M_K, L.M_K, rng)
        if d is None: return None
        return ('charg', [(L.M_K, d[0][1], d[0][2], d[0][3], d[0][4], +1),
                          (L.M_K, d[1][1], d[1][2], d[1][3], d[1][4], -1)], None)
    d = L.two_body(ms, M_K0, M_K0, rng)
    if d is None: return None
    is_short = rng.random() < 0.5
    drest = []; kaon_of = []
    for i in (0, 1):
        for (m, e, px, py, pz, ch) in L._K0_charged(M_K0, is_short, rng):
            b = L.boost([(m, e, px, py, pz)], d[i][1], d[i][2:5])[0]
            drest.append((b[0], b[1], b[2], b[3], b[4], ch)); kaon_of.append(i)
    return ('KS' if is_short else 'KL', drest, (d, is_short, kaon_of))


def _accept(dau, births):
    good = []
    for (px, py, pz, e, ch, q), (x, y, z) in zip(dau, births):
        if not _birth(px, py, pz, e, MINE, MINPT): continue
        if _track(x, y, z, px, py, pz, e, ch)[0]: good.append(q)
    if len(good) < 2: return 0.0
    return 1.0 if (any(q > 0 for q in good) and any(q < 0 for q in good)) else 0.0


def s_arrays(ms, s, rng, ndecay=8, ndecay_all=6):
    tag, e, px, py, pz, x, y, z = s
    pmod = math.sqrt(px*px+py*py+pz*pz)
    if pz <= 0 or pmod <= abs(pz): return None
    eta = 0.5*math.log((pmod+pz)/max(pmod-pz, 1e-300))
    if not (L.ETA_S_MIN < eta < L.ETA_S_MAX): return None
    x0, y0, z0 = x*1e-3, y*1e-3, z*1e-3
    dirx, diry = px/pz, py/pz
    nz = len(ZG)
    acc = {k: np.zeros(nz) for k in ('C', 'KS', 'KL')}   # prompt treatment, sums
    new = {k: np.zeros(nz) for k in ('C', 'KS', 'KL')}   # flight treatment, sums
    cnt = {'C': 0, 'KS': 0, 'KL': 0}
    nkk = 0; fl = []
    for _ in range(ndecay):
        kk = _split_kk(ms, rng)
        if kk is None: continue
        nkk += 1
        kind, drest, kinfo = kk
        lab5 = L.boost([(a[0], a[1], a[2], a[3], a[4]) for a in drest], e, (px, py, pz))
        dau = [(b[2], b[3], b[4], b[1], c[5], c[5]) for b, c in zip(lab5, drest)]
        if kind == 'charg':
            cnt['C'] += 1
            for k, zd in enumerate(ZG):
                bx = x0 + (zd-z0)*dirx; by = y0 + (zd-z0)*diry
                v = _accept(dau, [(bx, by, zd)]*len(dau))
                acc['C'][k] += v; new['C'][k] += v
        else:
            cnt[kind] += 1
            d4, is_short, kaon_of = kinfo
            ctau = CTAU_KS if is_short else CTAU_KL
            disp = []
            for i in (0, 1):
                k4 = d4[i]
                kl = L.boost([k4], e, (px, py, pz))[0]
                kp = math.sqrt(kl[2]**2 + kl[3]**2 + kl[4]**2)
                lam = (kp/k4[0])*ctau
                u = -math.log(max(rng.random(), 1e-12))
                disp.append((lam*u*kl[2]/kp, lam*u*kl[3]/kp, lam*u*kl[4]/kp))
                fl.append((1.0 if is_short else 0.0, lam, kp/k4[0], kl[4]/kp))
            dob = [disp[j] for j in kaon_of]
            for k, zd in enumerate(ZG):
                bx = x0 + (zd-z0)*dirx; by = y0 + (zd-z0)*diry
                acc[kind][k] += _accept(dau, [(bx, by, zd)]*len(dau))
                births = [(bx+dd[0], by+dd[1], zd+dd[2]) for dd in dob]
                new[kind][k] += _accept(dau, births)
    nall = 0
    Aa = np.zeros(nz)
    br, _ = L.channel_branching(ms)
    names = list(br.keys()); cum = np.cumsum([br[n] for n in names])
    for _ in range(ndecay_all):
        u = rng.random(); idx = int(np.searchsorted(cum, u)); idx = min(idx, len(names)-1)
        rest = L.s_decay_products(ms, names[idx], rng)
        if not rest: continue
        nall += 1
        lab5 = L.boost([(a[0], a[1], a[2], a[3], a[4]) for a in rest], e, (px, py, pz))
        dd = [(b[2], b[3], b[4], b[1], a[5], a[5]) for b, a in zip(lab5, rest)]
        for k, zd in enumerate(ZG):
            bx = x0 + (zd-z0)*dirx; by = y0 + (zd-z0)*diry
            Aa[k] += _accept(dd, [(bx, by, zd)]*len(dd))
    if nall: Aa /= nall
    res = {}
    for k in ('C', 'KS', 'KL'):
        if cnt[k]:
            res['Ac'+k] = acc[k]/cnt[k]
            res['An'+k] = new[k]/cnt[k]
        else:
            res['Ac'+k] = np.zeros(nz); res['An'+k] = np.zeros(nz)
        res['n'+k] = cnt[k]
    # total KK arrays (mixture over the three sub-channels)
    Ac = (acc['C']+acc['KS']+acc['KL'])/nkk if nkk else np.zeros(nz)
    An = (new['C']+new['KS']+new['KL'])/nkk if nkk else np.zeros(nz)
    return dict(Ac=Ac, An=An, Aa=Aa, eta=eta, pz=pz, z0=z0, nkk=nkk,
                flight=np.array(fl) if fl else np.zeros((0, 4)), **res)


# ------------------------- eta-conditioned S sampling through the real model ----------
_B = _H = None
def _load():
    global _B, _H
    _B = RS.load_sample(RS.SAMPLE_B)
    _H = RS.load_sample(RS.SAMPLE_H)


def _sel(ms, e, px, py, pz):
    if pz <= 0: return False
    pmod = math.sqrt(px*px+py*py+pz*pz)
    if pmod <= abs(pz): return False
    eta = 0.5*math.log((pmod+pz)/max(pmod-pz, 1e-300))
    return L.ETA_S_MIN < eta < L.ETA_S_MAX


def _disp(m, p3, px, py, pz, ctau, rng):
    d = (p3/m)*ctau*(-math.log(max(rng.random(), 1e-12)))*1e3
    return d*px/p3, d*py/p3, d*pz/p3


def sample_S_B(ms, nsel, rng):
    sel = np.isin(_B['pdg'], [521, -521, 511, -511, 531, -531]); idxs = np.where(sel)[0]
    brs = {bm: P.br_B_to_S(ms, bm) for bm in ('B+', 'B0', 'Bs')}
    pack = {}
    for bm, brd in brs.items():
        nm = [k for k, v in brd.items() if v > 0] or ['K']
        w = np.array([max(brd[k], 0.0) for k in nm]); w = w/w.sum() if w.sum() > 0 else np.ones(len(nm))/len(nm)
        pack[bm] = (nm, np.cumsum(w))
    ctau = {'B+': L.sw.HBAR_C/P.GAM_BP, 'B0': L.sw.HBAR_C/P.GAM_B0, 'Bs': L.sw.HBAR_C/P.GAM_BS}
    out = []; nd = 0
    while len(out) < nsel and nd < 400*nsel:
        nd += 1
        i = int(idxs[rng.randrange(len(idxs))])
        bm = {521:'B+', -521:'B+', 511:'B0', -511:'B0', 531:'Bs', -531:'Bs'}[int(_B['pdg'][i])]
        e = _B['e'][i]; px = _B['px'][i]; py = _B['py'][i]; pz = _B['pz'][i]
        m = math.sqrt(max(e*e-px*px-py*py-pz*pz, 1e-9)); p3 = math.sqrt(px*px+py*py+pz*pz)
        if p3 <= 0: continue
        nm, cum = pack[bm]
        j = min(int(np.searchsorted(cum, rng.random())), len(nm)-1)
        d = L.two_body(m, P.B_CHANNEL_MASS[nm[j]], ms, rng)
        if d is None: continue
        s4 = L.boost([d[1]], e, (px, py, pz))[0]
        if not _sel(ms, s4[1], s4[2], s4[3], s4[4]): continue
        dx, dy, dz = _disp(m, p3, px, py, pz, ctau[bm], rng)
        out.append(('B', s4[1], s4[2], s4[3], s4[4], _B['x'][i]+dx, _B['y'][i]+dy, _B['z'][i]+dz))
    return out, nd


def sample_S_Bquartic(ms, nsel, rng):
    sel = np.isin(_B['pdg'], [521, -521, 511, -511, 531, -531]); idxs = np.where(sel)[0]
    ctau = {'B+': L.sw.HBAR_C/P.GAM_BP, 'B0': L.sw.HBAR_C/P.GAM_B0, 'Bs': L.sw.HBAR_C/P.GAM_BS}
    out = []; nd = 0
    while len(out) < nsel and nd < 600*nsel:
        nd += 1
        i = int(idxs[rng.randrange(len(idxs))])
        bm = {521:'B+', -521:'B+', 511:'B0', -511:'B0', 531:'Bs', -531:'Bs'}[int(_B['pdg'][i])]
        qm = [q for q in P.quartic_channels(ms) if q[0] == bm]
        if not qm: continue
        e = _B['e'][i]; px = _B['px'][i]; py = _B['py'][i]; pz = _B['pz'][i]
        m = math.sqrt(max(e*e-px*px-py*py-pz*pz, 1e-9)); p3 = math.sqrt(px*px+py*py+pz*pz)
        if p3 <= 0: continue
        wq = sum(q[3] for q in qm); u = rng.random(); c = 0.0; selc = qm[-1]
        for q in qm:
            c += q[3]/wq
            if u <= c: selc = q; break
        d = L.two_body(m, ms, ms, rng) if selc[1] == '2b' else None
        if d is None and selc[1] != '2b':
            if m < 2*ms + selc[2]: continue
            d = L.nbody(m, [ms, ms, selc[2]], rng)
            if d is None: continue
            d = [d[0], d[1]]
        if d is None: continue
        dx, dy, dz = _disp(m, p3, px, py, pz, ctau[bm], rng)
        x0, y0, z0 = _B['x'][i]+dx, _B['y'][i]+dy, _B['z'][i]+dz
        for s4 in d:
            sl = L.boost([s4], e, (px, py, pz))[0]
            if not _sel(ms, sl[1], sl[2], sl[3], sl[4]): continue
            out.append(('Bq', sl[1], sl[2], sl[3], sl[4], x0, y0, z0))
    return out, nd


def sample_S_H(ms, nsel, rng):
    ev, first = np.unique(_H['ev'], return_index=True); idxs = np.sort(first)
    out = []; nd = 0
    while len(out) < nsel and nd < 60*nsel:
        nd += 1
        i = int(idxs[rng.randrange(len(idxs))])
        d = L.two_body(P.M_H, ms, ms, rng)
        if d is None: continue
        pair = L.boost(d, _H['e'][i], (_H['px'][i], _H['py'][i], _H['pz'][i]))
        for s4 in pair:
            if not _sel(ms, s4[1], s4[2], s4[3], s4[4]): continue
            out.append(('H', s4[1], s4[2], s4[3], s4[4], _H['x'][i], _H['y'][i], _H['z'][i]))
    return out, nd


def main():
    ms = float(sys.argv[1]); nsel = int(sys.argv[2]) if len(sys.argv) > 2 else 500
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 12345
    rng = random.Random(seed); t0 = time.time(); _load()
    scan = json.load(open(RUN+'/outputs/data/signal_scan.json'))
    thetas = np.array(scan['theta2'])
    rows = list(csv.DictReader(open(RUN+'/outputs/data/exclusion_boundary.csv')))
    bnd = {}
    for r in rows:
        if abs(float(r['lumi_fb'])-25.0) < 1e-9: bnd[round(float(r['mass_GeV']), 6)] = r
    bm_ = np.array(sorted(bnd))
    def interp_edge(col):
        y = np.array([float(bnd[m][col]) if bnd[m][col] not in ('', 'nan') else np.nan for m in bm_])
        ok = np.isfinite(y) & (y > 0)
        return float(np.exp(np.interp(math.log(ms), np.log(bm_[ok]), np.log(y[ok])))) if ok.sum() >= 2 else float('nan')
    edges = {'long': interp_edge('theta2_long_lived_edge'), 'short': interp_edge('theta2_short_lived_edge')}
    sc = scan['signal_per_fb_per_channel']; skeys = list(sc); smass = np.array([float(k) for k in skeys])
    o = np.argsort(smass); smass = smass[o]; skeys = [skeys[i] for i in o]
    def scan_chan(t2, c):
        vals = [float(np.exp(np.interp(math.log(t2), np.log(thetas), np.log(np.maximum(
                    np.array([r.get(c, 0.0) for r in sc[k]]), 1e-300))))) for k in skeys]
        return float(np.exp(np.interp(math.log(ms), np.log(smass), np.log(np.maximum(np.array(vals), 1e-300)))))
    sB = {nm: scan_chan(t2, 'B') for nm, t2 in edges.items() if np.isfinite(t2)}
    sH = {nm: scan_chan(t2, 'H') for nm, t2 in edges.items() if np.isfinite(t2)}
    print('m=%.4f edges %s' % (ms, {k: '%.3e' % v for k, v in edges.items()}))
    samples = {}
    S, nd = sample_S_B(ms, nsel, rng); samples['B'] = S; print('  Bmix %d/%d draws' % (len(S), nd))
    S, nd = sample_S_Bquartic(ms, max(100, nsel//4), rng); samples['Bq'] = S; print('  Bq %d/%d' % (len(S), nd))
    S, nd = sample_S_H(ms, nsel, rng); samples['H'] = S; print('  H %d/%d' % (len(S), nd))
    out = {}
    for name, S in samples.items():
        keys = ['Ac', 'An', 'Aa', 'AcC', 'AnC', 'AcKS', 'AnKS', 'AcKL', 'AnKL']
        vals = {k: [] for k in keys}; meta = []; fl = []
        for s in S:
            r = s_arrays(ms, s, rng)
            if r is None: continue
            for k in keys: vals[k].append(r[k])
            meta.append((r['eta'], r['pz'], r['z0'], r['nC'], r['nKS'], r['nKL']))
            if len(r['flight']): fl.append(r['flight'])
        out[name] = {k: np.array(v) for k, v in vals.items()}
        out[name]['meta'] = np.array(meta)
        out[name]['flight'] = np.concatenate(fl) if fl else np.zeros((0, 4))
        print('  %s: %d S (%.1fs)' % (name, len(vals['Ac']), time.time()-t0))
    np.savez_compressed(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'f04_m%.4f.npz' % ms),
                        mass=ms, ZG=ZG, edges=np.array([edges['long'], edges['short']]),
                        sB=np.array([sB.get('long', 0), sB.get('short', 0)]),
                        sH=np.array([sH.get('long', 0), sH.get('short', 0)]),
                        **{name+k: out[name][k] for name in ('B', 'Bq', 'H')
                           for k in ('Ac', 'An', 'Aa', 'AcC', 'AnC', 'AcKS', 'AnKS', 'AcKL', 'AnKL', 'meta', 'flight')})
    print('wrote f04_m%.4f.npz (%.1fs)' % (ms, time.time()-t0))

if __name__ == '__main__':
    main()
