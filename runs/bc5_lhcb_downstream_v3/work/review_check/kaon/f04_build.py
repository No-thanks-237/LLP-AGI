
# -*- coding: utf-8 -*-
"""F04: neutral-kaon flight in S->KK (downstream LHCb). Read-only audit script.

(a) = visible rate exactly as production code computes it (daughters born at the S vertex,
      on the S flight line, at each z-grid point of the acceptance table).
(b) = same final state but each neutral kaon's flight sampled (exponential, ctau_KS/ctau_KL)
      and its daughters transported from the true decay point; decays beyond the tracking
      volume fall out through the real kernel (BornAfterFront / OutsideSciFi).
"""
import sys, os, math, json, time, random
import numpy as np

RUN = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, RUN+'/work/analysis'); sys.path.insert(0, RUN+'/work/response'); sys.path.insert(0, RUN+'/work/widths')
import bc5_lib as L, bc5_production as P, run_scan as RS
import downstream_ctypes as dk

M_K0  = L.M_K0L                 # 0.497611
CTAU_KS = 0.026844              # m, PDG K_S c*tau
CTAU_KL = L.sw.HBAR_C / P.GAM_KL  # m, from the production K_L width (15.337 m)
ZG = np.array(L.ZGRID)

_track = dk.transport_track
_birth = dk.passes_birth_kinematics
MINE, MINPT = dk.BIRTH_CUTS['minE'], dk.BIRTH_CUTS['minPT']

def _passes(daughters, births):
    """daughters: [(px,py,pz,E,chg,q)], births: [(x,y,z)] same length. -> 1.0/0.0"""
    good = []
    for (px,py,pz,e,ch,q),(x,y,z) in zip(daughters,births):
        if not _birth(px,py,pz,e,MINE,MINPT):
            continue
        if _track(x,y,z,px,py,pz,e,ch)[0]:
            good.append(q)
    if len(good) < 2:
        return 0.0
    return 1.0 if (any(q>0 for q in good) and any(q<0 for q in good)) else 0.0

# ---------------------------------------------------------------------------------
# S four-momentum populations drawn through the real production code
# ---------------------------------------------------------------------------------
def sample_S_B(ms, n, rng):
    smp = RS.load_sample(RS.SAMPLE_B)
    sel = np.isin(smp['pdg'], [521,-521,511,-511,531,-531])
    idxs = np.where(sel)[0]
    out = []
    tries = 0
    while len(out) < n and tries < 40*n:
        tries += 1
        i = int(idxs[rng.randrange(len(idxs))])
        pdg = int(smp['pdg'][i])
        bm = {521:'B+', -521:'B+', 511:'B0', -511:'B0', 531:'Bs', -531:'Bs'}[pdg]
        mes = {'kind':'B','bmeson':bm,'x':smp['x'][i],'y':smp['y'][i],'z':smp['z'][i],
               'e':smp['e'][i],'px':smp['px'][i],'py':smp['py'][i],'pz':smp['pz'][i]}
        got = RS.sample_meson_decay_S(ms, mes, rng)
        for s in got or []:
            out.append(('B',)+tuple(s))
    return out

def sample_S_Bquartic(ms, n, rng):
    smp = RS.load_sample(RS.SAMPLE_B)
    sel = np.isin(smp['pdg'], [521,-521,511,-511,531,-531])
    idxs = np.where(sel)[0]
    out = []; tries = 0
    while len(out) < n and tries < 40*n:
        tries += 1
        i = int(idxs[rng.randrange(len(idxs))])
        pdg = int(smp['pdg'][i])
        bm = {521:'B+', -521:'B+', 511:'B0', -511:'B0', 531:'Bs', -531:'Bs'}[pdg]
        mes = {'kind':'B','bmeson':bm,'x':smp['x'][i],'y':smp['y'][i],'z':smp['z'][i],
               'e':smp['e'][i],'px':smp['px'][i],'py':smp['py'][i],'pz':smp['pz'][i]}
        qm = [q for q in P.quartic_channels(ms) if q[0]==bm]
        if not qm: continue
        wq = sum(q[3] for q in qm)
        u = rng.random(); c=0.0; selc=qm[-1]
        for q in qm:
            c += q[3]/wq
            if u<=c: selc=q; break
        kind = '2b' if selc[1]=='2b' else '3b'
        got = RS.sample_quartic_decay_S(ms, mes, kind, selc[2], rng)
        for s in got or []:
            out.append(('Bq',)+tuple(s))
    return out

def sample_S_H(ms, n, rng):
    smp = RS.load_sample(RS.SAMPLE_H)
    ev, first = np.unique(smp['ev'], return_index=True)
    idxs = np.sort(first)
    out = []
    while len(out) < n:
        i = int(idxs[rng.randrange(len(idxs))])
        d = L.two_body(P.M_H, ms, ms, rng)
        if d is None: continue
        pair = L.boost(d, smp['e'][i], (smp['px'][i], smp['py'][i], smp['pz'][i]))
        for s4 in pair:
            out.append(('H', s4[1], s4[2], s4[3], s4[4],
                        smp['x'][i], smp['y'][i], smp['z'][i]))
    return out

# ---------------------------------------------------------------------------------
# one S -> per-z acceptance arrays for the KK channel (current and flight-corrected)
# and for the BR-averaged all-channel response (current, as production code does it)
# ---------------------------------------------------------------------------------
def _split_kk(ms, rng):
    """Return (kind, daughters_rest, kaon_rest) replaying s_decay_products' KK branch.
    daughters_rest: [(m,e,px,py,pz,chg)], kaon_rest: None or [kaon4, kaon4] (mass,e,px,py,pz)"""
    if rng.random() < 0.5:
        d = L.two_body(ms, L.M_K, L.M_K, rng)
        if d is None: return None
        return ('charged', [(L.M_K,d[0][1],d[0][2],d[0][3],d[0][4],+1),
                            (L.M_K,d[1][1],d[1][2],d[1][3],d[1][4],-1)], None)
    d = L.two_body(ms, M_K0, M_K0, rng)
    if d is None: return None
    is_short = rng.random() < 0.5
    out = []; kaons = []
    for i in (0,1):
        kd = L._K0_charged(M_K0, is_short, rng)
        kaons.append((i, is_short, d[i]))
        for (m,e,px,py,pz,ch) in kd:
            b = L.boost([(m,e,px,py,pz)], d[i][1], d[i][2:5])[0]
            out.append((b[0],b[1],b[2],b[3],b[4],ch))
    return ('neutral', out, kaons)

def s_arrays(ms, s, rng, ndecay, allch=True):
    """For one S (tag,e,px,py,pz,x_mm,y_mm,z_mm) compute
       Acur_KK[k], Anew_KK[k], Aall[k] on ZG. Returns (arrays, nkk, nall)."""
    tag,e,px,py,pz,x,y,z = s
    pmod = math.sqrt(px*px+py*py+pz*pz)
    eta = 0.5*math.log((pmod+pz)/max(pmod-pz,1e-300))
    if pz<=0 or pmod<=abs(pz) or not (L.ETA_S_MIN < eta < L.ETA_S_MAX):
        return None
    x0,y0,z0 = x*1e-3, y*1e-3, z*1e-3
    dirx,diry = px/pz, py/pz
    nz = len(ZG)
    Acc_cur = np.zeros(nz); Acc_new = np.zeros(nz); Aall = np.zeros(nz)
    # ---- KK channel (forced) ----
    nkk = 0
    for _ in range(ndecay):
        kk = _split_kk(ms, rng)
        if kk is None:
            continue
        nkk += 1
        kind, drest, kaons = kk
        lab5 = L.boost([(a[0],a[1],a[2],a[3],a[4]) for a in drest], e, (px,py,pz))
        daughters = [(b[1],b[2],b[3],b[4],c[5],c[5]) for b,c in zip(lab5,drest)]
        # transport from S vertex (current)
        for k,zd in enumerate(ZG):
            bx = x0+(zd-z0)*dirx; by = y0+(zd-z0)*diry
            Acc_cur[k] += _passes(daughters, [(bx,by,zd)]*len(daughters))
        # true kaon flight (new): one displacement per neutral kaon
        offs = [np.zeros(3) for _ in daughters]
        if kind=='neutral':
            # map daughter index -> kaon index (order preserved by _split_kk)
            kaon_disp = []
            for (i, is_short, k4) in kaons:
                klab = L.boost([k4], e, (px,py,pz))[0]
                kp = math.sqrt(klab[2]**2+klab[3]**2+klab[4]**2)
                ctau = CTAU_KS if is_short else CTAU_KL
                lam = (kp/k4[0])*ctau*(-math.log(max(rng.random(),1e-12)))
                kaon_disp.append((lam*klab[2]/kp, lam*klab[3]/kp, lam*klab[4]/kp))
            # rebuild the daughter->kaon map by replaying _K0_charged counts is unreliable;
            # instead recompute per kaon with a dedicated second pass
        else:
            kaon_disp = None
        # For the neutral case we regenerate daughters per kaon with explicit mapping:
        if kind=='neutral':
            is_short = kaons[0][1]
            disc = []
            for i in (0,1):
                k4 = kaons[i][2]
                klab = L.boost([k4], e, (px,py,pz))[0]
                kp = math.sqrt(klab[2]**2+klab[3]**2+klab[4]**2)
                lam = (kp/k4[0])*CTAU_KS if is_short else (kp/k4[0])*CTAU_KL
                u = -math.log(max(rng.random(),1e-12))
                disp = (lam*u*klab[2]/kp, lam*u*klab[3]/kp, lam*u*klab[4]/kp)
                for (m,ee,ax,ay,az,ch) in drest_by_kaon[i]:
                    lab5b = L.boost([(m,ee,ax,ay,az)], e, (px,py,pz))[0]
                    disc.append((lab5b[1],lab5b[2],lab5b[3],lab5b[4],ch,ch,disp))
            for k,zd in enumerate(ZG):
                bx = x0+(zd-z0)*dirx; by = y0+(zd-z0)*diry
                births = [(bx+dd[6][0], by+dd[6][1], zd+dd[6][2]) for dd in disc]
                Acc_new[k] += _passes([(d[0],d[1],d[2],d[3],d[4],d[5]) for d in disc], births)
        else:
            # charged KK: identical to current
            Acc_new += Acc_cur
    # ---- all channels, current treatment (replicates acceptance_profile) ----
    nall = 0
    if allch:
        br,_ = L.channel_branching(ms)
        names = list(br.keys()); cum = np.cumsum([br[n] for n in names])
        for _ in range(ndecay):
            u = rng.random(); idx = int(np.searchsorted(cum,u)); idx = min(idx,len(names)-1)
            rest = L.s_decay_products(ms, names[idx], rng)
            if not rest: continue
            nall += 1
            chg = [a[5] for a in rest]
            lab5 = L.boost([(a[0],a[1],a[2],a[3],a[4]) for a in rest], e, (px,py,pz))
            dd = [(b[1],b[2],b[3],b[4],c,c) for b,c in zip(lab5,chg)]
            for k,zd in enumerate(ZG):
                bx = x0+(zd-z0)*dirx; by = y0+(zd-z0)*diry
                Aall[k] += _passes(dd, [(bx,by,zd)]*len(dd))
    if nkk>0: Acc_cur/=nkk; Acc_new/=nkk
    if nall>0: Aall/=nall
    return Acc_cur, Acc_new, Aall, nkk, nall

def zweight(ms, theta2, e, px, py, pz, z0):
    ctau = L.sw.HBAR_C/(L.sw.gamma_total(ms)*theta2)
    alpha = ms/(ctau*max(pz,1e-9))
    w = np.zeros(len(ZG)); m = ZG>z0
    w[m] = alpha*np.exp(-(ZG[m]-z0)*alpha)
    return w
