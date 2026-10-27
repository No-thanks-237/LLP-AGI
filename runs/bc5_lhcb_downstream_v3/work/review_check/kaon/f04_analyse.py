
# -*- coding: utf-8 -*-
"""F04 analysis: KK-channel visible-rate ratio (b)/(a), KK share, effect on the total signal.
Sub-channel decomposition of S->K0K0bar: 50% K+K-, 25% K_S K_S, 25% K_L K_L."""
import sys, os, math, json, glob
import numpy as np
RUN = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, RUN+'/work/analysis'); sys.path.insert(0, RUN+'/work/response'); sys.path.insert(0, RUN+'/work/widths')
import bc5_lib as L
HERE = os.path.dirname(os.path.abspath(__file__))
ZG = np.array(L.ZGRID)
FRONTZ = 7.7

scan = json.load(open(RUN+'/outputs/data/signal_scan.json'))
THE = np.array(scan['theta2']); SC = scan['signal_per_fb_per_channel']
_SK = list(SC); _SM = np.array([float(k) for k in _SK]); _o = np.argsort(_SM); _SM = _SM[_o]; _SK = [_SK[i] for i in _o]


def scan_chan(ms, t2, c):
    vals = [float(np.exp(np.interp(math.log(t2), np.log(THE), np.log(np.maximum(
                np.array([r.get(c, 0.0) for r in SC[k]]), 1e-300))))) for k in _SK]
    return float(np.exp(np.interp(math.log(ms), np.log(_SM), np.log(np.maximum(np.array(vals), 1e-300)))))


def zweight(ms, theta2, pz, z0):
    ctau = L.sw.HBAR_C/(L.sw.gamma_total(ms)*theta2)
    alpha = ms/(ctau*max(pz, 1e-9))
    w = np.zeros(len(ZG)); m = ZG > z0
    w[m] = alpha*np.exp(-(ZG[m]-z0)*alpha)
    return w


def boot(fn, n, nb=400, seed=3):
    rng = np.random.default_rng(seed)
    v = np.array([fn(rng.integers(0, n, n)) for _ in range(nb)])
    v = v[np.isfinite(v)]
    return float(v.std(ddof=1)) if len(v) > 1 else float('nan')


def analyse(path):
    d = np.load(path)
    ms = float(d['mass']); earr = d['edges']; sB = d['sB']; sH = d['sH']
    br, _ = L.channel_branching(ms); brkk = br['KK']
    out = {'mass': ms, 'br_KK': brkk,
           'edges': {'long': float(earr[0]), 'short': float(earr[1])},
           'sB': {'long': float(sB[0]), 'short': float(sB[1])},
           'sH': {'long': float(sH[0]), 'short': float(sH[1])}}
    for iedge, edge in enumerate(('long', 'short')):
        t2 = float(earr[iedge]); E = {}
        for name in ('B', 'Bq', 'H'):
            Ac = d[name+'Ac']; An = d[name+'An']; Aa = d[name+'Aa']; meta = d[name+'meta']
            if len(Ac) == 0: continue
            pz = meta[:, 1]; z0 = meta[:, 2]
            W = np.array([zweight(ms, t2, pz[i], z0[i]) for i in range(len(Ac))])
            pc = np.trapezoid(W*Ac, ZG, axis=1); pn = np.trapezoid(W*An, ZG, axis=1); pa = np.trapezoid(W*Aa, ZG, axis=1)
            def st(idx):
                c = pc[idx].mean(); nn = pn[idx].mean(); a = pa[idx].mean()
                return (nn/c if c > 0 else np.nan, brkk*c/a if a > 0 else np.nan, (nn-c)/a if a > 0 else np.nan)
            m = dict(n=int(len(pc)), pc=float(pc.mean()), pn=float(pn.mean()), pa=float(pa.mean()))
            m['r'], m['f'], m['dss'] = [float(x) for x in st(np.arange(len(pc)))]
            m['r_se'] = boot(lambda ii: st(ii)[0], len(pc))
            m['dss_se'] = boot(lambda ii: st(ii)[2], len(pc))
            m['f_se'] = boot(lambda ii: st(ii)[1], len(pc))
            order = np.argsort(-pc); tot = pc.sum()
            m['neff'] = float(tot*tot/np.sum(pc*pc)) if tot > 0 else float('nan')
            m['top10'] = float(pc[order[:10]].sum()/tot) if tot > 0 else float('nan')
            # sub-channel decomposition (same z weight): 50% C, 25% KS, 25% KL
            sub = {}
            for tag, wt in (('C', 0.5), ('KS', 0.25), ('KL', 0.25)):
                acc_ = d[name+'Ac'+tag]; new_ = d[name+'An'+tag]
                if acc_.shape[0] == 0: continue
                sub[tag] = dict(p_cur=float(np.trapezoid(W*acc_, ZG, axis=1).mean()),
                                p_new=float(np.trapezoid(W*new_, ZG, axis=1).mean()), w=wt)
            m['sub'] = sub
            E[name] = m
        wB = float(sB[iedge])*25.0; wH = float(sH[iedge])*25.0
        E['weights'] = dict(B=wB, H=wH)
        if 'B' in E and 'H' in E:
            E['total_dss'] = (wB*E['B']['dss'] + wH*E['H']['dss'])/(wB+wH)
            E['wBfrac'] = wB/(wB+wH)
        out[edge] = E
    # flight diagnostics (long edge only; theta2-independent kinematics)
    for srcn in ('B', 'H'):
        fl = d[srcn+'flight']
        if len(fl) == 0: continue
        isS = fl[:, 0] > 0.5
        for tag, sel in (('KS', isS), ('KL', ~isS)):
            a2 = fl[sel]
            if len(a2) == 0: continue
            lam = a2[:, 1]; bg = a2[:, 2]; cz = a2[:, 3]
            out['flight%s_%s' % (srcn, tag)] = dict(n=int(len(a2)), mean=float(lam.mean()), median=float(np.median(lam)),
                                      bg=float(bg.mean()),
                                      frac_beyond_front_from_2m=float((2.0+lam*cz > FRONTZ).mean()))
    fl = d['Bflight']
    if len(fl):
        isS = fl[:, 0] > 0.5
        for tag, sel in (('KS', isS), ('KL', ~isS)):
            a = fl[sel]
            if len(a) == 0: continue
            lam = a[:, 1]; bg = a[:, 2]; cz = a[:, 3]
            out['flight_'+tag] = dict(n=int(len(a)), mean=float(lam.mean()), median=float(np.median(lam)),
                                      bg=float(bg.mean()),
                                      frac_beyond_front_from_1m=float((1.0+lam*cz > FRONTZ).mean()),
                                      frac_beyond_front_from_2m=float((2.0+lam*cz > FRONTZ).mean()))
    return out


res = {}
for path in sorted(glob.glob(os.path.join(HERE, 'f04_m*.npz'))):
    a = analyse(path); res['%.4f' % a['mass']] = a
    print('\n===== m=%.4f  BR_KK=%.3f' % (a['mass'], a['br_KK']))
    for tag in ('B_KS','B_KL','H_KS','H_KL'):
        if 'flight'+tag in a:
            f = a['flight'+tag]
            print('  %s flight: n=%d mean=%.2f m median=%.2f m bg=%.1f  P(z_dec>frontZ|z0=2m)=%.3f' %
                  (tag, f['n'], f['mean'], f['median'], f['bg'], f['frac_beyond_front_from_2m']))
    for edge in ('long', 'short'):
        e = a[edge]
        print('  -- edge %-5s theta2=%.3e  sB25=%.2f sH25=%.2f (B %.2f)' %
              (edge, a['edges'][edge], a['sB'][edge]*25, a['sH'][edge]*25, e.get('wBfrac', np.nan)))
        for name in ('B', 'Bq', 'H'):
            if name not in e: continue
            m = e[name]
            print('     %-3s n=%4d pcur=%.3e r=%.4f+-%.4f f=%.3f+-%.3f dss=%+.4f+-%.4f neff=%.0f' %
                  (name, m['n'], m['pc'], m['r'], m['r_se'], m['f'], m['f_se'], m['dss'], m['dss_se'], m['neff']))
            if name in ('B', 'H') and m.get('sub'):
                s = m['sub']
                print('        sub: pC=%.4f KS_prompt=%.4f KS_true=%.4f  KL_prompt=%.4f KL_true=%.4f' % (
                    s['C']['p_cur'], s['KS']['p_cur'], s['KS']['p_new'], s['KL']['p_cur'], s['KL']['p_new']))
        if 'total_dss' in e:
            print('     TOTAL dss=%+.4f  (B %+.4f @w=%.2f, H %+.4f)' % (
                e['total_dss'], e['B']['dss'], e['wBfrac'], e['H']['dss']))
        # local log-log slope of s vs theta2 (for the boundary shift estimate)
        _key = min(SC.keys(), key=lambda k: abs(float(k)-a['mass']))
        s_all = np.array([r.get('B', 0)+r.get('K', 0)+r.get('H', 0) for r in SC[_key]])
        i = int(np.argmin(np.abs(np.log(THE/a['edges'][edge]))))
        lo = max(0, i-2); hi = min(len(THE)-1, i+2)
        slope = (math.log(max(s_all[hi], 1e-300))-math.log(max(s_all[lo], 1e-300)))/(math.log(THE[hi])-math.log(THE[lo]))
        print('     local dln s/dln theta2 = %.2f  -> boundary shift dtheta2/theta2 = %+.3f' % (
            slope, -e['total_dss']/slope if slope != 0 else np.nan))
json.dump(res, open(os.path.join(HERE, 'f04_results.json'), 'w'), indent=1, default=float)
print('\nwrote f04_results.json')