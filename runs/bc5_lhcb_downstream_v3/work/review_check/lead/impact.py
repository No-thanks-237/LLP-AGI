
# -*- coding: utf-8 -*-
"""Impact study: run the production scan at a subset of masses with a physics variant
   applied in-memory, and compare the exclusion crossing at 2.3 / L.  Production code is
   NOT modified; the variant is a monkeypatch that lives only in this process."""
import argparse, json, math, os, random, sys, time
import numpy as np

RUN = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3"
sys.path.insert(0, os.path.join(RUN, "work", "analysis"))
sys.path.insert(0, os.path.join(RUN, "work", "widths"))
import run_scan as R
import bc5_lib as L
import bc5_production as P
import scalar_widths as sw

OUTDIR = os.path.dirname(os.path.abspath(__file__))
LEVEL_S = 2.3
LUMI = 25.0


def _rotmat(rng):
    q = np.array([rng.gauss(0, 1) for _ in range(4)]); q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array([[1-2*(y*y+z*z), 2*(x*y-z*w),   2*(x*z+y*w)],
                     [2*(x*y+z*w),   1-2*(x*x+z*z), 2*(y*z-x*w)],
                     [2*(x*z-y*w),   2*(y*z+x*w),   1-2*(x*x+y*y)]])


def apply_variant(name):
    if name == 'base':
        return
    if name == 'rot':
        orig = L.s_decay_products
        def f(ms, channel, rng):
            d = orig(ms, channel, rng)
            if d and channel in L._LIB_MAP:
                M = _rotmat(rng)
                d = [(a[0], a[1], *(M @ np.array(a[2:5])), a[5]) for a in d]
            return d
        L.s_decay_products = f
        return
    if name == 'kLdead':
        orig = L._K0_charged
        def f2(m_k, is_short, rng):
            return orig(m_k, is_short, rng) if is_short else []
        L._K0_charged = f2
        return
    if name == 'multi4pi':
        orig = L.s_decay_products
        def f3(ms, channel, rng):
            if channel == 'multi':
                if ms < 4*L.M_PI:
                    return []
                return L.nbody(ms, [L.M_PI]*4, rng) and \
                       [(d[0], d[1], d[2], d[3], d[4], c)
                        for d, c in zip(L.nbody(ms, [L.M_PI]*4, rng), (+1, +1, -1, -1))]
            return orig(ms, channel, rng)
        L.s_decay_products = f3
        return
    if name == 'k0half':
        # F09-2: a PDG 311/-311 is a flavour state = 50% K_L + 50% K_S; only the K_L component
        # can give pi0 S (K_S is CP suppressed).  Halve the neutral-kaon mixing weight.
        origk = P.br_K_to_S
        def fk(ms, which='K+'):
            v = origk(ms, which)
            return 0.5*v if which == 'K0L' else v
        P.br_K_to_S = fk
        return
    if name == 'gamma':
        if 'gamma' not in L.CHANNELS:
            L.CHANNELS.append('gamma')
        return
    if name == 'bxi':
        _install_b_model(name)
        return
    if name in ('bincl', 'bsext'):
        _install_b_model(name)
        return
    if name == 'libmass':
        orig = L.HadronLib.sample
        def f4(self, ms, channel, rng):
            got = orig(self, ms, channel, rng)
            if not got:
                return got
            c = self.chan.get(channel)
            i = int(np.argmin(np.abs(c['mass_nodes'] - ms)))
            node = float(c['mass_nodes'][i])
            lam = ms / node
            if abs(lam - 1.0) < 1e-6:
                return got
            return [(a[0]*lam, a[1]*lam, a[2]*lam, a[3]*lam, a[4]*lam, a[5]) for a in got]
        L.HadronLib.sample = f4
        return
    raise SystemExit('unknown variant ' + name)


GF, MT, VTS, VTB = 1.1663787e-5, 173.0, 0.0404, 0.9991
XI_W = 3*math.sqrt(2)*GF*MT*MT*abs(VTS*VTB)/(16*math.pi**2)


M_PHI = 1.019461


def _install_b_model(name):
    """Winkler-only pieces.  'bxi': only the b->s loop factor becomes Winkler's CKM expression.
       'bincl': B+/B0 take eq:inclusive below 4.7 GeV and eq:BKphi above; B_s unchanged.
       'bsext': B_s keeps going above m_Bs-m_phi as B_s -> phi S (eq:BKphi with A_phi ~ A_K*);
                B+/B0 unchanged."""
    orig = P.br_B_to_S
    g_w = (P.M_B_QUARK/sw.V_HIGGS)*XI_W

    def inclusive(ms, mX, GX):
        return g_w*g_w*(mX*mX-ms*ms)**2/(32*math.pi*mX**3*GX)

    def two_body_br(ms, mX, GX, mY, A=None):
        if mX < mY + ms:
            return 0.0
        if A is None:
            fK = 0.33/(1.0-ms*ms/37.5)
            ME2 = 0.25*(mX*mX-mY*mY)**2/(P.M_B_QUARK-P.M_S_QUARK)**2*fK*fK
        else:
            lam = (mX*mX-(mY-ms)**2)*(mX*mX-(mY+ms)**2)/mX**4
            ME2 = 0.25*mX**4*lam/(P.M_B_QUARK+P.M_S_QUARK)**2*A*A
        lam_s = math.sqrt(max((mX*mX-(mY-ms)**2)*(mX*mX-(mY+ms)**2), 0.0))/(mX*mX)
        return g_w*g_w*ME2*lam_s/(16*math.pi*mX*GX)

    def A_Kstar(ms):
        return 1.36/(1.0-ms*ms/27.9) - 0.99/(1.0-ms*ms/36.8)

    if name == 'bxi':
        def br(ms, meson='B+'):
            out = orig(ms, meson)
            f = (XI_W/P.XI_SB)**2
            return {k: v*f for k, v in out.items()}
        P.br_B_to_S = br
        return

    if name == 'bincl':
        def br(ms, meson='B+'):
            if meson == 'Bs':
                return orig(ms, 'Bs')
            mX, GX = (P.M_BP, P.GAM_BP) if meson == 'B+' else (P.M_B0, P.GAM_B0)
            if ms < 4.70:
                return {'B_incl': inclusive(ms, mX, GX)}
            return {'B_2body': two_body_br(ms, mX, GX, P.M_KP)}
        P.br_B_to_S = br
        P.B_CHANNEL_MASS['B_incl'] = 0.0
        P.B_CHANNEL_MASS['B_2body'] = P.M_KP
        return

    if name == 'bsext':
        def br(ms, meson='B+'):
            if meson != 'Bs':
                return orig(ms, meson)
            m_sw = P.M_BS - M_PHI
            if ms < m_sw:
                return orig(ms, 'Bs')
            return {'B_2body': two_body_br(ms, P.M_BS, P.GAM_BS, M_PHI, A=A_Kstar(ms))}
        P.br_B_to_S = br
        P.B_CHANNEL_MASS['B_2body'] = M_PHI
        return
    raise SystemExit('unknown B model ' + name)


def crossings(th, s, level):
    out = []
    f = np.log10(np.maximum(s, 1e-300)) - math.log10(level)
    for i in range(len(th)-1):
        if f[i] == 0:
            out.append(th[i])
        elif f[i]*f[i+1] < 0:
            t = f[i]/(f[i]-f[i+1])
            out.append(10**(np.log10(th[i]) + t*(np.log10(th[i+1])-np.log10(th[i]))))
    return out


def run_mass(ms, nsamp, ndecay, seed):
    rng = random.Random(seed)
    tab = R.run_mass(ms, nsamp, rng, ndecay=ndecay, verbose=False)
    th = np.exp(np.linspace(math.log(R.THETA2_MIN), math.log(R.THETA2_MAX), 240))
    s, eps = [], []
    for t2 in th:
        d = R.signal_per_lumi(ms, tab, float(t2), 1.0, R.SIGMA_BB_PB, R.SIGMA_INEL_PB,
                              R.N_GEN_BBAR, R.N_GEN_MINBIAS, R.N_GEN_HIGGS)
        s.append(d['B'] + d['K'] + d['H'])
        eps.append((d['eps_B'], d['eps_K'], d['eps_H']))
    s = np.array(s)
    cr = crossings(th, s, LEVEL_S/LUMI)
    lo = min(cr) if cr else float('nan')
    hi = max(cr) if len(cr) >= 2 else float('nan')
    row = {'m': ms, 'lo': lo, 'hi': hi, 'n_cross': len(cr)}
    if np.isfinite(lo):
        i = int(np.argmin(np.abs(np.log10(th)-np.log10(lo))))
        row['eps_at_lo'] = eps[i]
    return row


def init_yield_corr():
    """Reproduce run_scan.main()'s runtime B-yield normalisation (it is 0 until main() runs)."""
    smpB = R.load_sample(R.SAMPLE_B)
    if smpB is not None and R.N_GEN_BBAR > 0:
        sel = np.isin(smpB['pdg'], [521, -521, 511, -511, 531, -531])
        per_event = float(sel.sum())/float(R.N_GEN_BBAR)
        R.B_YIELD_CORR = R.REF_B_MESON_YIELD/per_event if per_event > 0 else 1.0
    else:
        R.B_YIELD_CORR = 1.0
    print('  B_YIELD_CORR = %.4f' % R.B_YIELD_CORR, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--variant', required=True)
    ap.add_argument('--masses', required=True)
    ap.add_argument('--nsamp', type=int, default=1200)
    ap.add_argument('--ndecay', type=int, default=8)
    ap.add_argument('--seed', type=int, default=20261016)
    ap.add_argument('--tag', default=None)
    a = ap.parse_args()
    libdir = os.path.join(RUN, 'work', 'magnus', 'hadron_lib')
    assert L.load_hadron_lib(libdir), 'no hadron library'
    init_yield_corr()
    apply_variant(a.variant)
    masses = [float(x) for x in a.masses.split(',')]
    res = []
    t0 = time.time()
    for ms in masses:
        r = run_mass(round(ms, 6), a.nsamp, a.ndecay, a.seed)
        res.append(r)
        print('%-10s m=%-9.5f lo=%s hi=%s (%.0fs)' % (
            a.variant, ms, ('%.3e' % r['lo']) if np.isfinite(r['lo']) else 'None',
            ('%.3e' % r['hi']) if np.isfinite(r['hi']) else 'None', time.time()-t0), flush=True)
    tag = a.tag or a.variant
    with open(os.path.join(OUTDIR, 'impact_%s.json' % tag), 'w') as f:
        json.dump({'variant': a.variant, 'nsamp': a.nsamp, 'ndecay': a.ndecay,
                   'seed': a.seed, 'level_s_over_lumi': LEVEL_S/LUMI, 'rows': res}, f, indent=1)
    print('wrote impact_%s.json' % tag)


if __name__ == '__main__':
    main()