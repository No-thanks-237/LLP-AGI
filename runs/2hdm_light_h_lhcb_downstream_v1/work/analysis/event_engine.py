# -*- coding: utf-8 -*-
"""Event engine: probability algebra, normalisation and s(mH, param) for the LHCb-downstream
2HDM light-H scan, plus the BC5 v3.6 regression path.

DESIGN.md section 4, exactly:

    p_i(tanb)  = sum_ch BR_ch(mH,tanb) * p_window(mH, pz_i, z0_i, A_ch(i), ctau_H)   [1]
    o_hist     = w_hist * (1 - prod_{i in hist} (1 - p_i))                          [2]
    o_parent   = sum_hist o_hist            (histories of one parent are exclusive)  [3]
    o_event    = 1 - prod_{parents in event} (1 - o_parent)                         [4]
    eps_c      = yield_corr_c * 0.4 / N_gen_c * sum_events o_event                  [5]
    s          = 1e3 * L * sum_c sigma_c * eps_c          (L in fb^-1, sigma in pb) [6]

The 0.4 reconstruction efficiency is applied exactly ONCE per collision, inside eps.
The production branching ratio w_hist is applied exactly once per history, AFTER the OR
over that history's scalars ([2]) -- never per scalar.
"""
from __future__ import annotations
import json, math, os, sys, time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(RUN, 'work', 'response'))

import twohdm_lib as L
import models as M
import response_table as RT

PB_PER_FB = 1e3                 # 1 fb^-1 = 1e3 pb^-1
S90 = -math.log(0.1)            # 2.302585093: zero background, zero observed events, 90 % CL
LUMIS = (25.0, 300.0)


# ======================================================================================
# probability algebra
# ======================================================================================
def per_scalar_p(table, br, ctau_m, mH):
    """p_i = sum_ch BR_ch * p_window(A_ch) for every H candidate (DESIGN eq. [1])."""
    A = table['A']
    pz = table['pz']
    p = np.zeros(len(pz))
    if 'hist_lam' in table and table.get('model') != 'bc5':
        # v3 path: the parent decay position is integrated out analytically
        # (Rao-Blackwellisation of the exponential parent flight), which removes the
        # dominant variance source of the short-lived crossings.
        lam = np.asarray(table['hist_lam'], dtype=float)[table['scal_hist']]
        for ic, ch in enumerate(table['channels']):
            b = float(br.get(ch, 0.0))
            if b <= 0.0:
                continue
            p += b * L.p_window_vertex(mH, pz, A[:, ic, :], ctau_m, lam)
        return p
    z0 = table['z0']
    for ic, ch in enumerate(table['channels']):
        b = float(br.get(ch, 0.0))
        if b <= 0.0:
            continue
        p += b * L.p_window(mH, pz, z0, A[:, ic, :], ctau_m)
    return p


def _group_starts(idx):
    """Start offsets of runs of equal consecutive values (idx must be sorted)."""
    if len(idx) == 0:
        return np.zeros(0, dtype=np.int64)
    return np.flatnonzero(np.r_[True, idx[1:] != idx[:-1]])


def event_responses(table, p, hist_weights):
    """Return (event_ids, o_event) implementing eqs. [2]-[4]."""
    hist = table['scal_hist']
    n_hist = len(table['hist_event'])
    w = np.empty(n_hist)
    for i, k in enumerate(table['hist_wkey']):
        w[i] = hist_weights[k]
    # Horvitz-Thompson weight of the pz-stratified parent subsample (1.0 for uniform
    # sampling, so the BC5 regression path is unchanged).
    if 'hist_scale' in table:
        w *= np.asarray(table['hist_scale'], dtype=float)
    if len(hist) == 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0)
    # [2] OR over the scalars of one history, then apply the history weight ONCE
    st_h = _group_starts(hist)
    log1mp = np.log1p(-np.clip(p, 0.0, 1.0 - 1e-15))
    o_in = -np.expm1(np.add.reduceat(log1mp, st_h))
    o_hist = w * o_in
    # [3] histories of one parent are mutually exclusive -> add.  A parent is identified by
    # the PAIR (collision, parent), because parent indices are only unique inside a group
    # (review BUG C); the histories are sorted by (event, parent, ...) in the table builder.
    par = table['hist_parent']
    ev = table['hist_event']
    if len(par):
        st_p = np.flatnonzero(np.r_[True, (par[1:] != par[:-1]) | (ev[1:] != ev[:-1])])
    else:
        st_p = np.zeros(0, dtype=np.int64)
    p_par = np.add.reduceat(o_hist, st_p)
    ev_of_par = ev[st_p]
    # [4] parents of one collision -> OR
    st_e = _group_starts(ev_of_par)
    o_ev = -np.expm1(np.add.reduceat(np.log1p(-np.clip(p_par, 0.0, 1.0 - 1e-15)), st_e))
    return ev_of_par[st_e], o_ev


def group_efficiency(table, ev_ids, o_ev, n_gen, yield_corr, counts=None):
    """eps = yield_corr * (n_tot/n_smp) * 0.4/N_gen * sum_events o ; se from the event variance.

    The (n_tot/n_smp) factor scales the SAMPLED parents back to the full parent sample:
    a subsample of parents must not rescale the signal (BC5 channel_efficiency).

    Collisions whose response is exactly zero stay in the sample (they are listed in
    table['samp_events']), otherwise the MC variance would be underestimated.
    """
    ev_samp = table['samp_events']
    out = {}
    for g, name in enumerate(M.GROUP_NAMES):
        ev_g = ev_samp.get(name)
        if ev_g is None or len(ev_g) == 0 or n_gen.get(name, 0) <= 0:
            out[name] = (0.0, 0.0, 0)
            continue
        o = np.zeros(len(ev_g))
        if len(ev_ids):
            pos = np.searchsorted(ev_g, ev_ids)
            ok = (pos < len(ev_g))
            ok[ok] = ev_g[pos[ok]] == ev_ids[ok]
            o[pos[ok]] = o_ev[ok]
        n_gen_g = float(n_gen[name])
        n_tot, n_smp = (counts or {}).get(name, [1, 1])
        scale = (float(n_tot) / float(n_smp)) if n_smp else 1.0
        norm = yield_corr[name] * L.RECO_EFF * scale / n_gen_g
        eps = norm * float(o.sum())
        var = float(o.var()) if len(o) > 1 else 0.0
        se = norm * math.sqrt(max(var, 0.0) * len(o))
        out[name] = (eps, se, len(o))
    return out


# ======================================================================================
# one parameter point
# ======================================================================================
def evaluate(table, model, mH, param, lumi_fb=25.0, sigma=None, extra=None):
    """s and its MC error at one (mH, param) point.

    sigma (optional) lets the caller pass model.sigma_pb(mH,param) once for a whole scan
    (it is the expensive production call); otherwise it is evaluated here.
    """
    mH = float(mH)
    br, gamma_tot = model.decay_br(mH, param)
    ctau = model.ctau_m(mH, param)
    if sigma is None:
        sigma = model.sigma_pb(mH, param)
    w = model.history_weights(mH, param, table)
    p = per_scalar_p(table, br, ctau, mH)
    ev_ids, o_ev = event_responses(table, p, w)
    eff = group_efficiency(table, ev_ids, o_ev, table['n_gen'], table['yield_corr'],
                           counts=table.get('counts'))
    out = {'mH': mH, 'param': float(param), 'lumi_fb': float(lumi_fb),
           'Gamma_H_GeV': float(gamma_tot), 'ctau_m': float(ctau),
           'eps': {}, 'eps_se': {}, 's': {}, 's_se': {}, 'n_coll': {},
           'sigma_group': {k: float(sigma.get(k, 0.0)) for k in M.GROUP_NAMES},
           's_total': 0.0, 's_total_se': 0.0}
    if 'sigma_h_total' in sigma:
        out['sigma_h_pb'] = float(sigma['sigma_h_total'])
        out['BR_hHH'] = float(sigma.get('BR_hHH_detail', 0.0))
    var_tot = 0.0
    for name in M.GROUP_NAMES:
        eps, se_eps, ncol = eff[name]
        sg = float(sigma.get(name, 0.0))
        s_g = PB_PER_FB * lumi_fb * sg * eps
        se_g = PB_PER_FB * lumi_fb * sg * se_eps
        out['eps'][name] = eps
        out['eps_se'][name] = se_eps
        out['s'][name] = s_g
        out['s_se'][name] = se_g
        out['n_coll'][name] = int(ncol)
        out['s_total'] += s_g
        var_tot += se_g * se_g
    out['s_total_se'] = math.sqrt(var_tot)
    out['br_min'] = min(br.values()) if br else 0.0
    for ch in ('gaga', 'ee', 'mumu', 'tau', 'tautau', 'pipi', 'KK', 'multi', '4pi',
               'cc', 'ss', 'gg', 'bb', 'e', 'mu', 'gamma'):
        if ch in br:
            out['br_' + ch] = br[ch]
    out['br'] = br
    if extra:
        out.update(extra)
    return out


# ======================================================================================
# grids and crossings
# ======================================================================================
def evaluate_param_grid(table, model, mH, params, lumi_fb=25.0, sigma=None, verbose=False,
                        extra_meta=True):
    rows = []
    for i, par in enumerate(params):
        extra = model.extra_meta(mH, float(par)) if extra_meta else None
        rows.append(evaluate(table, model, mH, float(par), lumi_fb=lumi_fb, sigma=sigma,
                             extra=extra))
        if verbose and (i + 1) % 20 == 0:
            print('    param %d/%d s=%.3f' % (i + 1, len(params), rows[-1]['s_total']))
    return rows


def find_crossings(params, s_values, level=S90, lifetime_sign=-1):
    """All sign changes of s(level) in log(parameter), linearly interpolated in log(param).

    Edge labelling uses the model's lifetime direction:
        lifetime_sign = +1 if a LARGER parameter means a LONGER lifetime (2HDM: tanb),
        lifetime_sign = -1 if a larger parameter means a SHORTER lifetime (BC5: theta^2).
    The long-lived edge is the crossing where moving toward a longer lifetime takes s
    below the threshold, i.e. where sign(ds/dparam) = -lifetime_sign; the other one is the
    short-lived edge.  (Verified against the BC5 CSV: BC5 rising crossing = long-lived.)
    """
    lp = np.log(np.asarray(params, dtype=float))
    sv = np.asarray(s_values, dtype=float) - level
    out = []
    for i in range(len(lp) - 1):
        if not np.isfinite(sv[i]) or not np.isfinite(sv[i + 1]):
            continue
        if sv[i] == 0.0:
            out.append({'param': float(params[i]),
                        'edge': 'long_lived' if (i + 1 < len(sv) and sv[i + 1] > 0) else 'short_lived'})
            continue
        if sv[i] * sv[i + 1] < 0.0:
            # s(tanb) varies over many decades, so interpolate ln s (not s) against ln param:
            # a linear-in-s interpolation biases a steep crossing by tens of percent.
            s1, s2 = s_values[i], s_values[i + 1]
            if s1 > 0.0 and s2 > 0.0:
                l1, l2 = math.log(s1), math.log(s2)
                f = (l1 - math.log(level)) / (l1 - l2)
            else:
                f = sv[i] / (sv[i] - sv[i + 1])
            f = min(max(f, 0.0), 1.0)
            dsign = 1 if sv[i + 1] > sv[i] else -1
            # long-lived edge: moving toward a longer lifetime lowers s below threshold,
            # i.e. sign(ds/dparam) = -sign(dlifetime/dparam)
            edge = 'long_lived' if dsign * int(lifetime_sign) < 0 else 'short_lived'
            out.append({'param': float(np.exp(lp[i] + f * (lp[i + 1] - lp[i]))),
                        'edge': edge, 'ds_dparam_sign': dsign})
    return out


def refine_crossing(table, model, mH, lumi_fb, p_lo, p_hi, level=S90, n_iter=40,
                    max_step=1.15):
    """Bisect s(param) = level between two parameters (log space, local window).

    The sampled tanb grid is coarse relative to the very steep 2HDM s(tanb); the reported
    boundary should not depend on the grid.  A local bisection costs n_iter evaluations
    (milliseconds each) and gives the crossing to 1e-4 in ln(param).  The bracket is
    expanded by max_step until the two ends straddle the level.
    """
    def s_at(p):
        return evaluate(table, model, mH, float(p), lumi_fb=lumi_fb)['s_total']
    lo, hi = sorted((float(p_lo), float(p_hi)))
    slo, shi = s_at(lo), s_at(hi)
    for _ in range(25):
        if (slo - level) * (shi - level) <= 0.0:
            break
        lo = max(lo / max_step, 1e-30) if slo < level else lo
        hi = hi * max_step if shi < level else hi
        slo, shi = s_at(lo), s_at(hi)
    if (slo - level) * (shi - level) > 0.0:
        return None
    for _ in range(n_iter):
        mid = math.sqrt(lo * hi)
        sm = s_at(mid)
        if (slo - level) * (sm - level) <= 0.0:
            hi, shi = mid, sm
        else:
            lo, slo = mid, sm
        if hi / lo < 1.0002:
            break
    pm = math.sqrt(lo * hi)
    return {'param_refined': pm, 's_at_refined': s_at(pm), 'bracket': [lo, hi]}


def sigma_for_grid(model, mH, params):
    """sigma_pb dict for the first parameter (cross sections are param-dependent for the
    2HDM Higgs channel); evaluate() accepts it per point, so caller may recompute."""
    return model.sigma_pb(mH, float(params[0]))


# ======================================================================================
# CLI: one mass at a time, so the Lead can drive the full 2-D scan mass by mass
# ======================================================================================
def _param_grid(kind, lo, hi, n):
    if kind == 'log':
        return np.exp(np.linspace(math.log(lo), math.log(hi), n))
    return np.linspace(lo, hi, n)


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description='response table + s(mH,param) for one mass')
    ap.add_argument('--model', choices=('2hdm', 'bc5'), default='2hdm')
    ap.add_argument('--mass', type=float, required=True)
    ap.add_argument('--param-lo', type=float, default=None)
    ap.add_argument('--param-hi', type=float, default=None)
    ap.add_argument('--nparam', type=int, default=40)
    ap.add_argument('--param-grid', choices=('log', 'lin'), default='log')
    ap.add_argument('--lumi', default=None,
                    help='fb^-1 value, or "both" (default) for 25 and 300')
    ap.add_argument('--b-mode', default=None, choices=(None, 'foresee', 'paper_mb', 'exclusive'),
                    help='B -> X_s H rate prescription (production module DEFAULT_B_MODE)')
    ap.add_argument('--nsamp', type=str, default=None,
                    help='JSON {"B":1500,"inel":1500,"h":1500} parents per group')
    ap.add_argument('--ndecay', type=str, default=None,
                    help='JSON per-channel decay samples, e.g. {"*":10,"gg":6}')
    ap.add_argument('--seed', type=int, default=RT.DEFAULT_SEED)
    ap.add_argument('--samples', type=str, default=None)
    ap.add_argument('--no-cache', action='store_true')
    ap.add_argument('--force', action='store_true')
    ap.add_argument('--out', type=str, default=None, help='write JSON here')
    ap.add_argument('--quiet', action='store_true')
    args = ap.parse_args(argv)

    verbose = not args.quiet
    if args.model == '2hdm':
        model = M.TwoHDMModel(verbose=verbose, b_mode=args.b_mode)
        lo = 5.0 if args.param_lo is None else args.param_lo
        hi = 1e5 if args.param_hi is None else args.param_hi
    else:
        model = M.BC5Model(verbose=verbose)
        lo = 1e-20 if args.param_lo is None else args.param_lo
        hi = 1e-3 if args.param_hi is None else args.param_hi
    nsamp = json.loads(args.nsamp) if args.nsamp else None
    ndecay = json.loads(args.ndecay) if args.ndecay else None
    if nsamp:
        nsamp = {k: int(v) for k, v in nsamp.items()}
    samples = RT.Samples(directory=args.samples or RT.SAMPLES_DEFAULT, verbose=verbose)
    import twohdm_lib as _L
    if args.model == '2hdm':
        libdir = os.path.join(RUN, 'work', 'magnus', 'hadron_lib')
    else:
        libdir = os.path.join(os.path.abspath(os.path.join(RUN, '..', 'bc5_lhcb_downstream_v3')),
                              'work', 'magnus', 'hadron_lib')
    if _L.load_hadron_lib(libdir):
        if verbose:
            print('  hadron library: %s' % libdir)
    tab = RT.get_table(model, args.mass, samples, seed=args.seed, nsamp=nsamp, ndecay=ndecay,
                       cache=not args.no_cache, force=args.force, verbose=verbose)
    params = _param_grid(args.param_grid, lo, hi, args.nparam)
    if args.lumi is None or str(args.lumi).strip().lower() == 'both':
        lumis = list(LUMIS)
    else:
        lumis = [float(args.lumi)]
    res = {'model': model.name, 'mass_GeV': float(args.mass), 'param_name': model.param_name,
           'param_lifetime_sign': getattr(model, 'param_lifetime_sign', None),
           's90': S90, 'param_grid': [float(x) for x in params], 'lumis': lumis,
           'b_mode': getattr(model, 'b_mode', None),
           'table_meta': tab['meta'], 'table_channels': list(tab['channels']),
           'table_nsamp': dict(tab['nsamp']), 'table_ndecay': dict(tab['meta'].get('ndecay', {})),
           'table_counts': dict(tab.get('counts', {})),
           'table_n_scalar': int(len(tab['scal_hist'])), 'curves': {}}
    t_eval0 = time.time() if False else None
    for lumi in lumis:
        rows = evaluate_param_grid(tab, model, args.mass, params, lumi_fb=lumi)
        s_list = [r['s_total'] for r in rows]
        crossings = find_crossings(params, s_list,
                                   lifetime_sign=getattr(model, 'param_lifetime_sign', -1))
        # refine every grid-level crossing with a local bisection (grid-independent boundary)
        for c in crossings:
            j = int(np.argmin(np.abs(np.log(np.asarray(params, dtype=float) / c['param']))))
            lo_i = max(0, j - 1)
            hi_i = min(len(params) - 1, j + 1)
            got = refine_crossing(tab, model, args.mass, lumi,
                                  params[lo_i], params[hi_i])
            if got:
                c.update(got)
        # flag whether a crossing lies outside the BR(B->XsH) < 0.1 validity region
        valid = None
        if args.model == '2hdm' and hasattr(model, 'extra_meta'):
            try:
                valid = model.P.validity_boundaries(args.mass)['tanb_BR_gt_0.1']
            except Exception:
                valid = None
        for c in crossings:
            c['outside_validity_B'] = bool(valid is not None and c.get('param_refined', c['param']) < valid)
        res['curves'][str(lumi)] = {'rows': rows, 'crossings': crossings}
        print('  L=%.0f fb^-1: crossings %s' % (
            lumi, ', '.join('%s %.6g' % (c['edge'], c['param']) for c in crossings) or 'none'))
    timing = {'mass_GeV': float(args.mass), 'model': model.name,
              'table_build_s': float(tab['meta'].get('t_build_s', 0.0)),
              'table_response_s': float(tab['meta'].get('t_response_s', 0.0)),
              'n_scalar': int(len(tab['scal_hist'])),
              'n_scalar_nonzero_response': int(tab['meta'].get('n_scalar_nonzero_response', 0)),
              'n_hist': int(len(tab['hist_event'])),
              'cached': bool(tab['meta'].get('t_build_s', 1e9) < 1e-6),
              'curves': {str(l): {'nparam': len(res['curves'][str(l)]['rows']),
                                  'crossings': res['curves'][str(l)]['crossings']}
                         for l in lumis}}
    res['timing'] = timing
    print('TIMING ' + json.dumps(timing))
    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, 'w') as fh:
            json.dump(res, fh, indent=1)
        print('wrote %s' % args.out)
    return res


if __name__ == '__main__':
    main()
