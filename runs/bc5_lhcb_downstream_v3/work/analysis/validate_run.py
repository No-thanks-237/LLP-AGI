# -*- coding: utf-8 -*-
"""Post-hoc validation of the run: sample consistency, MC error, lifetime-integral convergence."""
from __future__ import annotations
import json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
import run_scan as R
import bc5_lib as L
import scalar_widths as sw
import bc5_production as P


def main():
    rep = {}
    # ---- 1. b bbar sample: forward acceptance versus the LHCb measurement ---------------
    smp = R.load_sample(R.SAMPLE_B)
    if smp is not None:
        eta = 0.5 * np.log((np.sqrt(smp['px']**2 + smp['py']**2 + smp['pz']**2) + smp['pz']) /
                           np.maximum(np.sqrt(smp['px']**2 + smp['py']**2 + smp['pz']**2) - smp['pz'], 1e-12))
        frac_acc = float(np.mean((eta > 2.0) & (eta < 5.0)))
        rep['bbar_sample'] = {
            'n_b_mesons_kept': int(len(eta)),
            'eta_min': float(eta.min()), 'eta_max': float(eta.max()),
            'fraction_of_kept_in_2_lt_eta_lt_5': frac_acc,
            'sigma_bb_total_pb_assumed': R.SIGMA_BB_PB,
            'lhcbl_forward_reference_ub': 79.8,
            'note': 'The LHCb measurement sigma(pp -> b bbar X; 2<eta<5) = 79.8 ub cannot be used '
                    'directly as the normalisation because the Pythia sample is already eta-cut at '
                    'generation (1.5 < eta < 5.5), so the 4pi b bbar cross section is the correct '
                    'normalisation. The forward/4pi ratio implied here is therefore a lower bound '
                    'on the true acceptance fraction and is quoted for information only.',
            'open_item': 'replace the assumed 500 ub 4pi b bbar cross section by a documented '
                         'FONLL/NNLO value and rerun without the generation-level eta cut',
        }
        # ---- 2. lifetime-integral convergence: extend the grid and compare ---------------
        # Quadrature convergence of the lifetime integral with the SAME A(z): evaluate A on a
        # 25-node grid and compare the trapezoidal integral on 25 nodes against every other node
        # (13 nodes, the production grid).  Using the same A isolates the quadrature error from
        # the Monte-Carlo noise of A itself.
        ms = 2.0
        zgrid25 = np.linspace(1.0, 2.5, 25)
        z13 = zgrid25[::2]
        import random
        rng = random.Random(1)
        sel = np.where(smp['pdg'] == 521)[0][:60]
        devs = []
        for i in sel:
            mes = {'kind': 'B', 'bmeson': 'B+', 'x': smp['x'][i], 'y': smp['y'][i], 'z': smp['z'][i],
                   'e': smp['e'][i], 'px': smp['px'][i], 'py': smp['py'][i], 'pz': smp['pz'][i]}
            ss = R.sample_meson_decay_S(ms, mes, rng)
            if not ss:
                continue
            s = ss[0]
            p4 = (s[0], s[1], s[2], s[3])
            A25 = L.acceptance_profile(ms, p4, (s[4], s[5], s[6]), rng, ndecay=40, zgrid=zgrid25)
            ctau = sw.HBAR_C / (sw.gamma_total(ms) * 1e-8)
            p13 = L.lifetime_probability(ms, p4, s[6] * 1e-3, A25[::2], ctau, zgrid=z13)
            p25 = L.lifetime_probability(ms, p4, s[6] * 1e-3, A25, ctau, zgrid=zgrid25)
            if p25 > 1e-6:
                devs.append(abs(p13 - p25) / p25)
        rep['lifetime_integral_convergence'] = {
            'n': len(devs), 'median_rel_dev_13vs25_nodes': float(np.median(devs)) if devs else None,
            'max_rel_dev': float(np.max(devs)) if devs else None,
            'note': 'same A(z) sampled on 25 nodes; the trapezoidal lifetime integral on the 13-node '
                    'production grid is compared with the 25-node result over 1 < z_dec < 2.5 m',
        }

    # ---- 3. widths cross-check ---------------------------------------------------------
    wpath = os.path.join(HERE, '..', 'widths', 'reference_comparison.json')
    if os.path.exists(wpath):
        rep['widths_reference_comparison'] = json.load(open(wpath))
    ppath = os.path.join(HERE, 'production_validation.json')
    if os.path.exists(ppath):
        rep['production_reference_comparison'] = json.load(open(ppath))
    hpath = os.path.join(HERE, 'hadronisation_validation.json')
    if os.path.exists(hpath):
        rep['hadronisation_library'] = json.load(open(hpath))

    out = os.path.join(RUN, 'outputs', 'data', 'run_validation.json')
    with open(out, 'w') as f:
        json.dump(rep, f, indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == '__main__':
    main()
