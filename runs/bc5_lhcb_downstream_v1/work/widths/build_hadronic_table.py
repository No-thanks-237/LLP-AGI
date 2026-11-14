# -*- coding: utf-8 -*-
"""
Build the hadronic (pi pi, K K) decay-rate table for the scalar from the digitised
Winkler figure (work/widths/reference_figures/decayrates_curves.csv).

Why this file exists
--------------------
The pi pi and K K rates around the f0(980) require the two-channel Muskhelishvili-Omnes
form factors of Winkler (arXiv:1809.01876), built from the Roy-Steiner pi pi / K K phase
shifts of Hoferichter et al.  The phase-shift tables are not available in this workspace, so
we take the *published* Winkler rates from his own Fig. 2 (decayrates.pdf), digitised with
sub-0.002 px calibration residuals.

Self-calibration
----------------
The same figure also contains the S -> mu+ mu- rate, for which the exact analytic result
Gamma = G_F m_S m_mu^2 beta^3/(4 sqrt2 pi) * theta^2 is known.  We use that curve to
measure the residual scale error of the digitisation and correct the pi pi / K K curves by
the same factor.  This turns the digitisation into a calibrated measurement.

Output: work/widths/hadronic_table.csv  (mass_GeV, gamma_pipi_GeV, gamma_KK_GeV) at theta^2 = 1
"""
import os
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import scalar_widths as sw

FIG = os.path.join(HERE, 'reference_figures', 'decayrates_curves.csv')
OUT = os.path.join(HERE, 'hadronic_table.csv')
REPORT = os.path.join(HERE, 'hadronic_table_report.json')

NGRID = 3000         # 5.7e-4 GeV spacing: resolves the K K threshold (which rises within 5e-4)
M_MIN, M_MAX = 0.2792, 2.0          # pi pi threshold .. the 2 GeV matching scale of [W]


def main():
    d = pd.read_csv(FIG)

    # ---- 1. calibration from the exactly known muon curve -------------------------------
    mu = d[d.curve_id == 'muons'].sort_values('mass_GeV')
    mu = mu[(mu.mass_GeV > 0.4) & (mu.mass_GeV < 7.2)]          # avoid threshold + frame
    exact = np.array([sw.gamma_ll(m, sw.M_MU) for m in mu.mass_GeV])
    ratio = mu.value.to_numpy() / exact
    # a residual stroke-line bias is a pure multiplicative offset on a log axis
    calib = float(np.median(ratio))
    print('muon calibration: median(fig/exact) = %.5f  (scatter %.2e)' % (calib, np.std(ratio)))

    # ---- 2. common mass grid -------------------------------------------------------------
    def curve(cid):
        s = d[d.curve_id == cid].sort_values('mass_GeV')
        # collapse non-monotonic duplicates by averaging over a fine bin
        s = s.groupby(np.round(s.mass_GeV, 5), as_index=False).value.mean()
        return s.mass_GeV.to_numpy(), s.value.to_numpy() / calib

    mp, gp = curve('pions')
    mk, gk = curve('kaons')

    grid = np.exp(np.linspace(np.log(M_MIN), np.log(M_MAX), NGRID))
    gp_g = np.interp(grid, mp, gp, left=0.0, right=float('nan'))
    gk_g = np.interp(grid, mk, gk, left=0.0, right=float('nan'))

    # ---- 3. low-mass completion -----------------------------------------------------------
    # Below ~0.33 GeV the pi pi line in the figure merges into the 5.012e-10 GeV frame floor,
    # so the digitised value there is a lower bound.  Near threshold the form factors are
    # reliably given by lowest-order ChPT ([W] Eq. formfactors), which is what we use there.
    FLOOR = 5.012e-10
    clipped = ~np.isfinite(gp_g) | (gp_g <= 1.5 * FLOOR)
    n_clipped = int(np.sum(clipped))
    chpt_pi = np.array([sw.gamma_pipi_chpt(m) for m in grid])
    gp_g = np.where(clipped, chpt_pi, gp_g)
    # kaons: below their threshold the channel is closed
    chpt_k = np.array([sw.gamma_KK_chpt(m) for m in grid])
    gk_g = np.where(np.isnan(gk_g), chpt_k, gk_g)
    gk_g = np.where(grid < 2 * sw.M_K_CHARGED, 0.0, gk_g)

    # total hadronic rate implied by the figure (for the record)
    multi = np.array([sw.gamma_multimeson(m) for m in grid])

    df = pd.DataFrame({'mass_GeV': grid, 'gamma_pipi_GeV': gp_g, 'gamma_KK_GeV': gk_g})
    df.to_csv(OUT, index=False, float_format='%.8e')
    print('wrote %s  (%d rows)' % (OUT, len(df)))

    # ---- 4. report ------------------------------------------------------------------------
    import json
    rep = {
        'source_figure': 'references/model/arXiv-1809.01876v1/decayrates.pdf (digitised)',
        'digitisation_file': 'reference_figures/decayrates_curves.csv',
        'muon_calibration_factor_median_fig_over_exact': calib,
        'muon_calibration_scatter': float(np.std(ratio)),
        'mass_grid_GeV': [float(grid[0]), float(grid[-1])],
        'clipped_points_replaced_by_chpt': n_clipped,
        'chpt_formula': 'Winkler Eq. (formfactors) lowest order: Gamma_pi^0=m_pi^2, Delta_pi^0=0, Theta_pi^0=s+2m_pi^2 (and kaon analogues)',
        'units': 'GeV at theta^2 = 1',
        'cross_check_vs_senscalc_example': {
            'm_S_GeV': 1.0,
            'winkler_figure_pipi_GeV': float(np.interp(1.0, grid, gp_g)),
            'winkler_figure_KK_GeV': float(np.interp(1.0, grid, gk_g)),
        },
        'note': 'SensCalc tabulated widths differ from the Winkler figure for the hadronic and '
                'gluon channels (see work/widths/reference_comparison.json); the Winkler figure is '
                'used as the production input because prompt.md fixes Winkler/Boarska as the reference.',
    }
    with open(REPORT, 'w') as f:
        json.dump(rep, f, indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == '__main__':
    main()
