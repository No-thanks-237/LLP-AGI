# Independent cross-check of the 2HDM scan aggregation (Lead request, 2026-10-25 22:4x)

Owner: production (read-only on outputs/data and work/analysis; writes only here).
Inputs: outputs/data/signal_scan.json (2640 points), per-mass JSONs
work/analysis/scan_2hdm{,_paper_mb,_exclusive}/m_*.json, work/analysis/cache/response_2hdm_m*_v1.npz,
outputs/data/exclusion_boundary.csv, outputs/data/decay_widths_grid.csv.

## Item 1 - s_total = 1e3*L*sum_g sigma_g*eps_g   -> EXACT

Four masses (nearest scan-grid nodes to 0.3/1.0/3.0/4.5: 0.289535, 1.0, 3.07067, 4.55151)
x tanb grid nodes 105.284357 and 290.733536 x L = 25, 300 (16 combinations):
max |s_recomputed - s_stored| / s_stored = 2.0e-16 (floating-point summation only).
Examples: (mH=1.0, tanb=105.284357): s(25)=3.4940307946e+04 recomputed identical;
s(300)=4.1928369535e+05 identical.

## Item 2 - component check and eps_h convention   -> EXACT / consistent

s_comp[g] = 1e3*L*sigma_group[g]*eps[g] for every group at all 16 points (max diff 0.0);
s_total - sum(s_comp) = 0.0 exactly.

Convention found (work/analysis/models.py lines 491-535, event_engine.py):
  sigma_group['h'] = sigma_h_total * BR(h->HH)     [exactly 1.000000000000000 at 4 points]
  history weight for h = sigma_process[g] / (sigma_h_total * h_frac[g])  (no BR(h->HH))
  eps_h = pure event-OR probability; s_h = 1e3*L*sigma_h_total*BR_hHH*eps_h
Independent verification: I recomputed eps_h from the stored response table
(work/analysis/cache/response_2hdm_m1.000000_v1.npz) with my own OR algebra, my
higgs_sector for the process weights, and the shared p_window kernel:
  eps_h(mine)/eps_h(stored) = 1.00000000 at (mH=1; tanb = 5, 10.7107, 105.284, 290.734),
  where BR(h->HH) spans 3.0e-2 ... 2.7e-7.
=> no double counting of BR(h->HH); the convention is the same as BC5's published eps_H.

Minor observation (0.17 %, below MC, no action): the Pythia h sample covers 6 of the 9
sigma processes (ggF, VBF, WH, ZH, bbH, ttH; h_frac sums to 1) while the weights are
normalised by the full sigma_h_total, so the total h-group weight is
1 - (tHj+tHW+tHb)/sigma_total = 0.99831 and s_h is 0.169 % low.

## Item 3 - production-rate modes vs the variant scans

Production module (frozen v1.0.0) at the scan tanb grid node 105.284357008:
  mH=1.0     : foresee 3.46986580e-3, paper_mb 3.31786152e-3, exclusive 1.88213072e-3
               paper/foresee = 0.956193, excl/foresee = 0.542422
  mH=3.07067 : foresee 1.63443160e-3, paper_mb 7.91083435e-4, exclusive 1.42908978e-3
               paper/foresee = 0.484011, excl/foresee = 0.874365
The scans use exactly these BRs: row['br_B_to_XsH'] - module = 0.0 (all four cases).
Variant/foresee ratios from the per-mass JSONs (L=25):
  mH=1.0     paper_mb s_B=0.95618890 (expected 0.95619304, -4.3e-6); exclusive s_B=0.54241942 (-4.3e-6)
  mH=3.07067 paper_mb s_B=0.48398182 (expected 0.48401134, -6.1e-5); exclusive s_B=0.87431167 (-6.1e-5)
  s_total ratios equal the s_B ratios to 1e-9 (at these points inel = 0 and h <= 5e-9, so
  s is pure B; sigma_group['B'] is the fixed bb cross section, so s_B ratio = rate ratio).
RESIDUAL FULLY EXPLAINED: the sample's B group is 50.32 % B0 (2981 B+ / 3019 B0 histories;
the shared table is identical, seed 20261025, same cache file and sample hashes).  In
mode='foresee' my module uses the per-meson mass in (1-mH^2/mB^2)^2, so
BR(B0)/BR(B+) = 0.9273585 at mH=1 (0.9273504 in paper_mb/exclusive); weighting the two
meson types by the sample composition predicts
   eps ratio paper_mb/foresee = 0.95618899 (mH=1), 0.48398320 (mH=3.07)
   eps ratio exclusive/foresee = 0.54241947 (mH=1), 0.87431415 (mH=3.07)
against the observed 0.95618890, 0.48398182 / 0.54241942, 0.87431167
=> agreement 1e-7 (mH=1) and 2.5e-6 (mH=3), i.e. the variant scans match the production
rate ratio including the B+/B0 mixture to <= 2.5e-6, which is 2000x below the per-point
MC errors (3.2 % at mH=1, 5.6 % at mH=3).

## Item 4 - long-lived edge, mH=1, 25 fb^-1

Boundary: long-lived crossing tanb = 357.03 +- 1.47 (0.41 %); grid neighbours 290.73 / 374.78.
Decay table: Gamma(1; 290.73) = 5.8322e-17 GeV, Gamma(1; 374.78) = 3.5097e-17 GeV;
Gamma ~ tanb^-2 here -> ctau(357.03) = 5.10 m (Gamma = 3.87e-17 GeV).
Representative H energy from B -> K H: E_H(B rest) = (mB^2+mH^2-mK^2)/(2mB) = 2.711 GeV
(p* = 2.52 GeV), i.e. E_H/E_B = 0.514.  From the response cache (6000 B-group H at mH=1):
beta*gamma = p/m distribution: p10 1.46, median 10.74, mean 54.8, p90 131.
=> decay length Lambda = beta*gamma*ctau = 7.5 m (p10) / 54.8 m (median) / 279 m (mean).
**NOT 1-2 m**: at the edge Lambda is ~34x the optimal 1.6 m decay length (median), i.e.
only ~2-3 % of the H decay inside the 1-2.5 m window.
The Lambda ~ 1.6 m condition is met at tanb = 63.4 (ctau = 0.160 m), where s(25 fb^-1)
= 1.58e6, i.e. 7e5 times above threshold.

Why the edge is where it is: s ~ BR_prod(tb^-2) x f_visible(tb^-4) x p_window(tb^-2) = tb^-8.
  * BR(B->XsH) ~ tb^-2 (verified);
  * visible (charged-track) fraction: BR(pi pi) = 1.469e-4 (tb=105) -> 9.149e-7 (tb=375),
    ratio 6.2e-3 = (375/105)^-4 (tb^-4), because all fermionic channels scale as xi_H^f^2
    ~ tb^-6 while Gamma_tot ~ tb^-2 -> at tb=357 BR(gamma gamma) = 1 - 1.3e-6;
  * p_window ~ 1/Lambda ~ tb^-2 in this long-lifetime regime.
Measured slope: s(290.73)/s(374.78) = 7.562 vs (374.78/290.73)^8 = 7.58 (p = 7.967).
Conclusion: the crossing at tanb=357 is real and the scan is self-consistent, but it is
NOT a "decay length ~ detector" edge - at the boundary Lambda ~ 55 m and the LHCb-downstream
charged-track signal dies because the H has become 100 % di-photon, not because the H
decays too late.  Suggested wording for the final report: the mH=1 long-lived edge is
driven by the collapse of the visible branching fraction into gamma gamma (composition),
with the 1/Lambda suppression as a sub-leading tb^-2 factor.

## Differences > 2 % (items 1-2) or > MC (item 3)?
None.  Item 1: 2e-16; item 2: exact + 0.17 % h-weight note; item 3: <= 6.1e-5 vs >= 3.2 % MC.
