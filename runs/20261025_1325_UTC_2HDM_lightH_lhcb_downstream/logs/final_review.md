# Final independent review of the 2HDM light-H scan (T5/model)

Read-only review of the completed scan. Inputs read: outputs/data/signal_scan.json (66 mass x 40
tanb x 2 lumi, b_mode=foresee), outputs/data/exclusion_boundary.csv, outputs/data/
exclusion_s_curves.csv, work/analysis/scan_2hdm/m_*.json, work/common/make_curves.py (to identify
the exact metric convention), work/analysis/event_engine.py, work/response/twohdm_lib.py,
work/magnus/samples/bbar_heavy.txt.gz (B kinematics).  Nothing outside work/model/audit/** and
logs/ was written.

Review scripts: work/model/audit/final_review_check12.py, final_review_check2.py,
final_review_check3.py, final_review_check34.py.

## 1. Crossing positions vs the s(tanb) grid  -- PASS

Method: for each (mH, L) take s(tanb) from signal_scan.json, find every sign change of
s - s90 (s90 = -ln 0.1 = 2.302585), interpolate linearly in (ln s, ln tanb).

| mH | L [fb^-1] | grid crossing tanb | dlns/dlntanb | type | CSV tanb | rel. diff |
|---|---|---|---|---|---|---|
| 0.289535 (0.3) | 25 | 179.3913 | -7.999 | long | 179.391 | 1.6e-06 |
| 0.289535 | 300 | 244.7423 | -8.000 | long | 244.742 | 1.4e-06 |
| 1.000000 | 25 | 6.0601 | +68.324 | short | 6.0601 | 4.0e-07 |
| 1.000000 | 25 | 357.0292 | -7.967 | long | 357.029 | 4.3e-07 |
| 1.000000 | 300 | 5.8437 | +68.324 | short | 5.8437 | 3.4e-08 |
| 1.000000 | 300 | 487.4913 | -7.988 | long | 487.491 | 5.3e-07 |
| 3.50112 (3.5) | 25 | 24.8203 | +12.384 | short | 24.8203 | 1.9e-06 |
| 3.50112 | 25 | 273.3656 | -5.453 | long | 273.366 | 1.5e-06 |
| 3.50112 | 300 | 21.5034 | +23.309 | short | 21.5034 | 1.5e-06 |
| 3.50112 | 300 | 410.7682 | -6.621 | long | 410.768 | 5.0e-07 |

n_crossings agrees everywhere, and the edge classification is verified by the sign of
dlns/dlntanb at each crossing: s falls with tanb at all long-lived edges (-5.4 ... -8.0) and
rises at all short-lived edges (+12 ... +68).  Maximum relative difference in tanb: 1.6e-06
(median 1.4e-06), i.e. the CSV is exactly the ln-ln interpolation of the published grid.

(The long-lived slopes are close to -8: Br(B->X_s H) ~ tanb^-6 from xi_H^f ~ 1/tanb^3 and the
lifetime acceptance ~ 1/(c tau) ~ tanb^-2 in the long-lifetime regime.)

## 2. MC-error metric -- PASS (exact reproduction)

make_curves.boundary_table() records, for each crossing, the LOWER bracket grid point i and
computes  rel_ds = ds[i]/s[i]  and  err = tanb_cross * rel_ds / |dlns/dlntanb|_centred(i-1,i+1).
Re-deriving the CSV from signal_scan.json with exactly this convention reproduces it:

| mH | L | rel_ds_at_long (CSV) | ds/s at bracket i | tanb err (CSV) | tanb*(ds/s)/|slope| | rel diff |
|---|---|---|---|---|---|---|
| 0.289535 | 25 | 0.034482 | 0.034482 | 0.77332 | 0.77332 | 3.4e-07 |
| 0.289535 | 300 | 0.034484 | 0.034484 | 1.05502 | 1.05502 | 1.7e-06 |
| 1.000000 | 25 | 0.032715 | 0.032715 | 1.46804 | 1.46804 | 3.3e-06 |
| 1.000000 | 300 | 0.032780 | 0.032780 | 2.00144 | 2.00144 | 2.3e-06 |
| 3.50112 | 25 | 0.037659 | 0.037659 | 2.00806 | 2.00806 | 2.1e-06 |
| 3.50112 | 300 | 0.032717 | 0.032717 | 2.11754 | 2.11754 | 2.2e-06 |

So rel_ds_at_long is the standard error of s divided by s at the crossing bracket point
(not ds/s90: that would be 2.303x smaller, e.g. 0.0142 instead of 0.0327 at mH=1), and the
quoted tanb error is tanb*(ds/s)/|dln s/dln tanb| as documented.  Max mismatch 3.3e-06.
nb the short-lived errors use the same formula but their rel_ds is 0.64-0.86, i.e. an
~65-86% MC error on s (see the flag in section 5).

## 3. Lifetime consistency of the long-lived edge  -- quantified, relation explained

c tau from the scan (interpolated in ln tanb), and beta*gamma_z = pz/mH of the H from the real
B sample (work/magnus/samples/bbar_heavy.txt.gz, 1.36M B+/B0 rows, median |p_B| = 26 GeV,
23% above 100 GeV, 9.2% above 300 GeV, max 4.7 TeV), B -> K H two-body decay with the correct
Lorentz boost, selecting 2 < eta_H < 5 (22-25% of B decays):

| mH | L | edge tanb | c tau [m] | beta gamma_z (p10/med/p90) | c tau beta gamma_z [m] | / window(1.5 m) |
|---|---|---|---|---|---|---|
| 1.000 | 25 | 357.029 | 5.10 | 8.2 / 30.7 / 130 | 157 (42 - 664) | 105 |
| 1.000 | 300 | 487.491 | 9.51 | 8.2 / 30.9 / 131 | 294 (78 - 1250) | 196 |
| 3.501 | 25 | 273.366 | 0.0697 | 2.9 / 10.2 / 41 | 0.71 (0.20 - 2.86) | 0.47 |
| 3.501 | 300 | 410.768 | 0.1573 | 2.9 / 10.2 / 41 | 1.60 (0.45 - 6.52) | 1.07 |
| 0.2895 | 25 | 179.391 | 53.1 | 30 / 110 / 456 | 5.8e3 | 3.9e3 |
| 0.021 | 25 | 21.875 | 2071 | 417 / 1515 / 6326 | 3.1e6 | 2.1e6 |

So "c tau beta gamma ~ window" holds only for the heaviest H: at mH=3.5 the long-lived edge sits
exactly at the decay-length scale (0.7 m at 25 fb^-1, 1.6 m at 300 fb^-1, i.e. 0.5-1.1x the
1-2.5 m window).  At mH=1 it is already 1e2x the window and at mH<=0.3 it is 1e3-1e6x.
The relation: s = 1e3 L sum_c sigma_c eps_c with eps_B ~ Br(B->X_s H) * <A p_window>,
Br ~ tanb^-6 (xi_H^f ~ 1/tanb^3) and p_window ~ 1/(c tau beta gamma) ~ tanb^-2 in the
long-lifetime regime, hence s ~ tanb^-8.  The acceptance is maximised where
c tau beta gamma ~ (z2-z1) = 1.5 m, which for mH=1 happens at tanb ~ 35 (and tanb ~ 393 for
mH=3.5).  When the maximum signal is enormous (s_max/s90 = 1.5e9 at mH=1, 1e11 at mH=0.29) the
crossing is pushed far down the tanb^-8 tail, i.e. to c tau beta gamma >> window; when the
maximum signal is only ~100x threshold (mH=3.5: s_max/s90 = 134) the edges lie near the
acceptance maximum, i.e. at the window scale.  Both behaviours are consistent with the same
estimator; the window scale is the position of the acceptance peak, not of the long-lived edge.
No inconsistency with the scan.

## 4. Low-mass end (mH = 0.021 GeV) -- PASS

From the FORESEE-normalised c tau of the scan (c tau ~ tanb^2 in the diphoton-dominated regime:
107.5 m at tanb=5, 180.4 at 6.45, 300 at 8.31, 498 at 10.7, 826 at 13.8) and beta gamma_z
(median 1495 for 2<eta<5):
 * c tau beta gamma_z at tanb=5 is 1.6e5 m, i.e. 1e5x the decay window; s(5) = 3.16e5
   (25 fb^-1) and 3.79e6 (300 fb^-1), both >> s90 = 2.30.
 * the decay length would be ~1 m only at tanb ~ 0.015 (extrapolating c tau beta gamma_z ~ tanb^2:
   tanb = 5*sqrt(1/1.6e5) = 0.0125; for a 2.5 m window 0.02), i.e. two orders of magnitude below
   the scan range 5...1e5.
 Therefore only one crossing exists in the scanned range - the long-lived edge at 21.88 (25 fb^-1)
 and 29.84 (300 fb^-1) - exactly as the boundary table reports (n_crossings = 1).  The short-lived
 edge for this mass is far below tanb = 5 and is correctly absent.

## 5. FLAG: the short-lived edges are not physical boundaries

The short-lived crossings (mH = 1: 6.060/5.844; mH = 3.5: 24.82/21.50) are produced by the
extreme high-momentum tail of the B sample, not by the bulk of the production:
 * at mH=1, tanb=6.06, Br(B->X_s H) = 1.004 (and 1.08 at 300 fb^-1): the inclusive formula is
   saturated, Br > 1, i.e. the point is outside the model's validity (DESIGN 3.1b flags
   Br > 0.1; validity.tanb_BR_gt_0.1 = 20.3 at mH=1).  For mH=3.5 the short-lived points have
   Br = 0.006-0.008, so they are formally inside the validity bound.
 * the required c tau beta gamma is ~1 mm or less (mH=1: c tau = 1.0e-4 m, beta gamma_z ~ 30),
   so p_window = exp(-1/(c tau beta gamma_z)) is nonzero only for the handful of B mesons with
   pz >~ 400 GeV (9.2% of the sample has |p| > 300 GeV, 1.7% > 1 TeV, max 4.7 TeV): for
   c tau = 1.01e-4 m, p_window = 3.2e-22 at pz=200 GeV, 4.6e-15 at 300, 1.8e-11 at 400,
   5.0e-05 at 1000 GeV.
 * consequently the MC error on s at those points is ds/s = 0.64-0.99 (the boundary table
   quotes them as tanb_short with 4-9% errors, but rel_ds_at_short is 65-99%).
 Recommendation: do not present the low-tanb (short-lived) side as a physical contour.  The
 physically meaningful boundary of this model/analysis is the long-lived edge (tanb_long); on the
 low-tanb side the model is either invalid (Br(B->X_s H) > 1 at mH <~ 1 GeV) or the "signal" is
 driven by the extreme momentum tail with O(1) MC uncertainty.  The delivered
 outputs/data/exclusion_boundary.csv does not carry the validity flag column that
 work/common/final_analysis.py computes (short_outside_validity); if the figure closes the
 contour at low tanb it should be drawn/shaded as "model invalid / not used".

## 6. Optional item: paper reach overlay -- done

LightH.pdf (the paper's right panel) contains the FASER (blue dashed) and FASER2 (red) reach
curves as vector drawings.  Digitised by rendering the page at 4x (fitz, 1728x2304) and selecting
the pure red / pure blue pixels, with the axes calibrated from the PDF text spans:
mH = 0.1, 1, 10 at px = 68.45, 313.35, 558.15 (244.85 pt/decade) and tanb = 10, 100, 1000, 1e4 at
py = 337.4, 254.4, 171.1, 87.6 (83.27 pt/decade).  Script: work/model/audit/
final_review_digitise.py, figure: work/model/audit/paper_reach_overlay.png.

Digitised paper reach (tanb at fixed mH):
  FASER  : 39 (mH=0.1), 133 (0.3), 444 (1.0), 561 (3.0)
  FASER2 : 459 (0.1), 1330 (0.3), 4324 (1.0), ~1e4 plateau up to mH~3.9 then closure
This work (long-lived edge): 179-245 (mH=0.29), 357/487 (mH=1), 273/411 (mH=3.5).
So the LHCb-downstream reach lies above the paper FASER reach (by ~1-4x in tanb) and below the
FASER2 reach (by ~1-2 orders of magnitude) - qualitatively consistent with a 25-300 fb^-1
downstream experiment smaller than FASER2.  The kinks at mH ~ 0.27-0.35 (pipi/KK thresholds) and
~0.5-1 GeV in our curve are visible and match the expected threshold structure.

## Verdicts

1. Crossings and edge classification: PASS (1.6e-06 max deviation, signs verified).
2. MC-error metric: PASS (exact, 3.3e-06 max; ds/s, not ds/s90; tanb error formula as documented).
3. Long-lived-edge lifetime: consistent; c tau beta gamma spans 0.7 m (mH=3.5) to 3e6 m (mH=0.02)
   and the relation to the window is explained by s ~ Br * p_window ~ tanb^-8.
4. Low-mass end: PASS (single long-lived crossing; short-lived edge at tanb ~ 0.015 << 5).
5. FLAG: the short-lived edges are tail-driven / model-invalid with 65-99% MC errors and should be
   excluded from, or explicitly flagged in, the final exclusion figure.
