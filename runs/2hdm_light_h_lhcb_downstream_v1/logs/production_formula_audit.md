# Production formula audit — Type-I 2HDM light H vs Li et al. arXiv:2212.06186

Task: **task-7** (production-audit).  Scope: independent, formula-by-formula audit of the
**production** module against the "Light_scalar" paper.
Owner of the code: production.  This audit edits no physics code and rebuilds no response table;
the Lead applies any fix.

| item | value |
|---|---|
| audited module | work/production/twohdm_production.py, sha256 `6381bb26049be4c241d80c9b0a43787115897edf3072f4ea969d6fcb1c868c79` |
| audited helper | work/production/b_channels.py, sha256 `6cf70445852fe027fe45d865dc4a1ce94fc9a11bfae26ffcbfa5f758f8651c4d` |
| paper | references/model/arXiv-2212.06186v1/FASER_light.tex (line numbers below refer to this file) |
| paper numerics | references/tools/FORESEE/Models/2HDM_H/2HDM_H.ipynb (cells 6-7) |
| Boiarska cross-check | references/model/arXiv-1904.10447v4/main.tex |
| audit script | work/review_check/prodformula_audit.py |
| machine-readable results | work/review_check/production_formula_audit.json |
| run date | 2026-10-26 |

**Verdict.** No physics formula error was found: every production identity reproduces the paper
equation to machine precision (worst relative difference 3.1e-16 over 75 App.-D numbers, 112
complex xi_H values, both B-rate modes, the K amplitude/BR, the semileptonic rates and the h
sector).  The material audit findings are **prescription** choices, not algebra errors:
one HIGH (the nominal B -> X_s H rate does not use the phase space of the printed formula), one
MEDIUM (CKM parameter set), and several LOW/NOTE items that cannot move the long-lived boundary.

---

## 1. The key question: which B -> X_s H prescription matches the paper formula literally?

**'paper_mb' (m_b = 4.18 GeV in the (1 - mH^2/m_ps^2)^2 phase space) is the literal match.**
The printed equation (lines 177-180) is

    Br(B->Xs phi)/Br(B->Xc e nu) = 12 pi^2 v^2/m_b^2 (1-m_phi^2/m_b^2)^2 / f(m_c^2/m_b^2) |xi^{bs}/V_cb|^2

and the only m_b the paper states is m_b = 4.18 GeV (line 299).  The code's nominal mode
'foresee' uses m_ps = m_B = 5.27934 GeV (and imposes the mesonic endpoint m_B - m_K = 4.78566 GeV).

The paper is **internally inconsistent** here: its own numerical model (FORESEE 2HDM_H notebook,
the model that produced the paper's reach figure) uses the m_B spectator phase space
BR = 9.685e-13 (mB/Gamma_B)(1-mH^2/mB^2)^2 F(tanb)^2, and the text claims "m_H reach up to m_B"
(line 757).  The code's 'foresee' mode reproduces that FORESEE rate to **1.65%** flat in mH
(check C11), while 'paper_mb' reproduces the printed equation to machine zero.

-> Report, as the literal paper-formula result, mode **'paper_mb'**; report 'foresee' explicitly
as the paper-analysis/FORESEE variant.  They agree at low mH but diverge catastrophically near the
m_b threshold (rate ratio 0.96 at mH = 1, 0.042 at mH = 4, 0 above 4.18).

---

## 2. Formula-by-formula audit table (confirmed matches)

All "max |rel diff|" values are from an independent re-implementation of the paper LaTeX in
work/review_check/prodformula_audit.py, evaluated on grids of mH and tanb (see JSON for the grids).

| # | quantity | paper equation | code | independent check | max rel diff | verdict |
|---|---|---|---|---|---|---|
| 1 | App. D functions g0,g1,g2,X1,X2 (Type-I) | lines 919-929 | 201-242 | fresh transcription, 5 tanb x 3 quarks x 5 functions | 3.1e-16 | MATCH |
| 2 | lambda_HH+H- (lambda v^2 = 0) | 905-915 | 245-248 | identity; feeds the -g0 (2v/mW^2)lambda term | 0 | MATCH |
| 3 | xi_H^{ij}, i,j in {d,s,b}, complex | 710-711 | 278-322 | fresh implementation, 7 mH x 4 tanb x 4 pairs | 0.0 | MATCH |
| 4 | B->XsH ratio and (1-mH^2/m_b^2)^2 phase space + f(mc^2/mb^2) | 177-180, 299 | 328-365 (mode='paper_mb') | fresh implementation | 0.0 | MATCH |
| 5 | B->XsH nominal ('foresee', m_B phase space) | not the printed formula | 338-365, 106 | same formula with m_B; vs FORESEE model | 0.0 (mode) / 1.65% (FORESEE) | DIFFERENT PRESCRIPTION -> F1 |
| 6 | B0 rate: BR(B0)=0.92735 BR(B+) | line 180 literal says 0.104 for both | 108-112, 346, 365 | B0/B+ measured 0.92736 | 8.7e-6 vs 0.92735; 7.8% vs literal text | DOCUMENTED CHOICE -> F3 |
| 7 | M(K->pi H), Eq. Matrix_4FH (coherent 4-quark + FCNC) | 193-205 | 499-544 | complex amplitude and BR, 40 points, K+ and K0L | 0.0 | MATCH |
| 8 | Br(pi/K -> H e nu), (1-2n_h/(33-2n_l))^2=(7/9)^2 | 220-224 | 557-583 | fresh implementation | 0.0 | MATCH |
| 9 | g_hHH, Eq. ghHH (lambda v^2 = 0) | 892-895 | 609-619 | fresh implementation | 0.0 | MATCH |
| 10 | Gamma(h->HH) and BR(h->HH) | 609-612 | 664-692 | Gamma = g^2 beta/(8 pi m_h) (the 4700 convention) | 2.1e-16 | MATCH to equations; prose inconsistent -> F5 |
| 11 | kappa_g, Eq. coup_hgg (t,b loops) | 149-151 + App. A | 622-636 | Type-I identity kappa_g = xi_h^f | 2.2e-16 | MATCH |
| 12 | kappa_gamma, Eq. coup_hgaga (W + t,b,c,tau) | 152-154 | 639-661 | fresh implementation of A1, A12 | 1.1e-16 | MATCH; H+- loop omitted -> F4 |
| 13 | 2HDM h production scaling per process (kappa framework), sigma_total -> 55.4527 pb | 517-519, 632-634 (couplings) | 695-708 | per process vs kappa^2, SM limit | 4.4e-16 | MATCH (framework is DESIGN, not a paper equation) |
| 14 | CKM matrix (unitary, PDG standard parametrisation) | paper gives no CKM values | 77-83, 135-170 | unitarity 1.1e-16, GIM 2e-18 | - | PARAMETER CHOICE -> F2 |
| 15 | Boiarska exclusive b->s shape (rate normalisation only) | Boiarska eq. productionBRmesons | b_channels 106-182 | fresh master-formula implementation; sum at mS=0 | 0.0 (3.349 vs 3.3 theta^2, +1.5%) | MATCH (shape only, never added to the inclusive rate) |

Limit/consistency checks (all pass):

* xi_H^{ij} -> 0 as tanb -> infinity: measured local log-slope **-1.00000002** (i.e. exactly
  1/tanb) between tanb = 1e3 and 1e4; |xi(1e5)|/|xi(1e3)| = 0.010.  (Note: the tree-level Yukawa
  modifier xi_H^f ~ 1/tanb^3 is a different object; DESIGN line 90 is correct.)
* Br(B->XsH) monotonic decreasing in tanb: 0 violations over 40 log-spaced points in 5..1e5 at
  mH = 0.1, 1, 2, 4.
* App. D analytic limits: g1(x->0) = -3/4 (residual -3.2e-14); g0(x->0) = -3 cot^2b/(16 x_H)
  (rel. 4.3e-10).
* Paper's own FORESEE bracket F(c) (top-quark bracket, independent functional form) vs exact
  bracket: ratio 1.00699 ... 1.00711 (0.7%, flat in tanb) - this is their fit quality.
* B -> X_s H exclusive/inclusive: our Boiarska theta^2 sum at mS=0 is 3.349 vs their 3.3 (+1.5%);
  inclusive/exclusive ~ 1.96 at mH = 0.1 (their factor-2 model ambiguity).
* K -> pi H vs Boiarska master formula: FCNC-only ratio **0.99667** (pure quark-mass constants
  m_s = 0.096/0.0047 vs 0.0924/0.00467); the four-quark term is +0.13% of the rate but -6.8%
  interference (destructive), so total/FCNC = 0.933 at mH = 0.
* h -> HH vs paper: |g_hHH| matches the paper prose approximation -mh^2 c^2/(4v) to 1% at
  mH -> 0 (the paper's sign is opposite - irrelevant for Br); Br(h->HH) = 0.969 * 19.7 c^4 at
  (mH=0.021, tanb=5), versus the paper prose 75 c^4 (factor 3.8 - paper inconsistency, F5).
* Higgs-sector scaling: ggF uses kappa_g^2, VBF/WH/ZH kappa_V^2, ttH/tHj/tHb/bbH kappa_f^2,
  tHW the DESIGN interference bracket (kt^2+kV^2+2 kt kV cos(phi))/4 with cos(phi)=1 and the
  [0,pi] range reported.  All reproduce the module exactly; tHW is 2.7e-4 of sigma_h, so its
  phi choice cannot move the boundary (< 1e-3 relative in s).
* kappa_gamma omitting the H+ loop (kappa_gamma would go from 1.000 to 1.167 at large
  tanb) changes Gamma_h_other by 0.08% only, because gamma-gamma is 0.227% of it.

---

## 3. B -> X_s H rate prescriptions at tanb = 100 (absolute BR, B+)

| mH [GeV] | foresee (nominal) | paper_mb (literal Eq) | exclusive (Boiarska) | paper_mb/foresee | exclusive/foresee |
|---|---|---|---|---|---|
| 0.1 | 4.13491e-3 | 4.13314e-3 | 2.11162e-3 | 0.9996 | 0.511 |
| 1.0 | 3.84628e-3 | 3.67779e-3 | 2.08631e-3 | 0.9562 | 0.542 |
| 2.07163 | 2.96168e-3 | 2.35479e-3 | 1.98081e-3 | 0.7951 | 0.669 |
| 3.0 | 1.89701e-3 | 9.72939e-4 | 1.64692e-3 | 0.5129 | 0.868 |
| 3.99192 | 7.58887e-4 | 3.20189e-5 | 5.64192e-4 | 0.0422 | 0.743 |
| 4.18 | 5.76030e-4 | 0 (closed) | 4.42199e-4 | 0 | 0.768 |

---

## 4. Impact on the 90% CL long-lived tanb boundary

Method: s = 1e3 L sigma eps is exactly linear in each production branching ratio
(s_new = s_B f_B + s_inel f_inel + s_h f_h).  The delivered signal_scan curves
(work/analysis/scan_2hdm, work/analysis/scan_2hdm_paper_mb, scan_2hdm_exclusive) are re-weighted
row by row with production-module rate ratios and the 90% CL crossing s = 2.302585 is re-solved.
No response table is rebuilt.  **Validation:** re-weighting the foresee rows to paper_mb/exclusive
reproduces the Lead's independent variant scans to **4.6e-5** relative; the crossing finder
agrees with the Lead's refined crossings to 2e-5 (mH=1) and <= 3e-3 (steep nodes near 4 GeV).

Absolute 90% CL long-lived tanb boundary (s = 2.302585, zero background):

| mH [GeV] | foresee (nominal) L=25 / 300 | paper_mb L=25 / 300 | exclusive L=25 / 300 |
|---|---|---|---|
| 1.0 | 357.03 / 487.49 | 355.03 / 484.77 | 330.64 / 451.52 |
| 1.81693 | 298.76 / 411.68 | 292.37 / 403.01 | 281.01 / 387.84 |
| 2.07163 | 297.24 / 412.41 | 288.25 / 400.33 | 281.48 / 391.45 |
| 3.7287 | 271.43 / 415.45 | 186.98 / 308.10 | 259.54 / 400.75 |
| 3.99192 | 298.91 / 457.01 | 143.12 / 260.87 | 282.33 / 435.67 |
| 4.18 | 303.65 / 468.06 | **no crossing** (rate = 0) | 288.87 / 448.11 |

Relative shift of the long-lived boundary (this deviation minus nominal, in %):

| deviation from the delivered nominal | mH = 1 (L=25 / 300) | mH = 2 (L=25 / 300) | mH = 4 (L=25 / 300) |
|---|---|---|---|
| **paper_mb** (literal paper m_b phase space) | **-0.56 / -0.56** | **-2.79 / -2.71** | **-52.6 / -43.4** |
| **exclusive** (Boiarska channel sum rate) | -7.39 / -7.38 | -5.47 / -5.27 | -5.58 / -4.70 |
| CKM: DESIGN nominal magnitudes (Vts=0.0415) instead of unitary (0.040122) | +0.85 / +0.85 | +0.89 / +0.87 | +1.22 / +1.09 |
| m_t = 172.5 instead of 173.1 GeV | -0.12 / -0.12 | -0.12 / -0.12 | -0.17 / -0.15 |
| m_b = 4.2 instead of 4.18 GeV | -0.20 / -0.20 | -0.21 / -0.20 | -0.28 / -0.26 |
| m_c = 1.3 instead of 1.28 GeV | +0.26 / +0.26 | +0.28 / +0.27 | +0.38 / +0.34 |
| B0 rate = B+ rate (literal paper line 180) | +0.48 / +0.48 | +0.51 / +0.50 | +0.69 / +0.62 |
| h -> HH width x4 (L = -1/2 convention) | 0 / 0 | 0 / 0 | 0 / 0 |
| Br(h->HH) = 75 c^4 (paper prose) | 0 / 0 | 0 / 0 | 0 / 0 |
| kappa_gamma with the H+- loop | 0 / 0 | 0 / 0 | 0 / 0 |
| any K/pi formula deviation (< 20%) | 0 / 0 | 0 / 0 | 0 / 0 |

mH = 2 is ln-ln interpolated between the scan nodes 1.81693 and 2.07163 (the grid deliberately
has no 2.0 node: the author decay program doubles the hadronic channels exactly at 2 GeV).  mH = 4
is a short (0.2%) extrapolation from the nodes 3.7287 and 3.99192; the next node, 4.18, is exactly
where the literal formula closes, so no bracketing node exists.  At the exact mH = 3.99192 node the
paper_mb shift is -52.1% / -42.9%, consistent with the mH = 4 entry (-52.6% / -43.4%).

**Why the K/pi/h entries are zero:** at *every* long-lived boundary the B channel carries all of s.
The measured fractions are s_B/s = 1 - 1.8e-11 (mH = 1), 1 - 1.4e-9 (mH = 2), 1 - 2.1e-7
(mH = 4.18); s_inel/s = 0.0 exactly (K -> pi H closes at 0.354 GeV, pi -> H e nu at 0.140 GeV,
K -> H e nu at 0.494 GeV, so the whole inelastic channel is closed at mH >= 0.5), and s_h/s is at
most 4e-7.  The local logarithmic slope at the boundary is dln s/dln tanb = -8.0 (mH = 1),
-7.4 (mH = 2.07), -5.2 (mH = 3.99) at L = 25.  Hence B-rate changes map onto tanb with
Delta tanb/tanb ~ Delta BR/BR / 8 at low mass, and the h/K/pi formulas cannot move this boundary
at all.  The short-lived (low-tanb, outside-validity) crossing does move: h -> HH x4 shifts it by
-0.5% (mH=1), -0.9% (mH=2.07) (the lifetime turn-on is extremely steep), and the paper_mb mode
shifts it by +34% at mH = 4 (indicative only; the region has Br(B->XsH) >> 1).

---

## 5. Findings, severity-ranked

### F1 (HIGH) - nominal B -> X_s H rate does not use the printed formula's phase space
Paper: lines 177-180 (m_b) and line 299 (m_b = 4.18 GeV).
Code: work/production/twohdm_production.py:106 (DEFAULT_B_MODE='foresee'), 338-365
(_b_ratio_and_total), 368-407.
The nominal mode replaces (1-mH^2/m_b^2)^2 by (1-mH^2/m_B^2)^2.  Consequences (section 4):
boundary shift -0.56% at mH=1, -2.8% at mH=2, **-52.6% (L=25) / -43.4% (L=300) at mH=4**;
above mH = 4.18 the literal formula has *no* long-lived crossing at all, while the nominal mode
still gives 304 (L=25) / 468 (L=300) at mH = 4.18.
The paper's own analysis (FORESEE) and its "reach up to m_B" sentence use m_B, so the code's
choice is not arbitrary - but it is **not** the printed formula, and the paper never reconciles
the two.  Also note the code adds a mesonic endpoint m_B - m_K that neither version of the paper has.
Recommendation: report mode='paper_mb' as the literal-formula result (and the headline for a
formula-reproduction study), keep 'foresee' labelled as the FORESEE/paper-figure variant; never
present the two as a "systematic band" without stating that at mH >= 4 the difference is a factor
~25 in rate.

### F2 (MEDIUM) - CKM parameter set
Paper gives no CKM values (only |Vcb| in the formula).  Code: 77-83, 135-170.
The exactly-unitary set has |Vts| = 0.040122, |Vtd| = 0.008331, |Vcb| = 0.040800; the DESIGN.md
nominal magnitudes are |Vts| = 0.0415, |Vtd| = 0.0086.  Effect: |xi_H^{bs}| +3.35% -> BR(B->XsH)
+6.95% -> long-lived boundary **+0.85% (mH=1), +0.89% (mH=2), +1.2% (mH=4)**.  K->pi H moves
+13.2% (no boundary effect, channel closed).  The unitary set is the more defensible choice; the
audit only asks that the choice and its ~1% boundary effect be stated.

### F3 (LOW) - B0 lifetime scaling vs the literal paper sentence
Paper line 180 says Br(B->Xc e nu) = 0.104 "for both B0 and B+-", which read literally gives
BR(B0->XsH) = BR(B+->XsH).  Code uses BR(B0) = 0.92735 BR(B+) (108-112), matching the Boiarska
caption ("multiply by 0.93 for B0") and FORESEE (per-meson mB/Gamma_B).  Literal-paper variant
would raise s_B by 3.9% (B0 fraction of the sample 0.5008) -> boundary +0.48% (mH=1) to +0.69%
(mH=4).  Recommendation: keep the code and cite Boiarska/FORESEE, not the sentence.

### F4 (LOW) - kappa_gamma has no charged-Higgs loop
Paper eq. coup_hgaga (152-154) plus the H+- term eq. coup_hgaga_Hpm (522-524); the h H+H- coupling
is not suppressed, so a charged-Higgs loop also exists for h -> gamma gamma.  Code: 639-661.
Including the term raises kappa_gamma from 1.000 to 1.167 at large tanb (A0(tau_Hpm) -> 1/6),
i.e. Gamma(h->gamma gamma) +36%; but Gamma_h_other changes by only **0.08%** (gamma-gamma is
0.227% of it) and the long-lived boundary shift is **< 1e-3 %** (h channel is < 4e-7 of s).
No action needed for this analysis.

### F5 (NOTE) - paper's own h -> HH prose is inconsistent; code follows the equations
The paper prose (line 625) quotes g_hHH ~ -mh^2 c^2/(4v) and Br(h->HH) ~ 75 c^4, but its own
App. C eq. ghHH gives |g_hHH| = mh^2 c^2/(4v) (sign opposite) and the 4700 (g/v)^2 relation
(609-612) gives Br = 19.7 c^4 - a factor 3.8 below the prose.  The code implements App. C and
Gamma = g^2 beta/(8 pi m_h), i.e. the 4733 coefficient; the printed denominator "8 pi m_h^2" is
dimensionally inconsistent (it must be 8 pi m_h to give 4700).  Using x4 or 75 c^4 changes the
short-lived crossing by <= 0.9% and the long-lived boundary by 0.  Recommendation: state the
convention in the write-up; do not quote 75 c^4.

### F6 (NOTE) - K0L isospin coefficient taken as 1
The paper gives no SU(2) coefficient for K0L -> pi0 H (line 206: "can be obtained similarly").
Code (491-544) uses the same amplitude with neutral masses: K0L/K+ = 4.32 vs the Boiarska table
ratio 4.12 (= 4.128 if the charged-K matrix element is used).  The 4.7% difference is in a channel
that is closed above 0.354 GeV and has s_inel/s = 0 at every long-lived boundary.  Documented
choice; zero impact here.

### F7 (NOTE) - provenance/documentation drift (no physics impact)
The module is MODEL_TAG '2HDM-lightH-production-1.0.1' but __version__ = '1.0.0', and
logs/production_validation.md still lists a "frozen" sha256 4cc687ad... whereas the shipped file,
the outputs manifest (outputs/sha256_manifest.txt) and every scan row all record
**6381bb26...**.  The changelog attributes the change to a mass-grid path refactor only
(physics unchanged), but the validation log's frozen hash should be updated so the audit trail is
self-consistent.

### Ordering sanity (not a deviation)
At large tanb the observed local slope dlns/dlntanb ~ -8 (not -2) because s_B is the product of
BR(B->XsH) ~ 1/tanb^2 and the detection probability, which falls steeply with the growing ctau;
this is why a 4.6% rate change moves the boundary by only 0.56% at mH = 1, while at mH = 4 the
rate change is a factor 24 and the boundary moves by -52%.

---

## 6. Method validation and limitations

* Re-weighting linearity: s_B/BR(B->XsH) is constant at fixed (mH, tanb) to numerical precision;
  the response layer weight is exactly the production BR.  Reweighting reproduces the Lead's
  independent paper_mb/exclusive scans to 4.6e-5 relative (24 node/luminosity combinations).
* Crossing finder: own ln s vs ln tanb interpolation agrees with the Lead's refined crossings to
  2e-5 (mH = 1) and <= 3e-3 at the steep 3.7-4.2 GeV nodes.  All shifts are computed with the
  same finder on both sides, so the method bias cancels; the validation above quantifies it.
* Limitations: (i) the B0-fraction 0.5008 used for the literal-paper B0 variant is taken from the
  bbar_heavy sample composition (B+ 677192, B0 679352 rows), not from a per-parent response split;
  (ii) mH = 2.0 and 4.0 are interpolated between scan nodes (there is no 2.0 node by decay-program
  design and no 4.0 node); (iii) detector response, decay widths/BRs and detector acceptance are
  outside this audit; (iv) no scan was re-run.

## 7. Reproduction

    cd runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream
    python work/review_check/prodformula_audit.py     # writes production_formula_audit.json (0.9 s)

The script imports the production module read-only; the m_t sensitivity test monkey-patches
P.M_QUARK['t'] in memory and restores it, and the CKM variant uses P.ckm_matrix_from_magnitudes().
