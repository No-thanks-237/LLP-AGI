# T5 — independent analytic audit of the 2HDM light-H implementation

Owner: model (T5). Write scope: work/model/audit/**, logs/model*.  No FeynRules/UFO/MG5 (user
decision 2026-10-25; the earlier work/model/{fr,ufo_2HDM_typeI,wolfram_apps,mg5} were deleted and
nothing in this report depends on them).

## 0. Provenance and files

| item | value |
|---|---|
| Host Wolfram | 15.0.0 for Microsoft Windows (May 15 2026), Professional licence, "C:\Program Files\Wolfram Research\WolframScript\wolframscript.exe" |
| Python | D:\py\python.exe, 3.14.3 |
| Paper text | inputs/references/FASER_light.tex (arXiv:2212.06186v1) |
| Authors' FORESEE model | references/tools/FORESEE/Models/2HDM_H/{2HDM_H.ipynb, model/br/*.txt, model/ctau.txt} (FORESEE git checkout in temp/) |
| Our decay module | work/decay/build/light_scalar_decay.exe (author Light_scalar_decay, batch driver) |
| Production module | work/production/twohdm_production.py |
| Response library | work/response/twohdm_lib.py (read-only comparison) |
| Boiarska master formula | inputs/bc5_reuse/bc5_production.py (read-only import) |

Audit code (all new, written from the paper + the FORESEE notebook, not from the production module):

* work/model/audit/audit_2hdm.py — independent implementation of xi_H^{ij} (App. D), B->XsH,
  K->pi H coherent amplitude, pi/K semileptonic, g_C/g_hHH^(3)/BR(h->HH), 9-process kappa scaling.
* work/model/audit/audit_checks.py — term-by-term xi_H break-down, limits, CKM-phase effect,
  B/K/pi rates, g_C/kappa tables.
* work/model/audit/audit_foresee_tables.py -> foresee_vs_exe.json — FORESEE br/ctau vs our exe.
* work/model/audit/audit_foresee_B.py -> foresee_B_compare.json — FORESEE B->XsH formula vs ours.
* work/model/audit/audit_compare_prod.py — production module vs this implementation.
* work/model/audit/audit_width_conventions.py — SM ground-truth test of the gamma-gamma / gg
  normalisations (resolves the "1/4" and "9/4" flags of the decay owner).
* work/model/audit/audit_or_toy.py — lifetime-integral / event-OR / 0.4 / 1e3 toy audit.
* work/model/audit/wl/wl_checks.wls -> wl/wl_checks.out — Wolfram symbolic/numeric checks
  (App. D limits, tanb scaling, g_C expansion, 4700 (g/v)^2, App. A form-factor conventions).

Commands: `python audit_foresee_tables.py`, `python audit_foresee_B.py` ... , and
`wolframscript.exe -file wl/wl_checks.wls` (log wl/wl_checks.out, stdout captured).

## 1. FORESEE decay tables vs our author-code executable  (PASS, 5-digit)

* FORESEE tables: 3 columns (m_H, tanbeta, BR), 61 log masses 0.01-10 GeV x 201 tanb (1-1e5)
  = 12261 points for pi_pi, k_k, 4pi, c_c, s_s, tau_tau, g_g, mu_mu0; the ga_ga, e_e, mu_mu,
  e_e0, ga_ga0, mu_mu0 tables use a 31 x 101 sub-grid (tanb up to 1e10); ctau.txt is 31 x 101.
* Our executable was run on all 15197 unique points in one batch (stdin "<mH> <tanb>", 13 output
  columns per DESIGN section 2).

| channel | N pts | median abs(ln ratio) | max abs(ln ratio) |
|---|---|---|---|
| pi_pi, k_k, 4pi, c_c, s_s, tau_tau, g_g | 12261 | 3.2e-6 | 4.7e-5 |
| ga_ga | 2604 | 2.6e-11 | 1.6e-2 (*) |
| e_e | 1550 | 8.1e-6 | 1.1e-4 |
| mu_mu | 1550 | 8.5e-6 | 1.2e-4 |
| ctau (hbar*c = 2.0e-16 as in FORESEE) | 2604 | 1.1e-5 | 4.7e-5 (excl. *) |

`*` The single structured outlier is (m_H = 0.01, tanbeta = 1e10): our Gamma_tot = 4.93e-38 GeV,
so the eight non-gamma-gamma channels hit the program's hard-coded 1e-40 floor and each
contributes BR = 1.996e-3, giving BR(gamma gamma) = 0.98403 vs the table's 1.0000. This is a
numerical-floor artefact at tanb = 1e10, far outside the analysis range (tanb <= 1e5); it is
also the only point where ctau deviates (1.6e-2).

The "threshold mismatches" (e_e 30, mu_mu 570 points) are the FORESEE tables storing exactly 0
below threshold while the code returns its 1e-40 floor (BR ~ 1e-13); no physics difference.
Without the hbar*c correction the ctau ratio is systematically ln(2.0/1.973269804) = 1.346e-2
(the 1.01355 factor quoted by the lead), i.e. our 1.973269804e-16 vs FORESEE's 2.0e-16.

VERDICT: our compiled author program reproduces the authors' own FORESEE BR/ctau tables to their
printing precision on 12261 x 10 grid points. The decay module is validated against the paper's
own generated numbers, independently of the paper figures.

## 2. FORESEE B -> X_s H formula vs ours and vs the production module

FORESEE (2HDM_H.ipynb cell 7; pid0 511/521 -> pid1 321, same expression for B+/B0):

    br = 9.68544e-13 * mB / Gamma_B * (1-(m/mB)^2)^2
         * ( -3/4/tb - 0.1415/tb^3 - (1.0386/tb - 0.1415/tb^3)*sqrt(1-1/tb^2)
             - 0.608446664/tb^3 )^2 ,      tb = tanbeta, mB = 5.27934, Gamma_B+ = 4.0184e-13 GeV

(9.68544e-13*mB/Gamma_B = 12.7247.)

### 2a. The flavour function
Our exact App.-D bracket (t quark; u,c add <1e-4) divided by the FORESEE fitted bracket =
0.99294 at every tanb from 5 to 1e5 and every m_H from 0.021 to 4 GeV (m_H-independent to 1e-7:
the g0*lambda_HH+H- term is only 3.4e-5 of the bracket, so the FORESEE formula's lack of m_H
dependence is justified). Limits: bracket*tanb -> -1.801316 (ours) vs -1.788600 (FORESEE) as
tanb -> infinity; both have the 1/tanb + 1/tanb^3 structure; m_H -> 0 is identical (the bracket
does not vanish, only xi -> 0 as 1/tanb).

### 2b. Absolute rates (br_FORESEE / br_ours, paper m_b convention)
m_H -> 0: 0.9194 = (12.7247/13.646)*0.99294^2 (fitted-vs-exact constant 0.9325, bracket 0.98593);
m_H = 1: 0.9624 (x1.0458 phase space); 2: 1.188; 3: 1.877; 4: 24.2.
The growth is entirely the phase-space convention (1-m^2/m_B^2)^2 vs (1-m^2/m_b^2)^2 =
1.00002, 1.0458, 1.2338, 1.9498, 25.55 at m_H = 0.02, 1, 2, 3, 4 GeV. The authors' B channel
therefore closes at m_B = 5.279 GeV, the DESIGN/paper formula at m_b = 4.18 GeV.

VERDICT: the physics content agrees (tanb dependence 0.7%, normalisation 8% from the fitted
coefficients); the two formulas differ by an m_b-vs-m_B phase-space convention that grows
rapidly above m_H ~ 2 GeV.

### 2c. IMPORTANT — the production module's default is the m_B convention
work/production/twohdm_production.py has B_MODES = ('foresee', 'paper_mb', 'exclusive') and
DEFAULT_B_MODE = 'foresee'.  Comparison with this independent implementation:

| mode | max abs(ln ratio) vs audit, over 8 m_H x 6 tanb |
|---|---|
| 'paper_mb' (DESIGN section 3.2, m_b) | 4.4e-16 (machine precision) |
| 'foresee' (default, m_B) | 3.24 (i.e. up to x25.5) |

So: **the two implementations are identical; the only difference is the phase-space mass
choice.** In the default 'foresee' mode Br(B->XsH) is larger than the DESIGN formula by
+4.6% (m_H=1), +23% (2), +95% (3), x25.5 (4), and the B channel is open up to 5.279 GeV.
This matters for the m_H >~ 3 GeV reach (the B channel is the dominant production there) and is
a Lead decision, not a bug: the authors' own FORESEE model uses m_B, DESIGN section 3.2 writes m_b.

## 3. Independent re-implementation: xi_H, K->pi H, semileptonic (production identical)

### 3a. xi_H^{bs} term-by-term at tanb = 100, m_H = 1 (DESIGN 3a)
prefactor -4*sqrt(2)*G_F/(16 pi^2) = -4.178254e-07; c = 0.01, s = 0.99995,
lambda_HH+H- = -29.26825 GeV, x_H+ = 55.71248.

| k | x_k | g0 | g1 | g2 | bracket | term |
|---|---|---|---|---|---|---|
| u | 7.49e-10 | -3.365494e-07 | -0.75 | +0.0351515 | -4.265052e-02 | +3.06e-17 -6.73e-17i |
| c | 2.535537e-04 | -3.365268e-07 | -0.75 | +0.0350967 | -4.259567e-02 | +1.158403e-09 |
| t | 4.637075 | -2.738287e-07 | -0.75001398 | +0.0105130 | -1.801324e-02 | -9.039007e-06 +1.756e-07i |

xi_H^{bs} = -9.037848e-06 + 1.756021e-07i, |xi| = 9.039554e-06.  Wolfram re-derivation gives
9.041784e-06 using |V| real (2.5e-04 higher): the entire difference is the CKM phase.
Independently of that, the production module returns the *same* complex number (bit-identical).

Limits: |xi_bs|*tanb = 9.039554e-04 (tanb=100), 9.039512e-04 (1e3), 9.039512e-04 (1e5)
=> xi ∝ 1/tanb to 4.6e-06 at tanb=1e4, as the paper requires (Wolfram confirms: g0,g2 -> 0,
g1 -> -3/4 at tanb -> infinity).  m_H dependence: |xi_bs| changes by <1e-8 between m_H=0.021 and
5 GeV (only the g0 term, 3.4e-5 relative).

### 3b. B -> X_s H (paper formula, m_b = 4.18)
Br(B+ -> XsH): 4.1376e-03 (m_H->0), 3.6778e-03 (1), 2.4602e-03 (2), 9.7294e-04 (3),
2.9385e-05 (4), 0 for m_H >= 4.18; Br(B0)/Br(B+) = 0.930.  Br > 0.1 for tanb <~ 19 at m_H=1
(tanb=20: 9.20e-02; 30: 4.09e-02; 100: 3.68e-03) -> model-validity flag for the scan.
Production module ('paper_mb') = these numbers to machine precision.

### 3c. K -> pi H coherent amplitude
|A_4q|/|A_FCNC| = 0.0345 and tanb-independent (both amplitudes ∝ 1/tanb); the phase of the
complex xi_ds makes the interference reduce the rate by ~7% w.r.t. |A_FCNC| alone.
Br(K+->pi+H) = 3.6588e-04 (tanb=5), 9.1281e-07 (100), 9.1281e-09 (1000) at m_H->0;
Br(K0L->pi0H) = 4.30x the K+ value (3.941e-06 at tanb=100) -> x0.5 flavour projection = 1.97e-6.
Cross-check against Boiarska's master formula (bc5_production.br_K_to_S rescaled by
(xi_ds^2HDM/3.3e-6)^2): ratio 0.9972 (0.28% = different m_s, m_d inputs: 0.096/0.0047 here vs
0.0924/0.00467 in the BC5 module).  FCNC-only value from the paper matches Boiarska exactly in
form (same 1/(32 pi v^2 Gamma) * p/m^2 master formula).
Production module = these numbers exactly (K+ and K0L, ratio 1.000000).

### 3d. Semileptonic
Br(pi->H e nu) = 1.96395e-13 (tanb=100, m_H->0), Br(K->H e nu) = 3.91173e-12 (n_h=n_l=3 ->
(7/9)^2 included).  Production module identical (ratio 1.000000). K0L semileptonic = 0 by DESIGN.

## 4. Probability algebra: lifetime integral, event OR, 0.4, 1e3  (PASS)

(a) A=1 window probability: p = exp(-alpha(z1-z0)) - exp(-alpha(z2-z0)), z1,z2 = 1.0,2.5 m,
    alpha = m_H/(c tau p_z) = 1/Lambda_z.  36 combinations alpha in [0.05,50] m^-1, z0 in
    [0,2.5] m: max |work/response/twohdm_lib._p_window - analytic| = 1.1e-16 (machine precision).
    Brute-force decay-length MC (4e5 decays/config) agrees within binomial noise (max 1.8e-3 =
    2.2 sigma).  Also verified p -> 0 for c tau -> 0 (prompt) and p -> (1.5 m)/c tau for c tau >> 1 m.
(a2) For a *peaked* A(z) (Gaussian sigma = 0.35 m) the 13-node piecewise-linear interpolation
    used by the response differs from the continuous quadrature by <= 1.8e-3 absolute
    (<=0.4% relative at p ~ 0.5).  This is the dominant approximation in the response layer and
    is acceptable for the required 5% MC-level precision.
(b) Event OR 1-prod(1-p_i): 4e5 events, 1 candidate / 1-3 candidates uniform p / 1-6 candidates
    log-uniform p in [1e-40, 0.1]: agreement 0.5, 0.6, 1.2 sigma.  h->HH sub-case
    o = BR_hHH*(1-(1-p1)(1-p2)) with a Bernoulli gate: 0.4 sigma.  Mutually exclusive history
    weights eps = sum_c w_c eps_c: direct per-event MC 0.015067 vs 0.014860 (+1.4%, MC noise).
(c) 0.4: twohdm_lib.RECO_EFF = 0.4, applied once; the estimator s = 1e3 L sum_c sigma_c eps_c
    gives, e.g., sigma_bb = 5.6e8 pb, L = 25 fb^-1, eps = 0.4 -> s = 5.6e12 (one factor 0.4;
    applying it twice would give 2.24e12).  1 fb^-1 = 1e3 pb^-1 verified numerically for
    L = 25 and 300 fb^-1 and sigma in {5.6e8, 7.95e7, 55.45} pb.

## 5. Wolfram symbolic checks (wl/wl_checks.wls, host Wolfram 15)

* App. D, x -> x_H: g0 -> -cot^2(beta)/(12 x_H), g1 -> -(9 + 5 cot^2(beta))/12,
  g2 -> cot(beta)/4 * [-1/2 - (5/3)cot^2 + 6/(x_H-1)^2 - 3(1+x_H)ln x_H/(x_H-1)^3].
* App. D, x -> 0: g0 -> -3cot^2/(16 x_H), g1 -> -3/4, g2 -> cot(beta)(2 + 3 ln x_H)/4;
  x_H -> infinity: g0 -> -3cot^2 x/(16 x_H) (O(x/x_H)), g1 -> -3/4 + O(1/x_H).
* tanb -> infinity at x_t: (g0, g1, g2) -> (0, -3/4, 0) exactly -> all xi_H^{ij} -> 0 like 1/tanb.
* g_C (App. C, lambda v^2 = 0), exact small-c expansion:
      g_C = 15625 c^2/984 + 15625 c^4/3936 - m_H^2/492 = (s/2v)[ mh^2 c^2/2 - m_H^2(1+c^2/2) ] + O(c^4)
  so (i) the coefficient of the c^2 term is numerically identical to mh^2/(4v) = 15.88, but with
  the *opposite sign* to the paper's quoted "-mh^2 c^2/(4v)", and (ii) the -m_H^2/492 term must
  be kept: it dominates unless m_H^2 << mh^2 c^2/2 (m_H < 0.88 GeV at tanb=100).  At
  (m_H=1, tanb=100): exact g_C = -4.4457e-04 GeV vs -mh^2 c^2/(4v) = -1.5879e-03 GeV (ratio
  0.28).  Gamma(h->HH) ∝ g_C^2, so the sign is immaterial; the magnitude must use the exact form.
* Br(h->HH) vs 4700 (g/v)^2: (m_H=0.021, tanb=100) exact BR = 1.9699e-07 vs 4700(g_C/v)^2 =
  1.9562e-07 (ratio 1.0070, the 0.7% is the beta factor); (m_H=1, tanb=100) BR = 1.5456e-08.
* App. A: A_0, A_1/2, A_1 of the paper = exactly 1/2 of the Djouadi forms used in the standard
  SM formulae (verified at tau = 0.01, 0.1, 0.5; limits A_1/2: 2/3 vs 4/3, A_1: -3.5 vs -7),
  confirming the paper footnote.

### 5b. Resolution of the decay owner's "1/4" and "9/4" flags (SM ground truth, m_phi = 125 GeV)
* gamma-gamma: with all xi = 1 and no H+- loop, the author-code formula (Gamma_gaga in
  CPevenDecay.cpp = paper Eq. hgaga with the App.-A form factors, prefactor alpha^2 m^3/(64 pi^3 v^2))
  gives 9.3139 keV vs the PDG SM value 9.28 keV (ratio 1.0037).  The apparent 1/4 is the (1/2)^2
  from A^paper = (1/2) A^Djouadi, exactly cancelled by 1/(64 pi^3 v^2) = 4/(128 sqrt2 pi^3):
  Gamma(code) = Gamma(standard).  **The 1/4 is a convention artefact, NOT an error.**
* gg: code form alpha_s^2 m^3/(32 pi^3 v^2) |sum xi_q A_12^paper|^2 = 0.1748 MeV and the paper form
  G_F alpha_s^2 m^3/(36 sqrt2 pi^3) |(3/2) sum|^2 = 0.1745 MeV — identical, because
  (9/4)*(32/72) = 1.  This is the origin of the "9/4".  However the *physically correct* LO
  normalisation requires xi_g = 2 sum A_12^paper (not 3/2), giving 0.3101 MeV = 0.92 x PDG
  0.3358 MeV (residual = NLO QCD).  So the author code's Gamma_gg is 9/16 = 0.5625 of the
  standard LO value: a real normalisation defect shared by the paper and its own figure (the
  figure cannot reveal it).  The same 3/2 convention enters the xi_g-dependent hadronic
  channels (pi pi, K K, 4 pi).
  Impact on this analysis: all of gg/pi pi/K K/4 pi scale with xi_H^f ∝ 1/tanb^3, so
  BR(gg) + sum BR(hadronic) < 1e-6 at tanb >= 30, where the LHCb-downstream reach lies; at low
  tanb the H decays promptly (c tau ~ mm) and contributes nothing.  Gamma_tot is gamma-gamma
  dominated at large tanb and the gamma-gamma normalisation is verified correct to 0.4%.
  => do NOT patch the author's physics (DESIGN section 1 forbids physics edits); document only.

## 6. Verdict summary

1. Decay module: PASS — reproduces the authors' FORESEE BR/ctau tables to 5-digit printing
   precision (max |ln ratio| 4.7e-5) on 12261 points x 10 channels; only artefact is the 1e-40
   channel floor at tanb = 1e10.
2. Analytic production layer: PASS — the production module is bit-identical to this independent
   implementation for xi_H (complex), B->XsH, K->pi H (K+ and K0L) and pi/K semileptonic.
3. OPEN DECISION for the Lead: production DEFAULT_B_MODE = 'foresee' uses m_B in the B->XsH
   phase space; DESIGN section 3.2 writes m_b.  Effect: +4.6% (m_H=1) ... +95% (3) and the B
   channel closes at 5.279 instead of 4.18 GeV.  Pick one explicitly and record it.
4. gamma-gamma normalisation: correct (0.4% vs SM).  gg/hadronic 9/16 normalisation defect:
   real, documented, no effect on the downstream curve.
5. K->pi H: the four-quark amplitude is 3.45% of the FCNC one (coherent, destructive ~7%),
   tanb-independent; consistent with Boiarska's master formula to 0.3%.
6. Probability algebra: exact (p_window to 1e-16 vs analytic; OR identity verified by MC); the
   13-node A(z) grid interpolation is the dominant response approximation (<=1.8e-3 absolute in p).
7. g_C: the "-mh^2 c^2/(4v)" approximation has the right magnitude but the wrong quoted sign and
   omits the -m_H^2/(2v) term; the exact App.-C formula used by DESIGN is correct.
