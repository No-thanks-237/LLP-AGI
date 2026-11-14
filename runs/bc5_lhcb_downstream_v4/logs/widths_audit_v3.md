# Width-model audit v3 — scalar_widths.py against Winkler, arXiv:1809.01876 (scalar.tex)

Read-only, formula-by-formula audit. No production file was touched; all numbers below were
computed by re-implementing the paper's printed formulae independently (in a scratch script in
%TEMP%, deleted afterwards) and by driving work/widths/scalar_widths.py itself.  Line numbers refer
to work/widths/scalar_widths.py (411 lines) and to inputs/physics/winkler1809_scalar.tex (490 lines).
The paper's only numerical result for the pi-pi / K-K form factors is the figure fig:decayrates
(decayrates.pdf); the 'paper' value for those two channels is therefore the digitised figure
(work/widths/hadronic_table.csv), i.e. the same input our code uses.

---------------------------------------------------------------------------------------------------
## 1. Executive summary

1. **Our model is exactly Winkler's for m_S < 2 GeV** (ratio ours/paper = 1.0000 for e, mu, tau,
   pi-pi, K-K, 4pi at every requested mass): the analytic formulae are transcribed correctly and
   the hadronic input is the paper's own figure.
2. **It is NOT Winkler's for m_S >= 2 GeV.**  Two genuine bugs in the strong sector and one scheme
   choice: (i) the 3-loop alpha_s beta coefficients are wrong (2-loop term 2x too big, 3-loop 2x
   too small), giving alpha_s(2 GeV) = 0.3376 instead of 0.299 (+12.9%) and Gamma_gg 27.3% too big;
   (ii) the quark-mass anomalous dimension has the 1-loop coefficient 1 instead of 2, so mbar_q runs
   far too slowly (m_s(60) 22% high, m_b(m_Z) 3.35 instead of 2.82 GeV); (iii) the quark channels use
   [Boiarska] running masses instead of Winkler's fixed m_s = 95 MeV, m_c = 1.3 GeV.
3. **The 2 GeV junction is not smooth in our code but IS smooth with the paper's own formulae.**
   Ours: dispersive hadronic sum at 2.000 GeV = **7.7608e-8 GeV**, perturbative = **8.9367e-8 GeV**,
   ratio **1.1515 (+15.2%)**; total Gamma_S jumps 9.195e-8 -> 1.0378e-7 (**+12.9%**).  The single
   responsible term is **Gamma_gg, +27.3% (absolute +1.443e-8 GeV)**; Gamma_ss is -5.4% (-0.13e-8).
   With the paper's own inputs (alpha_s = 0.299, m_s = 95 MeV) the same sum gives 7.6195e-8 vs
   7.7608e-8, ratio **0.982 (-1.8%)**; Winkler's own figure is smooth to **0.003%**.  So the paper's
   'transits smoothly' claim is reproduced, and **C = 5.1e-9 does not need retuning**.
4. **Digitisation is not the weak link.** The figure's total-hadron curve equals the sum of its
   pi-pi + K-K + 4pi parts to 0.003%, the muon-calibration scatter is 0.08%, and axis fits are good to
   0.0016 px (4e-5 decade).  Read-off precision near 2 GeV is ~1% per channel -> +/-0.8e-9 GeV on the
   hadronic sum, 15x too small to explain the 1.18e-8 GeV jump.
5. **beta_2pi = sqrt(1-(4 m_pi)^2/m^2) is correct as coded** (the 4pi-threshold form): it reproduces
   Winkler's plotted 4pi curve at 0.5% from 0.7 to 2 GeV, while the 2-pion reading is off by 35% at
   0.7 GeV and 3.5% at 2 GeV.  Do not 'fix' it (see section 8).
6. **gamma-gamma**: the channel is absent from Winkler; ours is the standard SM-Higgs expression.  Its
   alpha_EM is taken at m_Z (1/127.95) instead of alpha(0) (1/137.036): at the SM anchor
   Gamma(h->gamma gamma, 125.25) = 10.99 keV vs the reference 9.39 keV (9.58 keV with alpha(0)).  The
   earlier '2.000x Boiarska' is real but is an error in [B]'s printed prefactor, not in our code.  At
   BC5 masses gamma-gamma is <= 4e-5 of Gamma_S, so this changes nothing numerically.
7. Ranked fixes and exact code changes in section 11; the list of things verified correct (do not
   reopen) in section 12.

---------------------------------------------------------------------------------------------------
## 2. Formula-by-formula table

Reference masses: 1.0 GeV (dispersive region) and 2.0 GeV (junction).  'paper' = independent
transcription of the printed formula.  Ours = value returned by scalar_widths.py.

| item | paper eq / line | our code line(s) | ours (GeV) | paper (GeV) | ratio | verdict |
|---|---|---|---|---|---|---|
| Gamma_ll, l=e | eq:mudecay, L125 | gamma_ll L219-224 | 1.7138e-13 | 1.7138e-13 | 1.0000 | correct (G_F m m_l^2 beta^3/(4 sqrt2 pi)) |
| Gamma_ll, l=mu | eq:mudecay, L125 | gamma_ll L219-224 | 6.8417e-9 | 6.8417e-9 | 1.0000 | correct |
| Gamma_ll, l=tau | eq:mudecay, L125 | gamma_ll L219-224 | 0 (m<2m_tau) | 0 | - | correct threshold |
| Gamma_gammagamma | **not in the paper** | gamma_gammagamma L227-236 | 9.520e-14 | [B] eq:Sto2gamma = 4.76e-14 | 2.000 | channel absent from [W]; amplitude correct, prefactor 2x [B] (see sec. 9) |
| Gamma_gg | L282-290 | gamma_gg L242-253 | (m=2) 6.7255e-8 | 5.2821e-8 | 1.2733 | formula identical; **alpha_s input wrong** (sec. 4) |
| Gamma_gg anchor | L282-290 | L242-253 | 1.7276e-4 GeV @125.25 | 1.7276e-4 | 1.000000 | prefactor and loop function exact; LO, no K-factor like the paper |
| Gamma_ss | eq:branching L275-279 | gamma_qq L275-298 | 2.2112e-8 | 2.3374e-8 | 0.9460 | beta^3, m_K threshold correct; running vs fixed mass |
| Gamma_cc | eq:branching + L280 | gamma_qq L275-298 | 0 @2.0 | 0 @2.0 | - | threshold 2m_D = 3.7393 GeV correct (opens at 3.74) |
| Gamma_bb | 'analogous', L280 | gamma_qq L275-298 | 0 @2.0 | 0 @2.0 | - | threshold 2m_B = 10.5593 GeV correct |
| Gamma_4pi | L292-295 (C=5.1e-9) | gamma_multimeson L301-305 | 4.2312e-9 | 4.2312e-9 | 1.0000 | exact; beta is the 4m_pi-threshold form (sec. 8) |
| Gamma_pipi | L155 | gamma_pipi_KK L358-364 (table) | 3.4068e-7 | 3.4068e-7 | 1.0000 | identical (digitised figure is the input) |
| Gamma_KK | L156 | gamma_pipi_KK L358-364 (table) | 1.4962e-7 | 1.4962e-7 | 1.0000 | identical |
| ChPT form factors | eq:formfactors L189-193 | chpt_form_factors L311-316 | (m_pi^2,0,s+2m_pi^2; m_pi^2/2,m_K^2-m_pi^2/2,s+2m_K^2) | same | 1.0000 | correct (fallback only, used <0.33 GeV) |
| Gamma_pipi prefactor | L155: 3G_F/(16 sqrt2 pi m) beta | L326 | | | | correct |
| Gamma_KK prefactor | L156: G_F/(4 sqrt2 pi m) beta | L336 | | | | correct |
| widths() regime switch | L291 (perturbative for m>2 GeV, dispersive 'up to 2 GeV') | widths L381 (ms < 2.0) | - | - | - | one-sided: at exactly 2.0 ours is perturbative, paper says dispersive (no numerical effect if the junction is smooth) |
| gamma_total | sum of partial widths | L394-395 | - | - | - | trivially correct once partials are |
| ctau_m | - | L398-400 | HBAR_C/(Gamma theta^2) | - | - | correct (Gamma scales as theta^2) |

### 2.1 Ratio ours/paper over the requested mass scan

'paper' uses its own prescription: dispersive below 2 GeV, spectator at/above 2 GeV.
Below 2 GeV every ratio is exactly 1.0000 (same formulae, same digitised table).  At/above 2 GeV:

| m_S | e | mu | tau | pi-pi | K-K | 4pi | s sbar | c cbar | b bbar | g g | total | ctau |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.3 - 1.99 | 1.0000 | 1.0000 | - | 1.0000 | 1.0000 (>=1.0) | 1.0000 (>=0.8) | - | - | - | - | 1.0000 | 1.0000 |
| 2.0 | 1.0000 | 1.0000 | - | - | - | - | 0.9460 | - | - | 1.2744 | 1.1460 | 0.8726 |
| 2.5 | 1.0000 | 1.0000 | - | - | - | - | 0.8748 | - | - | 1.2090 | 1.1123 | 0.8991 |
| 3.0 | 1.0000 | 1.0000 | - | - | - | - | 0.8280 | - | - | 1.1700 | 1.0919 | 0.9158 |
| 5.0 | 1.0000 | 1.0000 | 1.0000 | - | - | - | 0.7307 | 0.6914 | - | 1.1025 | 0.8335 | 1.1997 |
| 10.0 | 1.0000 | 1.0000 | 1.0000 | - | - | - | 0.6421 | 0.6076 | - | 1.0576 | 0.7730 | 1.2937 |
| 20.0 | 1.0000 | 1.0000 | 1.0000 | - | - | - | 0.5796 | 0.5484 | 0.7623 | 1.0309 | 0.7566 | 1.3217 |
| 60.0 | 1.0000 | 1.0000 | 1.0000 | - | - | - | 0.5097 | 0.4823 | 0.6704 | 1.0062 | 0.6753 | 1.4808 |

Every entry is a one-parameter effect, verified to 4 digits:
  Gamma_ss/paper = (mbar_s(m_S)/0.095)^2 ; Gamma_cc/paper = (mbar_c(m_S)/1.30)^2 ;
  Gamma_bb/paper = (mbar_b(m_S)/4.18)^2 ; Gamma_gg/paper = (alpha_s^code(m_S)/alpha_s^3loop(m_S))^2.
e, mu, tau, pi-pi, K-K, 4pi are identical to the paper by construction.

---------------------------------------------------------------------------------------------------
## 3. The 2 GeV junction — the decisive numbers

All values in GeV at theta^2 = 1.  'paper formulas' = our digitised pi-pi/K-K/4pi on the dispersive
side and the printed formulae (alpha_s = standard 3-loop, m_s = 95 MeV) on the perturbative side.
cc and bb are identically zero at 2 GeV (2m_D = 3.739, 2m_B = 10.559).

| term | ours at 1.99 (dispersive) | ours at 2.000 | paper at 2.000 | ours/paper (2 GeV) |
|---|---|---|---|---|
| pi-pi | 3.7276e-9 | - | - | 1.0000 (same input) |
| K-K | 3.5310e-8 | - | - | 1.0000 (same input) |
| 4pi/eta-eta/rho-rho | 3.8577e-8 | - | - | 1.0000 (same formula) |
| **dispersive sum** | **7.7615e-8** (at 1.99) | 7.7608e-8 (at 2.000) | 7.7608e-8 | 1.0000 |
| s sbar | - | 2.2112e-8 | 2.3374e-8 | 0.9460 |
| c cbar | - | 0 | 0 | - |
| b bbar | - | 0 | 0 | - |
| g g | - | 6.7255e-8 | 5.2821e-8 | **1.2733** |
| **perturbative sum** | - | **8.9367e-8** | **7.6195e-8** | **1.1729** |
| leptons (e+mu+tau+gamma) | 1.4335e-8 | 1.4409e-8 | 1.4409e-8 | 1.0000 |
| **Gamma_S total** | **9.195e-8** (1.99) / 9.2018e-8 (1.999) | **1.0378e-7** | 9.0604e-8 | 1.1454 |

**Is the sum [pi-pi + K-K + 4pi] equal to [gg + ss + cc + bb] at m_S = 2.000 GeV, C = 5.1e-9?**
  * Ours: **no.** 7.7608e-8 vs 8.9367e-8, ratio **1.1515**, i.e. a **+15.15% upward step** in the
    hadronic width and **+12.9%** in Gamma_S (the leptons dilute it).
  * Paper's own formulae: **yes.** 7.7608e-8 vs 7.6195e-8, ratio **0.9818 (-1.8%)**.
  * Winkler's own figure (raw digitised, calibration cancels in the ratio): total hadrons at 2 GeV =
    7.81605e-8; ss(2)+gg(2) = 2.29124e-8 + 5.52502e-8 = 7.81626e-8 -> ratio **1.00003**.

**Which term is responsible:** only Gamma_gg.  Per-term discrepancies ours vs the paper's formulae:
  * g g: **+1.4434e-8 GeV (+27.3%)** — the whole jump;
  * s sbar: **-1.2619e-9 GeV (-5.4%)** — a small downward offset;
  * pi-pi / K-K / 4pi: +0 / +0 / +0 (identical input);
  * c/b: 0 (closed).
  Check: +1.4434e-8 - 0.1262e-8 = +1.3172e-8 = 8.9367e-8 - 7.6195e-8.  Since the paper's own
  perturbative side sits 1.8% below its dispersive side, the visible jump is 1.176e-8 (15.15%).
  So the answer to 'is the kink a term error or a paper inconsistency?' is unambiguous:
  **the paper's transition is smooth; our Gamma_gg is 27% too big at 2 GeV because alpha_s(2 GeV)
  is 12.9% too big.  Gamma_gg scales as alpha_s^2.**

### 3.1 C retuning (only if the alpha_s bug is NOT fixed — not recommended)

| scenario | required C [GeV^-2] | ratio to 5.1e-9 |
|---|---|---|
| current code, unchanged (alpha_s=0.3376, running m_s) | 6.631e-9 | +30.0% |
| alpha_s fixed, current running masses | 4.746e-9 | -6.9% |
| alpha_s fixed, Winkler fixed masses (m_s=95 MeV) | 4.910e-9 | -3.7% |
| Winkler's figure (self-consistent) | 5.1e-9 | 1.000 |

Side effects of the +30% retune on Gamma_S (the option that would hide the alpha_s bug):

| m_S | 1.0 | 1.3 | 1.5 | 1.8 | 1.99 | 2.0 |
|---|---|---|---|---|---|---|
| Delta Gamma_S / Gamma_S | +0.25% | +2.43% | +4.18% | +8.99% | +12.59% | (+11.3% if the dispersive branch were used) |

Recommendation: **fix alpha_s, keep C = 5.1e-9.**  The residual -1.8% (paper's masses) to -3.5%
(running masses) is smaller than the 1% digitisation plus alpha_s-input uncertainty and matches
the paper to the accuracy its own figure supports.

---------------------------------------------------------------------------------------------------
## 4. alpha_s(m_S) — the junction breaker (Lead item A, confirmed)

scalar_widths.py L109-111:
    da(aa) = -(b0/(2*pi))*a^2 - (b1/(4*pi^2))*a^3 - (b2/(64*pi^3))*a^4
The RGE in d(alpha)/d(ln mu) with b0 = 11-2nf/3, b1 = 102-38nf/3, b2 = 2857/2-5033nf/18+325nf^2/54 is
    d(alpha)/d(ln mu) = -(b0/(2*pi))*a^2 - (b1/(8*pi^2))*a^3 - (b2/(32*pi^3))*a^4 ,
so the 2-loop term is **2x too large** and the 3-loop term **2x too small**.  This was verified
independently: my RK4 integration with the code's own coefficient set reproduces the code's grid to
1e-4, and with the standard set gives the values below (nf thresholds at 1.27 and 4.18,
alpha_s(m_Z) = 0.1179; the same with 0.1176 in brackets).

| mu [GeV] | 1.0 | 1.27 | 1.5 | 1.777 | **2.0** | 3.0 | 4.18 | 5.0 | 10.0 | 20.0 | 60.0 | 91.19 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| our code | 0.6533 | 0.4738 | 0.41098 | 0.36386 | **0.33762** | 0.27271 | 0.23694 | 0.22320 | 0.18296 | 0.15558 | 0.12626 | 0.11790 |
| corrected 3-loop | 0.4608 | 0.3816 | 0.34604 | 0.31659 | **0.29907** (0.29697) | 0.25213 | 0.22396 | 0.21257 | 0.17791 | 0.15323 | 0.12587 | 0.11790 |
| code/corrected | 1.418 | 1.242 | 1.188 | 1.149 | **1.129** | 1.082 | 1.058 | 1.050 | 1.028 | 1.015 | 1.003 | 1.000 |

PDG/Bethke-type reference values: alpha_s(m_tau)=0.331, alpha_s(2 GeV)=0.301, alpha_s(5)=0.215,
alpha_s(10)=0.178 — the corrected column agrees within 1-2%; the code's column does not below ~4 GeV.
Bethke (hep-ph/0604172, as cited at L291 of the paper) uses exactly this 3-loop nf-threshold running,
so **alpha_s^Bethke(2 GeV) = 0.299 (+/-0.005)** and ours is **+12.9%**.  Consequence:
**Gamma_gg(2 GeV) is 1.273x too large** (alpha_s^2), and the discrepancy falls to 1.006x at 60 GeV.

alpha_s() is used in three places: gamma_gg L252 (dominant), the mbar_q running through _build_grids
L140-142 (via the same grid), and gamma_qq's optional NLO correction L291-296 (disabled).

---------------------------------------------------------------------------------------------------
## 5. The quark-mass running (second, independent bug)

scalar_widths.py L140-142:
    g1 = 101/12 - 5 nf/18
    x  = a/pi
    gam = x * (1 + g1*x)
with dm/dln(mu) = -gam*m.  The 1-loop MS-bar anomalous dimension is gamma_m = **2** (alpha_s/pi), not
1 x (alpha_s/pi): the standard 1-loop solution m ~ alpha_s^(4/beta0) requires gamma_m = 2 alpha/pi.
Equivalently the correct 2-loop parametrisation is gam = 2x (1 + (g1/2) x) = 2x + 7.306x^2 (nf=4), so
the code has the right x^2 coefficient but only half the 1-loop term.  Evidence, comparing mbar_q at
the 2 GeV reference (gamma_m prefactor, alpha_s set):

| channel | m_S | current | gamma_m=2 only | alpha_s fixed only | both fixed | [Boiarska] c(x) formula |
|---|---|---|---|---|---|---|
| s | 5 | 0.08121 | 0.07504 | 0.08232 | 0.07653 | 0.07577 |
| s | 10 | 0.07612 | 0.06728 | 0.07741 | 0.06895 | 0.06809 |
| s | 60 | 0.06782 | 0.05500 | 0.06911 | 0.05654 | 0.05568 |
| c | 10 | 1.01333 | 0.89563 | 1.03045 | 0.91780 | 0.90638 |
| b | 60 | 3.42243 | 2.96439 | 3.44472 | 2.99415 | 2.97532 |
| b | m_Z | 3.35236 | 2.85688 | 3.37430 | 2.88573 | 2.86653 |

Boiarska Eq. (running mass) gives m_b(m_Z) = 2.867 GeV; the PDG value is 2.82(9) GeV.  Our current
code gives 3.352 (+19%); with both fixes 2.886 (+2%).  So the mass-running error is confirmed by two
independent references.

Effect of fixing the running (gamma_m=2 AND corrected alpha_s), on Gamma_S:

| m_S | 3 | 5 | 10 | 20 | 60 |
|---|---|---|---|---|---|
| Delta Gamma_qq/Gamma_qq | -5.2% (s) | -11.2% (s,c) | -18.0% (s,c) | -23.5% (s,c), -15.8% (b) | -30.5% (s,c), -23.5% (b) |
| Delta Gamma_S/Gamma_S | -0.73% | -5.14% | -8.27% | -14.72% | -21.70% |

At 2 GeV there is **no** effect (the reference scale of the running is 2 GeV by construction), and
none below 2 GeV (gamma_qq is not called there).  This bug does not touch the junction, but it is a
real physics error for m_S >~ 4 GeV.

---------------------------------------------------------------------------------------------------
## 6. s s-bar channel: running vs Winkler's fixed m_s = 95 MeV (task item 2c)

Winkler's eq:branching (L275-279) plus the footnote at L274 fix m_s = 95 MeV, m_c = 1.3 GeV and
neglect u,d.  Our gamma_qq instead uses mbar_q(m_S) with mbar_s(2 GeV) = 0.0924 (docstring L47-49,
[B] L1469).  The beta factors are identical (beta^3 with the K/D/B threshold), so the ratio is purely
(mbar_q(m_S)/m_q^Winkler)^2:

| m_S | 2.0 | 3.0 | 5.0 | 10.0 | 20.0 | 60.0 |
|---|---|---|---|---|---|---|
| mbar_s [GeV] | 0.09240 | 0.08644 | 0.08121 | 0.07612 | 0.07232 | 0.06782 |
| (mbar_s/0.095)^2 = Gamma_ss ratio | 0.9460 | 0.8280 | 0.7307 | 0.6421 | 0.5796 | 0.5097 |
| mbar_c [GeV] (vs 1.30) | 1.2300 | 1.1507 | 1.0810 | 1.0133 | 0.9627 | 0.9028 |
| Gamma_cc ratio | (closed) | (closed) | 0.6914 | 0.6076 | 0.5484 | 0.4823 |
| Gamma_bb ratio (vs 4.18) | - | - | - | - | 0.7623 | 0.6704 |

So at 2 GeV the ss channel is only 5.4% below Winkler's ratio form (not 22% — the 22% in the previous
audit was against the *fixed-mass digitised curve*, which is the same statement integrated over mass).
This is a deliberate scheme choice ([B] uses running masses); it is not an arithmetic error, but the
code cannot simultaneously claim to reproduce Winkler's printed ratios.  Quantified decision:
  * Winkler-literal (fixed masses): Gamma_SS(2) = 2.3374e-8, junction -1.8%;
  * Boiarska running masses: Gamma_SS(2) = 2.2112e-8, junction -3.5%.
The junction closes better with Winkler's fixed masses, but both are 'smooth' at the few-% level.

---------------------------------------------------------------------------------------------------
## 7. Thresholds and beta powers (task item 4)

* Paper L280: 'the kinematic threshold is set by the lightest meson containing an s or c quark';
  eq:branching uses beta_K^3 and beta_D^3 (beta^3, not beta).
* Our gamma_qq L286-288: 'if ms <= 2.0*m_hadron: return 0', 'beta**3' — correct power and structure.
* Thresholds agree exactly: s -> K, 2*0.493677 = **0.98735 GeV**; c -> D, 2*1.86966 = **3.73932 GeV**;
  b -> B, 2*5.27966 = **10.55932 GeV** (Winkler's ratio form only spells out s and c; the b channel is
  'analogous').  The digitised cc curve starts at 3.722 GeV and gg/ss at 2.0, consistent with the
  same boundaries.
* m_K = 0.493677 is the charged kaon (lightest meson with an s quark), m_D = 1.86966 (D+-), m_B =
  5.27966 (B+-): correct choices.
* The one convention difference: L367 PERT_THRESHOLD = 2.0 and L381 'if ms < PERT_THRESHOLD' put
  exactly m_S = 2.0 on the perturbative branch, while the paper says the spectator model is valid for
  m_phi > 2 GeV and the dispersive result is used 'up to 2 GeV'.  Immaterial if the junction is smooth;
  listed in the fix table as optional.

---------------------------------------------------------------------------------------------------
## 8. The 4pi/eta-eta/rho-rho term and beta_2pi (Lead item C) — the code is RIGHT

The paper writes Gamma_{4pi,...} = C s_theta^2 m^3 beta_{2pi} (L292-295) but beta_{2pi} is not defined
in the text; the phrase 'channels opening above the 4pi threshold' and the plotted curve decide it.
Our code L305 uses sqrt(1-(4 m_pi)^2/m^2), i.e. the 4m_pi threshold.  Comparison with Winkler's own
plotted 4pi curve (digitised, rescaled by the muon-calibration 1.0095):

| m_S | 0.70 | 0.80 | 1.00 | 1.30 | 1.50 | 1.80 | 1.99 | 2.00 |
|---|---|---|---|---|---|---|---|---|
| C m^3 sqrt(1-16m_pi^2/m^2) (ours) | 1.0553e-9 | 1.8703e-9 | 4.2312e-9 | 1.0119e-8 | 1.5976e-8 | 2.8276e-8 | 3.8577e-8 | 3.9178e-8 |
| C m^3 sqrt(1-4m_pi^2/m^2) | 1.6042e-9 | 2.4471e-9 | 4.8973e-9 | 1.0943e-8 | 1.6912e-8 | 2.9383e-8 | 3.9794e-8 | 4.0401e-8 |
| figure (calibrated) | 1.0502e-9 | 1.8615e-9 | 4.2114e-9 | 1.0073e-8 | 1.5903e-8 | 2.8149e-8 | 3.8398e-8 | 3.8996e-8 |
| figure / ours | **0.9952** | 0.9953 | 0.9953 | 0.9954 | 0.9954 | 0.9955 | 0.9954 | 0.9953 |
| figure / 2-pion form | 0.6546 | 0.7607 | 0.8599 | 0.9204 | 0.9404 | 0.9580 | 0.9649 | 0.9652 |

The 4-m_pi form reproduces the plotted curve to **0.5% with a constant residual from 0.7 to 2 GeV**
(the 0.5% is the figure's own normalisation); the 2-pion reading is off by up to 35%.  Likewise the
curve only enters the visible frame at m = 0.615 GeV, just above 4 m_pi = 0.558 GeV, and C = 5.1e-9
is the constant that makes the paper's junction smooth.  **Conclusion: our beta is the one Winkler
actually used and calibrated C against.  Changing it to the 2-pion form would be a regression**
(it would add +2.08% at 0.5 GeV, +0.13% at 1, +0.66% at 1.3, +0.82% at 1.5, +1.18% at 2, +0.68% at
3 GeV, and would require C -> 4.92e-9 to keep the junction).  The notation in the text is ambiguous,
but the figure and the calibration are not.

---------------------------------------------------------------------------------------------------
## 9. Photon and gluon channels at the SM anchor (Lead item B)

Winkler's scalar.tex contains **no gamma-gamma equation** (searched: no Gamma_gammagamma, no photon
decay entry in fig:decayrates; the only photon channel in the workspace's reference sources is
Boiarska eq:Sto2gamma, L1425-1429 of boiarska1904_main.tex).  Our gamma_gammagamma is the standard
SM-Higgs expression (amplitude F_{1/2} for e,mu,tau,u,d,s,c,b,t and F_1 for W).

Absolute anchor, m = 125.25 GeV:
  * our code: |A|^2 = 44.052, Gamma = **10.987 keV** with ALPHA_EM_MZ = 1/127.95 (L37, L236);
  * same amplitude with alpha(0) = 1/137.036: **9.578 keV**;
  * Boiarska eq:Sto2gamma: 5.494 keV (alpha_mZ) / 4.789 keV (alpha(0)) -> ours/[B] = **2.0000**;
  * SM reference Gamma(h->gamma gamma) = 9.39 keV (NLO); LO with alpha(0) ~ 9.6 keV.
So (i) our normalisation/amplitude is right (reproduces the LO SM value); (ii) the 'exactly 2x
Boiarska' is a genuine factor 2 in [B]'s printed prefactor (or an absorbed 1/2 in their F_gamma
definition), not an error of ours; (iii) the only real defect on our side is alpha_EM at m_Z instead
of alpha(0), a 13.2% inflation of gamma-gamma.  Effect on BC5: gamma-gamma is 6.8e-7 of Gamma_S at
2 GeV and at most 4.0e-5 at 3 GeV, so the correction changes Gamma_S by **< 0.001%** at every mass.
Fix it for correctness, not for physics.

Gluon anchor, m = 125.25 GeV: our Gamma(h->gg) = **0.1728 MeV** at alpha_s(m_h) = 0.11229; the paper's
formula with the same alpha_s and quark list gives the identical 0.1728 MeV (ratio 1.000000); the
heavy-top-limit LO value alpha_s^2 m^3/(72 pi^3 v^2) = 0.1831 MeV; the SM NLO (PDG) value is 0.337 MeV
(K ~ 1.8).  Our normalisation is the correct LO one and, like Winkler, we include no K-factor
(gamma_gg's k_factor defaults to 1.0, L242).  The 5.6% difference from the heavy-top limit is the
finite-top-mass + light-quark loop content, not a normalisation error.

---------------------------------------------------------------------------------------------------
## 10. Digitisation of the pi-pi / K-K table (task item 2d) — better than assumed

Built by build_hadronic_table.py from the vector PDF (no image tracing): 3000 log-spaced points,
0.2792-2.0 GeV (L36-37); pi-pi/K-K rescaled by the muon-curve calibration 1.0094568 (L49-50, reported
scatter 0.080%); below 0.33 GeV the pi-pi curve merges into the frame floor and ChPT is substituted
(L70-78); the K-K threshold is forced at 2*M_K_CHARGED (L78).

Quality measures (all computed from the digitised CSV itself):
  * y-axis calibration: 42.1934 px/decade, label-fit residual max 0.0016 px -> 3.8e-5 decade
    (0.004% in rate); x-axis 161.73 px/decade, residual 0.0010 px.  Curve extraction is from vector
    path centre lines, so stroke width does not bias the value.
  * **Internal consistency at 2 GeV: the plotted total-hadron curve 7.81605e-8 equals
    pi-pi + K-K + 4pi = 3.80186e-9 + 3.49915e-8 + 3.93647e-8 = 7.81581e-8 to 0.003%.**
  * Across the 2 GeV boundary the plotted curves are continuous: hadrons 7.81538e-8 -> 7.81709e-8
    (1.0002); the raw ss(2)+gg(2) = 7.81626e-8 equals raw hadrons(2-)=7.81605e-8 (1.00003).
  * Read-off precision near the junction is therefore ~1% per channel (0.95% calibration systematic
    plus 0.1-0.5% plotting/interpolation).  Implications: +/-0.4e-9 on pi-pi+K-K, i.e. **+/-0.5% of
    the 7.7608e-8 hadronic sum, and +/-0.05% of Gamma_S.  The observed 1.18e-8 GeV jump is 15x larger
    and cannot come from the digitisation.**
  * Residual artefacts worth knowing (no scan mass sits there): the K-K curve drops 1.4% between
    m = 1.9926 and 2.0000, and the table's last row at exactly 2.000 GeV interpolates to 3.4664e-8
    instead of the 3.499e-8 plotted at 2.0000.  Effect at the junction: KK(2.0) is 0.9% low and
    KK(1.99) 0.8% high, i.e. <= 0.7e-9 GeV in the sum.
  * The table's upper end IS consistent with the junction: the 4pi analytic formula matches the
    plotted 4pi curve to 0.5%, and pi-pi+K-K at 2 GeV equals the figure's hadron total to 0.003%.

Nothing in build_hadronic_table.py uses an m_pi0 value: M_PI = M_PI_CHARGED = 0.13957039 is used for
the 2-pion phase space (the isoscalar pi-pi channel; charged-pion mass is the natural convention and
matches the paper's beta_pi).  M_PI0 = 0.1349768 is defined (L57) but never used — dead constant.
M_KAON = 0.495 (L61) enters only the ChPT fallback (used below 0.33 GeV), where it differs from the
physical K mass by 0.5%; no effect on the production table.

---------------------------------------------------------------------------------------------------
## 11. Ranked fix list (NOT implemented — for the Lead to apply)

'Delta Gamma_S' is the change in Gamma_S(theta^2=1) if that fix alone is applied to the current code.

### Fix 1 (rank 1) — alpha_s 3-loop beta coefficients.  scalar_widths.py:109-111
    old:  return -(b0/(2*math.pi))*aa**2 - (b1/(4*math.pi**2))*aa**3 \
                 - (b2/(64*math.pi**3))*aa**4
    new:  return -(b0/(2*math.pi))*aa**2 - (b1/(8*math.pi**2))*aa**3 \
                 - (b2/(32*math.pi**3))*aa**4
Impact: alpha_s(2) 0.33762 -> 0.29907; Gamma_gg -21.5% at 2 GeV; **removes the 2 GeV kink**.

| m_S | 0.5 | 1.0 | 1.3 | 1.5 | 2.0 | 3.0 | 5.0 |
|---|---|---|---|---|---|---|---|
| Delta Gamma_S/Gamma_S | 0 | 0 | 0 | 0 | **-13.95%** | **-11.31%** | -0.48% |

(with the mbar_q changing too: additional +0.2% at 3 GeV, +1.3% at 5 GeV, +1.6% at 10 GeV, cancelling
part of the gg drop only at high mass.)

### Fix 2 (rank 2) — quark-mass anomalous dimension.  scalar_widths.py:142
    old:  gam = x * (1.0 + g1 * x)
    new:  gam = 2.0 * x * (1.0 + 0.5 * g1 * x)
Impact: mbar_s(10) 0.07612 -> 0.06895; m_b(m_Z) 3.35 -> 2.89 (PDG 2.82); Gamma_qq -18% at 10 GeV.

| m_S | 0.5 | 1.0 | 1.3 | 1.5 | 2.0 | 3.0 | 5.0 |
|---|---|---|---|---|---|---|---|
| Delta Gamma_S/Gamma_S | 0 | 0 | 0 | 0 | 0 | -0.73% | **-5.14%** |

### Fix 3 (rank 3, a decision not a bug) — Winkler-literal quark masses.  scalar_widths.py:387-389
    old:  out['ss'] = gamma_qq(ms, MSB_S_2GEV, 2.0, M_K)
          out['cc'] = gamma_qq(ms, MSB_C_2GEV, 2.0, M_D)
          out['bb'] = gamma_qq(ms, MSB_B_MB,  MSB_B_MB, M_B)
    new:  W_S, W_C, W_B = 0.095, 1.30, 4.18          # Winkler L274
          out['ss'] = gamma_qq(ms, W_S, ms, M_K)     # mu_ref = ms -> no running
          out['cc'] = gamma_qq(ms, W_C, ms, M_D)
          out['bb'] = gamma_qq(ms, W_B, ms, M_B)
(keep the mbar_q call as-is; passing mu_ref = ms makes mbar_q return the input mass).  Use only if a
literal reproduction of Winkler's Fig. 3/eq:branching is wanted; the current running-mass choice is
[Boiarska]'s and is defensible on its own.
Impact (after Fix 1): Gamma_ss +5.7% at 2 GeV, +20.8% at 3, +36.8% at 5, +55.7% at 10, +72.5% at 20,
+96.2% at 60; Gamma_cc/bb similar.  Gamma_S: +1.22% at 2, +2.90% at 3, +20.5% at 5, +29.7% at 10.

| m_S | 0.5 | 1.0 | 1.3 | 1.5 | 2.0 | 3.0 | 5.0 |
|---|---|---|---|---|---|---|---|
| Delta Gamma_S/Gamma_S | 0 | 0 | 0 | 0 | +1.22% | +2.90% | +20.46% |

### Fix 4 (rank 4) — alpha_EM for the photon channel.  scalar_widths.py:37 (and L236)
    old:  ALPHA_EM_MZ = 1.0 / 127.95
    new:  ALPHA_EM_0  = 1.0 / 137.036   # on-shell photon coupling (use in gamma_gammagamma L236)
Impact: Gamma_gammagamma -13.2%; at the 125.25 GeV anchor 10.99 -> 9.58 keV (reference 9.39).

| m_S | 0.5 | 1.0 | 1.3 | 1.5 | 2.0 | 3.0 | 5.0 |
|---|---|---|---|---|---|---|---|
| Delta Gamma_S/Gamma_S | -0.000% | -0.000% | -0.000% | -0.000% | -0.000% | -0.001% | -0.000% |

### Fix 5 (optional) — put exactly 2.0 on the dispersive branch.  scalar_widths.py:381
    old:  if ms < PERT_THRESHOLD:
    new:  if ms <= PERT_THRESHOLD:
Paper's wording: spectator valid for m_phi > 2 GeV, dispersive used 'up to 2 GeV'.  Only matters at
exactly m_S = 2.0 (where the two branches differ by 3.5% after Fix 1).  No effect elsewhere.

### Explicitly NOT a fix
* beta_2pi (L305) — verified against the paper's own plotted curve (section 8).
* The 6-quark sum in gamma_gg (L239) — the complex branch of the paper's f(x) exists precisely so that
  quarks lighter than m_phi/2 are summed; u,d are negligible (<0.1%) but s contributes -4.6% to
  Gamma_gg at 2 GeV and <0.5% above 4 GeV.  Using only c,b,t would inflate the junction by +3.1%.
* The LO treatment (QUARK_NLO = False, L272) — the paper has no QCD K-factor for quarks or gluons, so
  LO is the correct comparison.  Boiarska's NLO factor is a different model (and at 2 GeV it would
  double the rate, breaking the junction by construction; see L257-271).
* ChPT form factors, thresholds, beta powers, lepton formula — correct (section 12).

---------------------------------------------------------------------------------------------------
## 12. Verified correct — do not reopen

1. **Gamma_ll (eq:mudecay, L125)**: identical formula, identical beta_l^3, identical G_F.
2. **Gamma_pipi / Gamma_KK prefactors and weights (L155-156)**: 3G_F/(16 sqrt2 pi m) beta_pi and
   G_F/(4 sqrt2 pi m) beta_K with 7/9, 7/9, 2/9 weights; verified at 1.0 GeV to 1e-7.
3. **ChPT form factors (eq:formfactors, L189-193)**: exactly the 6 entries and the s+2m^2 structure;
   the pi-pi/K-K ChPT fallback matches the printed equations (the fallback is only used below
   0.33 GeV, where it is reliable).
4. **The 7/9, 7/9, 2/9 weights and the T-matrix structure**: our code never re-derives the dispersive
   result; it uses Winkler's published figure, which is the only available form of that result.  No
   error can hide here beyond the digitisation, which is good to ~1%.
5. **Gamma_gg formula (L282-290)**: with the same alpha_s and the same quark list, the paper's
   |sum [x+(x-1)f]/x^2|^2/(32 pi^3) form and our |sum F_{1/2}|^2/(128 pi^3) form agree to 1.000000
   at every mass tested (125.25 and 2 GeV).  F_half(inf) -> -4/3 reproduces the heavy-top limit
   alpha_s^2 m^3/(72 pi^3 v^2) = 0.1831 MeV at 125 GeV (ours 0.1728 incl. finite top + light quarks).
   The complex branch for m_q < m_phi/2 is implemented correctly (log form, stable zeta).
6. **Gamma_4pi (L292-295)**: C = 5.1e-9, m^3, 4m_pi-threshold beta — matches Winkler's plotted curve
   to 0.5% over 0.7-2 GeV and is the constant that makes his junction smooth (section 8).
7. **Thresholds and beta powers for quarks**: beta^3, K/D/B meson masses, thresholds 0.98735 /
   3.73932 / 10.55932 GeV — exactly the paper's prescription (section 7).
8. **widths() structure, gamma_total, ctau_m**: branch switch, channel bookkeeping, theta^2 scaling and
   ctau = HBAR_C/(Gamma theta^2) are correct; the sums reproduce the partial widths to machine
   precision; the branching ratios sum to 1 by construction.
9. **The digitised hadronic table and its builder**: calibrated, internally consistent to 0.003% at
   2 GeV; no NaN/negative entries; thresholds handled; no m_pi0 contamination; M_PI0 unused.
10. **The 2 GeV matching concept**: keeping the dispersive result to 2 GeV and adding C m^3 beta_2pi
    is exactly Winkler's prescription; our implementation of the *concept* is right — only the alpha_s
    input breaks the smoothness.
11. **The 2x Boiarska gamma-gamma prefactor**: our gamma-gamma is the standard SM-Higgs normalisation
    (9.58 keV at 125.25 GeV with alpha(0), vs the 9.39 keV reference); the factor 2 is in [B]'s printed
    equation, not in our code (section 9).

---------------------------------------------------------------------------------------------------
## 13. Items in the code that are not in the paper, and paper items not implemented

Not in Winkler's scalar.tex (intentional additions, mostly from Boiarska):
  * the gamma-gamma channel entirely (no equation, no curve in [W]);
  * MS-bar running quark masses (eq:decay-quarks of [B]) instead of fixed m_s, m_c;
  * optional NLO QCD corrections (QUARK_NLO L272, Delta_QCD/Delta_t L291-297) — Winkler has none;
  * the k_factor argument of gamma_gg (L242, default 1.0 = LO, consistent with the paper);
  * the PERT_THRESHOLD switch (L367, L381), channel dictionary names, the digitised-table machinery
    (HADRONIC_TABLE L345, load_hadronic_table L348-355), the ChPT fallback and its series-stable
    F_half/F_one implementations (L177-213), HBAR_C/ctau_m, ALPHA_EM_MZ.
  * 'alpha_s(m_phi) from Bethke' is approximated by our own 3-loop run instead of Bethke's fitted
    parameterisation — acceptable in principle, but it is the locus of Fix 1.

In the paper but not implemented (by design, because the input data are not in the workspace):
  * the two-channel Muskhelishvili-Omnes construction itself (eq:unitary L200-217, T-matrix L206-216,
    eq:omnessolution L221-229, the Omega_ij from Moussallam, p = 0.73, q = 0.52 at L260-263).  We
    replace the *result* by the digitised figure — documented in build_hadronic_table.py and in
    work/widths/README.md.  Any residual uncertainty of that replacement is the ~1% digitisation error,
    not a formula error.

---------------------------------------------------------------------------------------------------
## 14. Caveats

* 'Paper' values for pi-pi/K-K are the digitised figure by necessity; the dispersive calculation
  itself was not re-derived (the Hoferichter/Moussallam phase-shift and Omega_ij inputs are not present
  in the workspace and reconstructing them is out of scope of a read-only audit).
* The absolute alpha_s reference is the standard 3-loop nf-threshold running with alpha_s(m_Z)=0.1179;
  Bethke's fitted parameterisation may differ by O(1%) at 2 GeV, i.e. O(2%) in Gamma_gg — far smaller
  than the 12.9% effect found.
* All Gamma_S numbers are theta^2 = 1; branching ratios and exclusions scale trivially.
* No production script, table or output was modified; the only file written is this report.
