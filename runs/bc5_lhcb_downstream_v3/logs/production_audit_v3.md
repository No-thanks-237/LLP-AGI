# Production audit v3 — every S-production formula and number vs arXiv:1904.10447v4 and arXiv:1809.01876

Scope: read-only audit of `work/analysis/bc5_production.py`, `work/analysis/run_scan.py`,
`work/analysis/bc5_lib.py`, `work/widths/scalar_widths.py` (only the parts entering production), the
parent samples and the cross-section inputs. **No file was modified and no scan was run.**
Every number below was produced by importing the modules read-only with `D:\py\python.exe`
(temp scripts in `%TEMP%`, deleted afterwards) or read from the .tex.

Sources used, with labels as printed in the paper: `eq:productionBRmesons` (L731), `tab:BR` (L263–285),
`tab:xi` (L659–670), `eq:BR-h-to-SS` (L355), `tab:b-production-quartic` (L394–420), `eq:bstossbranching`
(L1397), `eq:br-quartic-3-body` (L1402), `tab:meson-decay-constants` (L1383), Apps. C.1–C.5 (L1111–1346),
L682 (B_s excluded), L740 (K_S^0 → 0).

---

## 0. Executive summary

1. **The production formulas are essentially correct.** `eq:productionBRmesons` is implemented factor for
   factor; all 10 B→X_s S channels reproduce `tab:BR` at m_S→0 to −4 %…+6 % (sum 3.356 vs 3.3); the
   `tab:b-production-quartic` BR(0) values are reproduced from `eq:bstossbranching`/`eq:br-quartic-3-body`
   to 0.97–1.16, and the paper's own closing masses to <1 % for 8/10 three-body channels.
2. Every form factor in the code is **bit-identical** to an independent transcription of Apps. C.1–C.5
   (all 10 channels, all parameters, all q²-poles, checked at 5 q² values; max |Δ| = 0).
3. **α normalisation is right**: α(1 GeV) = 0.719525, dimension GeV, and it inverts the same convention as
   `eq:BR-h-to-SS` (α = 1 GeV ⇒ BR(h→SS) = 0.01951 vs the paper's "≈2.0·10⁻²").
4. The two non-paper heuristics are the only shape defects: the 3-body quartic power-fit is within
   +3…+7 % of the exact q²-dependent shape (worse only for K₂*(1430), K₁(1400), K₀*(1430)), and the kaon
   quartic exponential is 1.5–3.7× too steep — but the K channel is ≤0.08 % of s, so it does not matter.
5. The **real** problems are outside the formulas, in the parent-sample normalisation: the b-b̄ sample
   carries **1.320 B-hadrons per generated event against 1.818 expected** (the excited b-hadrons are kept
   stable by `tau0Max=0` and are not in the `BONLY` list) ⇒ the B channel is ~1.38× low; and B_s is fed
   through B→K form factors (a b→u transition) although the paper gives no B_s channel at all.
6. The Higgs LO subprocess mixture (ggF = 79.7 % of the sample vs 87.6 % of Σσ_H) biases the h-channel
   efficiency by only **0.06–1.3 %**, because ggF dominates both mixtures and the per-subprocess
   acceptances differ by only ±25 %.
7. D mesons are negligible by 9–10 orders of magnitude (D/B ≈ 5·10⁻⁹) and the two prompt.md §3.1/§3.3
   generation requirements are not met — but only the first (HardQCD b-b̄ instead of
   SoftQCD:inelastic) matters, and only through item 5 above.
8. **None of the findings explains the known 2–3.5× excess of the v3 contour over the reference curve.**
   The production sector is verified to ≤7 % on the mixing channels; the residual must be in the
   response/normalisation chain, not here.

### Top fixes by impact on the long-lived (small-θ²) 25 fb⁻¹ boundary

| rank | finding | factor | Δθ²/θ² (25 fb⁻¹ long-lived edge) | confidence |
|---|---|---|---|---|
| 1 | B-hadron yield per bb̄ event in the parent sample (1.320 vs 1.818) | B channel ×1.38 | **−19 %** at m_S≲0.2 (f_A=1); −5 % at 1.06; −3 % at 2.0 | medium-high |
| 2 | B_s fed with b→s (B→K) form factors although the paper has no B_s channel | B mixing ×0.898 | **+3…+5 %** (max at m_S≲0.2; +2…+3 % at 1–2 GeV) | high (paper), medium (correction size) |
| 3 | σ(pp→bb̄) = 500 μb assigned (FONLL ≈ 560 μb; paper's χ_bb̄=8.6·10⁻³ implies 684 μb with σ_inel) | ±11 % | **∓5.5 %** at low mass, ∓4.5 % at 2 GeV | medium |
| 4 | 3-body quartic power-fit shape (`_quartic_ratio`) instead of the q²-dependent Eq. (eq:br-quartic-3-body) | 0.97–1.07 on the quartic rate | **+1.2 %** at 1.06 GeV; ≤0.3 % elsewhere | high |
| 5 | `tab:BR` per-channel agreement ±6 %, net +1.7 % on ΣBR | ±5 % | **∓1…2.5 %** | high |
| 6 | Higgs LO subprocess mixture vs HXSWG | +0.06…1.3 % on ε_H | **≤1 %** | medium |
| 7 | kaon quartic exponential (1.5–3.7× too steep) | 1.5–3.7 on k_K | **0 %** (f_K ≤ 0.08 %) | high |
| 8 | K→πS BR 5.3 % below `tab:BR` | 0.947 | **0 %** (K mixing acceptance is zero) | high |

(Δθ²/θ² = −[δ_A f_A + δ_C(1−f_A)]/(1+f_A), where f_A = mixing-B share of s and f_B,quartic is measured
from outputs/data/signal_scan.json at the boundary θ², so these are the actual v3 numbers: f_A =
1.00/0.87/0.71/0.29/0.44/0.21 and f_B,quartic = 0/0.10/0.20/0.38/0.14/0.01 at m_S = 0.10/0.30/0.48/1.06/
1.98/3.17 GeV.)

---

## 1. Findings table

Locations are `file:line` of the v3 files. "Our value" is what the code produces, "paper value" what the
.tex prints. Only F1, F2, F4 change the answer at the >1 % level.

| # | item | code location | paper eq./table | our value | paper value | discrepancy | impact on s | impact on θ² boundary | proposed fix | confidence / how to verify |
|---|---|---|---|---|---|---|---|---|---|---|
| F1 | B-hadron yield per bb̄ event (parent sample) | work/magnus/samples/bbar_heavy.txt (1055619 rows, max id 799999 ⇒ n_gen 800000); normalisation.json; run_scan.py:54,461 | none (sample issue) | **1.320** B-hadrons per generated event (1.498 per event-with-B; 11.9 % of events have none) | 2·(0.404+0.404+0.101) = **1.818** ground-state B mesons per bb̄ event | **0.726** ⇒ B channel 1.38× low | B channel −27 % | **−19 %** at m_S≲0.2, −5 % at 1.06, −3 % at 2.0 | let the L=1 b-hadrons decay: drop `ParticleDecays:tau0Max = 0.0` (or add ids 513,523,533,5132,5232,5332,10413… to the output and run `moreDecays()`), then recount B±/B⁰/B_s per event | medium-high. Verify by regenerating/recounting with excited decays on. |
| F2 | B_s in the mixing channel | bc5_production.py:130–158 (br_B_to_S(ms,'Bs')); run_scan.py:250,259 counts pdg 531 | **no B_s row exists**: tab:BR has only B^± (L263–285) and L682 explicitly excludes B_s; the only B_s entry is the *quartic* 2-body B_s→SS (L405) | same 10 b→s channels with B→K form factors, m_X=5.36688, Γ_Bs | — | **B_s = 10.10 %** of the b-hadron rows, 10.2 % of the B-channel rate | −10 % of the mixing-B rate | **+3.1…+5.1 %** (drop or correct) | remove B_s from the mixing treatment (paper's own stance) or replace by B_s→φ S, B_s→η(′)S, B_s→f₀S — those need B_s→φ/η form factors that the paper does NOT provide, so they are not computable from this source; dropping B_s is the paper-consistent choice | high for the paper statement; the *size* of the physics correction is unknown (could be ×0.5–×2 of the deleted 10 %) |
| F3 | σ(pp→bb̄) normalisation | run_scan.py:56–57, normalisation.json B.sigma_pb = 5.0e8 pb | not in the paper (external input) | 500 μb | FONLL/NNLO ≈ 560 μb; paper's tab:meson-amounts χ_bb̄ = 8.6·10⁻³ at 13 TeV × σ_inel = 79.5 mb ⇒ 684 μb; × σ_pp^tot ≈ 100 mb ⇒ 860 μb | 0.89 (vs FONLL) / 0.73 (vs paper χ) | B channel ±11 % per 56 μb | **∓5.5 %** at low mass, ∓4.5 % at 2 GeV | quote a measured/FONLL value with its uncertainty (±15 %) and carry it as a systematic | medium; one-line literature check, no code verification possible |
| F4 | 3-body quartic shape heuristic | bc5_production.py:250–281 (_ps3_integral, _quartic_ratio) | eq:br-quartic-3-body (L1402) with M_XX'(q²) of Apps. C.1–C.5 | constant |M|² phase space raised to a fitted power through (0,BR₀) and (m_close,0.1BR₀) | should integrate |M(q²)|² over q²∈[4m_S²,(m_X−m_X')²] with the same form factors | exact/code = 1.007/1.032/1.074/1.029/0.775 (B-quartic total at 0.2/0.5/1.0/1.5/2.0 GeV); per channel at 1 GeV: K 1.055, K*892 1.064, K0*700 1.077, K1_1270 1.058, K*1410 1.019, K*1680 1.053, π 1.054, K0*1430 1.104, **K2*1430 0.658**, **K1_1400 0.524** | quartic-B rate +3…+7 % at 0.5–1 GeV | **+1.2 %** at 1.06 GeV; +0.24 % at 0.48; −0.3 % at 1.98; ≤0.05 % at 0.3 | replace _quartic_ratio by I_exact(m_S)/I_exact(0) with I_exact(m) = ∫dq² |M_XX'(q²)|² √(E₂*²−m²)√(E₃*²−m_X'²); _exact formula in §4.1 | high; reproduction of the paper's closing masses validates it |
| F5 | tab:BR per-channel agreement | bc5_production.py:81–158 | eq:productionBRmesons, tab:BR | π 0.960, K 1.004, K*892 1.016, K*1410 1.000, K*1680 1.024, K0*700 0.964, K0*1430 **1.060**, K1_1270 **1.047**, K1_1400 0.993, K2*1430 1.015; Σ = 3.35597 | 3.3 (eq:total-exclusive) | Σ +1.7 % | B mixing rate +1.7 % | **−0.85 %** at low mass (f_A=1) | none: within the paper's own form-factor uncertainties (Ball/Sun fits carry ±10–20 %) | high; production_validation.json reproduced independently |
| F6 | Higgs subprocess mixture | bundles/higgs_only2.tar.gz:parents.cc (HiggsSM:all=on); higgs_meta.json; run_scan.py:466 uses Σσ_H×0.01×ε_H | prompt.md §3.1/§3.3 wants MG5+Pythia8 | ggF (code 902) = 79.66 % of 199105 events; LO mixture; Σσ_Pythia = 32.22 pb but normalised to Σσ_H = 55.4527 pb | Σσ_H weights: ggF 87.6 %, VBF 6.8 %, WH 2.5 %, ZH 1.6 %, bbH 0.43 %, ttH 0.91 %, tHj 0.13 % | truth/sample efficiency = **1.0006/1.0062/1.0013/1.0116/1.0045** at m_S = 0.3/1/2/5/20 GeV | +0.06…+1.3 % on the H channel | **≤1 %** | reweight events by σ_i^HXSWG/σ_i^Pythia per info.code() (the code is already column 2 of higgs.txt) | medium; the mapping 902=ggF is fixed by σ(code 902) = 25.7 pb and run_report.md "ggF is 79 % of the generated sample" |
| F7 | Kaon quartic shape | bc5_production.py:301–321 (10**(-ms/mclose)) | should be eq:bstossbranching β for K→SS and eq:br-quartic-3-body for K→SSπ | 10^(−m/0.252) / 10^(−m/0.140) / 10^(−m/0.136) | exact: code/exact = 0.646/0.438/0.270 at m_S=0.05/0.10/0.20 (K→SS); 0.59/0.56/**1.45** at 0.05/0.10/0.15 (K0L→SSπ⁰) | 1.5–3.7× too steep | K quartic rate 1.5–3.7× low | **0 %** (K ≤ 0.08 % of s) | use _ps2_ratio(ms,m_K) for K→SS and the exact 3-body integral for K→SSπ, i.e. the same code paths as the B channels | high (formula), no impact |
| F8 | K→πS normalisation | bc5_production.py:161–173 | L269–270 tab:BR: 1.7e-3 (K^±), 7.0e-3 (K⁰_L) | 1.6099e-3, 6.9551e-3 | 1.7e-3, 7.0e-3 | 0.947, 0.994 | K mixing rate | **0 %** (K mixing acceptance = 0; only K quartic contributes) | none; the 5.3 % would need m̄_s = 0.09495 (vs 0.0924) or τ_K = 1.172 ns (vs 1.238 ns) ⇒ 2-significant-figure rounding | high |
| F9 | B_s quartic 2-body | bc5_production.py:227–240 row ('Bs','2b',4.0e-10,2.670) | tab:b-production-quartic L405, eq:bstossbranching | BR(0) = 4.079e-10, closing 2.66999 | 4.0e-10, 2.670 | 1.020 | exact | **0 %** | none — charge-allowed and correct | high |
| F10 | B⁰→SS quartic 2-body | bc5_production.py:236 row ('B0','2b',1.1e-11,2.624) | L413 | BR(0) = 1.275e-11 from eq:bstossbranching with f_B0=0.19 | 1.1e-11 | 1.16 | exact | 0 % | none: the code uses the paper's 1.1e-11; the 16 % is the paper's internal inconsistency | high |
| F11 | 3-body quartic BR(0) anchors | bc5_production.py:227–240 | tab:b-production-quartic | from eq:br-quartic-3-body: K 1.363e-10, K*892 9.08e-11, K*1410 1.904e-11, K0*700 1.181e-10, K0*1430 3.65e-11, K2* 2.48e-11, K*1680 1.030e-11, π 4.74e-12, K1_1270 1.219e-10, K1_1400 7.18e-13 | 1.4e-10, 9.1e-11, 1.9e-11, 1.2e-10, 3.5e-11, 2.5e-11, 9.9e-12, 4.7e-12, 1.2e-10, 7.3e-13 | 0.97–1.04 | — | 0 % | none | high |
| F12 | closing masses of the 3-body channels | bc5_production.py:227–240 (col 5) | tab:b-production-quartic col 3 | exact-shape closing masses: K 1.9981, K*892 1.6999, K*1410 1.3567, K0*700 1.6206, K1_1270 1.4751, K*1680 1.2524, π 2.1495, K1_1400 1.5322, **K0*1430 1.5021**, **K2*1430 1.2740** | 1.998, 1.701, 1.358, 1.621, 1.478, 1.240, 2.149, 1.545, **1.621**, **1.499** | <1 % for 8/10; the two exceptions are paper duplicates/typos | — | ≤0.3 % | keep the tabulated values (used only as the power-fit anchor); note the two suspect rows | high |
| F13 | K_S⁰→πS treatment | bc5_production.py:161–173 returns 0 for K0S | L740: "we assume that the corresponding branching ratio vanishes"; the source has a *commented-out* estimate 1.2e-5/θ² at closing mass 0.360 (L271) | 0 | 0 | — | — | 0 % | none | high |
| F14 | final-state masses of the 3-body modes | bc5_production.py:229–239 (0.700, 0.892, 1.270, 1.430, 1.410, 1.430, 1.680, 1.400) vs B_CHANNELS (0.824, 0.89555, 1.253, 1.425, 1.414, 1.4324, 1.718, 1.408) | the paper prints **no** resonance masses (only K^± 494, K⁰_S,L 498 MeV, L677) | two different mass sets in the same file | PDG: 0.8955, 1.253, 1.425, 1.414, 1.4324, 1.718, 0.824 | 0.892 vs 0.8955 (0.4 %), 1.270 vs 1.253 (1.4 %), 1.680 vs 1.718 (2.2 %), **0.700 vs 0.824 (15 %)** | 0.94–1.04 on the quartic BR(0) | ≤0.2 % | use one PDG mass set everywhere (B_CHANNELS values are better) | high (numbers), low (impact) |
| F15 | B⁰ counterpart prescription for the quartic channels | bc5_production.py:343–347 (same channel, own m_B0, ×0.93) | tab:b-production-quartic caption: "in the case of decays of B⁰ mesons, all the given branching ratios should be multiplied by a factor 0.93" | 0.93×BR(B⁺), shape recomputed with m_B0 | multiply by 0.93 (shape implicitly unchanged) | — | +2 % relative to a pure rescale | <0.1 % | none; the literal reading is implemented | high |
| F16 | B0_SCALE vs the exact width ratio | bc5_production.py:225 (0.93) and :143 (GAM_B0) | L283/L313/L418: Γ_B⁺/Γ_B⁰ ≈ 0.93 | mixing: 0.927350 (exact, from lifetimes); quartic: 0.93 | 0.93 | 0.3 % | <0.1 % | <0.1 % | optionally use Γ_B⁺/Γ_B⁰ also for the quartic B⁰ counterpart | high |
| F17 | kaon quartic multi-mode history | run_scan.py:301–305 (hid=1 for both K⁰_L modes) | — | both K⁰_L modes are put in the *same* history, so the group weight is br_s[start_h] = the first mode's BR only (4.40e-13) and the two mode samples are OR-ed | the two modes are exclusive: BRs should be summed as the history weight (4.41e-13) and one branch sampled | 0.3–1 % of the K⁰_L quartic rate | negligible | 0 % | sample one quartic kaon mode with normalised weights, summed BR as the history weight (as already done for the B quartic at run_scan.py:268–279) | high (code reading), no impact |
| F18 | D mesons absent | bc5_production.py (no D channel); prompt.md §3.1 requires K, D, B | tab:BR-3-body D→Seν: 5.2e-9/θ²/f(x); L391: D⁰→SS strongly suppressed; no D row in any table | 0 | 5.2e-9/θ² (leptonic only) | **D/B = (2χ_cc̄/2χ_bb̄)·(5.2e-9/3.3) = 3.37·1.58·10⁻⁹ = 5.3·10⁻⁹**; **K/D = (1.65/0.058)·(1.7e-3/5.2e-9) = 9.3·10⁶** | 5.3e-9 of the B channel | 0 % | document the check (the prompt is met *by demonstration*, not by implementation) | high |
| F19 | prompt §3.3: "hadron parents from Pythia8 SoftQCD:inelastic" | parents.cc + normalisation.json | prompt.md §3.3 | B sample: SoftQCD:inelastic=off; HardQCD:gg2bbbar/qqbar2bbbar=on; pTHatMin=0; K sample: SoftQCD:inelastic=on | SoftQCD:inelastic for all hadron parents | the bb̄ sample is a *hard-process* LO sample (Pythia's own σGen = 1.366 mb) renormalised by hand to σ_bb = 500 μb | folded into F1 | folded into F1 | either regenerate with SoftQCD:inelastic (physical bb̄ admixture for free) or keep it and fix F1 | high (documented deviation) |
| F20 | prompt §3.3: "Higgs hard production from MG5+Pythia8" | parents.cc: HiggsSM:all=on | prompt.md §3.3 | pure Pythia8 LO (Σσ = 32.22 pb, ggF 79.7 %) | MG5+Pythia8 (NLO) | −42 % on Σσ (but renormalised to Σσ_H); mixture 8 pp off | ≤1.3 % (F6) | ≤1 % | acceptable as a documented shape approximation; the mixture bias is measured at ≤1.3 % | medium |

---

## 2. Checked and already correct (do not re-open)

### 2.1 Master formula, eq:productionBRmesons (L731)
br_B_to_S (bc5_production.py:156): br = xi² m_Q² M²/(32π v² Γ_X)·p/m_X². The paper:
BR = (1/Γ_X) θ² |ξ_ij|² m_Qi² |M_XX'(m_S²)|²/(32π v²)·|p_S|/m_X². This is exactly the 2-body width
|M|²p/(8πm_X²) with |M| = (θ/2)(m_Qi/v)ξ M_XX' from eq:matrixelementhhS; the p/m_X² is the correct
phase-space normalisation (the paper's own L695–697 expression for |p_S| was reproduced algebraically).
Γ_X is the **total** width — the code divides by GAM_BP/GAM_B0/GAM_BS/GAM_KP/GAM_KL, i.e. partial/Γ_total.
Lifetimes: B⁺ 1.638 ps, B⁰ 1.519 ps, B_s 1.520 ps, K⁺ 1.2380e-8 s, K⁰_L 5.116e-8 s (PDG 2018 values,
which the paper does not print — verified as the standard ones). Gamma_BP/Gamma_B0 = **0.927350**.

### 2.2 ξ_ij and quark masses — tab:xi (L659–670)
XI_SB=3.6e-4, XI_DS=3.3e-6, XI_DB=7.9e-5 = exactly ξ_sb, ξ_ds, ξ_db of tab:xi.
M_B_QUARK=4.18, M_S_QUARK=0.0924, M_D_QUARK=0.00467 = the MS-bar values used in the paper's scalar-decay
appendix (m̄_s(2 GeV) = 0.0924, m̄_b(m_b) = 4.18). The b→d channel (π) correctly uses ξ_db and m_d−m_b in
the denominator; all b→s channels use ξ_sb and m_s−m_b (or m_b+m_s for the parity-odd ones). B_MQI is
uniformly m_b, i.e. m_{Q_i} (the initial quark), as in eq:productionBRmesons.

### 2.3 Form factors — Apps. C.1–C.5, bit-exact
An independent transcription of every formula and parameter of L1111–1346 was compared with
ff_ps, ff_sca, ff_vec892, ff_vec_hi, ff_psv, ff_ten, matrix_element at q² = 0, 0.25, 1, 4, 9 GeV² for all
10 channels: **maximum relative difference 0.00e+00** (identical floating-point values).
Parameters verified one by one:

| channel | form factor | code constants | table |
|---|---|---|---|
| K | f₀ = F₀/(1−q²/m_fit²) | 0.33 / 6.16 | tab:pseudoscalar-form-factor-parameters col "B^±,K^±": 0.33±0.04 / 6.16 |
| π | idem | 0.258 / 6.16 | col "B^±,π^±": 0.258±0.031 / 6.16 |
| K→π | f₀ = 0.96 const | 0.96 (m_fit=None) | col "K,π": 0.96 / ∞ |
| K₀*(700) | F₀/(1−a q²/m_B²+b(q²/m_B²)²) | 0.46, 1.6, 1.35, m_B=5.3 | tab:scalar-form-factor-parameters |
| K₀*(1430) | idem | 0.17, 4.4, 6.4 | idem |
| K*(892) | r₁/(1−q²/m_R²)+r₂/(1−q²/m_fit²) | r₁=1.364, r₂=−0.99, m_R=5.279, m_fit=√36.8 | tab:vector-form-factor-parameters1 (m_R = m_B⁺) |
| K*(1410) | (1−2m_V²/(m_B²+m_V²−q²))ξ_∥ + (m_V/m_B)ξ_⊥ | ξ_⊥=0.28, ξ_∥=0.22, m_B=5.3 | tab:vector-form-factor-parameters2 |
| K*(1680) | idem | 0.24, 0.18 | idem |
| K₁(1270) | [sinθ m_A V₀^A + cosθ m_B V₀^B]/m_{K₁} | F₀^A=0.22, a_A=2.4, b_A=1.78; F₀^B=−0.45, a_B=1.34, b_B=0.69; θ=−34°, m_A=1.31, m_B=1.34; m_{K₁}=1.253 | tab:pseudovector-form-factors-parameters1/2 |
| K₁(1400) | [cosθ m_A V₀^A − sinθ m_B V₀^B]/m_{K₁} | same, m_{K₁}=1.408 | idem |
| K₂*(1430) | F₀/[(1−q²/m_X²)(1−a_T q²/m_X²+b_T(q²/m_X²)²)] | 0.23, 1.23, 0.76, poles at m_X (5.279) | app:tensor L1344–1346 |

Algebraic expressions of the 10 matrix elements (signs/phase conventions included; only |M|² matters) all
match: the pseudoscalar M_XP = (m_X²−m_P²)/(m_Qj−m_Qi) f₀ (K, π), the scalar
M_XS̃ = i(m_X²−m_S̃²−q²)/(m_Qi+m_Qj) f₊ (K₀*), M_XV = −2m_X|p_V|/(m_Qi+m_Qj) A₀ (K*),
M_XA = 2m_X|p_A|/(m_Qj−m_Qi) V₀ (K₁), M_XT = −(1/(m_Qi+m_Qj))√(2/3)(m_X|p_T|²/m_T)·2A₀ (K₂*) —
including the abs() of the |m_s−m_b| denominator (a phase, drops in |M|²). The K₂* branch's
`if p<=0: return 0` is the correct kinematic closure.

### 2.4 Quartic normalisation, alpha_from_bh (bc5_production.py:184–190)
BR(h→SS) = α²|p_S|/(16π m_h² Γ_h) (L355) ⇒ Γ(h→SS) = α²β/(32π m_h) (the 1/2! for identical final-state
particles is implicitly included in the paper's 16π). alpha_from_bh inverts exactly this:
Γ_SS = b_h/(1−b_h)·Γ_h^SM (correct: BR = Γ_SS/(Γ_SM+Γ_SS)) and α = √(32π m_h Γ_SS/β).
**Dimension: GeV.** Numbers: α(1 GeV) = **0.719525**, α(1e−6) = 0.719479, α(2 GeV) = 0.719658.
The paper's eq:BR-h-to-SS normalisation is reproduced exactly: with α = 1 GeV and β = 1,
Γ_SS/Γ_h = 1/(32π m_h Γ_h) = 0.01951 vs the paper's "≈2.0·10⁻²" (L356). So α = 0.7195 ⇔ BR(h→SS) = 0.01
⇔ 0.7195²×0.01951 = 0.0101 ✓, and the paper's tables (α = 1 GeV ⇔ BR ≈ 0.02, L393) use the **same**
convention. Also verified: "α < 3.1 GeV" from BR_inv < 0.19 (L363) is the *linearised* inversion
(0.19×51.25 ⇒ α = 3.12); the exact inversion would give 3.47 — a 0.5 %-level convention difference,
irrelevant here.

### 2.5 Quartic 3-body and 2-body absolute normalisation (strong validation)
Implementing eq:br-quartic-3-body with matrix_element(ch,√q²,m_X,m_X') (the q²-dependent form factors) and
eq:bstossbranching with tab:meson-decay-constants reproduces the paper's table:
2-body: B_s→SS 4.079e-10 (paper 4.0e-10, 1.02); K⁰_L→SS 4.450e-13 (4.4e-13, 1.01); K⁰_S→SS 7.789e-16
(7.8e-16, 1.00); B⁰→SS 1.275e-11 (1.1e-11, 1.16). 3-body BR(0): F11, 0.97–1.04 for all ten channels.
**The paper's closing masses are reproduced to <1 % by this exact shape for 8 of 10 channels**
(K 1.9981/1.998, K*892 1.6999/1.701, K*1410 1.3567/1.358, K⁰*700 1.6206/1.621, K₁(1270) 1.4751/1.478,
K*(1680) 1.2524/1.240, π 2.1495/2.149, K₁(1400) 1.5322/1.545; r_exact(m_close) = 0.1000/0.0995/0.0995/
0.0998/0.0987/0.1058/0.1002/0.0958). Two exceptions, both suspected source artefacts: K₀*(1430)
(exact 1.5021 vs printed 1.621 — the identical 1.621 is printed for K₀*(700)) and K₂*(1430) (exact 1.2740
vs 1.499). This is the strongest evidence that the paper's curves are the q²-dependent
eq:br-quartic-3-body, not the constant-|M|² phase space.

### 2.6 Other verified items
* tab:b-production-quartic transcribed row-by-row (12 B rows + 5 K rows in QUARTIC/QUARTIC_KAON): every
  BR(0) and closing mass matches the .tex, including the 0.93 B⁰ prescription (F15).
* B_CHANNELS thresholds: pS returns exactly 0 at m_S = m_X−m_X' (checked at thr×0.999, ×0.9999, ×1.0,
  ×1.0001 for B⁺→K⁺S: 6.557e-2, 2.070e-2, 0, 0) and p_S(0) = (m_X²−m_X'²)/(2m_X) to all digits.
* m_S→0 limit of every channel is flat (drift <2e-13 between m_S = 1e-9 and 1e-6) — no numerical
  instability at the low-mass end of the scan.
* Σbr_B_to_S = 3.35597; excluding π (eq:total-exclusive sums over X_s) 3.34349, vs 3.3±0.8.
* ΣSIGMA_H = 55.4527 pb; the 9 values are the YR4 13 TeV numbers (ggF 48.58, VBF 3.782, WH 1.369,
  ZH 0.883, bbH 0.2384, ttH 0.5065, tHj 0.0740, tHW 0.0152, tHb 0.0046).
* σ_inel = 79.5 mb = 7.95e10 pb ✓ (Pythia 8.316 SoftQCD:inelastic value quoted in normalisation.json).
* n_gen vs the files: B 800000 (max event id 799999 ✓, 1055619 rows ✓), K 200000 (330171 rows ✓),
  H 200000 (max id 199999 ✓, 1235376 rows ✓, 199105 unique events).
* S widths untouched by production: bc5_production uses only sw.V_HIGGS; run_scan uses only sw.HBAR_C and
  sw.gamma_total (lifetime). Γ_μμ(1 GeV)/θ² = 6.842e-9 GeV = y_μ²m_Sβ³/(8π) (β³ = 0.934);
  Γ_S(1 GeV)/θ² = 6.488e-8 GeV ⇒ cτ = 30.4 m at θ² = 1e-10.

---

## 3. Parent-sample / normalisation section (the part that actually matters)

### 3.1 b-hadron yield (F1)
bbar_heavy.txt: 1055619 rows, 704535 unique event ids out of n_gen = 800000 (max id 799999).
Composition B⁺ 0.4503, B⁰ 0.4477, B_s 0.1009, B_c 0.0010. Therefore **1.3195 b-hadrons per generated
event**, of which 1.318 are B mesons. With the PDG b-hadron fractions (B⁺ 0.404, B⁰ 0.404, B_s 0.101,
Λ_b 0.085, B_c 0.001) two b quarks give **1.818 ground-state B mesons per bb̄ event**. The 27.4 % deficit
cannot be explained by b-baryons alone (9.1 % of b quarks): it needs ~18 % of b-hadrons to be excited
states, and indeed 11.9 % of events contain no B-hadron row at all. The cause is visible in the generator
(parents.cc): `ParticleDecays:limitTau0 = on` + `ParticleDecays:tau0Max = 0.0` keep every b-hadron
stable, and the output filter is BONLY = {521,511,531,541} — so the L=1 b-hadron resonances (B*, B**, B_s*,
Λ_b*) that in reality decay into the ground-state B mesons are never counted. The B_s/(B⁺+B⁰) ratio in the
sample is 0.112 vs PDG 0.125, i.e. B_s is *more* deficient, as expected since the B_s** fraction is larger —
an independent confirmation of the mechanism.
Consequence: the B channel is **1.818/1.320 = 1.38× low**; in s −27 %, in θ² −19 % at m_S≲0.2 GeV (f_A=1),
−5 % at 1.06 GeV, −3 % at 2.0 GeV. Verification is cheap: add the L=1 ids to the output (or call
`pythia.moreDecays()` after `next()`) and recount; a re-generation is not strictly needed.

### 3.2 Higgs subprocess mixture (F6)
higgs_only2.tar.gz:parents.cc uses `HiggsSM:all = on` with `25:mayDecay = off`, |η_h|<6, and writes
`pythia.info.code()` in column 2. Deduplicated per event the mixture is 902: 158607 (79.66 %), 907: 16146
(8.11 %), 905: 6914 (3.47 %), 906: 5956 (2.99 %), 901: 5543 (2.78 %), 904: 3760 (1.89 %), 908: 1594
(0.80 %), 909: 581 (0.29 %), 903: 4 events. Pythia's own σGen for this sample is **32.22 pb** (job.log
line 192), i.e. 58 % of Σσ_H = 55.4527 pb — so the sample *must* be renormalised to Σσ_H (which
run_scan.py:466 does) and the only residual question is the mixture. Measured per-code efficiencies
(J = p/θ² at θ² = 1e-13, the long-lived limit; 2000 events/code, 5–14 % stat. error) at m_S = 1 GeV:
901 2.82e6, 902 2.70e6, 905 2.76e6, 904 2.69e6, 906 2.43e6, 907 2.28e6, 908 1.86e6, 909 1.45e6.
Re-weighting with the Σσ_H fractions gives **truth/sample = 1.0006, 1.0062, 1.0013, 1.0116, 1.0045 at
m_S = 0.3, 1, 2, 5, 20 GeV** (two assignment hypotheses, plus a mapping-independent ggF-vs-rest bound
1.0006–1.0125). So the LO mixture is a **≤1.3 %** effect, not a 20 % one. (A hypothetical 10 % error in the
h-channel normalisation would move θ² by −2.6 % at 1.06 GeV and −6.4 % at 3.17 GeV, where the H share is
33 % and 78 %.)

### 3.3 The other cross-section inputs
σ(pp→bb̄) = 5.0e8 pb is an *assumption* (documented as such in normalisation.json); FONLL/NNLO at 13 TeV
is ≈ 560 μb and the paper's own χ_bb̄ = 8.6e-3 (tab:meson-amounts) times σ_inel = 79.5 mb implies 684 μb —
i.e. the assigned value is at the low end, +11…+37 %. Since the B channel is the whole signal below
0.2 GeV, this is a first-order systematic (F3). σ_inel and Σσ_H enter only through the K and H channels:
σ_inel's 79.5 mb is Pythia's own inelastic value, consistent with the sample; a 10 % error there is a 0 %
error on the answer because K is negligible.

### 3.4 B_s (F2) — what the paper literally says
L682: "The production of a scalar from mesons is possible through the flavor changing neutral current, so
the production from D_s, J/ψ, B_s, B_c and Υ mesons does not have any advantage with respect to the
production from D⁰, D⁺, B⁻ and B⁰, **while their amount at any experiment is significantly lower.
Therefore, we will discuss below only production from later mesons.**"
tab:BR (L263–285) contains no B_s row; eq:total-exclusive sums over X_s for B^±/B⁰ only. The only B_s
entry in the paper is the *quartic* 2-body B_s→SS (L405, 4.0e-10), which is charge-allowed and correctly
implemented (F9).
Our code applies the ten B→K form factors to B_s with ξ_sb. That is not the physical process:
B_s → K S requires b→u (|V_ub|²-suppressed), so the code **over-estimates** B_s→K S by many orders of
magnitude in reality; the physical b→s FCNC with an s̄ spectator gives B_s → φ S, B_s → η(′)S, B_s → f₀S,
whose rates need B_s→φ/η transition form factors that are **not in the paper and not in
boiarska1904_formfactors.json**. Since B_s is 10.10 % of the b-hadron rows and 10.2 % of the B-channel
rate, dropping it moves the long-lived boundary by +3.1…+5.1 %. Recommendation: drop B_s from the mixing
channel and state that the B_s→φ/η S channels are not computable from the source; do *not* silently keep
the b→u extrapolation (keeping it is "conservative" only if the true B_s→φ S rate is smaller, which is not
established).

### 3.5 D mesons (F18) — the arithmetic
tab:BR-3-body gives BR(D→Seν)/f(x)/θ² = 5.2·10⁻⁹; the inclusive B→X_s S sum is 3.3/θ². With the paper's
production fractions (χ_cc̄ = 2.9·10⁻², χ_bb̄ = 8.6·10⁻³ at 13 TeV) the number of D and B mesons per
collision are 2χ_cc̄ = 0.058 and 2χ_bb̄ = 0.0172, so **D/B = (0.058/0.0172)·(5.2e-9/3.3) = 5.3·10⁻⁹** and
**K/D = (1.65/0.058)·(1.7e-3/5.2e-9) = 9.3·10⁶** (K multiplicity per collision 1.65 measured in
minbias_light.txt: 330171 rows / 200000 events). The D⁰→SS quartic mode is "strongly suppressed" by the
paper's own words (L391) and no D row appears in tab:b-production-quartic. D mesons can therefore be
declared negligible by a documented check (~10⁻⁹ of the signal), which is what prompt.md §3.1 should be
taken to allow.

### 3.6 prompt.md §3.3 deviations (F19, F20)
* "强子母样本统一由 Pythia8 SoftQCD:inelastic=on 生成": the K sample is; the B sample is a HardQCD
  gg/qq̄→bb̄ LO sample with pTHatMin=0 renormalised by hand (Pythia's σGen for it is 1.366 mb, i.e. 2.7×
  the assigned 500 μb). It carries no minimum-bias overlay, so per-collision B multiplicities are only
  meaningful after the (wrong, F1) hand normalisation. Same defect family as F1/F3.
* "Higgs 硬产生使用 MG5+Pythia8": pure HiggsSM:all (LO, ggF 79.7 %, Σσ 32.2 pb). Renormalisation to Σσ_H
  is done; the mixture bias is measured at ≤1.3 % (F6), so the practical risk is the *shape* of the Higgs
  pT/y spectrum (ggF pT is softer at LO), which can only be quantified by running the response on an MG5
  sample — outside this audit.
* "按 Boiarska 文献逐道纳入运动学允许的 K、D、B 等介子过程及电荷共轭": B⁰/B̄⁰ and B^± are both counted
  (pdg ±521, ±511) with the correct 0.927 B⁰/B⁺ ratio ✓; K⁰_S is absent from the sample and its production
  is assumed zero per L740 ✓; charge-conjugate channels are covered ✓.

---

## 4. The two non-paper heuristics, quantified

### 4.1 _quartic_ratio (F4) — what the paper prescribes and the exact replacement
The paper prescribes **no** closing-mass interpolation. App. D.1 (L1363–1408) gives
M(X→X'SS) = (α/m_h²)(m_Qi/2v)ξ_ij M_XX'(q²) and Eq. eq:br-quartic-3-body (L1402):
Γ = |ξ|²m_Qi²α²/(512π³m_X³v²m_h⁴) ∫_{4m_S²}^{(m_X−m_X')²} |M_XX'(q²)|² √(E₂*²−m_S²)√(E₃*²−m_X'²) dq²
with E₂* = √q²/2, E₃* = (m_X²−q²−m_X'²)/(2√q²). The *exact* replacement (drop-in for _quartic_ratio; note
matrix_element(ch,√q²,m_X,m_X') already returns M at q² because it builds q2 = mS² and p = pS(m_X,m_X',√q²)):

    def _ratio_exact(ms, channel, mX, mXp, n=2000):
        """I(ms)/I(0), I = ∫ |M_XX'(q²)|² √(E2²−ms²) √(E3²−mXp²) dq²   [eq:br-quartic-3-body]"""
        lo, hi = 4*ms*ms, (mX-mXp)**2
        if hi <= lo: return 0.0
        tot = 0.0
        for i in range(n):                    # q² = lo+(hi−lo)sin²(πt/2) removes both √ endpoints
            t = (i+0.5)/n; s2 = sin(0.5*pi*t)**2
            q2 = lo+(hi-lo)*s2; dq2 = (hi-lo)*0.5*pi*sin(pi*t)/n
            sq = sqrt(q2); e2 = 0.5*sq; e3 = (mX*mX-q2-mXp*mXp)/(2*sq)
            a2, b2 = e2*e2-ms*ms, e3*e3-mXp*mXp
            if a2 > 0 and b2 > 0:
                M = matrix_element(channel, sq, mX, mXp)
                tot += M*M*sqrt(a2)*sqrt(b2)*dq2
        return tot                            # caller divides by the same integral at ms=0

(the absolute normalisation can stay anchored on tab:b-production-quartic BR(0), itself validated to
3–4 %: F11). Error of the current power-fit, per channel, ratio exact/code at m_S = 0.5/1.0/1.5/2.0 GeV:

| channel | 0.5 | 1.0 | 1.5 | 2.0 |
|---|---|---|---|---|
| K | 1.021 | 1.055 | 1.079 | 0.999 |
| K*(892) | 1.025 | 1.064 | 1.058 | 0.659 |
| K₀*(700) | 1.033 | 1.077 | 1.038 | 0.794 |
| K₀*(1430) | 1.066 | 1.104 | **0.618** | — |
| K₁(1270) | 1.026 | 1.058 | 0.976 | 0.047* |
| K₁(1400) | **0.745** | **0.524** | 0.911 | — |
| K₂*(1430) | **0.899** | **0.658** | **0.266** | — |
| K*(1410) | 1.008 | 1.019 | 0.960 | — |
| K*(1680) | 1.018 | 1.053 | 1.014 | — |
| π | 1.020 | 1.054 | 1.083 | 1.054 |

(*K₁(1270) at 2.0 GeV is at 99.8 % of its threshold, both shapes ≈ 0; the ratio is meaningless there.)
Weighted over the ten channels with their tabulated BR(0) the **total** B-quartic rate is code/exact =
1.007/1.032/1.074/1.029/0.775 at 0.2/0.5/1.0/1.5/2.0 GeV (the 2-body B_s/B⁰→SS rows, which are exact in
the code, are included and dilute the error). Because the quartic-B share of s is 0/10/20/38/14 % at
0.10/0.30/0.48/1.06/1.98 GeV, the propagated effect on the boundary is **+1.2 % at 1.06 GeV, +0.24 % at
0.48, −0.30 % at 1.98, ≤0.05 % at 0.3**.

### 4.2 Kaon quartic exponential (F7)
quartic_kaon_channels still uses br0·α(m_S)²·10^(−m_S/m_close). The exact shapes (2-body β and the
const-|M|² 3-body integral, which is *exact* here because f₀^{Kπ} = 0.96 has no q² dependence) give:

| m_S | K⁰_L→SS code/exact | K⁰_L→SSπ⁰ code/exact | K⁺→SSπ⁺ code/exact | total K-quartic code/exact |
|---|---|---|---|---|
| 0.02 | 0.836 | 0.765 | 0.760 | 1.20 |
| 0.05 | 0.646 | 0.590 | 0.583 | 1.55 |
| 0.10 | 0.438 | 0.555 | 0.556 | 2.28 |
| 0.15 | 0.318 | 1.446 | 1.697 | 3.13 |
| 0.20 | 0.270 | closed | closed | 3.70 |

i.e. up to 3.7× too steep at 0.2 GeV (and *too shallow* above the 3-body closure at 0.15–0.18 GeV). The
propagated effect is nevertheless **zero**: the K channel is ≤0.08 % of s at every boundary point
(outputs/data/signal_scan.json: s_K/s_tot = 0 at 0.1–0.19 GeV, 8e-4 at 0.046). Fix for correctness, not
for the answer — the same code path as the B channels does the job.

---

## 5. Source ambiguities (what the paper literally prints)

1. **The "0.93" caption appears in both tables** (tab:BR L283 and tab:b-production-quartic L418), with the
   same wording; L313 derives it as Γ_B⁺/Γ_B⁰ ≈ 0.93 ("the branching ratios BR(B⁰→X'⁰S) differ from
   BR(B⁰→X_{s/d}⁰S) only by the factor Γ_B⁺/Γ_B⁰", with the B⁰ mislabel). So it applies to the mixing
   channels as well, not only to the quartic table — and it is one factor, not two: the code's mixing
   channels get it exactly (Γ_B⁰ = ħ/1.519 ps ⇒ ratio 0.9274), the quartic B⁰ counterparts get the paper's
   rounded 0.93. **Recommendation: keep as is.**
2. **K→SS closing mass 0.252 GeV is above the kinematic closure** m_K/2 = 0.24881. The exact 2-body shape
   reaches BR/BR(0) = 0.1 at 0.24756, so the printed 0.252 cannot be the "factor-10" point. The code uses
   0.252 only as the exponential anchor (harmless once F7 is fixed). K⁰_S→SS carries the same 0.252,
   consistent with "same shape, different normalisation".
3. **K₀*(1430) closing mass 1.621 is a duplicate of K₀*(700)'s**; the exact 3-body shape gives 1.5021.
   K₂*(1430): printed 1.499, exact 1.274. Both are almost certainly transcription/plot-reading artefacts;
   the code inherits them as power-fit anchors (≤0.3 % effect after F4 is fixed).
4. **The 3-body BR(0) of B⁰→SS (1.1e-11) is 16 % below the paper's own eq:bstossbranching** (1.275e-11
   with f_B0 = 0.19, ξ_db = 7.9e-5, m_B0 = 5.27965, ħ/1.519 ps). No action: the tabulated value is used.
5. **eq:BR-h-to-SS is linearised** (BR ≈ 2.0·10⁻²α²/GeV²) while a benchmark at BR = 0.01 needs
   Γ_SS = 0.01/0.99·Γ_SM. The code uses the exact inversion (0.5 % effect on α). The paper's α < 3.1 GeV
   bound (L363) is likewise the linearised one (exact: 3.47).
6. **No B_s channel is given** (L682) although B_s→SS *is* given (L405). The code's B_s mixing
   implementation is an extrapolation not supported by the source (F2).
7. **No resonance masses, no Γ_X values, no v/m_h/G_F numbers** are printed (E.1–E.3 of the reference md);
   the code's PDG-2018 values are appropriate, but the two different resonance-mass sets inside
   bc5_production.py (F14) should be unified.
8. **prompt.md §4.2 s = 3 vs the code's s_cut = 2.3** (outputs/data/exclusion_report.json: "s_cut": 2.3,
   "cl": 0.8997): the extracted boundaries correspond to 90 % CL, not the 95 % asked for in the prompt.
   A factor 3/2.3 = 1.30 in the s threshold is NOT a small effect: θ² ∝ s for the θ²-driven channels
   (h→SS, quartic, i.e. 23 % lower boundary) and θ² ∝ √s for the θ⁴-driven mixing channels (12 % lower
   boundary); at m_S ≲ 0.2 GeV, where f_A = 1, the shift is the full 12 %, at 2–3 GeV (f_A ≈ 0.2) it is
   ≈ 21 %. Out of scope of the production audit but the single largest known normalisation difference —
   worth confirming whether the final curves are meant to be 90 % or 95 % CL.

---

## 6. What this audit does NOT explain
The v3 long-lived boundary sits 2.0–3.5× above the digitised reference curve in 0.3–3 GeV
(logs/audit_verification.md §7.9) and 0.28× at 62.5 GeV. None of the findings above produces more than a
19 % shift in θ² (F1), and the production sector as a whole is verified to ≤7 % on the mixing channels and
≤1.3 % on the h channel. The residual has to be sought in: the detector response / acceptance integration
(SensCalc table-vs-curve factor, z-window, daughter cuts), the S decay-width hadronic model, the reference
curve's own digitisation, and the fact that the whole low-mass B channel rests on one hand-normalised LO
bb̄ sample (F1+F3+F19 are the same underlying defect and together can move the low-mass boundary by up to
~25 %).
