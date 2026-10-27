# BC5 model audit — F05 / F06 / F07 / F08 / F13

Auditor: audit-model.  Run dir RUN = E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3.
All scripts and machine-readable results are in `work/review_check/model/` (f05_*.py,
f06_*.py, f07_quartic.py, f08_bs.py, f13_reference.py; f0*_results.json).  Production files
were only read.

## Verdict table

| id | claim (GPT review) | verdict | decisive number | how measured |
|---|---|---|---|---|
| F05 | the merged multi-meson width is split 1/3-1/3-1/3 with no basis and the eta is not decayed | **CORRECT** | visible rate moves by **+13 % (physical 4pi isospin) to +20 % (all-4pi)** at m_S=1.99 GeV; <2 % at 1.0-1.3 GeV | fixed-kinematics controlled test (N=6000/composition, real kernel) + full B/K/H pipeline (nsamp=2000) at the 25 fb crossing |
| F06 | run_scan keeps the FIRST Pythia event-record copy of the Higgs, not the final one | **CORRECT** | eps_H(last)/eps_H(first) = 0.960+-0.045 (5 GeV), 0.945+-0.041 (30), 0.972+-0.018 (62) at the crossing lifetimes -> boundary theta^2 too low by 3-5 %; invariant-mass choice changes eps by <=0.05 % | event-record statuses from a local Pythia 8.316 listing + 3000 Higgs parents, ndecay=8, real kernel |
| F07 | the X->X'SS total rate uses the q^2-dependent matrix element but the generator uses unweighted nbody | **CORRECT** | affected (quartic) channels are **41 % / 5.7 % / 0 %** of the signal at the 25 fb crossing for m_S=1/2/4 GeV (53 % / 11 % / 0 % at 300 fb) — NOT tiny at 1 GeV | run_mass with the real samples; quartic histories separated by hid |
| F08 | B_s: zero-mass recoil + hard cutoff at 4.347 GeV is not a complete model | **PARTLY** | B_s is 1.7 % / 1.0 % / 0.05 % / 0 % of the signal at 3.5/4.0/4.3/4.4 GeV; cutting it moves the 25 fb crossing by +0.046 % (+0.004 % at 300 fb) -> invisible in the curve | production rates x parent counts x scan channel fractions; eq.(BKphi) evaluated for the continuous phi-S rate |
| F13 | the validation model maps single-charge-state reference columns onto isoscalar widths, and validate_widths checks the running-mass gamma_qq | **CORRECT** | leftover of the reference table at 0.5 GeV = 0.29670 = BR(pi+pi-)/2 exactly; Gamma_ours(KK)/Gamma_ref(K+K-) = 2.00-2.46 and Gamma_ours(pipi)/Gamma_ref(pi+pi-) = 1.36-1.63 (1.5 expected); ss running/fixed = 0.90 ... 0.58 over 2.2-7 GeV | reference CSV + widths(); digitised Winkler table; validate_widths' own functions |

---

## F05 — 'multi' composition (bc5_lib.py:305-320)

Code: the total width of the merged channel is the Winkler toy width
Gamma = C m^3 sqrt(1-(4m_pi)^2/m^2), C = 5.1e-9 GeV^-2 (scalar_widths.py:340-344);
bc5_lib.py:305-312 splits it into three equal pieces (2pi+2pi-, pi+pi-2pi0, eta eta pi+pi-),
then the eta is written with charge 0 and never decayed (so eta->pi+pi-pi0, 22.9 %, is lost),
and for m_S < 2m_eta+2m_pi = 1.3749 GeV the third branch silently falls back to 4 charged
pions (313-316).  BR(multi) from channel_branching: 0.008 / 0.081 / 0.173 / 0.420 at
m_S = 1.0 / 1.3 / 1.6 / 1.99 GeV — the review's "42 % at 1.99 GeV" is right.

Controlled test (f05_multi_composition.py + f05_addendum.py): fixed scalar four-momentum
(|p| = 60 GeV, eta_S = 3.0 and 4.0, vertex (0,0,1.5 m)), real bc5_lib._passes_pair transport
kernel, N = 6000 decays per composition; the base acceptance is the real channel mixture
(A_base = 0.72 at eta=3, 0.40 at eta=4 for m=1.99; the production value is ~0.44, so the
eta=4 point is the production-like one).  Delta R/R = BR_multi (A_var - A_cur)/A_base.
The quoted SE of each A is <=0.7 % (see f05_results.json), i.e. the acceptance differences
are resolved at >10 sigma.

| m_S | BR_multi | A_cur | A(4 charged pi) | dR/R | A(rho rho) | dR/R | dR/R(eta decays) | dR/R(isospin 4pi) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.00 | 0.008 | 0.305 | 0.385 | +0.2 % | closed | - | -0.0 % | -0.0 % |
| 1.30 | 0.081 | 0.344 | 0.424 | +1.6 % | closed | - | -0.0 % | -0.4 % |
| 1.60 | 0.173 | 0.228 | 0.475 | +10.9 % | 0.519 | +12.8 % | 0.0 % | +6.4 % |
| 1.99 | 0.420 | 0.266 | 0.536 | **+28.5 %** | 0.559 | **+30.9 %** | +0.9 % | **+17.3 %** |

(eta_S = 4.0, production-like; at eta_S = 3.0 the same entries are +14.7/+16.5/+2.8/+9.2 %.
"isospin 4pi" = 2/3 (2pi+2pi-) + 1/3 (pi+pi-2pi0), the standard isoscalar 4-pi composition.
At m=1.0/1.3 the current code's fallback already IS that composition, hence no effect.)

Full-pipeline cross-check (f05_production.py, real B/K/H parent samples, nsamp=2000,
ndecay=8, at the 25 fb scan crossing theta^2): m_S=1.99: s=2.073+-0.088 (current),
2.496+-0.106 (all-4pi, +20 %), 2.339+-0.097 (isospin 4pi, +13 %).  m_S=1.6:
+4.4 % / -0.8 % (the effect scales with BR_multi, so it is small there).
Conclusion: at the hadronic edge the arbitrary composition is a **10-20 % systematic on the
visible rate** (up to +31 % for a rho-rho composition in the controlled test).  This is a
real model uncertainty, not a cosmetic one.

Fix: keep Gamma_multi fixed (do not retune C), but define the exclusive composition with
isospin weights and decay the eta; e.g. replace the equal split by
  x<2/3 -> 2pi+2pi-, x<1 -> pi+pi-2pi0  (4pi isoscalar)
and add eta eta pi pi / rho rho branches only with a stated weight, decaying each eta with
its PDG charged modes.  Better still use the Pythia hadronisation library for multi (as is
already done for gg/ss/cc/bb, bc5_lib.py:256-259) and quote the difference to the analytic
composition as the systematic.  Do this together with F03 (the 4-body generator used here is
the non-uniform nbody of F03) and F01.

---

## F06 — which Pythia Higgs copy (run_scan.py:380-405; generator)

What is written: the generator loops over the whole event record and writes **every**
p.id()==25 entry with |eta|<6, in record order, columns
event, pythia.info.code(), px, py, pz, E, xProd, yProd, zProd, tProd
(bundles/higgs_only2.tar.gz parents.cc:96-99; identical archived copy
work/physics_review_20261024/higgs_parents_archived.cc:85-89).  There is no status, no
index, no mother/daughter information, and column 2 is the *subprocess* code (902=ggF, ...),
not the status.  run_scan.py:386-387 therefore does
np.unique(ev, return_index=True) -> the **first** row of each event.

Which copy is physical: a complete Pythia 8.316 event listing (saved as
f06_pythia_status_excerpt.txt, from E:/AAAwork/LLP-AGI/temp/pythia_job/out4_3.txt) shows the
same particle in the record as
  status -22 (hard process, pT = 0, E = m_H) -> status -44 (parton-shower copies, momentum
  updated by ISR) -> status **-62 (last, final pre-decay copy)**.
So the correct copy is the final status-62 copy, i.e. the **last** row per event; the
hard-process copy is not even present in higgs.txt: with pT = 0 its pseudo-rapidity is
infinite, so the |eta|<6 cut always rejects it (verified: 0 of 1,235,376 rows have
pT < 1e-4 and only 20 have |E-125.25|<0.01 GeV).  The first written row is thus the first
shower copy (-44).

Effect (f06_higgs.py; 3000 Higgs parents, ndecay=8, identical parent subset for both
variants, real kernel; ratio eps(last)/eps(first), bootstrap SE over the 3000 parents):

| m_S | ctau at 25 fb crossing | ratio at ctau=10 m | ratio at crossing lifetime |
|---:|---:|---:|---:|
| 5.0 | 15.1 m | 0.960 +- 0.045 | ~0.96 (boundary theta^2 too low by ~4 %) |
| 30.0 | 94.7 m | 0.945 +- 0.041 | ~0.945 (~5.5 %) |
| 62.0 | 680 m | 0.972 +- 0.018 | ~0.97 (~3 %) |

(The full ctau scan: 0.974/0.964/0.960/0.960 at 5 GeV; 0.984/0.943/0.945/0.946 at 30 GeV;
1.044/0.991/0.972/0.972 at 62 GeV for ctau=0.1/1/10/100 m.)  Since the H channel dominates
the signal above ~5 GeV, the first-copy choice makes the Higgs efficiency 3-5 % too high,
i.e. the crossing theta^2 3-5 % too low — a few-percent effect, not the 33 % of the review's
F01 table (that 0.667 was the rotation effect, not first/last).  The review's own
first/last diagnostic (0.9445 at 62.5 GeV from 10k parents) agrees with ours.

Invariant mass: the rows carry the forced width (mean 125.2500 GeV, BW with Gamma=4.1 MeV,
0.1-99.9 % quantiles 125.096-125.407).  Replacing P.M_H by the row's invariant mass in
two_body changes eps_H by at most 5e-4 (ratios 1.0000-1.0005).  It is therefore not
numerically required, but it is the consistent thing to do; the bigger issue is the copy.

Fix: in the generator write, in addition to the four-momentum, the status and the
mother/daughter indices (or at least status) and do NOT pre-filter with |eta|<6 on the
Higgs (the eta_S cut is applied later, per scalar).  In run_scan select the last occurrence
with status 62 per event (the final pre-decay copy) instead of return_index=first, and use
that row's invariant mass in L.two_body.  If the sample cannot be regenerated, "last row per
event" is the correct proxy for status 62 (document it) — not the first row.

---

## F07 — quartic channels: rate integral vs generator

Confirmed inconsistency:
* bc5_production.py:301-334 integrates |M_XX'(q^2)|^2 explicitly:
  m = matrix_element(channel, sqrt(q2), mX, mXp); mod2 = m*m inside the q^2 integral;
* but quartic_channels() (bc5_production.py:423-451) returns only
  (mes, kind, mfin, BR) — the channel label is dropped at line 444 — and
  run_scan.sample_quartic_decay_S (run_scan.py:223-258) calls
  L.nbody(m, [ms, ms, spec], rng) (line 234): uniform q^2, no matrix element, and the
  generator does not even receive the channel label.  The total rate and the differential
  decay therefore describe different models.  The size of the shape error is not small:
  for the 3-body modes |M(q^2)|^2 varies over the allowed q^2 by factors
  M^2(hi)/M^2(lo) = 2.1-5.1 (B+K), 3.6-8.7 (B+pi), 0.15-0.31 (K0*(700), K0*(1430)) at
  m_S = 1-2 GeV; for the p-wave channels the matrix element vanishes at the endpoint.

How much of the signal is affected (f07_quartic.py: real samples, nsamp=1200, ndecay=8,
quartic histories isolated by hid, same pipeline as run_scan; s25 at the scan crossing):

| m_S | quartic / total @25 fb | quartic / total @300 fb | quartic / B-channel @25 fb |
|---:|---:|---:|---:|
| 1.0 | **41.2 %** | **52.9 %** | 62.7 % (93.9 % at 300 fb) |
| 2.0 | 5.7 % | 11.0 % | 9.1 % (38.1 % at 300 fb) |
| 4.0 | 0 % | 0 % | 0 % |

So the review's parenthetical "they should be tiny" is wrong at m_S ~ 1 GeV: the quartic
channels are theta^2-independent and form the floor of the B channel exactly where the
long-lived boundary is at low mass; at 1 GeV they are 40-53 % of the signal.  Whether that
is worth fixing: yes for m_S <~ 1.5-2 GeV (it moves the low-mass long-lived edge), no for
m_S >~ 3 GeV.

Fix: propagate the channel label through quartic_channels() -> sample_quartic_decay_S() and
importance-sample q^2 with weight |M_XX'(q^2)|^2 * sqrt(E2*^2-ms^2) * sqrt(E3*^2-mfin^2)
(rejection or per-event weight), then generate the SS angles isotropically in the q rest
frame and boost.  Do it together with the F03 phase-space rewrite; the integral and the
generator must share one matrix-element function.

---

## F08 — B_s inclusive rate and the 4.347 GeV cutoff

Code: br_Bs_inclusive (bc5_production.py:135-168) uses the zero-mass spectator formula with
a zero-mass recoil (B_CHANNEL_MASS['Bs_incl'] = 0.0, line 132, used in run_scan.py:190-191)
and returns 0 for m_S >= m_Bs - m_phi = 4.347419 GeV (line 163).  The exclusive B_s -> phi S
rate from Winkler eq.(BKphi) with the K* form factor and m_Bs in place of m_B
(f08_bs.py: bs_phi_S) is

| m_S [GeV] | 0.5 | 1.0 | 2.0 | 3.0 | 4.0 | 4.30 | 4.34 | 4.3474 | 4.40 |
|---|---|---|---|---|---|---|---|---|---|
| BR_phiS/theta^2 | 0.455 | 0.456 | 0.453 | 0.421 | 0.195 | 0.0179 | 1.2e-3 | 1.6e-7 | 0 |
| BR_inclusive | 4.53 | 4.29 | 3.41 | 2.18 | 0.910 | 0.590 | 0.551 | 0.544 | 0 |
| ratio | 0.10 | 0.11 | 0.13 | 0.19 | 0.21 | 0.030 | 0.0022 | 3e-5 | - |

So (a) the phi-S channel closes at exactly 4.347419 GeV: the "continuous B_s -> phi S" rate
at 4.4 GeV is **exactly zero**, not a finite endpoint value; (b) just below the cutoff the
inclusive formula gives 0.544 while the only open exclusive channel gives 1.2e-3 (factor
~450) — the inclusive formula is invalid near the endpoint (Winkler says so himself), and
the "cliff" is therefore an artefact of extending it to the phi threshold.  The channels
that stay open above 4.347 are eta' S (threshold m_Bs-m_eta' = 4.409 GeV), f0(980) S
(~4.377) and eta S (4.819), for which no B_s form factors are provided.

Impact on the result: B_s is 10.07 % of the B-meson rows.  Its share of the total signal at
the 25 fb crossing is 1.67 % (3.5 GeV), 1.00 % (4.0), 0.05 % (4.3) and 0 (>=4.4); at the
scan point 4.3446 GeV the B channel is only 0.43 % of the signal (H dominates), so removing
B_s shifts the long-lived crossing by +0.046 % (25 fb) / +0.004 % (300 fb).  The cliff is
invisible in the crossing curve.  The zero-mass recoil is a similar-size approximation:
at m_S = 3.5 GeV it makes the scalar momentum 1.542 GeV instead of 1.294 GeV for a ~1 GeV
recoil (16 %), i.e. <=0.3 % of the total via the 1.7 % B_s share.

Verdict PARTLY: the criticism is correct as a statement about model completeness, but it is
not worth a form-factor calculation now.  Fix (cheap): document the B_s treatment as a
flagged approximation with a +-(0.5-1.7 %) uncertainty on the total for 3.5-4.3 GeV; use the
exclusive phi-S rate of eq.(BKphi) below 4.347 (much smaller than the inclusive value above
~4 GeV) or, at minimum, stop trusting the inclusive formula above ~4.0 GeV; do not
extrapolate above 4.347 and state that eta/eta'/f0 S are not implemented.

---

## F13 — reference-BR validation model (reference_br_model.py, validate_widths.py)

(1) Isospin factors are real and the digitised Winkler curve is the isoscalar rate.
At m_S = 0.5 GeV the reference table's columns sum to 0.70330, so the missing fraction is
0.29670; BR(pi+pi-) = 0.59340, and 0.29670 = BR(pi+pi-)/2 **exactly** — the missing part is
the pi0pi0 channel, i.e. the table's pi+pi- column is the single charge state while our
widths() 'pipi' is the isoscalar total (1.5x).  Likewise Gamma_ours(KK)/Gamma_ref(K+K-) =
2.00-2.46 for the kaon channel.  The 1.5 factor also works directly on the digitised curve:
ours(pipi)/[1.5*BR_ref(pi+pi-)*Gamma_ref_tot] = 1.03 at 0.5 GeV (2.27e-8 vs 2.21e-8);
it deviates (0.30-0.57) only in the f0(980) region 1.0-1.6 GeV, where the digitised Winkler
curve and the SensCalc table differ in shape (the known ~1.5x hadronic discrepancy).
reference_br_model.py:38-41 maps pi+pi- -> 'pipi' and K+K- -> 'KK' without the 1.5/2 factors,
and s_decay_products then applies 2/3 and 1/2 (bc5_lib.py:276-283, 284-304), so the
validation model's visible pi-pi / K-K rate is 2/3 and 1/2 of the reference value; the
missing pieces are dumped into 'multi' (which has a different, higher charged acceptance),
so the errors partly cancel: current/corrected visible rate = 1.157 (1.0 GeV), 1.010 (1.3),
0.958 (1.6), 0.982 (1.99) with fixed kinematics (eta=3, A(pipi)=0.53-0.56, A(KK)=0.80-0.87,
A(multi)=0.65-0.91).  i.e. <=16 % at 1 GeV, <=4 % above 1.6 GeV.  Note also that the
reference table itself only sums to 0.70/0.64/0.47/0.42 at 0.5/1.0/1.5/1.99 GeV, so 36-58 %
of Gamma_tot is assigned to 'multi' in the validation model — a larger fudge than the
isospin factor.

Fix in reference_br_model.widths():
    out['pipi'] = 1.5 * br['pipi'] * g_tot    # isoscalar pi pi
    out['KK']   = 2.0 * br['KK']   * g_tot    # isoscalar K K
    rest = g_tot - sum(out.values()) ...      # only the genuinely missing channels
and state explicitly that the table is incomplete below 2 GeV.

(2) validate_widths.py:47-51 tests the wrong function.  The ss/cc checks call
sw.gamma_qq(m, MSB_S_2GEV, 2.0, M_K) and sw.gamma_qq(m, MSB_C_2GEV, 2.0, M_D), i.e. the
MS-bar masses run from 2 GeV to m_S, while the production widths() path
(scalar_widths.py:434-443) uses the fixed Winkler masses W_M_S=0.095, W_M_C=1.30 with
mu_ref = m_S (no running).  Ratio running/fixed: ss = 0.902 (2.2 GeV), 0.849 (2.5),
0.785 (3.0), 0.702 (4.0), 0.649 (5.0), 0.583 (7.0); cc = 0.664 (4.0), 0.614 (5.0),
0.552 (7.0).  So the ss/cc leg of the published validation tests a function up to 42 %
different from the production one; the old reference_comparison.json cannot be quoted as
validation of the production path.
Fix: compare against sw.widths(m)['ss'] / ['cc'] (or
sw.gamma_qq(m, sw.W_M_S, m, sw.M_K)) and keep the running-mass curve, if at all, only as a
separate "Boiarska/SensCalc model" comparison; regenerate reference_comparison.json.

---

## Summary of recommended order

1. F05 + F07 + F03 together (they all touch the S-decay / quartic generators): implement a
   correct uniform n-body sampler, propagate the quartic channel label and |M(q^2)|^2, and
   define the multi composition with isospin weights + eta decays (or the Pythia library).
   These change the low-mass (m_S <~ 2 GeV) long-lived edge by 10-30 %.
2. F06 (3-5 % in the Higgs-dominated region): write the particle status in the generator and
   select status 62 (last copy); use the row's invariant mass.
3. F13: apply the 1.5/2 isospin factors and fix the validate_widths ss/cc comparison before
   quoting the reference-BR scan as a validation.
4. F08: keep the B_s channel as a flagged approximation (+-0.5-1.7 % above 3.5 GeV); stop
   extending the inclusive formula to 4.347 GeV; no form-factor work needed now.
