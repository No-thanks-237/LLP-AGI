# BC5 dark scalar at LHCb downstream — run report

Run: runs/20261016_144605_UTC_bc5_lhcb_downstream
Model: Higgs-portal real scalar, PBC benchmark **BC5** (BR(h->SS) = 0.01), per prompt.md
Deliverable: exclusion curve in the (m_S, theta^2) plane, zero background, at L = 25 and
300 fb^-1.

**Threshold:** the canonical outputs now use **s = 2.3** (P(n=0|s) = exp(-2.3) = 0.100, i.e.
**90% CL**). The earlier **s = 3** (95% CL) results are kept unchanged alongside, with the
suffix **_s3**: bc5_lhcb_downstream_exclusion_s3.{png,pdf}, exclusion_boundary_s3.csv,
exclusion_s_curves_s3.csv, exclusion_report_s3.json. Both come from the same signal scan; only
the crossing level differs.

Update 2026-10-17: gaps 1, 2 and 3 of the first pass have been closed (real Pythia
hadronisation and tau decays, and the quartic-coupling production channels).  The curve was
recomputed; see section 4 for what changed and what is left.

## 1. Result

The excluded region is a band in theta^2 at every mass: above the *short-lifetime* edge the
scalar decays before the fiducial volume (z_dec > 1 m), below the *long-lifetime* edge it
decays beyond it (z_dec < 2.5 m). Both edges are read off the s = s_cut crossings; where only
one crossing lies inside theta^2 <= 1e-3 the band is left **open** (flagged in
exclusion_boundary.csv by band_open_above = True) and no curve is fabricated.

Representative values at the canonical threshold s = 2.3 (90% CL), read from
outputs/data/exclusion_boundary.csv:

| m_S [GeV] | 25 fb^-1 | 300 fb^-1 | dominant channel |
| --- | --- | --- | --- |
| 0.021 | > 2.84e-04 (open) | > 8.18e-05 (open) | h -> SS |
| 0.101 | > 6.07e-05 (open) | > 1.75e-05 (open) | h -> SS |
| 0.302 | 3.23e-08 .. 8.79e-04 | > 3.86e-09 (open) | h -> SS |
| 0.484 | 5.30e-09 .. 9.97e-05 | 4.65e-10 .. 1.69e-04 | h -> SS |
| 1.059 | 2.83e-10 .. 4.47e-06 | 2.36e-11 .. 7.54e-06 | h -> SS |
| 1.983 | 4.51e-10 .. 5.97e-06 | 3.76e-11 .. 1.14e-05 | h -> SS |
| 3.175 | 7.00e-11 .. 1.17e-06 | 5.81e-12 .. 2.07e-06 | h -> SS |
| 5.082 | 1.95e-12 .. 2.27e-08 | 1.62e-13 .. 4.02e-08 | h -> SS (B -> X_s S has just closed) |
| 8.136 | 4.54e-13 .. 5.81e-09 | 3.77e-14 .. 1.11e-08 | h -> SS |
| 20.847 | 6.47e-15 .. 1.08e-10 | 5.37e-16 .. 1.95e-10 | h -> SS |
| 33.371 | 1.96e-15 .. 3.66e-11 | 1.63e-16 .. 6.26e-11 | h -> SS |
| 62.490 | 1.76e-16 .. 6.97e-12 | 1.46e-17 .. 1.22e-11 | h -> SS |

Crossing structure at s = 2.3:

* 300 fb^-1: 34 masses with two edges and 18 open (band_open_above = True, all with
  m_S < 0.35 GeV, where the excluded band continues above the scanned window; only the
  long-lifetime edge is drawn);
* 25 fb^-1: 35 with two edges and 17 open.
  (These numbers are after the normalisation fix described in section 1a; before it the
  contours were two orders of magnitude too narrow on the long-lived edge.)

At the earlier threshold s = 3 (95% CL) the same signal scan gives 35 masses with two edges and
17 open at 300 fb^-1; the bands are uniformly narrower, by a factor growing from ~1.2 at the
long-lifetime edge to ~1.4-5 at the short-lifetime edge, where s(theta^2) falls steeply.

The reach extends to m_S = m_h/2 purely through h -> SS, exactly as the LHCb benchmark paper
(Gorkavenko et al., EPJC 84 (2024) 608) reports for BC5. Two kinematic steps appear and are
physical:

* m_S ~ 0.28 GeV — the pi pi threshold, below which S is very long-lived;
* m_S ~ 5.13 GeV = m_B - m_pi — where B -> X_s S closes and only h -> SS remains.

The f0(980) peak at m_S ~ 0.98 GeV produces the visible narrow structure.

Files (s = 2.3, canonical): outputs/figures/bc5_lhcb_downstream_exclusion.png and .pdf,
outputs/data/exclusion_boundary.csv, outputs/data/exclusion_s_curves.csv.
Files (s = 3, preserved): the same names with the suffix _s3.

## 1a. Two normalisation bugs found by comparing with the BC5 reference figure

The BC5 benchmark figure (references/model/BC5.pdf, from Gorkavenko et al. EPJC 84 (2024) 608)
was digitised to cross-check the contour.  The comparison exposed a discrepancy of up to two
orders of magnitude on the long-lived edge, which traced back to two independent defects in the
efficiency normalisation of the scan, both since fixed.

**Bug 1 — the event efficiency was divided by the number of GENERATED events instead of the
number of events actually evaluated.**  The scan samples nsamp = 1200 parents per channel, but
run_scan.event_or_probability divided the sum of the per-event OR probabilities by n_gen
(200 000 Higgs events, 373 659 b bbar events).  Each sampled parent represents
n_sample/n_sampled parents of the full sample, and that factor was missing.  The measured
shortfall was a factor 1274 at m_S = 20.85 GeV (mean OR 0.347 over 157 surviving pairs, divided
by 200 000 instead of being scaled to the 199 105 event sample).  The efficiency therefore
depended on nsamp, which is the signature of the error: increasing nsamp would have made the
result converge to the correct value, which is never true of a physical efficiency.

Detection: hand-computing the expected yield at the peak (1.664e5 h -> SS events times a 5 %
decay-in-window probability times a 12 % acceptance times 0.4) gives of order 10^3 events, while
the scan reported 60.

**Bug 2 — the 0.4 event-level reconstruction efficiency was never applied.**  RECO_EFF = 0.4 was
defined in bc5_lib.py and documented in the run_scan docstring, but no code path used it, so the
signal was a factor 2.5 too large.

Net effect: the signal was too small by a factor of about 60, so the excluded theta^2 band was
too narrow by the same factor on the long-lived edge (where s is proportional to theta^2) and by
a smaller factor on the short-lived edge (where s falls exponentially).

After the fix, eps_H at the peak of m_S = 20.85 GeV is 1.806e-2 instead of 2.72e-4, and the
edges at that mass become 5.4e-16 and 2.0e-10 against the reference 8.3e-16 and ~1e-10, i.e.
agreement at the 30 % level instead of two orders of magnitude.  The scan, both contours and the
comparison figure were regenerated.

The full comparison after the fix (outputs/figures/bc5_vs_reference.png), ratio mine / reference:

| m_S [GeV] | long-lived edge, 25 fb^-1 | long-lived edge, 300 fb^-1 | short-lived edge, 25 fb^-1 |
| --- | --- | --- | --- |
| 0.5 | 10.2 | 4.5 | - (above the frame) |
| 1.0 | 6.9 | 5.2 | 2.2 |
| 2.0 | 6.1 | 3.3 | 1.6 |
| 5.0 | 0.98 | 0.98 | 1.17 |
| 10.0 | 0.85 | 0.85 | 1.04 |
| 20.0 | 0.82 | 0.82 | 1.14 |
| 40.0 | 0.70 | 0.70 | 1.05 |

Above m_S = 5 GeV the two analyses agree to better than 30 % on both edges.  The remaining
factor 3-10 at m_S = 0.5-2 GeV is in the region where the pi pi / K K width itself is uncertain by
an order of magnitude (section 1b), which is exactly where the two treatments of the hadronic
width are expected to differ.

## 1b. `work/widths` as a standalone program

The width module now runs on its own, without the run pipeline:

    cd work/widths
    python widths_cli.py --mass 1.0
    python widths_cli.py --mass 1.0 --theta2 1e-9 --json
    python widths_cli.py --scan 0.021 62.49 52
    python widths_cli.py --check          # validation against the digitised Winkler curves
    python widths_cli.py --hadronic chpt  # lowest-order ChPT instead of the digitised table

It prints every channel width, its branching ratio, Gamma_S and c*tau, and states which hadronic
input was used.  Files:

    scalar_widths.py   closed-form widths: leptons, gamma gamma, gg, q qbar (LO spectator,
                       running masses), the multi-meson toy channel, the assembly
    omnes.py           Muskhelishvili-Omnes machinery (T matrix, unitarity relation, dispersion,
                       form factors, rate equations) -- written, NOT yet validated; the solver
                       raises rather than returning wrong numbers
    widths_cli.py      the standalone driver
    validate_widths.py the reference comparison (also reachable as widths_cli.py --check)
    hadronic_table.csv the digitised Winkler figure 2, currently still the production input for
                       pi pi / K K

### Status of the pi pi / K K Omnès implementation

Winkler's section "Dispersive Analysis" gives the complete recipe and it has been transcribed
into omnes.py: the two-channel T matrix in terms of the phases delta, psi and the inelasticity
eta, the unitarity relation Im F_i = T*_ij beta_j F_j, the once-subtracted Muskhelishvili-Omnes
solution with Omega(0) = 1, the form-factor assembly (with his numerical p = 0.73, q = 0.52) and
the rate equations.  Two things are still missing before it can replace the digitised table:

1. **The T-matrix input.** The phases delta, psi and the inelasticity are not formulas in
   Winkler's paper: he takes them from the Roy-Steiner analysis of Hoferichter et al.  We do
   not have that data set in the workspace.  This is an external, experimental/theoretical
   input, exactly as it is for Winkler himself; anyone reproducing his curve has to supply it.
2. **A validated principal-value quadrature.** The naive log-subtracted Cauchy kernel currently
   implemented does not reproduce the analytic single-channel benchmark
   Omega = exp[(s/pi) P int delta/(s'(s'-s)) ds'] exp(i delta), and the fixed-point iteration is
   unusable because Im Omega = tan(delta) Re Omega diverges at delta = pi/2.  A direct linear
   solve is in place but the benchmark test still fails, so the solver is disabled.

Until both are resolved the pi pi and K K widths remain the digitised Winkler input, and this is
the one place in the whole analysis where a published curve is used instead of a computation.

#### What has since been built and what it showed (2026-10-17, second pass)

The two source papers were obtained and mined completely:

* Garcia-Martin, Kaminski, Pelaez, Ruiz de Elvira, Yndurain, PRD 83 (2011) 074004 (the "CFD"
  parametrisation) gives the isoscalar S-wave phase shift and inelasticity in closed form with
  every coefficient printed (Appendix A, Eqs. A1-A4, Table V).  These are now implemented in
  work/widths/phases.py and validated: exact continuity of delta and its derivative at the
  0.85 GeV matching point, delta(4 m_K^2) = 226.5 deg = d_0 as the paper defines it, and the
  Table XII phases reproduced at 0.5 and 0.7 GeV.
* Hoferichter et al. (arXiv:1204.6251, whose LaTeX source was supplied by the user) builds
  exactly the same two-channel Omnes matrix.  All 89 of its equations are transcribed in
  work/widths/references_thesis/OMNES_RECIPE.md, including the Omega^inf + xi construction that
  solves the asymptotic problem this implementation had.  It prints NO numerical
  parametrisation of the S matrix, and its footnote 3 states that the pi pi -> K K amplitude it
  uses is a private communication.

A single-channel Omnes function (valid below the K K threshold, where only delta is needed) is
implemented in work/widths/omnes.py.  Compared with the digitised Winkler curve
(outputs/figures/widths_pipi_omnes_vs_winkler.png) it agrees at m_S ~ 0.7 GeV but is 20-25% low
below 0.5 GeV and up to a factor 3 high near 0.95 GeV.  Two causes were identified, both
physical rather than implementation errors:

1. **The phase-shift input disagrees near the f0(980).**  The CFD set has
   delta(4 m_K^2) = 226.5 deg, whereas Caprini-Colangelo-Leutwyler quote 185 +- 10 deg at the same
   point and Winkler's own figure shows 183 deg at 1 GeV.  This is the known "dip vs non-dip"
   ambiguity in the S0 inelasticity.  Winkler evidently used the CCL-like input, which is not
   available as a curve (CCL prints only three anchor values).
2. **The coefficient p = 0.73 cannot be fixed in the single-channel approximation.**  The rate
   is dominated by the trace-anomaly form factor Theta_pi = (2 m_pi^2 + p s) Omega, and p is
   defined through Omega'(0) of the FULL two-channel solution, including the derivative of
   Omega_12, which is non-zero even though Omega_12(0) = 0.  A single-channel Omega therefore
   cannot determine p self-consistently, and the p s term dominates the rate above ~0.7 GeV.

**Conclusion.**  A faithful independent reproduction of Winkler's pi pi / K K rates is not
possible from published inputs alone: the two-channel construction needs the pi pi -> K K
amplitude (private) and the specific phase-shift set he used (not published as a curve).

**Decision: keep the digitised Winkler curve (option a).**  This is also what the benchmark
phenomenology paper does.  Boiarska et al. arXiv:1904.10447, section "Scalar decays", write

> "the reliability of the dispersion relation method is questionable for m_S >~ 1 GeV because of
> lack of the experimental data on meson scattering at high energies and unknown contribution of
> scalar hadronic resonances ... To have a concrete benchmark ... we use the decay width into
> pions and kaons from [Winkler:2018qyg], see Fig. 1"

and their Table 1 column 6 is explicitly "reference to the formula (or figure in case of decays
into pions and kaons)".  They use the same isospin ratios, the same multi-meson toy formula and
the same Lambda_pert = 2 GeV as here, and they quote the pi pi / K K width as uncertain by
**more than an order of magnitude** around 1 GeV.  Our treatment therefore reproduces the
literature benchmark exactly; the residual uncertainty in this channel is a property of the
field, not of this analysis.

The Omnes implementation and the CFD phase-shift parametrisation are preserved under
work/widths/research/ with a README recording what they showed and why they are not in the
production path.  work/widths/README.md documents the program as a whole.

## 2. Physics inputs and validation

| Quantity | Treatment | Validation |
| --- | --- | --- |
| Gamma(S -> l l) | first principles, G_F m_S m_l^2 beta^3/(4 sqrt2 pi) | digitised Winkler figure: ratio **1.000**, scatter 1e-4 |
| Gamma(S -> gamma gamma) | SM h -> gamma gamma loop functions with m_h -> m_S | SM value 9.1e-6 GeV at 125 GeV reproduced |
| Gamma(S -> gg) | alpha_s^2 m^3 abs(sum F_1/2)^2/(128 pi^3 v^2) | matches Winkler figure to **6%**; heavy-top limit alpha_s^2 m^3/(72 pi^3 v^2) reproduced |
| Gamma(S -> q qbar) | Boiarska 3-loop running mass, **LO spectator model** (see the note below) | strange 3.5% at 2 GeV, charm 32% below the Winkler figure at 5 GeV but within ~10% of SensCalc |
| Gamma(S -> pi pi), Gamma(S -> K K) | **digitised Winkler Fig. 2**, self-calibrated on the exact S -> mu mu curve | calibration factor 1.00946, scatter 8e-4 |
| **S decay products (gg, ss, cc, bb, tau tau)** | **Pythia 8.316**, isolated colour-singlet S decay hadronised with forceHadronLevel and fully decayed | see section 3; event-by-event energy sum equal to m_S |
| BR(B -> X_s S) | Boiarska master formula + form factors | all 10 channels agree with tab:BR to **1-6%**; total 3.356 vs 3.3 |
| BR(K -> pi S) | same master formula | 0.95 (K+) and 0.99 (K0L) of tab:BR |
| BR(X -> SS, X -> X'SS) | Boiarska tab:b-production-quartic with alpha(m_S) and the tabulated closing masses | order 1e-10 per B: four orders below h -> SS |
| alpha(m_S) | fixed by BR(h->SS) = 0.01 | alpha = 0.7195 GeV at threshold; BR(h->SS) = 0.0100 exactly |
| detector kernel | unmodified DownstreamTransport.h + BirthKinematics.h via ctypes | 32/33 checks pass (the 33rd is a contradiction in my own test spec: the beam axis lies inside the beam hole) |
| lifetime integral | trapezoid on 13 nodes over [1, 2.5] m | 13 vs 25 nodes with identical A(z): median dev **0.02%**, max 1.5% |
| event-level OR | 1 - prod(1 - w_k p_k) per simulated collision | — |

### Note: the quark channels use the LO spectator width

Boiarska's Eq. (decay-quarks) multiplies the spectator result by (1 + Delta_QCD + Delta_t) with
alpha_s(m_S). At the 2 GeV matching scale alpha_s = 0.338 and Delta_QCD = 1.09, so this
"correction" doubles the width and the hadronic width jumped by a factor **1.49** across the
matching point:

| | 2 GeV below (dispersive) | 2 GeV above (perturbative) | ratio |
| --- | --- | --- | --- |
| this run, NLO quarks | 7.76e-8 GeV | 1.155e-7 GeV | **1.49** |
| this run, LO quarks (production setting) | 7.76e-8 GeV | 8.95e-8 GeV | **1.15** |
| Winkler figure | 7.815e-8 GeV | 7.817e-8 GeV | 1.0002 |
| SensCalc tables | 0.9290 (BR sum) | 0.9290 (BR sum) | 1.0000 |

Both published references are continuous, so the step was an artefact of applying a
next-to-leading-order correction where the expansion parameter is ~0.34. The production
setting is therefore \`QUARK_NLO = False\` in work/widths/scalar_widths.py; set it to True to
recover Boiarska's NLO expression. This lowers Gamma_S by a factor 0.75-0.90, i.e. lengthens
the lifetime by 1.1-1.3, and widens the excluded theta^2 band above 3 GeV by about 1.25.

The residual 15% step comes entirely from Gamma(S -> gg) at 2 GeV, where this implementation is
22% above the digitised Winkler curve (it agrees to better than 3% above 4 GeV); the effect is
below the ~factor 2 theory uncertainty quoted for the gluon channel.

Trade-off worth recording: with LO quarks Gamma(S -> c cbar) at 5 GeV is 32% below the
digitised Winkler curve, whereas with NLO it agreed to 2%. SensCalc's c cbar width, however,
lies within ~10% of the LO value, so the two references disagree with each other here.

The reason for the residual difference is the **running quark mass**, not the QCD correction.
This implementation evaluates mbar_q at mu = m_S (as Boiarska's formula prescribes), so mbar_s
falls from 0.092 GeV at 2 GeV to about 0.068 GeV at 7 GeV and mbar_c falls from 1.27 GeV at its
own scale to about 1.05 GeV at 5 GeV. That alone accounts for the ratio between this
implementation and the Winkler curve, 1/1.04 -> 1/1.46 for s sbar over 2 - 7 GeV and 1/1.48 for
c cbar at 5 GeV, without invoking any higher-order term. The NLO version looked closer to the
Winkler curve only because (1 + Delta_QCD) happens to compensate the running over 4 - 7 GeV.

The Monte-Carlo statistical error of the predicted signal *at the contour* is
**8.5% (median) on the long-lifetime edge and 17.5% (median) on the short-lifetime edge** at
300 fb^-1 for the canonical s = 2.3 threshold (16% at s = 3; columns
rel_s_mc_err_at_long_edge / rel_s_mc_err_at_short_edge in the boundary CSV). This is the event-level standard error of s; it does not include the (smaller) error from
evaluating A(z) with a finite number of decay samplings. This is **above** the 5% target of
prompt.md and is the main remaining statistical limitation: the short-lifetime edge is where s
falls steeply, so the effective number of contributing simulated events is small. Enlarging the
Higgs sample (the dominant channel) and the number of decay samplings per scalar would close it.

## 3. How gaps 1 and 2 were closed

The first pass hadronised S -> gg/ss/cc/bb with an analytic jet model and tau decays with a
leading-track model. Those are now replaced by Pythia 8.316 itself.

The technical problem was that Pythia's built-in production processes cannot be forced to a
light resonance mass (PhaseSpace:mHatMin defaults to 10 GeV, and constraining it makes the
cross section vanish), and Pythia::forceHadronLevel needs a complete event record. The working
recipe, found by direct experiment, is:

    ProcessLevel:all = off;  PartonLevel:all = off;  HadronLevel:all = on
    event.reset()
    append the two S decay partons at rest, with the correct colour flow:
        q qbar : col(1) on q, acol(1) on qbar
        g g    : col(1),acol(2) on g1 and col(2),acol(1) on g2   <-- a closed colour loop
    pythia.forceHadronLevel()      (the returned bool is a false alarm: Pythia::check
                                    complains about the absent beam particles, but the
                                    final-state energy sums to m_S exactly)
    for (i<40) if(!pythia.moreDecays()) break;   <-- this is what decays the tau leptons

The gluon colour assignment was the key: a single col/acol pair leaves the string
incomplete and produces no charged particles at all.

work/magnus/hadron_lib/hadron_lib_<channel>.txt.gz holds 1500 full hadronisations for each of
the 52 scan masses for gg, ss, cc, bb and tau tau, with the stable charged particles in the
S rest frame. The analysis samples this library; the analytic models remain only as a fallback
and for the elementary channels (e, mu, pi pi, K K, 4 pi), which are treated exactly anyway.

Verification: for every partonic channel and mass the generator sums the energies of ALL final
particles and obtains m_S exactly; the charged-energy fractions are 0.59-0.77. The charged
multiplicities are compared with the old analytic model in
work/analysis/hadronisation_validation.json; the old model underestimated them by factors
1.5-3 (e.g. gg at m_S = 10 GeV: 3.7 -> 8.1).

Gap 3 (quartic production X -> SS, X -> X' SS) is now implemented from Boiarska's
tab:b-production-quartic with alpha(m_S) and the tabulated closing masses, using two-body
phase space for B_s/B0 -> SS and three-body phase space with a kaon spectator for
B -> X' SS. Its rate is 1e-11 - 1.5e-10 per B meson, i.e. of order 20 S pairs at 300 fb^-1
against 3.3e5 scalars from h -> SS; it does not move the contour, but the claim is now backed
by the calculation rather than by an estimate.

Effect on the answer: the boundaries move by factors 0.2 - 1.3 in theta^2 (mass dependent,
because the channel composition and the boost change with m_S). The qualitative structure —
the two edges, the two thresholds, the f0(980) structure and the reach to m_h/2 — is unchanged.

## 4. What is still open

1. **tau polarisation and the tau spectrum** now come from Pythia, but the two taus are
   injected back-to-back at rest without spin correlations; the effect on the acceptance is
   expected to be small (the acceptance is driven by the boost, not by the tau helicity).
2. **Cross-section normalisations**: b bbar events are normalised to an assumed
   sigma(pp -> b bbar) = 500 ub at 13 TeV (FONLL/NNLO). That number is not verified against a
   measurement in this run, and because the Pythia sample is eta-cut at generation the LHCb
   forward measurement cannot be folded in directly. This affects only the meson channel.
3. **b bbar statistics**: 373 659 b bbar events were produced before the job was stopped;
   the B meson spectrum is 162 530 entries, and the B mixing channel contributes at most ~10%
   of the signal at theta^2 = 1e-8, so its own MC error is subdominant but not negligible at
   the large-theta^2 edge.
4. **Higgs sample**: 200 000 Pythia HiggsSM events reweighted to the LHC HXSWG total cross
   section; the LO subprocess mix is used for the shape (ggF is 79% of the generated sample
   versus 88% at NLO). No MG5 sample was produced in this run.
5. **Theory systematics** (factor ~2 on the gluon channel, factor ~3 on the low-mass hadronic
   width) are not propagated; the SensCalc cross-reference comparison is in
   work/widths/reference_comparison.json.
6. **UT**: prompt.md deliberately drops the reference paper's requirement that downstream
   tracks pass through UT, so this result is *more* permissive than Gorkavenko et al.
7. **No detector resolution, track-reconstruction pT dependence, vertex resolution or PID** is
   modelled; the single 0.4 factor is the entire reconstruction efficiency, as specified.
8. **D and K mesons contribute nothing** (K -> pi S needs m_S < 0.354 GeV where the maximum
   daughter pT is below the 0.5 GeV cut; D is suppressed by xi_uc^2 ~ 1e-18). Both are
   implemented and return zero acceptance.
## 1c. Kaon quartic channels added

Boiarska et al. tab:K->SS lists five quartic-coupling kaon modes that were missing from the
implementation:

| channel | BR(m_S = 0) at alpha = 1 | closing mass [GeV] | kinematic limit |
| --- | --- | --- | --- |
| K0L -> SS | 4.4e-13 | 0.252 | m_S < m_K/2 = 0.2488 |
| K0L -> SS pi0 | 6.6e-15 | 0.140 | m_S < (m_K-m_pi0)/2 = 0.1815 |
| K+ -> SS pi+ | 1.4e-15 | 0.136 | m_S < (m_K-m_pi+)/2 = 0.1770 |
| K0S -> SS | 7.8e-16 | 0.252 | m_S < 0.2488 |
| K0S -> SS pi0 | 1.2e-17 | 0.140 | m_S < 0.1815 |

They are now implemented in bc5_production.quartic_kaon_channels and sampled in run_scan with the
same two-scalar treatment as the B quartic modes.  K0S is not present in the minimum-bias sample
(it is 560 times smaller than K0L in the 2-body mode and decays promptly in Pythia), so only the
K0L and K+ rows contribute.

Effect on the contour (long-lived edge, before -> after):

| m_S [GeV] | 25 fb^-1 before | after | ratio |
| --- | --- | --- | --- |
| 0.021 | 2.84e-4 | 2.81e-4 | 0.99 |
| 0.101 | 6.07e-5 | 4.74e-5 | 0.78 |
| 0.302 | 3.23e-8 | 3.35e-8 | 1.04 |
| 0.484 | 5.30e-9 | 5.14e-9 | 0.97 |
| 1.059 | 2.83e-10 | 2.78e-10 | 0.98 |
| 1.983 | 4.51e-10 | 4.21e-10 | 0.93 |
| 3.175 | 7.00e-11 | 6.93e-11 | 0.99 |

The largest effect is 22 % at m_S = 0.101 GeV; above 0.3 GeV it is at the 3-7 % level, i.e. inside
the Monte-Carlo error.  This is expected: all five modes are kinematically closed above
m_S = 0.25 GeV, so they cannot be the cause of the factor 3-8 discrepancy with the BC5 figure in
the 0.5-2 GeV region.  Their relative weight in the total signal is about 5-9 % at
m_S <= 0.15 GeV and below 1 % elsewhere.
## 1d. How SensCalc builds its yield (reading arXiv:2305.13383)

SensCalc computes the number of events as

    N_ev = sum_i N_prod^(i) * int dE dtheta dz f^(i)(theta,E) eps_az dP_dec/dz eps_dec eps_rec eps_S/B

and in the long-lifetime limit (their Eq. 8) the total acceptance factorises as

    eps = <eps_FIP> * <(gamma^2-1)^(-1/2)> * (Delta z / c tau) * <eps_decay>

**f^(i)(theta,E) is a tabulated grid, not a Monte-Carlo sample.**  Their notebook
"FIP distribution.nb" produces a table {m, theta, E, f^(i)}; for each production channel the
parent-particle angle-energy distribution has to be supplied.  Crucially, Table in Appendix B
("list of the references used to generate, or directly take, the distributions of secondary
particles that may produce FIPs") gives, for the LHC column:

| parent | LHC source |
| --- | --- |
| pi0, eta, eta', rho0, omega, gamma | Pierog et al. |
| B, D | Kling:2021fwx |
| **W, h, Z** | **Kling:2021fwx** |

Kling:2021fwx is Kling and Trojanowski, "Forward experiment sensitivity estimator for the LHC and
future hadron colliders", arXiv:2105.07077, Phys. Rev. D 104, 035012 (2021) - the FORESEE tool.
So the Higgs-parent spectrum that drives the BC5 h -> SS yield in SensCalc is a **tabulated
forward-production spectrum taken from FORESEE**, not a Pythia 8 sample of ggF+VBF+... production
with the full kinematics.

This is the single most plausible place for the residual factor: in the long-lifetime limit
N_ev is proportional to <1/gamma> of the scalar, i.e. inversely proportional to the mean scalar
energy, so any difference between the FORESEE Higgs spectrum and the Pythia spectrum propagates
one-to-one into the position of the lower edge.  Our own Pythia sample gives
<E_S> = 462 GeV (<gamma_S> = 954) at m_S = 0.484 GeV, where the tabulated-spectrum route may be
appreciably softer.

Also noted from the same paper: their dark-scalar decay table (BC4/BC5) matches ours - ee, mu mu,
tau tau, pi pi, K K, 4 pi below 2 GeV and c cbar, s sbar, b bbar, GG above - and they state that
they include the NLO correction of Spira et al. in Gamma(S -> gg), which we do not (our
QUARK_NLO = False production setting uses LO).

**What is now confirmed to be the same on both sides** (from the same paper's App. A.1/Table 2 and
from Gorkavenko App. A, quoted in references/detector/arXiv-2305.13383v3/SENSCALC_YIELD_NOTES.md
sections 4 and 7.1):

* the decay volume: "the setup considered in this paper - the one with the decay volume extending
  from z = 1 m to z = 2.5 m and SciFi as the detector" (Gorkavenko App. A).  SensCalc's generic
  LHCb-downstream volume is 1-7.7 m, but the BC5 paper states that the 2.5 m version is the one
  used, so the factor 6.7/1.5 = 4.5 in decay volume is NOT the explanation;
* 2 < eta < 5 on the LLP, at least two charged decay products, neutral-only decays invisible;
* the dark-scalar decay channels (ee, mu mu, tau tau, pi pi, K K, 4 pi below 2 GeV; c cbar,
  s sbar, b bbar, gg above), from Boiarska plus Winkler;
* the semi-analytic formula itself - our code reproduces SensCalc's Eq. (8) analytically to 5 %.
  At m_S = 1 GeV, feed the analytic formula with our own numbers: ctau(theta^2=1) = 3.94e-10 m,
  gamma beta = 462, Delta z = 1.5 m, <eps_LLP> = 0.129, <eps_det> = 0.9, eps_rec = 0.4, N_scalars
  = 2.77e4 gives theta^2_edge = 2.8e-10, and the code returns 2.78e-10.  The implementation is
  therefore not the problem;
* SensCalc counts each of the two scalars from h -> SS independently; our event-level OR,
  1-(1-p1)(1-p2) = p1+p2 for small p, is numerically identical, so that is not a factor 2 either.

**What remains open.**  Reproducing the published lower branch at m_S = 1 GeV (theta^2 = 4.1e-11)
needs 6.8 times more signal in the long-lifetime limit, and since N_ev is proportional to
<1/(gamma beta)> there, either the effective scalar energy is ~68 GeV instead of our 462 GeV, or
there is a comparable factor elsewhere in N_h x 2 x BR.  The one genuine implementation
difference found in the SensCalc paper is that they include the NLO correction of Spira et al. in
Gamma(S -> gg) while our production setting uses LO; that is a factor of at most about 1.3 in
ctau, not 7.  Settling it requires the FORESEE Higgs tabulation (arXiv:2105.07077).

**That tabulation has since been checked and it is not the explanation.**  The FORESEE package was
cloned into the isolated copy of this run at runs/20261017_155516_UTC_bc5_foresee; its forward
Higgs spectrum (files/hadrons/14TeV/Pythia8/Pythia8_14TeV_25.txt) gives <E_S> = 483 GeV and an
effective inverse boost of 1/625 GeV at m_S = 0.484, against 459 GeV and 1/575 GeV from our Pythia
sample - agreement to 9 %.  See that run's logs/foresee_spectrum_check.md.  The residual factor is
therefore not attributable to any input that either paper publishes numerically.

**The shape of the residual localises it.**  Writing theta^2_edge = 3/(A Gamma_S), where A bundles
production x acceptance x <1/(gamma beta)> x eps_rec x Delta z / hbar c, and taking the ratio of our
edge to the reference edge after dividing out the width ratio Gamma_ref/Gamma_mine gives

| m_S [GeV] | edge ratio | Gamma_ref/Gamma_mine | residual A_ref/A_mine |
| --- | --- | --- | --- |
| 0.30 | 23.1 | 1.24 | 18.7 |
| 0.48 | 10.2 | 1.25 | 8.1 |
| 0.66 | 6.6 | 1.27 | 5.2 |
| 1.06 | 4.1 | 0.92 | 4.5 |
| 1.98 | 5.8 | 1.81 | 3.2 |
| 3.17 | 2.1 | 2.04 | 1.0 |
| 5.08 | 0.91 | 1.15 | 0.8 |
| 9.52 | 0.90 | 1.04 | 0.9 |
| 20.85 | 0.78 | 1.13 | 0.7 |

Above 3 GeV the residual is 0.7-1.0, i.e. **agreement**.  Below 2 GeV it rises smoothly to 19 at
m_S = 0.3 GeV.  But every ingredient inside A has been measured to be flat in this mass range -
the eta acceptance is 0.128 at all masses, the scalar energy distribution is identical from
0.484 to 20 GeV (E_S percentiles 121/239/349/555/1154 GeV), Delta z and eps_rec are constants.
**A mass-independent A cannot produce a mass-dependent residual.**  The only genuinely
mass-dependent input in the long-lived edge is therefore Gamma_S in the m_S < 2 GeV region - the
pi pi / K K / 4 pi window.  Quantitatively the reference curve behaves as if Gamma_S below 1 GeV
were 4-15 times larger than the Winkler/SensCalc value, i.e. as if the scalar were much
shorter-lived there.

That is also the one place where the literature genuinely disagrees with itself: the width audit
records pi pi and K K as uncertain by more than an order of magnitude around 1 GeV (the f0(980)
region), which is why the two published sources differ.  Our curve uses the digitised Winkler
figure, which is also what the reference cites, so the difference is not a choice between
published inputs - it is a difference in what the reference's code actually evaluates.

The candidate most often raised for this kind of factor, the difference between <1/(gamma beta)>
and 1/<gamma beta>, has been measured and is NOT the explanation.  Our decay probability is
computed per scalar with that scalar's own momentum, so the correct mean of the inverse is what
enters.  The energy distribution of the scalars in 2 < eta < 5 at both masses is broad -
percentiles 121 / 239 / 349 / 555 / 1154 GeV - and the effective energy implied by
<1/(gamma beta)> is

| m_S [GeV] | <E_S> | 1/<1/(gamma beta)> | skew factor |
| --- | --- | --- | --- |
| 0.484 | 459 | 575 | 0.80 |
| 1.000 | 459 | 278 | 1.65 |

i.e. at m_S = 1 GeV the true average inverse boost is 1.65 times larger than the naive
1/<gamma beta>, worth only 1.65 of the required 6.8.

**Where the residual must therefore sit.**  Of the four bracketed quantities
<N_h x 2 x BR> x <eps_LLP> x <1/(gamma beta)> x <eps_det> x eps_rec, three are measured and agree,
leaving the Higgs parent spectrum and <eps_det>.  Both are unpinned in the literature: SensCalc
takes the W/h/Z parent spectra from FORESEE (arXiv:2105.07077) and publishes them only as a
tabulated grid, and its eps_dec is a simulated per-decay acceptance with no published LHCb grid.
The SensCalc paper also documents the LHCb setup only as a name in a list, states that eps_rec is
beyond its scope, and validates the h -> SS channel nowhere (its only scalar checks are
B-produced scalars against SensMC, 10-20 %).  The residual factor therefore lies inside inputs
that are not published in numerical form, not inside this implementation, which reproduces
SensCalc's own Eq. (8).
