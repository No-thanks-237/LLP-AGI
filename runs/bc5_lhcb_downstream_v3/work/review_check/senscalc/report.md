# SensCalc audit (task-4): what the reference actually does for the BC5 scalar curve

Author: teammate senscalc-audit. Date: 2026-10-25. Write scope: this directory only.
Nothing outside runs/bc5_lhcb_downstream_v3/work/review_check/senscalc/ was modified.

## 0. What was read and how

Local copy audited: E:/AAAwork/LLP/code/SensCalc-main/ (git snapshot of the public
main branch, downloaded 2025; notebooks are plain text). No Wolfram kernel was
launched. Binary .m/.mx tables were not decoded; where a binary table matters I used
the text export of the *same 2023 tables* that is archived in
runs/bc5_lhcb_downstream_v3/work/physics_review_20261025/ (that archive also carries the
2023 git tree JSON, sha 6664a923869c6b21a3aba3d4e8429ec3ae9d4edd, i.e. the version whose
date matches BC5.pdf).

Notebook text was extracted with this directory's extract_nb3.py (a WL-box parser that
renders RowBox/SuperscriptBox/FractionBox back to input form, including 10^-x, which is
otherwise easy to misread). All quoted expressions below are verbatim from the extracted
code; the file:line numbers are lines of the *raw* .nb file.

Numbers were re-derived by running code (scripts listed in section 10), never taken from
the review text.

## 1. Verdict table

| # | Question | Answer | Evidence |
|---|---|---|---|
| 1a | B -> X_s S production: inclusive spectator formula or exclusive sum? | **Exclusive sum by default** (Boiarska tab, BtoXsStotal.dat). The inclusive spectator table exists as a one-line alternative and is NOT selected. No 4.7 GeV B->K S switch in either option. | 3. Scalar sensitivity.nb:506, 529-546 |
| 1b | K -> pi S | Tabulated Boiarska-type tables KplustoPiS.dat / K0LtoPiS.dat, normalised per parent lifetime; only defined for facility SPS, so it never enters the LHCb-downstream scalar curve. | 3. Scalar sensitivity.nb:679-690; prod-pheno-Scalar.nb:66-68 |
| 1c | quartic / double-scalar channels | Tabulated BR x phase space from Boiarska (BplustoXsSStotal.dat, BstoSS.dat) + an analytic h -> SS rate. 3-body |M|^2 is approximated by B -> K S S with a fake pole mass sqrt(38) GeV. | 3. Scalar sensitivity.nb:562-592, 641-654; prod-pheno-Scalar.nb:341,354,368-395 |
| 2 | Higgs production | Mother spectrum = spectra/SM particles/DoubleDistr_LHC_h.txt (theta 1e-5..pi, E 125.2..7366 GeV, no header). sigma(pp->h)=55 pb at LHC (552 at FCC-hh), a single number, no per-process split, no Pythia/MG provenance in the file or manual. b_h = BR(h->SS) = 2e-2 alpha^2 sqrt(1-4m^2/125^2), evaluated at the user's BR (0.01 for the shipped curve); the *tabulated* quartic tables use 0.15 and are rescaled by b_h/0.15. | experiments.nb:222; 3. Scalar sensitivity.nb:426-447, 641-654, 4223, 6793-6794, 7487 |
| 3 | epsilon_rec = 0.4 and 0.87 | 0.4 is present and applied exactly once, as a flat factor in NevInt, only for experiments whose name contains "LHCb-downstream". **0.87 is not applied anywhere**; no S/B retention factor of any form exists in the scalar path. | 3. Scalar sensitivity.nb:6630-6645, 6793-6794; grep of the whole repo finds no 0.87 outside plot grids |
| 4 | widths / BR table | BRs from phenomenology/Scalar/decay widths/BrRatioScalar.mx (4 descriptions; default = "1809.01876", i.e. Winkler). Lifetime from ctauS.mx, used as ctau(m,theta2)=ctau_table(m)/theta2 with log-log interpolation. Implied Gamma_tot at 5.08/13.02/33.37/53.42/62.49 GeV = 9.09e-6 / 1.425e-4 / 8.82e-4 / 1.390e-3 / 1.606e-3 GeV (from the 2023 text export of the same table). | llp-decays.nb:175, 203, 236, 13341; 3. Scalar sensitivity.nb:1505-1518 |
| 5 | acceptance | Acceptances/LHCb-downstream-D/Acceptance_LHCb-downstream-D_for_Scalar_decayDescription_Default.m; applied as a *differential* A(z) (azimuthal fraction x decay-product acceptance) inside the z integral, log-interpolated; geometry: 1.0 < z < 2.5 m, 2 < eta < 5 (annular cylinder), detector plane at 9.4 m, dipole 3-7 m with 4 T m; photons/K_L/n not detectable. | experiments.nb:16744-17000; 3. Scalar sensitivity.nb:1518, 1595-1630, 3400-3430, 4598-4640 |
| 6 | luminosity / N_pot / criterion | Shipped curve: Npot = 1.75e15 (= 25 fb^-1 x 70 mb), NevMin = 2.3, eps = 0.4. The shipped Nevents tables are at Npot = 2.16e16 (= 300 fb^-1 x 72 mb) and are rescaled by Npot/NpotDefault. I closed the shipped json with these numbers to 0.1%. | 3. Scalar sensitivity.nb:6621, 7487-7488, 7625; json file name; section 6 below |
| 7 | provenance | The local copy is NOT the 2023 version that produced BC5. It is a 2025 rewrite (files renamed, notebooks re-saved with Mathematica 14.2, manual dated 2025-04-09, Bremsstrahlung "AP-*" variants from 2409.11096 in the shipped curves). However all six production-probability .dat tables are byte-identical to the 2023 commit, and the h spectrum is the same table. | section 7 below |

## 2. Q1 — production prescription (the central question)

### 2.1 The switch that selects inclusive vs exclusive

3. Scalar sensitivity.nb:506 hard-codes

    SelectedBrChoice="Exclusive";

(the dialog that used to offer {Exclusive, Inclusive} is commented out at :510-515). The B
mixing weight is (:529-546, verbatim):

    PmotherToLLP[mLLP_,theta2_,BrhToSS_,"B-mixing",Facility]=
        theta2*(fbtob0[Facility]+fbtobplus[Facility])*
        intProd[If[SelectedBrChoice=="Exclusive","BtoXsStotal.dat","BtoXsSinclusive.dat"],
                mSM["Bcharged"],mSM["PiCharged"]];

with

    intProd[filename_,m1_,m2_]:=If[mLLP<m1-m2,Evaluate[Interpolation[Import[
        FileNameJoin[{directory["Production weights"],filename}],"Table"],
        InterpolationOrder->1][mLLP]],0.]

The dialog text (:510-515) labels the two options
"Sum over all exclusive decays B -> X_s/d + S" (Exclusive) and
"All production is approximated by the decays b -> s + S" (Inclusive), both citing 1904.10447.

**=> The reference curve, run with the shipped code, uses the exclusive sum.** This is the
opposite of our v3 production code, which uses Winkler's inclusive spectator formula
(bc5_production.py:113-199).

### 2.2 What the two tables actually contain (verified numerically)

Both are theta^2-normalised BRs on a 0->5.5 GeV grid (551 rows, step 0.01):
phenomenology/Scalar/Production probabilities/BtoXsStotal.dat and BtoXsSinclusive.dat.

* Exclusive table: BR/theta^2(0) = 3.30959; last non-zero mass 5.13 GeV;
  a sharp drop 0.0948 (4.78 GeV) -> 0.0133 (4.79 GeV) = closing of B -> K S at
  m_B - m_K = 4.786 GeV; residual 0.0133 -> 0.0052 between 4.8 and 5.1 GeV = B -> pi S
  closing at m_B - m_pi = 5.139 GeV. So BtoXsStotal is a *sum of exclusive two-body
  channels*, as advertised.
* "Inclusive" table: BR/theta^2(0) = 6.6291; it falls 0.0514 (4.0 GeV) -> 6.5e-3
  (4.10 GeV) -> 0 (4.12 GeV). The endpoint is the *parton-level* b -> s + S endpoint
  m_b - m_s = 4.18 - 0.0924 = 4.09 GeV. So this table is the pure spectator formula
  (Winkler eq:inclusive) with **no** B -> K S switch at 4.7 GeV: it simply ends at 4.1 GeV.
  (In the 2023 notebook this branch carries a commented-out
  BrBplusToXsSinclusive[4.7/5.3*mS] "Max[...]" patch, reference_2023_scalars.nb:349/369 -
  evidence that the 4.7 GeV switch idea was around, but it was never enabled.)
* Ratio inclusive/exclusive: 2.003 (m=0), 1.95 (0.5), 1.80 (1.0), 1.26 (2.0), 0.61 (3.0),
  0.057 (4.0) - i.e. the two prescriptions differ by O(2) exactly in the region where the
  B channel matters.

Cross-check against our model (run of bc5_production.py):
* our inclusive B+ : 5.398 at m=0, B0: 5.006  (vs SensCalc-inclusive 6.629: ours 23% lower;
  vs SensCalc-exclusive 3.310: ours 1.63x higher);
* our B -> K S at 4.70/4.75: 0.3093/0.2084 vs SensCalc-exclusive 0.2824/0.1961 (10%/6%
  agreement) - so the B -> K S piece is essentially the same formula in both codes;
* our inclusive formula has no parton-level endpoint: it gives 0.979 at m_S=4.0 where the
  SensCalc inclusive table gives 0.051. That is a real difference in our eq:inclusive
  implementation (it uses (m_B^2-m_S^2)^2 with no b-quark phase-space factor), worth
  flagging to the production/widths audit, but it is irrelevant to the reference: the
  reference does not use that table.

### 2.3 K -> pi S

3. Scalar sensitivity.nb:679-690:

    PmotherToLLP[mLLP_,theta2_,BrhToSS_,"K-mixing","SPS"]=
        theta2*(fstokplus["SPS"]*intProd["KplustoPiS.dat",mSM["Kcharged"],mSM["PiCharged"]]+
                fstok0l["SPS"]*intProd["K0LtoPiS.dat",mSM["KL"],mSM["Pi0"]]);

Two comments matter:
* the channel is only defined for Facility "SPS"; the LHC/FCC-hh/facility loop (:474-668)
  never defines K-mixing for LHC, and the shipped LHCb-downstream sensitivity json lists
  7 production modes without "K-mixing". So K -> pi S does not enter the BC5 curve at all.
* the tables themselves are Boiarska-type: at m=0 K+ = 1.4261876687078384e-3,
  K0L = 6.198330036596367e-3 (ratio 4.346, i.e. the K_L/K+ lifetime ratio 4.13 x matrix
  element). Our independent formula gives 1.610e-3 / 6.955e-3 (ratio 4.320): same
  prescription, ours 12-13% higher.

### 2.4 Quartic (double-scalar) channels

3. Scalar sensitivity.nb:562-592:

    PmotherToLLP[...,"B-quartic",Facility]=2*alphavalBrhToSS[mLLP,BrhToSS]^2*
        (fbtob0[Facility]+fbtobplus[Facility])*intProd["BplustoXsSStotal.dat",m_B^2,m_pi^2];
    PmotherToLLP[...,"Bs-quartic",Facility]=2*alphavalBrhToSS[mLLP,BrhToSS]^2*
        fbtobs[Facility]*intProd["BstoSS.dat",m_Bs^2,0.];

(the factor-2 in the comment at :546-552 is "due to the production of scalars in pairs").
The tabulated "BR(0)" files are Boiarska's Table 4 numbers; the 3-body matrix element used
by the EventCalc sampler is the B -> K S S approximation with a fake pole
(prod-pheno-Scalar.nb:341-395):

    MBpole=Sqrt[38.];
    (*Approximating the matrix element of the 3-body B decays by the matrix element of
      the process B->KSS*)
    Msquared3Body[LLPname,"B-quartic",mLLP_,E1_,E3_]=...=
        (1/(1-ScalarProduct[k-k1,k-k1]/MBpole^2)^2)//ExpandScalarProduct ...

(1/(1-q^2/38))^2 is the same 6.16 GeV pole that our ff_ps uses, so the two codes agree on
the shape of the 3-body matrix element.

## 3. Q2 — Higgs production

### 3.1 Spectrum file

* File: spectra/SM particles/DoubleDistr_LHC_h.txt (330,892 B, no header, no provenance
  comment). Columns = polar angle theta [rad] (94 values, log-spaced 1e-5 -> pi),
  energy E [GeV] (60 values, 125.21 -> 7366.4), double-differential weight.
* It is imported by the generic mother-particle machinery:
  codes/LLP distribution/mother-sampler.nb:52,74,126 (pattern "DoubleDistr_*_*.txt",
  "spectra/SM particles") and interpolated log-log in (theta,E) by
  MotherInterpolation / DoubleDistrMotherInt (mother-sampler.nb:3468-3560). The manual
  (section 4.1, p.11) states the mother distributions are *user-provided* inputs,
  "normalized by one", and lists no generator.
* The table is in fact not normalised to 1: integral f dtheta dE = 20.3 (checked). That is
  harmless because the mother spectrum is only used to sample the LLP angle-energy
  distribution, which is renormalised to 1 afterwards (manual section 4.1: "The
  distribution is normalized to one").
* I could not find any Pythia8/MadGraph/NNLO statement for this file anywhere in the repo,
  manual or README. What can be said from the data (section 3.3) is what its shape is.

### 3.2 Rates, processes, b_h

* sigma(pp->h) is a single number per facility, no per-process breakdown and no
  reweighting: experiments.nb:222

      \[Sigma]ppTohInpbFacility:=Association[{"FermilabBD"->0.,"SPS"->0.,"LHC"->55,
        "FCC-hh"->552.,"Serpukhov"->0.,"ESS"->0.,"SLAC-E137"->0.,"SuperKEKB"->0.}]

  and it enters as a fraction of the inelastic cross section (experiments.nb:145):

      PppTohExperiment[exp]=(\[Sigma]ppTohInpbFacility[...]/
          Max[10^-90.+\[Sigma]ppInpbFacility[...]])*AtargetExperiment[exp]^0.29

  with sigma_pp,inelastic = 72e9 pb (experiments.nb:63) and Atarget = 1 for pp.
  So the reference h yield is Npot * 55/72e9; no ggF/VBF/ttH mixture is modelled, in
  contrast to our 9-process mixture (bc5_production.py:508-522).
* b_h: 3. Scalar sensitivity.nb:426-438

      BrhtoSS[mLLP_,\[Alpha]_]=2*10^-2*\[Alpha]^2*Sqrt[1-4*mLLP^2/125^2]

  (m_h = 125 fixed in the formula) and :641-654

      PmotherToLLP[...,"h-quartic",Facility]=
          2*2*10^-2*\[Alpha]valBrhToSS[mLLP,BrhToSS]^2*Sqrt[1-4*mLLP^2/125^2]

  i.e. P = 2 x BR(h->SS) (two scalars per h decay), with alpha traded for the requested
  BR by alphavalBrhToSS (:448-462). Two different BR values appear:
  - :4223 and :6052   BrhToSSreference=0.15  -> the value used when the Nevents tables
    were tabulated;
  - :7487-7488        input1=NpotDefault, input2=2.3  and BrhToSSval=0.01 -> the value
    used for the shipped sensitivity curve (file name "...BrhToSS=0.01...").
  NevInt rescales quartic channels by BrhToSS/BrhToSSreference (:6793-6794), i.e. by 1/15
  for the shipped curve. I verified this end-to-end (section 6).

### 3.3 Is the h spectrum Pythia8-like? (measured, not assumed)

I compared the table with our own Pythia8 13 TeV Higgs sample
(work/magnus/samples/higgs.txt.gz, status 62, N=200000; median E=252 GeV, median pT=39 GeV):

* SensCalc dN/deta at eta = 0/1/2/3/4/5: 0.047/0.079/0.113/0.132/0.097/0.045;
  fraction 2<eta<5 = 0.285.
* Ours: dN/deta(2)=0.139, (3)=0.056, (4)=0.003; fraction 2<eta<5 = 0.090.
  (Our Higgs is much more central - as expected for Pythia8 ggF.)
* However, the quantity that actually enters the long-lived rate,
  K(m_S) = < 1_{2<eta_S<5} * m_S/p_z,S > for h -> S S built isotropically, comes out
  almost identical, because the S direction is smeared by the decay:

  | m_S [GeV] | K (SensCalc mother table) | K (our Pythia8 sample) | ratio |
  |---:|---:|---:|---:|
  | 5.08 | 2.363e-3 | 2.321e-3 | 1.018 |
  | 13.02 | 5.622e-3 | 5.367e-3 | 1.048 |
  | 33.37 | 1.164e-2 | 1.102e-2 | 1.057 |
  | 53.42 | 1.463e-2 | 1.479e-2 | 0.989 |

* Cross-check with the archived 2023 reference: reference_2023_Higgs.txt (2820 rows = 47
  theta x 60 E, one hemisphere, m_h=125.1) has exactly the same shape as the local 2025
  file, only normalised over one hemisphere (its dN/deta values are 2.0x the local
  full-sphere ones). The 2023 and 2025 reference h spectra are the same generator/table;
  only m_h changed 125.1 -> 125.21 and the theta range was symmetrised to 0..pi.

Caveat: the prior diagnosis work/physics_review_20261025 measured, with the *precomputed
scalar* distribution DoubleDistr_Scalar_h-quartic_LHC.m (read with Wolfram), a ratio
current/reference = 0.76-0.83 at these masses, i.e. a 17-24% difference that my
mother-spectrum calculation does not reproduce. That must come from the decay/derivation
step inside DoubleDistr_Scalar_h-quartic_LHC.m (or its log interpolation of the scalar
distribution) and cannot be checked without Wolfram; it is the one remaining spectral
difference worth resolving before blaming the mother spectrum.

## 4. Q3 — efficiencies

* epsilon_rec = 0.4, exactly once, as a flat factor:
  3. Scalar sensitivity.nb:6630-6645

      \[Epsilon]reco[mLLP_]=If[StringContainsQ[GivenExperimentForSensitivityComputation,
          "LHCb-downstream"],0.4, If[StringContainsQ[...,"CHARM-lepton"],
          If[mLLP<0.105*2,0.51,0.85],1]];

  and :6793-6794 inside NevInt:

      NevInt[...]=If[StringContainsQ[prod,"quartic"],BrhToSS/BrhToSSreference,1]*
                  \[Epsilon]reco[mLLP]*
                  If[within grid, 10^(Interpolation[...][Log10[mLLP],Log10[theta2]]),0];

  It is applied to every production channel (mixing and quartic) and is not
  mass-dependent for LHCb.
* **0.87 is not applied at all.** A repo-wide grep for "0.87" finds only plot-grid
  numbers; there is no epsilon_S/B, no retention factor, no equivalent anywhere in the
  scalar path (production, acceptance, sensitivity). The visible branching ratio is not a
  separate 0.87-like number either: it is read from the acceptance table as
  BrVis[mLLP] = dataAcceptances[[5]] (3. Scalar sensitivity.nb:1518).
* Consequence for us: if BC5.pdf used a SensCalc run, it cannot contain a 0.87 factor.
  Adding it in our v3 makes us weaker by 0.87 in signal and moves the long-lived edge by
  ~+15% in theta^2 (the archived efficiency_087_mass_dependence.json measures exactly
  1.15 at high mass; the local slope there is 0.996-0.999, consistent with N ~ theta^2).

## 5. Q4 — widths and BRs

* BR table: phenomenology/Scalar/decay widths/BrRatioScalar.mx (10.6 MB Wolfram binary
  dump), imported at codes/llp-decays.nb:13341 and used at 3. Scalar sensitivity.nb via
  BrRatiosDPdata. It is a list of **4 decay descriptions**; the first/default one is
  labelled (:236, :244-264)

      "Description following 1809.01876"
      "Description of decays into PiPi and KK following 2407.13587 (central/lower/upper)
       plus added fake 4Pi decay as in 1809.01876 ..."

  so the default scalar decay model is Winkler 1809.01876 (same as our
  widths/scalar_widths.py "strict Winkler"); the 2407.13587 variants are additional
  uncertainty descriptions only.
  The 2023 version of the same table (git blob f9d658b613, archived as
  reference_2023_BrRatios.m; its text export reference_2023_BR.csv has 1683 rows on a
  0.03-63.0 GeV grid step 0.01) gives BRs that sum to 1 and agree with our reference BR
  model to display precision (prior diagnosis section 4).
* Lifetime: phenomenology/Scalar/decay widths/ctauS.mx (imported at llp-decays.nb:175),
  used as (llp-decays.nb:203-207)

      c\[Tau]LLP[LLPsel,descr,mLLP_,theta2_]=
          (1/theta2)*10^Interpolation[Log10[ctauListXX[[i]][[2]]],
                                      InterpolationOrder->1][Log10[mLLP]]

  i.e. ctau(m,theta2) = ctau_table(m)/theta2, log-log interpolation; the table itself is
  ctau at theta^2 = 1 for each decay description. The 2023 text version of this table
  (reference_2023_ctau.dat, git blob 8f247610d7) covers 0.002-69.902 GeV with
  non-uniform steps 0.001/0.01/0.1. From it:
  Gamma_tot = hbar c / ctau at theta^2 = 1:

  | m_S [GeV] | ctau [m] | Gamma_tot [GeV] |
  |---:|---:|---:|
  | 5.081 | 2.172e-11 | 9.085e-6 |
  | 13.002 | 1.385e-12 | 1.425e-4 |
  | 33.402 | 2.237e-13 | 8.821e-4 |
  | 53.402 | 1.420e-13 | 1.390e-3 |
  | 62.502 | 1.229e-13 | 1.606e-3 |

  Consistency: the BR table and the lifetime table are generated offline from the same
  model; the code never recomputes Gamma_tot from the BRs at run time, so there is no
  internal cross-check. Comparing the 2023 lifetime table with our current widths model
  (scalar_widths.py) shows our Gamma is 10-50% larger at 30-60 GeV (our ctau at 33/53/62
  GeV is 1.72e-13/9.8e-14/8.2e-14 m), i.e. the reference lifetime is longer at high mass.
  I could not re-read the *2025* .mx table (binary), but its default description is the
  same 1809.01876 and the prior diagnosis measured the 2023-vs-current BR/lifetime
  agreement at ~0.3%, so this difference is presumably in our own width model rather than
  a change of the reference table - it should be checked by the widths audit.
* Does the implied lifetime affect the high-mass edge? In the long-lived limit the signal
  scales as 1/ctau, so a 10-50% Gamma difference is a 10-50% rate difference, i.e. a
  2.5-11% effect on the theta^2 edge - not the factor ~1.6 needed to explain the
  reference/current discrepancy.

## 6. Q5 — acceptance

* Table: Acceptances/LHCb-downstream-D/
  Acceptance_LHCb-downstream-D_for_Scalar_decayDescription_Default.m (42 MB binary dump).
  Imported at 3. Scalar sensitivity.nb:1595-1630 as
  {FullAcceptanceData0, ..., decayDescr, ctauLLP[m,theta2]} = dataAcceptances[[1,3,-2,-1]],
  BrVis = dataAcceptances[[5]], CrossSectionData = dataAcceptances[[4]].
* A(z) is applied as a **differential** factor inside the decay-position integral, not as a
  single efficiency. The integrand (3. Scalar sensitivity.nb:4598-4622) is

      IntegrandLLP = DoubleDistrLLPint[m,theta,E,prod] * PdecayDensity[m,coupling,theta,E,z]
                     * AzimuthalAcceptanceInt[theta,z] * DecayAcceptanceInt[m,theta,E,z]

  and IntegralLLP (ibid:4623-4640) integrates it over theta, log E and z in
  [zminTheta(theta), zmaxTheta(theta)] (the z window where the point is inside the decay
  volume). The fast production mode tabulates the same integrand on grids of 30 theta and
  30 z points (ibid:3400-3430, OutGridThetafinal/OutGridzfinal) and integrates the table.
  Interpolation is InterpolationOrder->1 on logarithms of the tabulated acceptance and of
  the data (BlockEpsilonDecay / BlockEpsilonAzimuthal, ibid:1588-1630), i.e. log-log
  linear interpolation.
* Detector parameters (codes/experiments.nb, LHCb-downstream-D entries of explistlhcb;
  block starts at :16744):
  - zToDecayVolumeExperiment = 1 m (:16780), dzDecayVolume = 1.5 m (:16914) ->
    decay volume 1.0 < z < 2.5 m;
  - GeometryExperiment = "Annular cylinder" (:16853), Rinner = z tan(theta(5)),
    Router = z tan(theta(2)) (:16971, :17000) with theta(eta)=2 atan(exp(-eta))
    (experiments.nb:5) => **2 < eta < 5**;
  - detector plane at 9.4 m (:17281 xLongMaxDetectorExperiment = 9.4 for the
    LHCb-downstream non-TDfull experiments; :17183-17185 for zToDetectorCenter);
    dipole magnet 3-7 m (:17338 zMagnetMin, :17352 zMagnetMax, :17370
    MagneticFieldExperiment = 4/(7-3) T m mean field -> 1 T);
  - ECALoption = False and detectable PDGs exclude gamma, K_L, n, nbar for this experiment,
    i.e. photons and K_L are not counted as visible.
  So the reference decay volume and eta cut are the same as our card (1 < z_dec < 2.5 m,
  2 < eta_S < 5). Only the magnet span differs from the paper's 3.5-7.5 m (prior
  diagnosis: 0.1-0.7% effect).

## 7. Q6 — luminosity, N_pot and the statistical criterion (closed)

* NpotDefault is parsed from the shipped Nevents file name (= 2.16e16 pp collisions) at
  3. Scalar sensitivity.nb:6621; the dialog takes input1=NpotDefault and input2=2.3
  (:7487-7488); the event yield is then scaled by Npot/NpotDefault
  (SensitivityExtractor, :7600-7640, "FractionBox[Npot,NpotDefault] * NevTot >= NevMin").
* The shipped sensitivity json
  Sensitivity domains/Scalar/LHCb-downstream-D/
  Sensitivity_Scalar_at_LHCb-downstream-D_BrhToSS=0.01_Nev=2.3_Npot=1.75e15.json
  therefore corresponds to **Npot = 1.75e15 collisions = 25 fb^-1 x 70 mb** (1.75e15 =
  2.5e4 pb^-1 x 7.0e10 pb), with **NevMin = 2.3** and BR(h->SS) = 0.01.
  The shipped Nevents tables are at 2.16e16 = 300 fb^-1 x 72 mb.
* I verified the whole chain by re-evaluating the shipped tables at 40 points of the json
  boundary (bilinear log-log interpolation, sum of the 7 production modes, eps=0.4,
  quartic rescaled by 0.01/0.15). Representative results (N = sum of tabulated Nevents at
  2.16e16; the json contour requires N*0.4*1.75e15/2.16e16 = 2.3):

  | m_S | mode that dominates | N (2.16e16) | N*0.4 | N*0.4*(1.75/21.6) |
  |---:|---|---:|---:|---:|
  | 0.021 | B-mixing | 70.98 | 28.39 | 2.30 |
  | 0.911 | B-mixing | 71.04 | 28.42 | 2.30 |
  | 1.540 | B-mixing | 71.00 | 28.40 | 2.30 |
  | 2.359 | B-mixing | 70.98 | 28.39 | 2.30 |
  | 13.68 | h-quartic | 1065.0 | 425.9 (x1/15 = 28.39) | 2.30 |
  | 33.39 | h-quartic | 1064.6 | 425.8 (x1/15 = 28.39) | 2.30 |
  | 53.41 | h-quartic | 1064.6 | 425.8 (x1/15 = 28.39) | 2.30 |

  (points that do not close are those outside the validity window of the tables, where
  the tabulated N is clipped to 0 or to the plateau; the 40 sampled on-contour points are
  consistent to <1%.) This proves (a) the json is internally consistent, (b) 0.4 and the
  2.3 criterion are applied exactly as described, (c) the shipped curve is the 25 fb^-1
  one, (d) the quartic tables carry BrhToSSreference=0.15.

## 8. Q7 — provenance: is the local copy the one that produced BC5?

* BC5.pdf was created 2023-10-30 (PDF metadata, per the prior diagnosis, and it matches
  Fig. 7 of the paper).
* The 2023 SensCalc tree that matches that date is commit
  6664a923869c6b21a3aba3d4e8429ec3ae9d4edd (repo maksymovchynnikov/SensCalc; the tree JSON
  is archived in run work/physics_review_20261025/reference_2023_tree.json). Its notebook
  names are different from the local ones: "2. FIP distribution.nb" (2,930,048 B),
  "3. Scalars sensitivity.nb" (376,284 B) versus the local "2. LLP distribution.nb"
  (472,892 B) and "3. Scalar sensitivity.nb" (440,530 B). The 2023 scalar notebook is
  archived in that folder and is what I quoted for the dropdown above.
* The local copy is a later snapshot: README references 2503.01760 (March 2025), the
  manual PDF was produced 2025-04-09, notebooks were last saved with Mathematica 14.2
  (Dec 2024) and the shipped sensitivity json lists "Bremsstrahlung-AP-Central/Lower/Upper"
  variants that come from 2409.11096 (Sept 2024). So **the local copy cannot be the
  code that produced the 2023 BC5 curve**, and its shipped json is a 2024/2025
  recomputation, not the original run.
* What did *not* change (verified with git blob hashes: sha1 of "blob <len>\0<content>"):
  the six production-probability tables are byte-identical between the 2023 commit and the
  local copy:
  BtoXsStotal.dat 4ccb2e9f77, BtoXsSinclusive.dat 9fc7590a50, BplustoXsSStotal.dat
  44edd88a94, BstoSS.dat 61b1f2dbcb, KplustoPiS.dat 06bceca5b1, K0LtoPiS.dat 12b2683541.
  The Higgs mother spectrum is the same table (2023 file = one hemisphere, m_h=125.1;
  local = symmetrised, m_h=125.21). sigma(pp->h)=55 pb is the same. So the *inputs* that
  dominate the high-mass end are unchanged between 2023 and now.
* What did change: the BR and lifetime tables were regenerated (2023 BrRatios.m
  183,844 B / ctauS.dat 176,859 B vs local BrRatioScalar.mx 10.6 MB / ctauS.mx 531 KB with
  4 decay descriptions), the LLP-distribution and sensitivity notebooks were rewritten,
  and the two-scalar (Di-decays) machinery from 2503.01760 appeared.
* Conclusion: the physics *prescription* in the local copy is the 2023 one for the parts
  that matter here (exclusive B production default, Winkler decays default, 55 pb h cross
  section, same h spectrum, 0.4 only, 2.3 events), but the local code is a later rewrite
  and the shipped json is not the original BC5 run. The prior diagnosis already showed
  that the local json is itself 0.55-0.62x the PDF along the high-mass long-lived edge,
  so "reference = local SensCalc json" is not a valid identification and the original
  run configuration cannot be reconstructed from this copy alone.

## 9. Implications for the high-mass difference (short)

What is now excluded as the cause of the high-mass difference, quoting numbers:

1. The h transverse-momentum/eta spectrum: the reference table and our Pythia8 sample
   differ strongly in dN/deta (fraction 2<eta<5: 0.285 vs 0.090) but the acceptance-weighted
   K factor agrees to 2-6% (1.02/1.05/1.06/0.99 at 5/13/33/53 GeV) because the h->SS decay
   smears the direction. Only the precomputed *scalar* spectrum suggests a 17-24%
   difference (prior diagnosis), which needs a Wolfram read to settle.
2. sigma(pp->h) (55 pb) and BR(h->SS) (0.01 through 2e-2 alpha^2 beta): the same numbers in
   both codes; they cannot produce a >20% effect.
3. The 0.87 efficiency: not present in SensCalc at all.
4. The high-mass lifetime/width: 2023-vs-now agreement quoted as 0.3% by the prior
   diagnosis; our own width model differs by 10-50% (our Gamma larger), which is worth
   checking but is a <=~10% effect on the edge, not the observed factor.
5. The B/K production channels: they are exactly zero on the high-mass long-lived edge
   (all B/quartic channels close below ~5.2 GeV; the h-quartic channel is 100% of the
   signal there), so the inclusive-vs-exclusive choice (section 2) cannot explain the
   high-mass offset - it matters at m_S < 4.7 GeV, where it is a factor 1.6.

What is a genuine model choice rather than a bug, and where it bites:
* B -> X_s S: reference = exclusive sum (Boiarska), ours = Winkler inclusive. At m=0 the
  rates are 3.31 vs 5.40 theta^2 (ours 1.63x larger); the crossing is at ~3.2 GeV; at and
  above 4.7 GeV both are effectively the same B -> K S formula. So a low/mid-mass
  difference of up to ~1.6 is expected from this choice and is NOT a bug in either code.
  (Separately, our eq:inclusive has no parton-level endpoint at 4.1 GeV - a real
  implementation difference, but one the reference never uses.)
* N_pot: the shipped reference curve is 25 fb^-1 (1.75e15 collisions). If BC5.pdf is that
  curve, comparing it with our 300 fb^-1 curve introduces a factor 12.3 in signal, i.e.
  a factor ~12 in theta^2 on the long-lived quartic branch. In the high-mass long-lived
  limit the reference is N ~ theta^2, so any comparison must be done at the same Npot.
* The two-scalar event treatment: the 2023 reference produced 2 scalars per h decay (factor
  2 in the h-quartic weight), which is what we also do; the new Di-decays OR machinery is a
  2024/2025 addition and is not the 2023 prescription.

## 10. Scripts (all in this directory, reproducible)

extract_nb3.py (WL notebook -> readable input form, keeps line references), inspect_prod.py,
sumcheck.py, shape.py, cmp_tables.py, hspec.py, hspec2.py, hspec3.py, hspec4.py,
cmp_hspec.py, kfactor2.py, cmp_sm.py, json_pts.py, check_npot.py, check_npot2.py, slice.py,
widths_test.py, ref2023.py, ref2023b.py, tree3.py, tree4.py, shacmp.py, parse_json2.py,
probe_mx.py, dump_probe.py.

All claims above with numbers were produced by one of these (D:/py/python.exe, numpy).
