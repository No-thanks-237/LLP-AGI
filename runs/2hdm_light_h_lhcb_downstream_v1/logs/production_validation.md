# Production module validation — Type-I 2HDM light H (DESIGN.md section 3)

Owner: production. Scope: work/production/**, logs/production*.
Status: **module complete, all DESIGN items a-f validated, unit tests 12/12.**
Version: **2HDM-lightH-production-1.0.0** (module MODEL_TAG / __version__).
Machine-readable numbers: work/production/validation.json (written by
work/production/validate_production.py; full stdout kept in
work/production/_validation_stdout.txt).

## 1. Deliverables and reproduction

| file | role |
|---|---|
| work/production/twohdm_production.py | sections 3.1-3.5; pure model physics, numpy complex CKM, no MC |
| work/production/b_channels.py | Boiarska b->s exclusive channel list / masses / kinematic shape |
| work/production/tests/test_production.py | 12 self-contained unit tests |
| work/production/validate_production.py | validation driver for items a-f + Lead requests |
| work/production/validation.json | machine-readable validation numbers |
| work/production/production_reference.csv | BR / sigma table on the run mass grid (reference) |

Commands (from the run directory):

    python work/production/tests/test_production.py     # 12/12 PASS
    python work/production/validate_production.py       # writes validation.json
    python -c "import sys; sys.path.insert(0,'work/production'); import twohdm_production as P; print(P.evaluate(1.0,100.0))"

Timing: evaluate_grid on the full 67 x 40 = 2680-point scan grid: 0.48 s (pure python,
no per-point numpy overhead).  Caches: g0,g1,g2 per tanb; Boiarska shape per mH.

Provenance (sha256, 2026-10-25):

    inputs/references/FASER_light.tex      ac5caa0268ca1be82dd767d2a21fefeab15e61d5c71d2d7dbfc203ee2a7e1a8f
    inputs/bc5_reuse/bc5_production.py     11fb6f925e7323746162f044f2ece55ad63293d3329ccec792f7ee9e12ec7367
    work/common/mass_grid.txt              3edeb6d4a09a69ac18f2c1bd8fd62e23a9ad0107495addc2579298752cf3c21f
    DESIGN.md                              e8df2cb127ec4ee7a82f002d80b4bbd653a3f29ae6c6de64b98a46fed71c2c7d
    (external, read-only) references/tools/FORESEE/Models/2HDM_H/2HDM_H.ipynb  38f64d32b5e16756ed933198af6a72dbbe41d2e752abe15ab8268d60d7af3741

The App. D functions were independently re-evaluated with the Wolfram Language MCP tool
(mcp__wolfram__WolframLanguageEvaluator, 30-digit precision); see item a below.
The B-rate was independently compared with the paper's own FORESEE 2HDM_H model
(references/tools/FORESEE/Models/2HDM_H/2HDM_H.ipynb); see section 4.

## 2. Conventions actually implemented (all from DESIGN.md section 1)

    m_h = 125, v = 246, m_A = m_Hpm = 600, lambda v^2 = 0
    c = cos(beta-alpha) = 1/tanb, s = sin(beta-alpha) = +sqrt(1-c^2)
    xi_H^V = c, xi_H^f = c - s/tanb = 1/(tanb^3 (1+s)), xi_h^V = s, xi_h^f = s + c/tanb
    lambda_HH+H- = -(2 m_Hpm^2 - m_H^2) c / v

xi_H^{ij} (paper eq. xi_H_2HDM, lower row = H), implemented literally:

    xi_H^{ij} = -(4 sqrt2 G_F)/(16 pi^2) sum_k V*_ki m_k^2
                [ g1(x_k,x_H) c - g2(x_k,x_H) s - g0(x_k,x_H) (2v/m_W^2) lambda_HH+H- ] V_kj

* the H row is + g1 cos(b-a) + g2 (-sin(b-a)) - g0 (2v/mW^2) lambda_HH+H-; we implement
  g1*c - g2*s - g0*(2v/mW^2)*lambda.  The g1*c and -g2*s terms are both negative at
  c > 0, and lambda_HH+H- < 0 flips the sign of the g0 term relative to g1*c.
* the three bracket terms are printed per loop quark below (item a).
* m_k: m_u = 2.2e-3, m_c = 1.28, m_t = 173.1 GeV; m_W = 80.385 GeV; x_H = m_Hpm^2/m_W^2 = 55.71.
* CKM: standard parametrisation, exactly unitary (see section 5 for the one documented
  difference from the DESIGN.md nominal |V| list).

## 3. Validation items a-f

### a. xi_H^{bs} hand calculation, mH = 1, tanb = 100  (PASS)

    prefactor -(4 sqrt2 G_F)/(16 pi^2) = -4.178254034734e-07
    lambda_HH+H- = -(2*600^2 - 1)*0.01/246 = -2.9268252033e+01 GeV
    2v/mW^2 * lambda                       = -2.2284959888e+00

per loop quark (bracket_k = g1_k c - g2_k s - g0_k*(2v/mW^2)*lambda):

| k | g1_k*c | -g2_k*s | -g0_k*k | bracket | V*_ks V_kb | m_k^2*ckm*br |
|---|---|---|---|---|---|---|
| u | -7.500000e-03 | -3.514977e-02 | -7.499990e-07 | -4.2650522886e-02 | 3.547e-04 -7.800e-04i | -7.321736e-11 +1.610060e-10i |
| c | -7.500000e-03 | -3.509492e-02 | -7.499487e-07 | -4.2595665762e-02 | 3.973e-02 +1.298e-06i | -2.772458e-03 -9.061040e-08i |
| t | -7.500140e-03 | -1.051249e-02 | -6.102262e-07 | -1.8013241758e-02 | -4.008e-02 +7.787e-04i | +2.163345e+01 -4.202763e-01i |

    sum_k = 21.630681759 - 0.420276348i
    xi_H^{bs} = -9.0378483335e-06 + 1.7560213449e-07i ,  |xi_H^{bs}| = 9.0395541156e-06

Independent checks:

* full numpy re-implementation of the paper formulas from scratch: xi_H^{bs} identical
  to machine zero (delta = 0.0e+00).
* Wolfram Language re-evaluation of App. D (g0,g1,g2,X1,X2 for k=u,c,t at tanb=100):
  maximum relative difference 2.0e-12 over all 15 numbers.
* Lead's independent bracket calculation: our top-quark bracket is -0.36091937 at
  tanb=5 and -0.01801324 at tanb=100, against the Lead's -0.3609 and -0.0180
  (agreement to the 4 quoted digits).  This fixes sign and normalisation of
  (g1 c - g2 s) and of the -g0 term path-by-path.

Per-term breakdown at tanb=5 (mH=1), top quark: g1*c = -1.511182e-01,
-g2*s = -2.049194e-01, -g0*k = -4.881809e-03, bracket = -3.609194e-01.
At tanb=100: -7.500140e-03, -1.051249e-02, -6.102262e-07, bracket = -1.801324e-02.

### b. Br(B -> X_s H) at mH -> 0, tanb = 100  (PASS)

Hand check at (mH=0.1, tanb=100), every factor printed:

    12 pi^2            = 118.43525281
    v^2/m_b^2          = 3463.51960807
    (1-mH^2/m_b^2)^2   = 0.99885567
    1/f(mc^2/mb^2)     = 1.97590341   (f = 0.50609761)
    |xi_H^{bs}/Vcb|^2  = 4.90885303e-08
    product            = 3.97417620e-02
    x Br(B->Xc e nu)   = 4.13314325e-03   (0.104)
    module total       = 4.13490964e-03   (mode 'foresee'; the 0.04 % difference is the
                                           m_B vs m_b phase-space factor)

Boiarska relation (their eq:exclusive-to-inclusive, exclusive/inclusive = 0.5):

    mH=0.1  : inclusive/exclusive = 1.958   (expected 2, -2.1 %)
    mH=0.021: 1.959 ; mH=1: 1.844 ; mH=2: 1.524 ; mH=3: 1.152 ; mH=4: 1.346
    Boiarska total exclusive at mS=0 = 3.3 theta^2 (their eq:total-exclusive);
    our Boiarska-normalised exclusive sum = 3.307 theta^2 (+0.2 %)

The inclusive (nominal) B rate and the exclusive-channel sum agree at the factor-2 level
exactly as Boiarska's own comparison shows; this is the model ambiguity that motivates
the rate-mode systematic (section 4).

Validity flags (tanb below which Br(B->XsH) > 0.1 / > 1; monotonic in tanb, checked):

    mH=0.021: 20.34 / 6.44     mH=0.5: 20.16 / 6.38     mH=1: 19.61 / 6.21
    mH=2:     17.43 / none     mH=3:   13.78 / none     mH=4:  8.67 / none

evaluate()['br_B_to_XsH_validity'] is 'OK' / 'SATURATED_BR_GT_0.1' /
'UNPHYSICAL_BR_GE_1' and ['outside_validity_B'] bool; the scan must record it.
The region tanb < ~20 (mH below ~3 GeV) is outside the validity of the
inclusive/effective description; those points are still computed but must be flagged.

### c. K -> pi H at mH -> 0, tanb = 100  (PASS)

Amplitude (paper Eq. Matrix_4FH):

    M = G_F^{1/2} 2^{1/4} xi_H^W [7 lambda (mK^2+m_pi^2-mH^2)/18 - 7 A_K mK^2/9]
        + (xi_H^{ds}/(2v)) m_s (mK^2-m_pi^2)/(m_s-m_d) f0 ,   A_K = 0, xi_H^W = c
    xi_H^{ds} = -7.4629867252e-08 - 3.1964607095e-08i ,  |xi_H^{ds}| = 8.118715e-08

At mH = 0 (p* = 0.227109 GeV):

    |FCNC|^2                 = 1.395091e-21   -> BR_FCNC-only = 9.733970e-07
    |4-quark|^2              = 1.660697e-24   -> BR_4q-only   = 1.158718e-09 (0.12 % of FCNC)
    interference 2Re(M4q Mfc*) = -8.849158e-23  (-6.76 % of |M|^2, destructive)
    |M|^2 total              = 1.308260e-21   -> BR(K+->pi+H) = 9.128125e-07

FCNC-only / [Boiarska master formula rescaled by (xi_H^{ds}/xi_Boiarska^{ds})^2] = 0.9972.
The four-quark term is negligible in rate but NOT in coherence: it shifts the total by
-6.2 % through the interference, so it must be kept (as DESIGN requires).
K0L -> pi0 H at mH=0: 3.943562e-06 per K0L; K0L/K+ = 4.320 (Boiarska table ratio 4.12).
K0S -> pi0 H = 0 (CP suppressed).  The response layer applies the 1/2 PDG-311 flavour
projection; the module deliberately does NOT include it.

### d. Br(pi -> H e nu), Br(K -> H e nu) at mH -> 0, tanb = 100  (PASS)

    Br(pi+ -> H e+ nu) = 1.963947e-13     (at tanb=1, c=1: 1.963947e-09)
    Br(K+  -> H e+ nu) = 3.911730e-12     (at tanb=1, c=1: 3.911730e-08)

Boiarska table tab:BR-3-body: BR(K->S e nu)/f(x)/theta^2 = 4.1e-8; our theta=1
normalisation is 3.91e-8 (-4.6 %, their table/rounding; the DESIGN formula is followed
exactly, including (1-m_mu^2/m_X^2)^2 where Boiarska print m_h).
K0L -> H e nu = 0 (not a physical two-body leptonic channel).  pi threshold at mH=0.13957.

### e. g_C, g3, Gamma(h->HH), BR(h->HH)  (PASS, one flagged paper inconsistency)

    g_C  = s/(2v)[(mH^2 - mh^2)(2 cot2b s c - c^2 + s^2) - mH^2]   (lambda v^2 = 0)
    g3   = -2 g_C
    Gamma(h->HH) = |g3|^2/(32 pi mh) sqrt(1-4mH^2/mh^2)
    BR(h->HH)    = Gamma(h->HH)/(Gamma_h^other + Gamma(h->HH))

Lead's requested points:

| mH | tanb | g_C | g3 | Gamma_hHH | Gamma_h_other | BR(h->HH) |
|---|---|---|---|---|---|---|
| 1 | 5  | +6.390935e-01 | -1.278187e+00 | 1.2999e-04 | 4.153201e-03 | 3.0350e-02 |
| 1 | 10 | +1.571491e-01 | -3.142983e-01 | 7.8599e-06 | 4.090811e-03 | 1.9177e-03 |

(the paper's rough expectations were ~3 % and ~0.2 % for BR(h->HH): reproduced).

Against the paper text approximation -mh^2 c^2/(4v):

    mH=0.021, tanb=100: g_C = +1.587050e-03 vs 1.587907e-03  -> |ratio| = 0.99946
    mH=0.1            : +1.567621e-03 vs 1.587907e-03        -> |ratio| = 0.98722
    mH=1              : -4.445741e-04 vs 1.587907e-03        -> |ratio| = 0.2800

The mH=1 deviation is not an error: the exact App. C bracket carries the extra term
-(mH^2 c^2/2), so the paper approximation only holds for mH^2 << mh^2 c^2; the
convergence as mH -> 0 is explicit above.  The paper's sign (-mh^2 c^2/(4v)) differs
from the App. C formula (+mh^2 c^2/(4v) at mH->0); BR uses g^2 so this is immaterial,
but it is flagged.

Against 4700 (g/v)^2 (paper eq:brhphiphi):

    4700 reproduces v^2/(8 pi mh Gamma_h^SM) = 4732.9  (0.7 %)
    the correct comparison includes the self-consistent denominator and the modified
    Gamma_h_other:
      BR = 4700(g/v)^2 * [Gamma_h^SM/Gamma_h_other] / (1 + Gamma_hHH/Gamma_h_other)
    reproduced to <0.4 % at every sampled point.  (The raw ratio BR/[4700(g/v)^2] is
    0.9567 at (mH=1, tanb=5) — that is exactly the 1/(1+BR/(1-BR)) factor — and 1.0069
    at tanb=100.)

FLAG: the paper text also claims Br(h->HH) ~ 75 c^4; the same paper's exact App. C
g_hHH plus their own 4700 relation give 19.7 c^4.  We implement the App. C formula and
the self-consistent BR (not 75 c^4, not a fixed 1 %).  This is a factor 3.8 in the
h->HH rate w.r.t. the paper's prose; recorded as a paper inconsistency.

### f. kappa -> 1 and sigma_total -> 55.4527 pb  (PASS)

| tanb | kappa_f - 1 | kappa_V - 1 | kappa_g - 1 | kappa_gamma - 1 | sigma_total [pb] |
|---|---|---|---|---|---|
| 10 | +4.987e-03 | -5.013e-03 | +4.987e-03 | -7.707e-03 | 55.886382 |
| 100 | +5.000e-05 | -5.000e-05 | +5.000e-05 | -7.695e-05 | 55.457037 |
| 1000 | +5.000e-07 | -5.000e-07 | +5.000e-07 | -7.695e-07 | 55.452743 |
| 1e5 | +5.000e-11 | -5.000e-11 | +5.000e-11 | -7.695e-11 | 55.452700 |

SIGMA_H sum = 55.452700 pb; sigma_total(1e5) - sum = 4.3e-09 pb.  kappa_g == kappa_f
exactly in Type-I (all fermion modifiers equal) — implemented as the explicit loop
ratio of DESIGN 3.5 and unit-tested.

Per-process at tanb=5 (kappa_f=1.019796, kappa_V=0.979796, kappa_g=1.019796), pb:

    ggF 50.5224 | VBF 3.6307 | WH 1.3142 | ZH 0.8477 | bbH 0.2479 | ttH 0.5268
    tHj 0.0770 | tHW 0.0152 | tHb 0.0048 | total 57.1867

each equals sigma_SM * kappa_process^2 exactly (ratio 1.000000), except tHW whose
interference bracket is 0.961 of kt^2 for cos(phi)=1 (see below).
At tanb=1000 all nine equal the SM values and the total is 55.4527 pb.

tHW bracket (DESIGN 3.5): sigma = sigma_SM [kt^2 + kV^2 + 2 kt kV cos(phi)]/4.
Baseline phi=0 is the maximum; the phi in [0,pi] bracket is
[(kt-kV)^2/4, (kt+kV)^2/4].  At tanb=5: base 1.5194e-02 pb, range
[6.08e-06, 1.5194e-02], weight in the total 0.0266 %; at tanb=100/1000: base
1.5200e-02 pb, weight 0.0274 %.  Negligible as expected (2.7e-4 of the total), and not
applied as a global kappa.

## 4. B-rate normalisations (Lead request) and the FORESEE cross-check

The paper's own FORESEE 2HDM_H model was found locally
(references/tools/FORESEE/Models/2HDM_H/2HDM_H.ipynb).  Its B rate is

    BR = 9.68544e-13 * (mB/Gamma_B) * (1-mH^2/mB^2)^2 * F(tanb)^2
    F(c) = -3/(4c) - 0.1415/c^3 - (1.0386/c - 0.1415/c^3) sqrt(1-1/c^2) - 0.608446664/c^3

and F is their fitted version of our bracket (g1 c - g2 s - g0 k).  Note the mB
(not m_b) in the phase space, which is how their reach figure gets to mH ~ 5.3 GeV.

(i) bracket check: our exact top-quark bracket / F(c) over tanb = 5 ... 1e5:

    1.00699, 1.00708, 1.00710, 1.00711, 1.00711, 1.00711, 1.00711, 1.00711, 1.00711

    -> 0.7 % agreement, flat in tanb (as the Lead expected).

(ii) BR check at 54 (mH, tanb) nodes with mH < 4, tanb in {5,10,30,100,300,1000}:

    our(mode='foresee')/FORESEE = 1.0162 ... 1.0165, mean 1.0164  -> 1.6 % agreement.
    B0: ours = 3.56688e-03 vs FORESEE 3.50925e-03, ratio 1.0164 (identical to B+ since
    both now carry BR(B0) = Gamma_B+/Gamma_B0 * BR(B+); before the B0 lifetime fix this
    ratio was 1.096).

So the nominal 'foresee' mode reproduces the paper's own B physics to ~1.6 %.

Rate modes (module DEFAULT_B_MODE = 'foresee'), tanb = 100:

| mH | foresee | paper_mb | exclusive | foresee/paper_mb | foresee/excl |
|---|---|---|---|---|---|
| 0.1 | 4.13491e-03 | 4.13314e-03 | 2.11162e-03 | 1.000 | 1.958 |
| 1.0 | 3.84628e-03 | 3.67779e-03 | 2.08631e-03 | 1.046 | 1.844 |
| 3.0 | 1.89701e-03 | 9.72939e-04 | 1.64692e-03 | 1.950 | 1.152 |
| 4.0 | 7.50698e-04 | 2.93849e-05 | 5.57866e-04 | 25.55 | 1.346 |
| 4.5 | 3.09408e-04 | 0 | 2.67380e-04 | - | 1.157 |

High-mass behaviour (tanb=100):

    mH=4.17: foresee 5.8532e-04  paper_mb 9.4503e-08  exclusive 4.4887e-04
             open channels: K, K*892, K0*700
    mH=4.20: 5.5761e-04 / 0 / 4.2873e-04
    mH=4.50: 3.0941e-04 / 0 / 2.6738e-04   (only K open)
    mH=4.78: 1.3440e-04 / 0 / 4.6991e-05
    mH=4.79: 0 / 0 / 0                     (above the mesonic endpoint)
    mH=5.00: 0 / 0 / 0

B0 lifetime scaling (all modes): Gamma(B->XsH) is meson independent while the total
widths differ, so BR(B0->XsH) = BR(B+->XsH) * Gamma_B+/Gamma_B0 with
Gamma_B+/Gamma_B0 = tau_B0/tau_B+ = 1.519/1.638 = 0.92735 (module constant
B0_LIFETIME_SCALE; same convention as Boiarska's caption and FORESEE's mB/Gamma_B).
Measured in the module: B0/B+ = 0.92736 (foresee; the 1e-5 residual is the per-meson
m_B in the phase space), 0.92735 (paper_mb), 0.92735 (exclusive).  The normalised
'channels' fractions are unchanged and sum(channels) == total for both mesons.
This also removes the previous 7 % B0 discrepancy with FORESEE (now ratio 1.0164 for B0,
same as for B+ instead of the earlier 1.096).

Mesonic endpoint (Lead ruling): the lightest strangeness-carrying recoil is a kaon, so
no B -> X_s H channel is open above m_B - m_K = 4.785663 (B+) / 4.785973 (B0).  In
mode='foresee' the module now returns total = 0.0 at/above that mass and sets
'beyond_mesonic_endpoint' = True (FORESEE's m_B phase space would otherwise extrapolate
unphysically to m_B).  Other modes keep their natural thresholds (paper_mb at 4.18,
exclusive per channel).  Measured (tanb=100, foresee): mH=4.70 -> 1.7805e-04,
4.78 -> 1.3440e-04, 4.7857 -> 0 (flag True), 4.79/5.0 -> 0.  Below the endpoint the
K channel is open (B->K H), so mH=4.5 still has 3.0941e-04 (the endpoint is 4.786, not
4.5).  'channel_shape_extrapolated' remains as a diagnostic for the open-shape fallback.
Between 4.18 and 4.786 the paper's m_b formula closes unphysically early (B->K H is
open until 4.786), which is why 'foresee' is the nominal mode.

## 5. CKM choice (Lead ruling: exactly unitary)

Standard parametrisation with (s12,s23,s13) = (|Vus|,|Vcb|,|Vub|) = (0.2243, 0.0408,
0.00382) and delta = gamma = 1.144 rad:

    0.974513 0.224298 0.003820
    0.224176 0.973694 0.040800
    0.008331 0.040122 0.999160       max |V^dag V - 1| = 1.1e-16
    GIM: sum_k V*_ki V_kj = 0 for i != j to 2e-18

The DESIGN.md nominal list (rounded, from different measurements) is not unitary: it
gives |Vtd Vts| = 3.569e-04 versus 3.3425e-04 for the unitary matrix (+6.8 %), hence
BR(K+->pi+H) = 1.0334e-06 instead of 9.1281e-07 (+13.2 %).  Per the Lead's ruling we use
the unitary matrix; ckm_matrix_from_magnitudes() reproduces the nominal list for the
record.  Impact on the final tanb boundary at low mH is ~1-2 %.

## 6. Other differences from the paper / approximations (documented, not silent)

1. B -> X_s H kinematic shape: Boiarska's exclusive-channel mixture is used as the
   recoil shape (renormalised to 1, b->s channels only; the b->d 'pi' entry dropped),
   while the absolute rate is the inclusive/mode formula.  Never summed on top.
2. B_s is not included (Boiarska and SensCalc/BC5 also exclude it).
2b. B+ / B0 lifetime: Gamma(B->XsH) is meson independent while the total widths differ,
   so BR(B0->XsH) = BR(B+->XsH) * Gamma_B+/Gamma_B0 = 0.92735 * BR(B+->XsH) in ALL
   modes and for the exclusive channels (BR = (1/Gamma_B) x partial width).  It cancels
   in the normalised 'channels' mixture, where sum(channels) == total for each meson.
3. kappa_Zgamma: Z gamma is W-loop dominated; we use kappa_V for it (0.15 % of the SM
   width, <0.02 % of Gamma_h_other) — documented approximation.
4. h -> HH uses the App. C formula and the self-consistent BR; the paper's prose value
   75 c^4 and the sign of its g_hHH approximation are not reproduced (flagged above).
5. The K0L rate is returned per K0L; the 1/2 PDG-311 flavour projection is applied by
   the response layer, not here.  The K0 isospin coefficient is not given in the paper;
   we use the K+-normalised relation (factor 1) as DESIGN 3.3 instructs.
6. eta/eta', Upsilon and D production are not included (the paper drops them; D is 11
   orders of magnitude smaller).
7. Quark masses for the K amplitude: m_s = 0.096, m_d = 0.0047 (DESIGN values; identical
   to the author program's Parameters.h lines 29-30).
8. m_t = 173.1 GeV (MSbar) is used in the xi_H loop; the difference to 172.5 GeV is a
   ~0.7 % effect on xi_H (m_t^2 weighting) and is absorbed by the B-rate comparison
   with FORESEE (which fits its own coefficients).

## 7. Unit tests (12/12 PASS)

test_couplings, test_ckm, test_g_functions_vs_wolfram, test_xi_limits, test_lambda_HH,
test_br_B_nominal, test_br_B_exclusive_variant, test_K_pi, test_semileptonic,
test_higgs_sector, test_grid_sanity, test_b_channels.
They cover: coupling identities (c - s/tb == 1/(tb^3(1+s)) for tanb up to 1e5), exact
CKM unitarity + GIM, the Wolfram g0/g1/g2 numbers, xi_H ~ 1/tanb, the lambda_HH formula,
sum(channels)==total in every B mode, 1/tanb^2 scalings, thresholds (pi / K / h->HH),
kappa_g == kappa_f, sigma_total -> 55.4527 pb, the 4700(g/v)^2 relation within 2 %, the
BR_gt_0.1 validity flags, and the b->s-only channel list.
## 8. Version, changelog, frozen artefacts

    MODEL_TAG = '2HDM-lightH-production-1.0.0'   (module-level, also __version__)

Changelog (all same day, DESIGN.md section 3):
 v1.0.0 frozen:
 * B-rate modes added per Lead: DEFAULT_B_MODE = 'foresee' (paper's FORESEE 2HDM_H
   treatment, B-meson phase space) + 'paper_mb' + 'exclusive' systematics.
 * BR(B0->XsH) = BR(B+->XsH) * Gamma_B+/Gamma_B0 = 0.92735 * BR(B+->XsH) in all modes
   (Gamma(B->XsH) meson independent; tau_B+ = 1.638 ps, tau_B0 = 1.519 ps).
 * mode='foresee' returns 0 at/above the mesonic endpoint m_B - m_K (4.785663 B+ /
   4.785973 B0) with 'beyond_mesonic_endpoint' = True.
 * exactly-unitary CKM adopted per Lead ruling; FORESEE cross-check added.
 * tHW interference bracket, validity flags (Br(B->XsH) > 0.1 / > 1) and
   tanb_valid()/validity_boundaries() added per Lead.
Frozen sha256 (recompute with sha256sum after any edit):
    work/production/twohdm_production.py   4cc687ad2bac7334763d893dde0edfb81310992638f1b158bf9cf4c0025a0ad6  [superseded, see below]
    work/production/validation.json        b94b4160ebb1b4d8f0444491f941bc49e369b480495bdc4948e469e6d6d6b5b6
    work/production/b_channels.py          6cf70445852fe027fe45d865dc4a1ce94fc9a11bfae26ffcbfa5f758f8651c4d

2026-10-26 (review round) provenance update:
 * the shipped module is MODEL_TAG = '2HDM-lightH-production-1.0.1' (see the file).  The only
   change since 1.0.0 is the work/common -> work/analysis/grids path refactor; the physics is
   identical (the independent formula audit of logs/production_formula_audit.md reproduces
   every production identity to 3.1e-16).
 * therefore the current frozen sha256 of twohdm_production.py is
       6381bb26049be4c241d80c9b0a43787115897edf3072f4ea969d6fcb1c868c79
   (this is the value recorded in outputs/sha256_manifest.txt and in every scan row);
   4cc687ad... above is the superseded 1.0.0 value.  b_channels.py / validation.json unchanged.
 * audit finding F7 also notes MODEL_TAG 1.0.1 vs __version__ = '1.0.0'.  Left as-is on
   purpose: editing __version__ would change the module sha256 and break the provenance link
   to the already-produced scan outputs.

Reference table work/production/production_reference.csv: 67-point mass grid x 8 tanb
values (5,10,30,100,300,1e3,1e4,1e5), columns br_B_foresee / br_B_paper_mb /
br_B_exclusive / K,pi BRs / BR_hHH / Gamma_h_other / sigma_total / validity flag.
