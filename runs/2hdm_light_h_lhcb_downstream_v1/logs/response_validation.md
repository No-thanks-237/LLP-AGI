# T4 — detector response + event engine: validation report

Owner: response.  Design contract: DESIGN.md section 4.  Machine-readable summary:
logs/response_validation.json; BC5 regression: work/analysis/bc5_regression.json.

## 1. Provenance

| item | sha256 |
|---|---|
| work/response/downstream_kernel.dll (copied VERBATIM from inputs/bc5_reuse) | 5bb46a1df5e0af72c2ddfd85cf9810d8d5e7760dc2da40e443ccf0f9b6569a08 |
| work/response/downstream_ctypes.py (verbatim) | 259cb6349fe586fd7510ee67d4ff245e4caa66d598a9081be035b075809ad765 |
| work/response/test_kernel.py (verbatim) | 1427ba479235ae1a56639177e6e20f68f36a36a44703fc3a4f605bd9f58bb213 |
| work/response/twohdm_lib.py (adapted from bc5_lib.py) | c27c3b8062daa2562a25962dd1cebeea6f9c9f95c2e439dd5d3b9f223d50bb59 |
| work/analysis/models.py | 3f930c960e1ef13141d0976cf112ed152c7372ca73e30789312c3763bc0fc5cf |
| work/analysis/response_table.py | (recorded at run time in every cache) |
| work/analysis/event_engine.py | 509d42edb809d2301ff393e82002dfe95343495fd3aff845368eefe344bd71d2 |
| work/analysis/bc5_regression.py | d6a3e9a606b76a26ff6da95c476720966cb07c7983cb4df7ada20c1ea5242d4e |

Every response-table cache and every scan JSON records `code_sha256` = sha256 of
twohdm_lib.py + models.py + response_table.py (in that order) together with the DLL sha256, the
detector configuration, the z grid, the sampling seed, nsamp and ndecay.

Detector configuration (identical to BC5 v3.6, copied into every cache's metadata):
By = 1 T for 3 < z < 7 m, SciFi frontZ 7.7 m / backZ 9.4 m, halfX 3.24 m, halfY 2.415 m,
beam hole 0.09 m, birth cuts E > 5 GeV and pT > 0.5 GeV, 1.0 < z_dec < 2.5 m,
2 < eta_H < 5, event reconstruction efficiency 0.4 (applied once per collision), 13-node z grid.

Input samples (copied into work/magnus/ by the samples owner; sha256 verified here):
bbar_heavy.txt.gz ac8e1a05...37f54b (n_gen 800000), higgs.txt.gz 7a62b1f0...e65dc0 (n_gen 200000),
minbias_pik.txt.gz cc157eee...9de41f (n_gen 500000; read as a 100000-collision event subsample,
section 7), hadron_lib_{gg,ss,cc,tautau}.txt.gz from work/magnus/hadron_lib/.

## 2. Transport kernel (cross-check i)

`python work/response/test_kernel.py` -> **32/33 PASS**.  The single failure is `a2`, the
documented specification conflict (a straight track at x=y=0 is reported BeamHole by the detector
geometry, which is the correct behaviour of requirement f), identical to BC5.  The suite
re-derives the exact uniform-field arc with an independent fine-step integrator, so the DLL is
checked against something other than itself.

## 3. Unit tests (cross-checks ii and iii)

`python work/analysis/test_engine.py` -> **14/14 PASS**:

1. kernel suite 32/33 (above);
2. transport fast path (the per-track call inside twohdm_lib) equals the public
   `downstream_ctypes.transport_track` bitwise on the same tracks;
3. p_window vs analytic A = 1, z0 before / inside the window: agreement to 1e-12.  The integral
   runs over the 13-node domain [1.0, 2.5] m only; for a birth point z0 < 1 m the survival factor
   exp(-alpha (1 - z0)) is exact.  Short-lifetime limit (alpha*h >> 1): agreement to 1e-10 (the old
   BC5 trapezoid was wrong by factors ~1.5 there);
4. p_window with a random piecewise-linear A(z) vs a 200001-point numerical quadrature: worst
   relative deviation **3.7e-09**;
5. per-channel decay kinematics IDENTICAL to bc5_lib for the same random stream: all channels
   (ee, mumu, tautau, pipi, KK, 4pi, cc, ss, gg) with the analytic models AND with the Pythia
   hadron library loaded.  This is the requirement that the per-channel response must not change
   BC5's decay treatment;
6. event-OR algebra vs a brute-force per-event loop: o_event exact to 1e-12, including (a) two
   exclusive histories of one parent (they ADD), (b) a two-scalar history (OR inside the history),
   (c) two parents in one collision (OR), and the check that equal parent indices in different
   collisions are NOT merged (review BUG C);
7. eps normalisation: eps = yield_corr * (n_tot/n_smp) * 0.4 / N_gen * sum_events o, with the 0.4
   applied exactly once; verified against the brute-force value to 1e-12;
8. response-table invariants: scalars sorted by history, histories sorted by (event, parent),
   sampled collisions unique and sorted.

## 4. Probability algebra (DESIGN section 4, implemented in event_engine.py)

    p_i(tanb) = sum_ch BR_ch p_window(mH, pz_i, z0_i, A_ch(i), c_tau)
    o_hist    = w_hist * (1 - prod_{i in hist}(1 - p_i))
    o_parent  = sum_hist o_hist              (histories of one parent are exclusive)
    o_event   = 1 - prod_{parents}(1 - o_parent)
    eps_c     = yield_corr_c * (n_tot,c / n_smp,c) * 0.4 / N_gen,c * sum_events o_event
    s         = 1e3 * L * sum_c sigma_c eps_c      (L in fb^-1, sigma in pb, 1 fb^-1 = 1e3 pb^-1)

The parent-sampling factor (n_tot/n_smp) is essential: the parent samples are subsampled and
without it every group is too small by that factor (B ~ 900, h -> HH ~ 133).  It was missing in the
first version and is now recorded per group in the table (`counts`).  The 0.4 is applied once,
inside eps.  MC error per group: se = norm * sqrt(n_coll * Var(o)) over the SAMPLED collisions,
including collisions with zero response.

## 5. BC5 v3.6 regression (cross-check iv, the end-to-end test)

Same engine + BC5 model adapter (BC5 v3.6 scalar_widths + BC5 production + theta^2 mixing), same
BC5 samples and BC5 hadron library, compared with the published
`bc5_lhcb_downstream_v3/outputs/data/exclusion_boundary.csv` (25 fb^-1 long-lived edge).
16 masses on the BC5 grid, nsamp B/inel/h = 2500/1500/2500, every crossing refined by a local
bisection so the comparison is grid-independent.

| mS [GeV] | published theta2_long | mine (refined) | ratio mine/published | delta_s/s at edge | curve ratio median |
|---|---|---|---|---|---|
| 0.483517 | 5.637e-10 | 5.3198e-10 | **0.9437** | 4.4% | 1.0763 |
| 0.661654 | 3.1948e-10 | 3.202e-10 | **1.0022** | 4.5% | 0.9956 |
| 0.905421 | 1.0733e-10 | 1.0437e-10 | **0.9724** | 3.7% | 1.0341 |
| 0.9884 | 4.527e-11 | 4.6274e-11 | **1.0222** | 3.5% | 0.9763 |
| 1.0 | 4.814e-11 | 4.9434e-11 | **1.0269** | 3.4% | 0.9709 |
| 1.2 | 1.3349e-10 | 1.3199e-10 | **0.9887** | 3.9% | 0.9522 |
| 1.44937 | 1.2641e-10 | 1.3224e-10 | **1.0461** | 3.8% | 0.9229 |
| 1.983347 | 1.3772e-10 | 1.3927e-10 | **1.0113** | 3.6% | 0.9948 |
| 2.2 | 1.144e-10 | 1.1705e-10 | **1.0232** | 3.5% | 1.0161 |
| 2.6834 | 7.6832e-11 | 7.7485e-11 | **1.0085** | 3.6% | 0.9988 |
| 3.174881 | 5.6152e-11 | 5.8643e-11 | **1.0444** | 3.9% | 0.9390 |
| 3.7579 | 2.6361e-11 | 2.6539e-11 | **1.0068** | 4.3% | 0.9918 |
| 4.34457 | 3.8572e-12 | 3.8403e-12 | **0.9956** | 4.6% | 1.0060 |
| 5.082255 | 1.37e-12 | 1.4737e-12 | **1.0757** | 4.8% | 0.9310 |
| 6.954656 | 4.114e-13 | 4.2212e-13 | **1.0261** | 4.5% | 0.9758 |
| 10.5834 | 1.3408e-13 | 1.4621e-13 | **1.0904** | 4.7% | 0.9187 |

**Summary (25 fb^-1 long-lived edge, 16 masses): median ratio 1.0167, median |deviation| 2.46%, max |deviation| 9.04% (13/16 within 5%), curve-ratio median 0.9841.**

Noise budget: BC5 quotes a boundary-position MC error of 3.4% median / 5.0% max for the published curve itself, and its event-level delta_s/s at the long-lived edge is 4.3-4.8% here (mine 3.4-4.8%).  A direct check at the noisiest point (mS = 10.58 GeV, 100% h -> HH) with three independent seeds gives boundary ratios 1.090 / 0.929 / 1.016 and curve medians 0.923 / 1.076 / 0.986, i.e. the +/-9% outliers are the engine's own MC scatter at the h-dominated masses, not a systematic.  The median ratio 1.017 and the curve median 0.984 show no bias.

## 6. 2HDM integration (production + decay modules)

Engine: `python work/analysis/event_engine.py --model 2hdm --mass <m> --nparam 40 --param-lo 5 --param-hi 1e5 --lumi both --b-mode foresee --out <json>`.
Imports work/production/twohdm_production.py (frozen API) and work/decay/decay_widths.py (author
Light_scalar_decay).  The response table contains only mass-dependent geometry/MC, so it is reused
for every tanb and for all three B-rate modes ('foresee' nominal, 'paper_mb', 'exclusive').

Timing at nsamp = 1500/1500/1500 (one process; the parent samples load once):
mH = 0.3 GeV: table 9.4 s (response 5.1 s, 13458 scalars), all 3 groups active;
mH = 1.0 GeV: table 14.0 s (response 9.8 s, 4500 scalars), inel group empty;
mH = 3.07 GeV: table 6.4 s at nsamp 500 (hadron-library channels ss/gg active).
The 66-mass scan is therefore ~15-25 s per mass; a cached re-evaluation of the full 40-tanb x
2-lumi grid costs 5.8 s including the refined crossings.  Per-process memory ~250 MB (samples
stored float32; the 100k-collision minbias subsample is a 70 MB npz).

Boundary at mH = 1.0 GeV (refined): 25 fb^-1 long-lived edge tanb = 357.86 (valid region),
short-lived edge tanb = 6.047 (flagged outside_validity_B: BR(B->X_s H) > 0.1 below tanb = 19.6);
300 fb^-1: 488.58 / 5.727.  At the long-lived edge delta_s/s = 6.5% at nsamp 1500 and the boundary
position error is 0.8% because dlns/dlntanb = -8; with nsamp B/h = 2500-3000 delta_s/s ~4.5%.

PHYSICS SHAPE (flagged to the Lead): the 2HDM signal is a tanb WINDOW.  Above the upper edge the H
is simultaneously long-lived (c tau ~ tanb^6) and gamma-gamma dominated (BR_gaga -> 1), so the
visible rate collapses; the lower edge sits where the H decays before the tracking volume.  At
mH = 1 GeV the lower edge is inside the region where BR(B->X_s H) > 1 (effective description
invalid) and must not be quoted as a limit.

## 7. Known limitations / decisions

* minbias pi/K sample: 22.8 M rows / 525 MB, read ONCE and reduced to a uniform random **event
  subsample of 100000 of 500000 collisions** (all their rows, float32 npz cache, seed 20261025,
  cache name carries the source-file tag).  N_gen = 100000 is used for the normalisation, so the
  estimator stays unbiased (a subsample of collisions is statistically equivalent for an
  event-level OR); the subsample is recorded in the table metadata.  Sanity: 45.68 rows/event in
  the cache (samples owner: 45.64).
* pi/K parent masses: PDG masses (0.13957039 / 0.493677 / 0.497611), not sqrt(E^2-p^2), because
  the %.8g text quantises E and makes the derived mass of a 0.14-0.50 GeV meson noisy above
  100 GeV.  B parents also use the PDG mass; h -> HH uses the per-event Higgs invariant mass
  (BC5 behaviour; 0.2% effect vs mh = 125).  The BC5 adapter keeps BC5's sqrt(E^2-p^2) for
  apples-to-apples regression.
* Channel treatment is BC5's exactly: ee, mumu two-body; tautau = Pythia library when present (else
  the BC5 analytic tau model); pipi = 1/3 invisible + 2/3 pi+pi-; KK = 1/2 K+K- + 1/2 K0K0bar with
  real K_S/K_L flight; 4pi = isospin 4/9, 4/9, 1/9; cc/ss/gg = Pythia library on the exact mass
  grid (analytic jet fallback only when the library is absent); gaga = no charged daughters, A = 0
  by construction.
* Closed channels are masked out (author-program mass windows in twohdm_lib.OPEN_2HDM, checked
  against the decay module's widths at 3 tanb values; a channel outside its window with a width
  above 1e-35 GeV raises).  Closed channels have A = 0 identically, which is physically exact.
* A_ch is estimated with ndecay decays per channel: 10 for the analytic channels, 6 for the
  hadronisation channels.  The MC error of A_ch enters s only through the average over thousands of
  candidates (A errors are independent between candidates), so it is negligible; the dominant MC
  error is the parent sampling.  The per-channel A values themselves carry the statistical error
  of ndecay samples.
* K0L: the production module returns the per-K0L branching ratio; the 1/2 PDG-311 flavour
  projection is applied in the response layer (history_weights).  No K0L -> H e nu (contract).
* B yield correction: derived at run time from f_Bs (truth 1.918 B+,B0,Bs per bb event):
  corr = 1.918*(1-f_Bs)/(N_B+,B0/N_gen) = 1.0172 for the delivered sample, applied once.  The BC5
  regression keeps BC5's own 1.0172 for apples-to-apples comparison.
* h -> HH: BR(h->HH) enters s through sigma_h = sigma_total * BR(h->HH); the per-event history
  weight is the process reweighting sigma_i/(sigma_tot f_i), so eps_h is a pure event-OR
  probability (the same convention as eps_B and BC5's eps_H).
* T4 review fixes (Lead code review 2026-10-25): scalar/history sort order (A), event-id offset
  before sorting (B), parent ids unique per group and parents grouped by (event, parent) (C),
  sorted unique sampled-event lists (D), model-aware long-lived/short-lived edge labels (E), and
  the missing parent-sampling normalisation.

## 8. Reproduction

    python work/response/test_kernel.py
    python work/analysis/test_engine.py
    python work/analysis/bc5_regression.py --nsamp 2500,1500,2500 --stride 2
    python work/analysis/event_engine.py --model 2hdm --mass 1.0 --nparam 40 \
        --param-lo 5 --param-hi 1e5 --lumi both --b-mode foresee --out work/analysis/scan_2hdm/m_1.0.json
