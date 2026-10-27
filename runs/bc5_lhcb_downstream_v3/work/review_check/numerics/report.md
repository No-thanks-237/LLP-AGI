# Audit of F03 / F10 / F11 / F12 - BC5 LHCb-downstream v3 (model tag v3.3-2026-10-21)

Auditor: audit-numerics (task-1). Production code was read only; every number below was
re-derived here. Scripts + corrected modules live in this directory:
`f03_current_tests.py`, `f03_current_3body.py`, `f03_correct_tests.py`, `f03_batch_tests.py`,
`correct_ps.py`, `batch_ps.py`, `f10_setup.py`, `f10_analysis.py`, `f10b.py`, `f10c.py`,
`f10d.py`, `f11_f12.py`, `f11b.py`, `f12_effect.py`.

## Verdict table

| claim | verdict | decisive number | how measured |
|---|---|---|---|
| **F03** `nbody()` is not uniform n-body phase space (massive 3-body, any 4-body) | **CORRECT** | massless 4-body `<E_i/M>` = 0.33346, 0.22315, 0.22176, 0.22164 +/- 0.00066 (exact 0.25); massive 3-body `<s23>` = 0.32667 +/- 0.00044 vs exact 0.352813 (ratio 0.926); s23 shape chi2 = 8486/19 | ran the production `bc5_lib.nbody` directly, 30k/150k events; exact Dalitz density and analytic dPhi3 marginal |
| **F10** 13-point lifetime trapezoid has a discrete bias | **PARTLY CORRECT** | A=1 benchmark reproduces the review table exactly (1.12694, 1.01443, 0.55473, 0.93280, 1.08737); but realistic code/exact = 1.00000000 at the long-lived edge (all masses) and 1.0250/1.0479 (m=1), 1.0220/1.0404 (m=3.56) at the short-lived edge; the z0-inside-window worst cases do not occur | `run_scan.run_mass(nsamp=300, ndecay=8)` for m=0.2802, 1.0, 3.5567; real (z0,pz,A) rows; exact piecewise-linear-A integral |
| **F11** adaptive/QA pass looks for crossings at level=2.3 although the stored signal is per fb | **CORRECT** | code (level=2.3, i.e. L=1 fb^-1): 0/67 > 5 %, median 4.82 %, max 4.99 %; L=25: **7/67**, max 5.84 %; L=300: **10/67**, max 6.36 %; short-lived edge: 46/46 > 5 %, median 27.1 % / 38.0 % | recomputed crossings of `signal_scan.json` at 2.3/L for L=25 and 300 |
| **F12** `gamma` missing from `CHANNELS`, remainder goes to `gg` | **CORRECT** | BR(gamma) = 1.371 % at 0.1 GeV, **8.919 % at m=0.18869 GeV**, 1.20 % at 0.2118, 2e-5 at 0.2802, <=4e-5 above 0.3 GeV; sum(BR over CHANNELS) = 0.9863 at 0.1 GeV; visible-rate effect = **exactly 0** for m < ~0.25 GeV (gg fallback is empty there) and <= 3e-5 elsewhere | `sw.widths()` vs `L.CHANNELS`, `channel_branching`; per-decay acceptance MC with the real kernel |

---

## F03 - `bc5_lib.py:71-120` nbody(): production sampler is NOT uniform phase space

**Mechanism.** Lines 86-90 of the production file state the factorization correctly, but the
rejection at lines 99-109 only weights the *first-level* |p*|. The recursive call (line 116)
returns a *normalised* (unit-probability) sub-decay, so the q^2 marginal is
`p*(q)` instead of the correct `p*(q) * Phi_(n-1)(q)`; the (n-1)-body phase-space **volume**
`Phi_(n-1)(q)` is missing. For massless 3-body `Phi_2` is constant, which is why the existing
massless-3-body test passes and the bug was missed. For massive 3-body
`Phi_2(q) = lambda^(1/2)(q^2,m2^2,m3^2)/q^2` is not constant; for 4-body even massless
`Phi_3(q) ~ q^2` is not constant.

**(a) massless 4-body, production code (30,000 events):**
`<E_i/M>` = 0.33346 +/- 0.00068, 0.22315 +/- 0.00066, 0.22176 +/- 0.00066, 0.22164 +/- 0.00066.
Correct: 0.25 for all. Ordering-dependent: `<E1 E2>/M^2` = 0.0690 vs `<E2 E3>/M^2` = 0.0451.
(Total four-momentum is conserved to 1e-16, as in the review.)

**(b) massive 3-body, independent exact integral** `dPhi3/ds ~ sqrt(lam(M^2,m1^2,s)) * sqrt(lam(s,m2^2,m3^2))/s`
(Simpson, 40001 points, exact `<s23>` = 0.352813 GeV^2 for 3 pions, M = 1 GeV - reproduces the
review's 0.352813):

| case | exact <s23> | production nbody | ratio |
|---|---|---|---|
| 3 pions, M=1 | 0.352813 | 0.326668 +/- 0.000438 | 0.926 |
| 3 x 0.5 GeV, M=2 | 1.583333 | 1.489141 +/- 0.001030 | 0.941 |
| 0.4/0.3/0.2 GeV, M=1 | 0.303566 | 0.293753 +/- 0.000091 | 0.968 |

The full s23 distribution is wrong: chi2/ndf against the exact density (20 equal-probability
bins) = **8486 / 19** for production vs 13.0/19 for the corrected sampler below.

**(c) permutation test.** 4 identical massless particles: production gives per-index means
0.333/0.223/0.222/0.222 (not 1/4 each) and unequal two-particle correlators - see (a).
A correct sampler must be permutation symmetric; the new one gives 0.25024, 0.24961, 0.24993,
0.25022 +/- 0.00025 and all pair products 0.0583-0.0584.

### Corrected sampler (validated; `correct_ps.py`, vectorised in `batch_ps.py`)

Exact, unweighted, no matrix element needed: decompose dPhi_n as a binary tree of two-body
splits. Every internal node contributes a factor `lambda^(1/2)(s_node; s_L, s_R)/s_node <= 1`,
so rejection sampling inside a constant box is exact (and generalises to any n):

* n=2: two-body (unchanged).
* n=3: sample s23 in [(m2+m3)^2, (M-m1)^2] with density
  `g(s) ~ sqrt(lam(M^2,m1^2,s)) * sqrt(lam(s,m2^2,m3^2)) / s` by rejection in a cached
  envelope (this includes the *second* two-body factor the production code misses).
* n=4: tree (12)|(34): sample s12, s34 uniformly in their allowed boxes and accept with
  `sqrt(lam(s12,m1^2,m2^2))/s12 * sqrt(lam(s34,m3^2,m4^2))/s34 * sqrt(lam(M^2,s12,s34))/M^2`
  (each factor <= 1); then two two-body decays and two boosts.

Validation of the new sampler (all in `f03_correct_tests.py` / `f03_batch_tests.py`):

| test | result | expected |
|---|---|---|
| massless 4-body `<E_i/M>` (200k) | 0.25024, 0.24961, 0.24993, 0.25022 +/- 0.00025 | 0.25 |
| massless 4-body pair products | 0.0583-0.0584 all equal | permutation symmetric |
| 4-momentum conservation | max abs(sum E - M) = 1.6e-12, max abs(sum p) = 1.6e-12 | 0 |
| 3 pions M=1 `<s23>` (150k) | 0.352782 +/- 0.000262 (ratio 0.99991); chi2/ndf = 13.0/19 | 0.352813 |
| massless 3-body `<x>`, `<x^2>` | 0.33387 +/- 0.00037, 0.16713 +/- 0.00031 | 1/3, 1/6 |
| 4-body massive M=2 (0.5,0.4,0.3,0.2): s12, s34 marginals | chi2/ndf = 21.2/19, 15.5/19; exact <s12> 1.221812 vs 1.222216 +/- 0.000607 | pass |
| n=4 equal-mass 0.3, M=2 means | 0.50090, 0.49998, 0.49931, 0.49981 +/- 0.00050 | equal |

**Cost** (same machine, `random.Random` / `numpy.default_rng`):

| sampler | n=3 pions M=1 | n=4 massless M=1 | n=4 massive (0.5,0.4,0.3,0.2) M=2 |
|---|---|---|---|
| production `nbody` (per call) | 24-33 us | 49-53 us | 51-58 us |
| new `correct_ps.nbody` (per call) | 8.7-11.3 us (2.9-3.3x faster) | 29-34 us (1.6x faster) | 54-56 us (~equal) |
| new `batch_ps.nbody_batch` (**per event**) | **0.63 us** | **3.4 us** | **6.4 us** |

So a per-call Python implementation cannot reach 5 us/call (the current one is 25-58 us); the
vectorised batch sampler does for n=3 and massless n=4 and is 50-80x cheaper per event than the
current code (1e7 decays: 6.3 s / 34 s / 64 s vs ~500 s). **Recommended fix:** replace
`nbody()` by `correct_ps.nbody`, and add `batch_ps.nbody_batch` for the scan loops (the
`multi` 4pi channel and K_L/tau 3-body channels are the users). Do not keep the "64 tries then
return None" branch - it silently drops phase-space points.

---

## F10 - `run_scan.py:419-429` lifetime integral with the 13-point ZGRID

**Analytic claim is correct.** With A=1 the code must give
`exp(-(max(1,z0)-z0)/L) - exp(-(2.5-z0)/L)`; re-deriving the review's table with the code's own
`run_scan.per_scalar_p` reproduces it *exactly*:

| z0 (m) | L=1/alpha (m) | code | exact A=1 | code/exact |
|---|---|---|---|---|
| 0 | 0.1 | 5.116e-5 | 4.540e-5 | **1.12694** |
| 0 | 0.3 | 3.5945e-2 | 3.5434e-2 | **1.01443** |
| 1.01 | 0.1 | 0.55473 | 1.00000 | **0.55473** |
| 1.01 | 1.0 | 0.72257 | 0.77463 | **0.93280** |
| 1.249 | 1.0 | 0.77615 | 0.71378 | **1.08737** |

The trapezoid error is O(h^2): with z0=0, A=1, the code/exact ratio is 1+1.3e-5 (alpha=0.1),
1+3.3e-4 (0.5), 1+1.3e-3 (1), 1+5.2e-3 (2), 1+3.2e-2 (5); doubling the node count divides the
error by 4.

**Realistic bias, however, is much smaller than the benchmark suggests.** In this run the
parent samples all sit at z=0 mm (`run_scan.load_sample` of bbar/minbias/higgs: all z=0.0000
mm); the only displacement is the parent B/K flight in `sample_meson_decay_S`. For
`run_mass(nsamp=300, ndecay=8)`:

| m_S | rows | z0 max (z0<1 m) | z0 inside 1..2.5 m | z0 > 2.5 m |
|---|---|---|---|---|
| 0.2802 | 1183 | 0.935 m (277 rows) | 17 rows (K only) but A(z) == 0 there | 889 rows (p=0) |
| 1.0 | 274 | 0.261 m | 0 | 0 |
| 3.5567 | 138 | 0.058 m | 0 | 0 |

Because z0 < 1 m for every row that has non-zero acceptance, the integral always starts at the
grid node z=1 and the review's large "birth point inside the window" errors (0.55-1.13) are not
realised. Signal-weighted `sum(p_code)/sum(p_exact PWL)` at the *actual* crossings from
`signal_scan.json`:

| m_S | long-lived edge 25 fb | long-lived edge 300 fb | short-lived edge 25 fb | short-lived edge 300 fb |
|---|---|---|---|---|
| 0.2802 | 1.00000000 | 1.00000000 | (band open; no edge) | (band open) |
| 1.0 | 1.00000001 | 1.00000000 | **1.025034** | **1.047865** |
| 3.5567 | 0.99999995 | 0.99999999 | **1.021985** | **1.040382** |

At the long-lived edge the bias is < 1e-7: for small alpha the trapezoid of f*A and the exact
integral of f*A_lin agree to O((alpha h)^2) because both are exact for constant A. At the
short-lived edge the code **overestimates p by 2.2-4.8 %** (the first grid cell is integrated by
the trapezoid as 0.5*alpha*h instead of (1-exp(-alpha*h))/alpha; with alpha*h up to 80 the
relative error reaches a few %). Effect on the final result: the short-lived boundary
theta^2 is shifted by at most ~2-5 %; the long-lived boundary (the lower edge, the
constraining one) is unaffected.

### Exact fix (validated): integrate the exponential exactly, A linear on each cell

`A(z)` is known only on ZGRID, so take it piecewise linear. On a cell [z1,z2] with
d=z2-z1, x=alpha*d, f1=exp(-alpha*(z1-z0)):

    w0 = -expm1(-x)/alpha                      # int_0^d exp(-alpha u) du
    w1 = (1-(1+x)*exp(-x))/alpha^2             # int_0^d u exp(-alpha u) du
    I_cell = alpha*f1*(A_k*w0 + s_k*((z1-z_k)*w0 + w1)),   s_k=(A_(k+1)-A_k)/h

plus a partial first cell [max(z0,1), z_(k0+1)] with z1=z0 for every row whose z0 lies in the
window. This formula equals the A=1 analytic answer to **2.2e-16** over a grid of 12 z0 values
(incl. z0 inside and outside the window) x 10 alpha values (incl. alpha=0, 1e-6, 500), is
exact for any node count, and its A-interpolation error is O(h^2) (verified on a smooth A).
Kept in this directory as `lifetime_exact_pwl()` in `f10b.py`/`f10d.py`; drop-in replacement
for `per_scalar_p`:

    def per_scalar_p_exact(ms, pz, z0, A, theta2):
        ctau = sw.HBAR_C / (sw.gamma_total(ms) * theta2)
        alpha = np.where(pz > 0, ms / np.maximum(ctau * np.maximum(pz, 1e-9), 1e-30), 0.0)
        return np.array([lifetime_exact_pwl(alpha[i], z0[i], A[i]) for i in range(len(pz))])

(loops over scalars; ~200k rows per mass point is cheap, ~0.2 s).

---

## F11 - `run_sensitivity.py:81-116, 212, 254`: level = 2.3 on a per-fb signal

`signal_per_fb_per_channel` is per fb^-1 (`make_curves.py:92-97` multiplies by 25/300 before
crossing at S95=2.3), but `crossing_errors(scan, args.s)` searches crossings of the **per-fb**
curve at level=2.3, i.e. it optimises/validates the **1 fb^-1** contour, not the 25/300 fb^-1
contours that are actually drawn. Reproduced from `outputs/data/signal_scan.json`:

| level searched | long-lived-edge rel MC error: median | max | > 5 % |
|---|---|---|---|
| 2.3 (code: L = 1 fb^-1) | 4.824 % | 4.985 % | 0/67 |
| 2.3/25 (L = 25 fb^-1) | 4.779 % | 5.837 % | **7/67** |
| 2.3/300 (L = 300 fb^-1) | 4.779 % | 6.356 % | **10/67** |

Masses above 5 % at the long-lived edge and the parent nsamp that would be needed
(`nsamp_new = nsamp * (rel/0.05)^2`):

* 25 fb^-1 (7): 3.71396 (5.61 %, 1400 -> 1760), 3.75790 (5.84 %, 1300 -> 1772),
  4.34457 (5.31 %, 1700 -> 1920), 5.08225 (5.14 %, 1900 -> 2007), 8.13552 (5.08 %, 2000 -> 2064),
  13.02309 (5.03 %, 1900 -> 1926), 24.38668 (5.09 %, 1900 -> 1968).
* 300 fb^-1 (10): the seven above plus 3.17488 (5.37 %, 1800 -> 2074),
  3.55670 (5.62 %, 1500 -> 1892), 11.13280 (5.00 %, 2100 -> 2102);
  3.71396 (6.07 % -> 2066), 3.75790 (6.36 % -> 2101).
  All needed nsamp are <= ~2100, below the 25000 cap - this edge is close.

The much larger problem the code never checks: at the **short-lived** edge (all 46 masses with
two edges) rel = 27.1 % median (max 98.9 %) at 25 fb^-1 and 38.0 % median (max 96.4 %) at
300 fb^-1, i.e. **46/46 masses exceed 5 %**. Reaching 5 % there needs nsamp ~5.5e4 (median,
25 fb) / 1.2e5 (median, 300 fb) and up to 1.2e6 - 44/46 (25 fb) and 46/46 (300 fb) masses
exceed the `--nsamp-max 25000` cap.

**Fix:** pass the luminosity, e.g. `level = args.s / L` and call
`crossing_errors(scan, args.s/L)` for each L in {25, 300} (the same L list as
`make_curves.LUMI`), and validate **both** the long- and short-lived edge (the boundary
*position* error from crossing s +/- ds is a different, smaller quantity on the steep wall -
the code already prints it; it should be the acceptance criterion at the short-lived edge).

---

## F12 - `bc5_lib.py:180` CHANNELS vs `scalar_widths.widths()`

Confirmed exactly as claimed: `widths()` returns {'e','mu','tau','gamma','pipi','KK','multi'}
for m<=2 GeV and {'e','mu','tau','gamma','ss','cc','bb','gg'} above, while `CHANNELS` lacks
`gamma`. `channel_branching()` (bc5_lib.py:425-433) divides by the *full* total, so
`sum(BR over CHANNELS) = 1 - BR(gamma)`; in `acceptance_profile` (470-473) the cumulative array
ends at `1-BR(gamma)` and `idx = min(searchsorted(cum,u), 9)` sends the whole gamma probability
to the last channel, `'gg'`.

BR(gamma) and the visible-rate consequence:

| m_S (GeV) | BR(gamma) | sum(BR over CHANNELS) | charged daughters of the 'gg' fallback | visible-rate change |
|---|---|---|---|---|
| 0.1 | 1.371 % | 0.986294 | none (empty list) | **0** |
| 0.18869 (max) | 8.919 % | 0.91081 | none (nearest lib node 0.3021 is >20 % away; analytic e_jet = 0.094 < 0.35) | **0** |
| 0.2118 | 1.196 % | 0.98804 | none (same reason) | **0** |
| 0.25821 | 4e-5 | ~1 | yes (lib node 0.3021) | < 5e-5 |
| 0.2802 | 2e-5 | ~1 | yes | < 2e-5 |
| >= 0.3 | <= 4e-5 | >= 0.99996 | yes | <= 3e-5 |

Per-decay acceptance MC with the real kernel (`f12_effect.py`, eta_S=3, pz=50 GeV, z=1.5 m,
600 decays/channel) gives A_mix = 0.70/0.97/0.65/0.54 at m = 0.189/0.212/0.5/0.988 GeV and
dA/A = BR(gamma)*P_gg/A_mix <= 3e-5 for m >= 0.28 GeV, exactly 0 below. So the **semantic bug
is real** (a channel with zero charged daughters is silently merged into 'gg'), but its effect
on the visible rate is *zero below ~0.25 GeV and <=3e-5 above* - the review's "except below
~0.2 GeV" qualifier is not realised: at the masses where BR(gamma) is large the gg fallback is
empty, so the visible rate is unchanged. **Fix:** add `'gamma'` to `CHANNELS` and let
`s_decay_products` return `[]` for it (it already does, line 331); or build the cumulative
array from `widths()` instead of `CHANNELS` so nothing can leak into the last channel.

---

## Bottom line for the final result

* **F03 is the only claim that changes physics numbers materially**: the production sampler
  distorts massive 3-body and all 4-body Dalitz distributions (s23 shape chi2 8486/19; 4-body
  energy means 0.33/0.22/0.22/0.22 instead of 0.25). It feeds the `multi` 4pi channel,
  K_L -> pi pi pi0 / pi e nu, tau -> 3pi and the quartic X -> X'SS 3-body generator. The
  corrected sampler in this directory passes all the same tests and is not slower (batch:
  ~50x cheaper per event).
* **F10 is real but essentially irrelevant at the long-lived edge** (bias < 1e-7) and only
  2-5 % at the short-lived edge; still worth the 20 lines because the exact formula is free.
* **F11 is real**: the QA was certifying the 1 fb^-1 contour; at 25/300 fb^-1 7/67 and 10/67
  masses miss the 5 % target, and the short-lived edge (46/46) is far from it - but that is a
  *bookkeeping/acceptance* failure, not a central-value failure.
* **F12 is real but numerically harmless** (< 1e-4).
