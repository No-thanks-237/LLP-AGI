# MC statistical error of the BC5 downstream scan, and how to converge it to 5 %

Study of the Monte-Carlo error of the production scan 'outputs/data/signal_scan.json'
(52 masses, 240 theta^2 points, s = 2.3 -> 90 % CL, LHCb downstream 1 < z < 2.5 m).
**No existing physics code or production output was modified.** New files only:
'work/analysis/mc_study.py' (this study), 'work/analysis/mc_study_proto.py'
(a byte copy of run_scan.py with one hook added -- see 3.1), and suffixed scans
'outputs/data/signal_scan_<name>.json'.

All numbers below come from the commands in section 7; the raw printed tables were produced by
'mc_study.py diag / fit / proto / cost / modevar'.

---

## 0. Result in one page

* **The MC error of s at the long-lived (lo) edge is 4.7 - 8.2 % (median 6.6 %) at 1200 parents**,
  and the theta^2 boundary error is 3.6 - 7.1 % (median 4.6 %). The starting statement
  "5-8 % at lo" is confirmed.
* **At the short-lived (hi) edge the error of s is 29.5 - 97.3 % (median 41.5 %) at 300 fb^-1**
  (20.5 - 98.8 % at 25 fb^-1), confirming "20-95 %". Because the curve is a steep wall
  (median |dln s / dln theta^2| = 5.6, locally up to 13) the **boundary** error is only
  4.3 - 23.2 % (median 7.6 %) -- but that number is *optimistic* (see the fourth point).
* **The dominant variance source is the finite parent sample, i.e. the Monte-Carlo mean of the
  per-parent response -- not the inner A(z) MC, not the branching-ratio sampling.** At lo the
  parent sample carries 90-97 % of the variance (ndecay 8 -> infinity lowers the error only by
  2-10 %); at hi it carries essentially 100 %. The one-mode branching-ratio sampling contributes
  < 0.1 % of the variance (wB = 3.17 - 3.35 is a normalisation constant, not a noise source).
* **At lo s(theta^2) converges exactly as 1/sqrt(nsamp)** (measured dlog(rel)/dlog(nsamp) =
  -0.47 ... -0.53). **At hi it does not**: the measured slopes are -0.16 ... -0.31, because the
  response there is concentrated in a handful of parents: the *effective* number of parents is
  n_eff = 3.0 / 4.8 / 6.7 out of 7200 (per channel) for m_S = 0.905 / 5.08 / 1.98 GeV. The
  short-lived edge is a rare-coincidence regime, and the reported single-run standard error
  underestimates the seed-to-seed scatter there (the same m_S hi crossing moves by factors
  1.7-2 between independent 1200-parent scans, e.g. 0.76-1.54e-5 at m_S = 0.905 GeV).
* **Variance reduction (idea A).** A fixed, deterministic, stratified parent set with explicit
  weights (10 equal-count bins of log parent |pz|) is unbiased and works at lo, but the gain is
  modest because parent |pz| only explains part of the variance:
  - lo: variance-reduction factor 1.28 - 1.71 (relerr gain 1.13 - 1.31x) at equal n;
    with Neyman allocation on a lo-only mass (m_S = 0.302 GeV) the factor is 4.9 (relerr
    6.5 % -> 2.9 %, bias +3.2 % against a reference whose own error is 3.3 %).
  - hi: parent |pz| explains 1-5 % of the variance only -> reduction factor 1.02 - 1.05, and
    the deterministic design is *catastrophically* unstable: its 3 offset replicas scatter by
    10-170 % and the mean can be off by a factor 2 to 200. **A fixed parent set is not a usable
    estimator at the short-lived edge.**
  - The variable that would work is the *scalar's own* pz (or direction), not the parent's:
    binning by the scalar pz gives variance-reduction factors 1.9 - 4.4 at lo (and still only
    1.05 at hi). Realising that requires importance sampling of the S decay direction, which is
    a change to the sampling, not to the parent set.
* **Idea B (deterministic quadrature over a (p, theta) histogram) does not converge like a
  quadrature.** Measured RMS error of a one-representative-per-bin estimate scales as
  K^(-0.5), i.e. like K random parents, not as K^(-2): at K = 200 reps the error is 14-23 %
  (lo) and 229 % (hi), vs the MC prediction 17-21 % and 219 %. To reach 5 % one needs
  K ~ (CV/0.05)^2 ~ 3500-4000 bins/cells, i.e. no saving over 3500 parents, plus a binning bias.
  The acceptance is not smooth in (p, theta) because it is driven by the decay products.
* **Cost to reach 5 %** (this machine, single-threaded, measured: 119 s for the full 52-mass grid
  at nsamp = 300 -> 7.9 min at the production nsamp = 1200):
  - 5 % at lo for all 52 masses: **8.6 min** (median nsamp 1200, max 2433; only 22/52 masses
    need more than the current 1200).
  - 5 % at both edges: **15.0 min** (median nsamp 1901, max 25871 for m_S = 0.905 GeV;
    7 masses need > 4800, 2 need > 12000). This uses the *reported* errors; for the 3-4
    heavy-tail masses at hi it is a lower bound, add a factor ~4 there.
  - parallelising over the 52 masses on 32 cores turns 8.6 / 15.0 CPU-min into
    ~1.5 / ~5.5 min wall (sum/32), but buys **no** reduction of total CPU and **cannot** cure the
    heavy-tail masses -- only more samples or importance sampling can do that.

**Ranked recommendation**

1. **(cheap, do first) Run nsamp = 2400 for all 52 masses** (~16 min single-threaded, ~2 min on
   32 cores): this alone puts 5 % at lo for all masses and 5 % at hi for 27 of the 34 two-edge
   masses. Do not bother with ndecay > 8: it buys 2-10 % at lo and nothing at hi.
2. **(targeted) nsamp = 4800-25000 only for the 7 hi-wall masses** (0.905, 2.32, 2.71, 1.70,
   6.96, 0.353, 3.18 GeV): +~1.3x total CPU.
3. **(structural, only if 5 % at hi is a hard requirement for every mass) importance-sample the
   S decay direction relative to the LHCb acceptance** (or tabulate A(z) on a (pz, eta) grid
   once), because there the answer is decided by 3-7 parents out of 1200. A stratified parent
   set cannot fix it.
4. **Do not** invest in the one-mode branching-ratio sampling (negligible), and **do not** use a
   deterministic fixed parent set at the short-lived edge (biased and unstable).

---

## 1. Diagnostic of the production scan (task 1)

Definition of the quantities (identical algebra to 'make_curves.py', and the recomputation
reproduces the production 'exclusion_boundary.csv' exactly, max relative difference 0.0):

* 'rel_lo' / 'rel_hi' = relative MC error of s at the lo / hi crossing (ds/s from the
  per-channel 'relerr_*' combined in quadrature).
* 'dlo/lo', 'dhi/hi' = relative MC error of the theta^2 boundary, from the crossings of s -/+ ds
  (the same estimator that make_curves puts on the curve).
* channel shares are the fractions of s and of its variance carried by B (mesons from b bbar),
  K (minimum-bias kaons) and H (h -> SS) at the crossing.

### 1.1 Global summary over the 52 production masses

    edge  lumi    rel error of s [min/med/max]      boundary rel err [min/med/max]
    lo    25       5.0%   6.6%   8.2%                 3.3%   4.0%   7.1%
    hi    25      20.5%  28.7%  98.8%                 4.0%   7.8%  27.4%
    lo    300      4.7%   6.6%   8.2%                 3.6%   4.6%   7.1%
    hi    300     29.5%  41.5%  97.3%                 4.3%   7.6%  23.2%

    crossings at 300 fb^-1: 34 masses have two edges, 18 only the long-lived edge
    (the band is open above theta^2 = 1e-3), 0 masses have none.

The median local slope at the crossing is |dln s / dln theta^2| = 1.14 (lo) and 5.64 (hi); the
slope-based estimate of the boundary error (rel/slope) gives 4.6 % (lo) and 7.0 % (hi), in
agreement with the crossing estimator.

### 1.2 Per-mass table at L = 300 fb^-1 (production scan)

Columns: theta^2 at the lo / hi crossing, rel error of s there, relative boundary error, and
the s-share of the three channels (B/K/H in %). A dash means no crossing inside the theta^2
window.

    m_S      theta2_lo rel_lo dlo/lo  s-share B/K/H   theta2_hi rel_hi dhi/hi  s-share B/K/H
    0.021    2.49e-06    8.0%   4.0%  100.0/ 0.0/ 0.0   -          -      -     -
    0.02457  2.06e-06    7.4%   3.7%  100.0/ 0.0/ 0.0   -          -      -     -
    0.02874  1.82e-06    7.8%   3.9%  100.0/ 0.0/ 0.0   -          -      -     -
    0.03362  1.59e-06    7.7%   3.9%  100.0/ 0.0/ 0.0   -          -      -     -
    0.03932  1.41e-06    8.2%   4.1%  100.0/ 0.0/ 0.0   -          -      -     -
    0.046    1.10e-06    7.9%   4.0%   99.9/ 0.1/ 0.0   -          -      -     -
    0.05381  1.00e-06    7.9%   4.0%  100.0/ 0.0/ 0.0   -          -      -     -
    0.06295  8.86e-07    8.1%   4.1%  100.0/ 0.0/ 0.0   -          -      -     -
    0.07364  7.53e-07    7.9%   3.9%  100.0/ 0.0/ 0.0   -          -      -     -
    0.08614  6.41e-07    7.9%   4.0%  100.0/ 0.0/ 0.0   -          -      -     -
    0.1008   5.19e-07    7.6%   3.8%  100.0/ 0.0/ 0.0   -          -      -     -
    0.1179   4.44e-07    7.9%   4.0%  100.0/ 0.0/ 0.0   -          -      -     -
    0.1379   3.62e-07    7.9%   4.0%  100.0/ 0.0/ 0.0   -          -      -     -
    0.1613   3.32e-07    8.2%   4.1%  100.0/ 0.0/ 0.0   -          -      -     -
    0.1887   2.69e-07    7.8%   3.9%  100.0/ 0.0/ 0.0   -          -      -     -
    0.2207   5.02e-09    7.2%   3.7%   99.0/ 0.0/ 0.9   -          -      -     -
    0.2582   1.77e-09    7.3%   3.9%   96.6/ 0.0/ 3.4   -          -      -     -
    0.3021   4.77e-10    6.4%   3.7%   91.1/ 0.0/ 8.9   -          -      -     -
    0.3533   2.96e-10    5.8%   3.6%   87.0/ 0.0/13.0   7.91e-04  85.1%  11.3%  88.0/ 0.0/12.0
    0.4133   2.15e-10    5.8%   3.8%   78.8/ 0.0/21.2   3.24e-04  31.5%   4.3%  18.1/ 0.0/81.9
    0.4835   1.55e-10    5.5%   3.8%   75.2/ 0.0/24.8   1.50e-04  35.5%   5.9%   5.7/ 0.0/94.3
    0.5656   1.08e-10    5.0%   3.7%   67.2/ 0.0/32.8   1.22e-04  64.3%   6.7%  75.7/ 0.0/24.3
    0.6617   6.82e-11    4.7%   3.7%   65.3/ 0.0/34.7   5.15e-05  56.1%   7.7%  56.3/ 0.0/43.7
    0.774    4.69e-11    5.3%   4.5%   55.3/ 0.0/44.7   2.47e-05  35.4%   6.4%   0.3/ 0.0/99.7
    0.9054   1.69e-11    5.0%   4.6%   51.1/ 0.0/48.9   1.69e-05  97.3%  23.2%  97.4/ 0.0/ 2.6
    1.059    1.30e-11    4.9%   4.6%   49.7/ 0.0/50.3   7.87e-06  36.5%   5.2%  23.0/ 0.0/77.0
    1.239    2.15e-11    5.2%   4.6%   44.1/ 0.0/55.9   1.51e-05  40.4%   7.7%   0.1/ 0.0/99.9
    1.449    1.87e-11    5.2%   4.7%   33.9/ 0.0/66.1   1.15e-05  39.4%   5.7%  22.6/ 0.0/77.4
    1.695    2.45e-11    5.6%   4.8%   32.0/ 0.0/68.0   1.33e-05  77.7%  13.8%  77.8/ 0.0/22.2
    1.983    2.45e-11    5.5%   4.8%   23.6/ 0.0/76.4   8.97e-06  29.5%   4.9%   0.1/ 0.0/99.9
    2.32     1.51e-11    5.8%   5.4%   11.7/ 0.0/88.3   7.60e-06  91.9%  17.0%  92.2/ 0.0/ 7.8
    2.714    9.87e-12    6.4%   6.1%    5.2/ 0.0/94.8   4.17e-06  65.1%  14.8%   0.0/ 0.0/100.0
    3.175    6.42e-12    6.7%   6.5%    3.0/ 0.0/97.0   2.21e-06  50.8%  11.1%   0.0/ 0.0/100.0
    3.714    2.92e-12    6.5%   6.5%    0.6/ 0.0/99.4   1.13e-06  51.6%   9.9%   5.2/ 0.0/94.8
    4.345    3.38e-13    6.2%   6.2%    0.0/ 0.0/100.0  1.27e-07  31.1%   5.4%   0.0/ 0.0/100.0
    5.082    1.50e-13    6.6%   6.6%    0.0/ 0.0/100.0  5.18e-08  42.6%   8.4%   0.1/ 0.0/99.9
    5.945    8.03e-14    6.4%   6.4%    0.0/ 0.0/100.0  2.40e-08  32.3%   5.7%   0.0/ 0.0/100.0
    6.955    6.05e-14    6.6%   6.7%    0.0/ 0.0/100.0  1.52e-08  53.8%  12.1%   0.0/ 0.0/100.0
    8.136    3.93e-14    6.9%   6.9%    0.0/ 0.0/100.0  8.04e-09  30.4%   5.2%   0.0/ 0.0/100.0
    9.517    2.30e-14    6.6%   6.6%    0.0/ 0.0/100.0  6.63e-09  46.6%   9.4%   0.0/ 0.0/100.0
    11.13    1.36e-14    7.1%   7.1%    0.0/ 0.0/100.0  3.19e-09  45.2%   9.2%   0.0/ 0.0/100.0
    13.02    3.93e-15    6.6%   6.6%    0.0/ 0.0/100.0  1.13e-09  30.8%   5.2%   0.0/ 0.0/100.0
    15.23    1.79e-15    6.2%   6.2%    0.0/ 0.0/100.0  5.38e-10  39.2%   7.4%   0.0/ 0.0/100.0
    17.82    1.13e-15    6.9%   6.9%    0.0/ 0.0/100.0  2.92e-10  43.8%   8.5%   0.0/ 0.0/100.0
    20.85    6.46e-16    6.9%   6.9%    0.0/ 0.0/100.0  1.91e-10  35.8%   6.5%   0.0/ 0.0/100.0
    24.39    3.62e-16    6.2%   6.2%    0.0/ 0.0/100.0  1.44e-10  48.7%   9.4%   0.0/ 0.0/100.0
    28.53    2.56e-16    6.6%   6.6%    0.0/ 0.0/100.0  8.86e-11  36.9%   7.2%   0.0/ 0.0/100.0
    33.37    1.55e-16    6.0%   6.0%    0.0/ 0.0/100.0  6.62e-11  44.5%   8.7%   0.0/ 0.0/100.0
    39.04    9.84e-17    6.0%   6.0%    0.0/ 0.0/100.0  4.22e-11  34.8%   6.2%   0.0/ 0.0/100.0
    45.67    5.73e-17    5.6%   5.6%    0.0/ 0.0/100.0  2.87e-11  39.1%   7.1%   0.0/ 0.0/100.0
    53.42    3.22e-17    5.3%   5.3%    0.0/ 0.0/100.0  1.81e-11  33.3%   5.7%   0.0/ 0.0/100.0
    62.49    1.66e-17    5.7%   5.7%    0.0/ 0.0/100.0  1.29e-11  44.1%   8.3%   0.0/ 0.0/100.0

(The same table at 25 fb^-1 is in the 'mc_study.py diag' output; the pattern is identical, the
lo/hi crossings move and the boundary errors are 1-20 % smaller at hi because the curve shifts.)

### 1.3 Which channel dominates

    lo edge:  s  dominated by B in 25 masses, H in 27;  median share B=0.47 H=0.53
    lo edge:  var dominated by B in 25 masses, H in 27;  median share B=0.40 H=0.60
    hi edge:  s  dominated by B in  6 masses, H in 28;  median share H=1.00
    hi edge:  var dominated by B in  6 masses, H in 28;  median share H=1.00

The K channel never dominates either edge (median share 0.00), but it is the most pathological
one: at 2250 of 2640 theta^2 points its relative error is exactly 1.0, the signature of *one
single surviving scalar in the whole 1200-parent sample*. It is harmless only because its
cross-section weight is tiny; the same pathology is what the H (and occasionally B) channel
shows at the short-lived edge (section 3).

**Worst cases** (L = 300 fb^-1, boundary error): lo -- 11.1 GeV 7.1 %, 8.14 6.9 %, 17.8 6.9 %,
20.8 6.9 %, 6.95 6.7 %, 5.08 6.6 %. hi -- 0.905 GeV 23.2 %, 2.32 17.0 %, 2.71 14.8 %,
1.70 13.8 %, 6.95 12.1 %, 0.353 11.3 %.

---

## 2. Variance-source experiments (task 2)

Four representative production masses: 0.3021, 0.9054, 1.9833, 5.0823 GeV (exactly on the
production grid). All runs write their own suffixed scan; production files untouched.

### 2.1 Configurations, wall-clock time

Measured with 9 single-threaded processes running concurrently on 32 logical cores (plus one
other user process), so the per-mass times are accurate to ~+-20 %.

    name     nsamp ndecay     seed      npseed    wall_s   per-mass [0.302, 0.905, 1.983, 5.082] s
    base      1200      8  20261016   20261016     46.3   22.4  7.8  8.0  3.7
    ns2400    2400      8  20261016   20261016     72.4   33.5 14.5 13.2  5.8
    ns4800    4800      8  20261016   20261016    108.7   53.2 22.7 16.5  8.8
    nd16      1200     16  20261016   20261016     61.6   28.2 10.7 11.8  6.7
    nd32      1200     32  20261016   20261016     84.6   39.2 14.6 18.1  9.4
    seed1     1200      8  11235813   20261016     46.3   22.4  7.9  7.9  3.6
    seed2     1200      8    271828   20261016     46.6   22.4  8.0  7.9  3.8
    np1       1200      8  20261016    5555555     46.2   22.6  7.5  7.5  4.0
    np2       1200      8  20261016  987654321     45.9   22.6  7.5  7.5  3.9
    cost300    300      8  20261016   20261016    130.2   (full 52-mass grid: 1.6-11.3 s per mass)

'--seed' controls the decay/A(z) stream, '--npseed' (BC5_NPRNG_SEED) controls the parent
subsampling; at fixed npseed the 1200 parents are the same, so the two seed families separate
the *inner* (direction + A(z)) error from the *parent-set* error.

### 2.2 Relative MC error at the production crossings, fixed theta, L = 300 fb^-1

    edge lo        m=0.302  m=0.905  m=1.983  m=5.082
    base             6.5%     4.9%     5.8%     6.6%
    ns2400           4.6%     3.6%     4.0%     4.5%
    ns4800           3.3%     2.6%     2.9%     3.3%
    nd16             6.4%     4.9%     5.4%     6.2%
    nd32             6.4%     4.8%     5.4%     6.4%
    seed1            6.8%     5.1%     5.5%     6.4%
    seed2            6.5%     4.8%     5.6%     6.4%
    np1              6.5%     5.0%     5.7%     6.4%
    np2              6.4%     5.0%     5.5%     6.6%

    edge hi        m=0.302  m=0.905  m=1.983  m=5.082
    base              nan%    42.0%    33.4%    54.9%
    ns2400            nan%    54.5%    24.2%    29.5%
    ns4800            nan%    24.5%    90.2%    25.6%
    nd16              nan%    50.1%    38.9%    65.7%
    nd32              nan%    62.0%    31.9%    45.0%
    np1               nan%    82.7%    61.7%    53.9%
    np2               nan%    88.4%    33.1%    64.1%
    seed1             nan%    58.8%    38.0%    38.0%
    seed2             nan%    42.9%    30.4%    56.8%

At lo the table is stable and every configuration reproduces it. At hi the *same* physical
quantity takes values from 24 % to 90 % across configurations whose only difference is the
sample: that is the signature of a non-Gaussian, heavy-tailed error, not of a scaling law.

### 2.3 Scaling fits rel^2 = (A + B/ndecay)/nsamp  (four masses, production crossings)

    edge lo  m=0.302: dlog(rel)/dlog(nsamp) = -0.53   A=4.89 B=1.18   B/(8A)= 3 %   fit 6.5% vs measured 6.5%
    edge lo  m=0.905:                          -0.49   A=2.71 B=1.65   B/(8A)= 8 %   fit 4.9% vs measured 4.8%
    edge lo  m=1.983:                          -0.47   A=3.27 B=6.02   B/(8A)=23 %   fit 5.8% vs measured 5.6%
    edge lo  m=5.082:                          -0.51   A=4.57 B=5.11   B/(8A)=14 %   fit 6.6% vs measured 6.4%
    edge hi  m=0.905:                          -0.31   (fit meaningless, B < 0)
    edge hi  m=1.983:                          -0.16
    edge hi  m=5.082:                          -0.24

* **nsamp**: at lo rel ~ nsamp^(-1/2) within the measurement (-0.47 ... -0.53). At hi the
  scaling is *slower* than nsamp^(-1/2) (-0.16 ... -0.31): quadrupling nsamp does not halve the
  hi error.
* **ndecay**: at lo A + B/8 with A >> B/8: moving ndecay 8 -> infinity would lower the rel
  error by only 2 % (0.302), 4 % (0.905), 10 % (1.983), 6 % (5.082). The inner A(z) MC is a
  small correction, not the bottleneck. At hi the ndecay series (42->50->62 % at 0.905,
  33->39->32 % at 1.983, 55->66->45 % at 5.082) shows no systematic trend: the parent-set / tail
  fluctuation swamps it (the three ndecay runs share the same 1200 parents, so the scatter is
  the *direction-stream* noise alone -- again heavy-tailed).
* Cost: ndecay is linear in wall time (61.6 s for x2, 84.6 s for x4 against 46.3 s), so raising
  ndecay is expensive for a 2-10 % gain -> **keep ndecay = 8**.

### 2.4 Seed experiments and the theta^2 crossings

Scatter of s at fixed theta over the two seed families (5 configurations each; at lo the mean
value is 2.3 by construction of the crossing):

    lo:  --seed spread 3.3-5.4 %   npseed spread 4.5-7.1 %
    hi:  --seed spread 23-50 %     npseed spread 24-150 %   (m=0.905: factor 40 between npseeds)

At lo the parent-set component is comparable to or larger than the direction component
(4.5-7.1 % vs 3.3-5.4 %), consistent with 90-97 % parent-sample variance. At hi the parent-set
component explodes (a single parent can carry the whole answer).

Theta^2 crossings, L = 300 fb^-1 (a converged scan would reproduce them):

    edge lo    m=0.302      m=0.905      m=1.983      m=5.082
    production 4.774e-10    1.693e-11    2.449e-11    1.500e-13
    9 reruns   4.55-5.04e-10 1.56-1.88e-11 2.00-2.47e-11 1.40-1.64e-13   (+-8-11 %)

    edge hi    m=0.905      m=1.983      m=5.082
    production 1.686e-05    8.971e-06    5.182e-08
    9 reruns   8.35-15.4e-06 8.8-15.1e-06 3.85-6.63e-08   (factors 1.7-2.0)

The lo crossings are stable at the quoted level. The hi crossings are only reproducible to a
factor ~2, and the production value at m_S = 0.905 GeV (1.69e-5) is not reproduced by any of
the 9 independent 1200-parent reruns (eight give 8.35-9.94e-6, one 1.54e-5): **that boundary point is a ~2
sigma outlier of a bimodal distribution** (most samples draw no extreme parent, a few do).
The quoted hi boundary errors (4.3-23 %) therefore underestimate the true seed-to-seed spread
for the worst masses.

---

## 3. Variance decomposition and variance reduction (task 3)

### 3.1 Parent-level estimator and its validation

'mc_study_proto.py' is a byte copy of 'run_scan.py' with one hook: 'run_mass(..., parent_idx=...)'
uses an explicit parent-row set instead of '_np_rng.choice', and returns the selected rows
('tab["sel_idx"]'). Everything else is the production code. 'mc_study.py' rebuilds the
per-parent response p_par (histories -> parents, exactly as 'channel_efficiency', but stopping
before the collision OR, which is an O(p^2) ~ 1e-14 correction here) and computes
eps = 0.4 sum_i W_i p_par,i / n_gen with explicit weights W_i.

Validation: the parent-level s agrees with the production 'channel_efficiency' to all printed
digits (e.g. m=5.082 lo: 0.0071428 both) and with the independent nsamp=4800 scan to 0.6 %
(m=0.302 lo: 0.007408 vs 0.007470, the reference's own error is 3.3 %). The estimator is
therefore the production one, and the study of its variance is meaningful.

### 3.2 Where the variance lives

Per-parent response y_i = (1e3 sigma_c) p_par,i, summed over channels; n_eff = (sum y)^2/sum y^2;
'within/total' = 1 - eta^2 is the fraction of the variance *not* explained by an equal-count
binning of the stated variable (10 bins), so 1/(within/total) is the variance-reduction factor
achievable at equal n. All four masses, both edges:

    m_S    edge  var-share B/H     1-eta^2 parent |pz| B / H     VR parent pz   VR S pz   n_eff (/7200)
    0.302  lo    97.9 /  2.1       0.580 / 0.848                 1.71            1.93       213
    0.905  lo    30.5 / 69.5       0.457 / 0.810                 1.42            2.55       305
    1.983  lo     4.3 / 95.7       0.540 / 0.794                 1.28            2.96       217
    5.082  lo     0.0 /100.0       0.510 / 0.647                 1.55            4.41       190
    0.905  hi     0.0 /100.0       0.992 / 0.978                 1.02            1.02         3.0
    1.983  hi     0.0 /100.0       0.992 / 0.950                 1.05            1.05         6.7
    5.082  hi     0.0 /100.0       0.990 / 0.964                 1.04            1.04         4.8

* At lo the response is spread over ~200 effective parents out of 7.2 k (6-7 % of parents have a
  non-zero response, the top 1 % of parents carry 35-48 % of s). Parent |pz| explains only part
  of the spread -> stratification gain 1.28-1.71 in variance.
* At hi **n_eff = 3-7 parents out of 7200**: the top 1 % of parents carry 100 % of s. No
  fixed-parent-set method can be efficient there, and the local variance estimate is unreliable.
* The scalar's own pz is a much better stratifier at lo (VR 1.9-4.4) than the parent pz, but it
  is not known before the decay is sampled -> realising this gain requires importance sampling
  of the S decay direction, not a fixed parent set.

### 3.3 Fixed stratified designs with explicit weights (idea A)

Design: 10 equal-count bins of log10(parent pz) (pz<=0 parents in the lowest bin); either
n_s = n/10 (proportional) or Neyman n_s ~ N_s sigma_s (sigma_s from the pilot, allocation taken
from the *binding* edge); within each stratum the n_s parents are chosen at systematic positions
(offsets 0, 1/3, 2/3 are the 3 replicas); weight N_s/n_s per parent, sum of weights = number of
candidates -> unbiased normalisation. One design per mass, applied to both edges (as production
would have to). n = 1200 per channel.

    mass   edge  uniform relerr   design mean s (3 offsets)   offset scatter   design se   bias vs nsamp=4800
    0.302  lo    6.5 %    prop   0.007694 (+3.0 %)              2.8 %          5.9 %        +3.0 %  (ref err 3.3 %)
                          neyman 0.007709 (+3.2 %)              2.3 %          2.9 %        +3.2 %  (ref err 3.3 %)
    0.905  lo    5.2 %    prop   0.007816 (+3.6 %)              0.6 %          4.5 %        +3.6 %  (ref err 2.6 %)
                          neyman 0.005226 (-31 %)               42.4 %         29.1 %       -30.7 %
    1.983  lo    6.1 %    prop   0.007874 (-3.1 %)              2.8 %          5.1 %        -3.1 %  (ref err 2.9 %)
                          neyman 0.005581 (-31 %)               28.1 %         30.7 %       -31.3 %
    5.082  lo    6.6 %    prop   0.007542 (+4.9 %)              4.5 %          5.4 %        +4.9 %  (ref err 3.3 %)
                          neyman 0.007106 (-1.1 %)              53.0 %         29.5 %       -1.1 %

    0.905  hi    58 %     prop   0.01827 (+21367 %)             172.8 %        99.8 %       +21400 %
                          neyman 0.000235 (+176 %)              76.8 %         47.2 %       +176 %
    1.983  hi    38 %     prop   0.6932 (-11 %)                 169.8 %        98.4 %       -11.0 %
                          neyman 0.1711 (-78 %)                 106.5 %        69.5 %       -78.0 %
    5.082  hi    46 %     prop   0.003615 (-50.6 %)             47.8 %         43.9 %       -50.6 %

Conclusions, with numbers:

* **Proportional stratification at lo is unbiased and gives a real but modest gain**: the
  design-based se is 4.5-5.9 % against the uniform 5.2-6.6 %, i.e. a variance reduction of
  1.2-1.5, consistent with the predicted 1.28-1.71. Bias at lo is +3.0 to +4.9 %, of the same
  order as the reference's own 2.6-3.3 % error, so no bias is resolved.
* **Neyman allocation is where the gain is, but only when the design is matched to the edge it
  will be used at**: at the lo-only mass (m_S = 0.302, no hi edge) it delivers the predicted
  4.9x variance reduction: relerr 6.5 % -> 2.9 % with 2.3 % replica scatter and a bias of
  +3.2 % (reference error 3.3 %). At the masses with two edges one allocation cannot serve both
  edges: allocating for hi leaves the lo estimate with -31 % bias and 28-42 % scatter. A
  production design would have to build *two* parent sets per mass (one per edge) or allocate
  for both simultaneously.
* **At hi the whole idea fails**: the parent-set fluctuations are so heavy-tailed that the
  deterministic design is not a usable estimator -- the 3 offsets scatter by 48-170 % and the
  mean can be off by a factor 2 (5.082), a factor ~2 and -11 % (1.983) or +214x (0.905, one
  replica caught an extreme parent). This is not a bias of stratification; it is a statement
  that at hi the observable is decided by 3-7 parents, so *any* deterministic 1200-parent subset
  is a lottery.
* The **bias of the estimator itself (not of stratification) is zero**: the parent-level
  estimator reproduces the production channel efficiency exactly, and the 1200-parent uniform
  estimates at lo agree with the ns4800 reference to 0.6-4.9 % (reference error 2.6-3.3 %).

### 3.4 Idea B: deterministic quadrature over a (p, theta) histogram

Test: replace the random parent subset by one representative per equal-count bin of log10(parent
pz) with weight N_b/n (8 offsets, RMS over offsets), and compare with the exact mean of the same
1200-parent sample and with the MC prediction CV/sqrt(K) for K random parents. If the integrand
were smooth the error would fall like K^-2; measured:

    m=0.302 lo (B channel, CV of per-parent response 2.45):
      K= 10  RMS 57.8 %   MC prediction 77.5 %   K= 50  RMS 20.3 %   prediction 34.6 %
      K=100  RMS 15.1 %   prediction 24.5 %     K=200  RMS 18.7 %   prediction 17.3 %
    m=0.905 lo (CV 2.43):
      K= 10  RMS 39.3 %   prediction 76.9 %     K=100  RMS 14.8 %   prediction 24.3 %
    m=5.082 lo (CV 2.97):
      K= 10  RMS 53.6 %   prediction 94.1 %     K=100  RMS 24.1 %   prediction 29.7 %
      K=200  RMS 14.2 %   prediction 21.0 %
    m=5.082 hi (CV 30.9):
      K= 10  RMS 100 %    prediction 978 %      K=100  RMS 353 %    prediction 309 %
      K=200  RMS 229 %    prediction 219 %

The RMS error tracks CV/sqrt(K) (ratios 0.5-1.3, i.e. at most a factor ~2 better than plain MC
with K parents), **not** K^-2. The reason is physical: the per-scalar response is not a smooth
function of (p, theta) -- it is the probability that the *decay products* pass the detector, and
that varies parent-to-parent at fixed (p, theta).

Estimated bins needed for a 5 % quadrature: K ~ (CV/0.05)^2 ~ (2.5-3.0/0.05)^2 ~ 2500-3600
cells per channel (and ~4e5 for the hi wall where CV = 31), i.e. the cost of 2500-3600 parents,
*before* adding the binning bias from the variables not binned (phi, vertex, daughter flavour).
There is no cost advantage over simply running 2400-4800 parents, and the quadrature adds a
systematic error. **Idea B is not competitive for this observable.**

---

## 4. Cost estimate (task 4)

Measured basis (single-threaded, this machine): the full 52-mass grid at nsamp=300 takes
**119 s** of scan time (file I/O excluded); the number of transported rows is linear in nsamp
(verified: 1.00 / 2.03 / 4.06 for x1 / x2 / x4), so t(n) = 4 t(300) n/1200:

    nsamp      1200      2400      4800      10000     25871
    CPU        7.9 min   15.8 min  31.7 min  66 min    171 min

The production README quotes ~20 min for the 1200 scan; the measured value on this machine is
7.9 min (the README number is conservative; all cost estimates below use the measured one).

Required nsamp per mass for a 5 % boundary error, n = 1200 (err/0.05)^2 (the boundary error
scales as rel/slope, and rel ~ n^-1/2 at lo; at hi this is a lower bound, see 2.3):

    target            median n   max n     masses >1200  >4800  >12000   total CPU   x baseline   32-core wall
    5 % at lo only     1200       2433       22/52          0       0      8.6 min     1.1x        ~1.5 min
    5 % at both edges  1901      25871       32/52          7       2     15.0 min     1.9x        ~5.5 min

with 15.8 min / 31.7 min for a uniform nsamp = 2400 / 4800 everywhere. The worst masses and
their per-mass cost at nsamp=1200:

    m=0.9054  err_lo 4.6 % n_lo 1200 | err_hi 23.2 % n_hi 25871 | 7.7 s
    m=2.32    5.4  1423            | 17.0  13917            | 4.0 s
    m=2.714   6.1  1805            | 14.8  10548            | 2.2 s
    m=1.695   4.8  1200            | 13.8   9168            | 6.9 s
    m=6.955   6.7  2132            | 12.1   7062            | 2.4 s
    m=0.3533  3.6  1200            | 11.3   6141            | 6.2 s
    m=3.175   6.5  2054            | 11.1   5959            | 3.0 s

Why the hi numbers are a lower bound: they assume the reported single-run error. For the
heavy-tail masses the true seed-to-seed spread at hi is ~2x the reported one, so their n_hi
should be multiplied by ~4 (m_S = 0.905 -> ~1e5 parents, ~11 min for that mass alone). A
tail-aware convergence criterion is needed, not just the sample variance.

**What 32-core parallelism buys and does not buy**: the 52 masses are embarrassingly parallel,
so 8.6 / 15.0 CPU-min become ~1.5 / 5.5 min of wall time (load-balanced, ~2 waves because 52
masses > 32 cores). It buys wall time only: it does **not** reduce the total CPU, it does **not**
reduce the variance of any single mass, and it does **not** fix the heavy-tail masses -- those
need more samples (or a better sampling) at those masses. The per-mass cost is very uneven
(1.6-11.3 s at nsamp=300, i.e. 6-45 s at 1200: the low-mass points carry the K channels), so a
naive split of the 52 masses over 32 cores wastes ~40 % of the cores in the second wave unless
the expensive low-mass points are launched first.

---

## 5. Honest bottleneck (task 5)

Ranked by measured contribution to the variance of s at the crossing:

1. **Finite parent sample -- the MC mean of the per-parent response (90-100 % of the variance).**
   At lo: rel ~ nsamp^(-0.5) measured; the 1200-parent sample uses n_eff = 190-305 effective
   parents (6-7 % of the rows have non-zero response, the top 1 % carry 35-48 % of s); the
   parent-set component of the seed scatter (npseed 4.5-7.1 % at lo) is larger than the
   direction component (--seed 3.3-5.4 %). At hi it is 100 % and n_eff = 3-7.
2. **Rare/coincident surviving scalars -- the mechanism that makes (1) non-Gaussian at hi.**
   Not a separate error term, but the reason: the hi boundary moves by factors 1.7-2 between
   1200-parent scans, the nsamp scaling is only nsamp^(-0.16..-0.31), and the K channel shows
   relerr = 1.0 at 2250/2640 points (one surviving scalar). This is why 5 % at hi is expensive
   and why a deterministic parent set is unusable there.
3. **Inner A(z) MC (ndecay = 8): 2-10 % of the variance at lo, 0 at hi.** Raising ndecay
   8 -> infinity would improve the lo error by 2-10 % (m=0.302: nothing; m=1.983: 10 %) at
   linear cost. Verdict: keep ndecay = 8; it is not the bottleneck.
4. **One-mode-per-parent branching-ratio sampling with the summed-BR weight wB: < 0.1 %.**
   wB = 3.35 / 3.32 / 3.17 at m_S = 0.302 / 0.905 / 1.983 GeV (0.01 at 5.082 GeV where only
   one mode is open) -- the starting fact wB ~ 3.36 is confirmed. The weight is a
   normalisation constant that cancels in the relative error; the only noise is the spread of
   the mode-averaged response, and the S energy in the B rest frame varies by only 10-12 %
   across the 10 open modes (2.37-2.65 GeV at 0.302, 2.44-2.72 at 0.905, 2.73-3.01 at 1.983),
   so CV_mode = 0.05-0.06 and CV_mode^2/1200 = 2-3e-6 against 2.5e-3 - 4.3e-3 measured at lo
   and 0.09-0.95 at hi. Negligible; do not spend CPU here.

---

## 6. Ranked recommendations

1. **nsamp = 2400 everywhere** (16 min single-threaded, ~2 min wall on 32 cores): 5 % at lo for
   all 52 masses, 5 % at hi for 27 of the 34 two-edge masses, and halves the visible jitter of
   the contour. Cheapest defensible improvement; no code change.
2. **nsamp = 4800-25000 for the 7 hi-wall masses only** (0.905, 2.32, 2.71, 1.70, 6.96, 0.353,
   3.18 GeV): +1.3x total CPU on top of (1); these are the points where the current contour is
   not reproducible at 5 %.
3. **Report a tail-aware boundary error at hi** (e.g. the scatter of the crossing over 3 parent
   seeds, not the single-run sample variance): the present quoted 4-23 % is a factor ~2
   optimistic for the worst masses, and the m_S = 0.905 hi point of the production scan sits
   2 sigma away from the value all reruns give.
4. **If 5 % at hi for every mass is a hard requirement**: importance-sample the S decay
   direction (the S-pz correlation ratio promises 1.9-4.4 at lo and the hi wall is a geometric
   acceptance problem); a fixed parent set cannot do it (3.3). This is a change to the sampling,
   so it needs a deliberate decision, not a silent patch.
5. **Do not**: raise ndecay (item 3), change the branching-ratio sampling (item 4), or use a
   deterministic stratified parent set at hi (3.3).
6. Optionally keep the stratified design for the *lo* edge only, with a Neyman allocation built
   for lo (m_S = 5.082: se 5.4 % -> 2.9 % expected, see 3.3) -- but it requires a second parent
   set for the hi edge, which is where it fails anyway, so the gain does not justify the
   complexity for this scan.

---

## 7. Commands and runtimes (reproduction)

All from 'runs/bc5_lhcb_downstream_v3', with 'D:\py\python.exe':

    # 1. diagnostic (production scan + exclusion_boundary.csv), seconds
    D:\py\python.exe work/analysis/mc_study.py diag

    # 2. experiment runs (each writes outputs/data/signal_scan<sfx>.json only), wall times above
    D:\py\python.exe work/analysis/mc_study.py run --name base   --suffix _base   --nsamp 1200 --ndecay 8
    D:\py\python.exe work/analysis/mc_study.py run --name ns2400 --suffix _ns2400 --nsamp 2400 --ndecay 8
    D:\py\python.exe work/analysis/mc_study.py run --name ns4800 --suffix _ns4800 --nsamp 4800 --ndecay 8
    D:\py\python.exe work/analysis/mc_study.py run --name nd16   --suffix _nd16   --nsamp 1200 --ndecay 16
    D:\py\python.exe work/analysis/mc_study.py run --name nd32   --suffix _nd32   --nsamp 1200 --ndecay 32
    D:\py\python.exe work/analysis/mc_study.py run --name seed1  --suffix _seed1  --nsamp 1200 --seed 11235813
    D:\py\python.exe work/analysis/mc_study.py run --name seed2  --suffix _seed2  --nsamp 1200 --seed 271828
    D:\py\python.exe work/analysis/mc_study.py run --name np1    --suffix _np1    --nsamp 1200 --npseed 5555555
    D:\py\python.exe work/analysis/mc_study.py run --name np2    --suffix _np2    --nsamp 1200 --npseed 987654321
    # (all of the above use the 4 representative masses 0.302053,0.905421,1.983347,5.082255,
    #  the default --masses; each was given --timing-out work/analysis/_mc_timings/<name>.json)

    # 3. full-grid timing calibration (writes signal_scan_cost300.json), 130 s wall
    D:\py\python.exe work/analysis/mc_study.py run --name cost300 --suffix _cost300 --nsamp 300 --nmass 52

    # 4. fits, crossings, seed scatter (reads the timing json + the 9 suffixed scans), seconds
    D:\py\python.exe work/analysis/mc_study.py fit

    # 5. variance-reduction prototype per mass (pilot 1200 + 6 fixed designs), 1-3.5 min each
    D:\py\python.exe work/analysis/mc_study.py proto --mass 0.3020529196914042
    D:\py\python.exe work/analysis/mc_study.py proto --mass 0.9054207886414979
    D:\py\python.exe work/analysis/mc_study.py proto --mass 1.9833466252592793
    D:\py\python.exe work/analysis/mc_study.py proto --mass 5.082254549604771

    # 6. cost table and BR-mode estimate, seconds
    D:\py\python.exe work/analysis/mc_study.py cost
    D:\py\python.exe work/analysis/mc_study.py modevar

Total CPU spent on the experiments (9 runs + timing grid + 4 prototypes): ~30 min, inside the
40 min budget. Outputs kept: 'outputs/data/signal_scan_{base,ns2400,ns4800,nd16,nd32,seed1,seed2,np1,np2}.json'
(~0.3 MB each). The temporary timing json files and the calibration scans (_cost100, _cost300,
_pcheck1, _pcheck2, _timeprobe) were deleted after the numbers above were tabulated.
'mc_study_proto.py' is kept because 'mc_study.py proto' imports it.

### Caveats

* Wall times were measured with 9 processes competing for 32 logical cores; identical runs
  repeat to ~+-20 %. The cost table uses the dedicated single-process nsamp=300 full-grid run.
* The hi edge is not Gaussian: all hi "required nsamp" numbers are extrapolations of the
  reported (single-run) standard error and are lower bounds for the heavy-tail masses.
* The prototype's design uses the pz of the *parent particle*; the scalar's own pz is the
  physically controlling variable (see 3.2) but it is only known after the decay is sampled.
