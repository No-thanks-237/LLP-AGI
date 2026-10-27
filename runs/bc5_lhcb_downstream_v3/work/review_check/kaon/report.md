# F04 + F09 audit (task-2): neutral-kaon flight, kaon parent sample / mixing / history

Read-only audit of run `runs/bc5_lhcb_downstream_v3` (model tag v3.3-2026-10-21).
All numbers re-derived here; scripts in this directory (`f04_mc.py`, `f04_analyse.py`,
`f04_summary.py`, `f09_cut_check.py`, `f09_loss.py`, `f09_verify.py`), results in
`f04_results.json` / `f04_m*.npz`. No production file was modified.

## Verdict table

| claim | verdict | decisive number |
|---|---|---|
| F04 neutral-kaon daughters placed at the S vertex | **CORRECT** | KK visible rate drops to r_KK = 0.88/0.86/0.83/0.82/0.81 at m_S=1.0-2.0; total s at the 25 fb^-1 long-lived crossing drops by -3.9/-10.5/-13.1/-10.1/-8.0 %; the crossing moves by +3.1/+7.0/+8.6/+6.5/+5.1 % in theta^2 |
| F09-1 K sample pre-cut 1.5<eta_K<5.5 | **CORRECT** (fact) / negligible impact | verified eta in [1.50003, 5.49999]; lost fraction of the K->pi S channel 7.7 % (m_S=0.1), 4.1 % (0.2), 1.2 % (0.3), 0 at 0.5 and 1.0 (channel closed); K channel is <=1.2e-4 of s => Delta s/s <= 1e-5 |
| F09-2 PDG 311/-311 treated as full K_L | **CORRECT** / negligible impact | code uses GAM_KL for every 311 (run_scan.py:357, bc5_production.py:207-219); a K0 flavour state has |K_L|^2=1/2 => factor 2. K0L is ~48 % of the K channel => K channel -24 % => Delta s/s <= 3e-5 |
| F09-3 K->SS and K->SSpi share hid=1 | **CORRECT** (1.783x) / negligible impact | production channel_efficiency: 1.513e-13 vs correct exclusive sum 8.485e-14 = 1.783; K quartic is <=2.3e-5 of s (only at theta^2->1e-20) => Delta s/s <= 1e-5 |

The K channel (B/K/H scan components) is 8e-5 ... 1.2e-4 of the total signal at the
long-lived crossing everywhere below 0.36 GeV and exactly 0 above, in agreement with the Lead.

## F04: how measured

For each m_S, S four-momenta were drawn through the real production code: B->X_s S mixing
(real channel mixture and br_B_to_S, B sample), B->X_s S quartic 3-body (125 S per mass),
and h->SS (real Higgs sample, both scalars). Only S with 2<eta_S<5 are kept (this selection
factor cancels in every ratio quoted). For every S the code's acceptance table A(z) is
computed on the production z grid (bc5_lib.ZGRID, 13 points) with 8 decays per point.

* (a) current code: all charged daughters born at the S decay vertex, on the S flight line
  `(x0+(z-z0)*px/pz, ...)`, exactly as `bc5_lib.acceptance_profile` (my reimplementation
  reproduces the production `_passes_pair` response and the same random stream as
  `s_decay_products(...,'KK',...)`).
* (b) flight-corrected: the S->KK branch is replayed as in the code (50 % K+K-, 50 % K0K0bar
  -> 50/50 K_SK_S / K_LK_L). For every neutral kaon the lab 4-momentum is boosted from the
  S rest frame and a decay length l = (p/m) c tau * (-ln U) is drawn with c tau(K_S)=2.684 cm
  and c tau(K_L)=15.337 m (from the production GAM_KL). Daughters are boosted from the kaon
  rest frame and transported from the true decay point with `downstream_ctypes.transport_track`
  plus the standard birth cuts E>5 GeV, pT>0.5 GeV.
* A decay point with z>frontZ=7.7 m (or outside the envelope) gives `BornAfterFront` /
  `OutsideSciFi` from the real kernel and contributes no pair: it is discarded as not
  reconstructable in this downstream model. Nothing else is special-cased.
* The z-integral uses the production weight w(z)=alpha exp(-(z-z0)alpha),
  alpha=m_S/(c tau_S pz), exactly as `run_scan.per_scalar_p`; the S lifetime is taken at the
  measured 25 fb^-1 long-lived crossing (and separately at the short-lived crossing).
* 500 S per source per mass; errors from bootstrap over S (400 resamplings). The combined
  total uses the scan's own B and H signals at the crossing as weights (B share 65-73 %).

### Results at the 25 fb^-1 long-lived crossing

| m_S | r_KK (b)/(a) combined | r_KK B-source | r_KK H-source | f_KK(B) | f_KK(H) | Delta s/s total | dtheta^2/theta^2 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.0 | 0.884 +- 0.010 | 0.9387 +- 0.0056 | 0.7061 +- 0.0102 | 0.354 | 0.290 | -3.9 % | +3.1 % |
| 1.2 | 0.861 +- 0.010 | 0.9087 +- 0.0077 | 0.7128 +- 0.0118 | 0.755 | 0.744 | -10.5 % | +7.0 % |
| 1.5 | 0.826 +- 0.010 | 0.8822 +- 0.0093 | 0.7002 +- 0.0108 | 0.772 | 0.726 | -13.1 % | +8.6 % |
| 1.8 | 0.816 +- 0.010 | 0.8708 +- 0.0092 | 0.6966 +- 0.0091 | 0.573 | 0.506 | -10.1 % | +6.5 % |
| 2.0 | 0.809 +- 0.010 | 0.8620 +- 0.0097 | 0.6967 +- 0.0102 | 0.446 | 0.373 | -8.0 % | +5.1 % |

f_KK is the KK-channel share of that source's visible signal; Delta s/s = sum_s w_s f_KK,s (r_KK,s - 1).
The theta^2 shift uses the measured local dln s/dln theta^2 = 1.27-1.58 at the crossing.
The Lead's crude upper bound (delete K_L daughters only) gave +9.6 % at 1.5 GeV; the full
treatment gives **+8.6 %**, i.e. inside that bound (the difference is the K_S flight loss in
the Higgs channel, see below).

At the short-lived crossing the same ratio is 0.4-0.75 but there the current signal is carried
by a handful of extreme high-pz S (effective sample 1-16), so those numbers are only indicative.

### Why r_KK depends on the source (this is the physically important part)

| | mean lab flight of K_S | beta*gamma | P(decay beyond front face, z0=2 m) |
|---|---:|---:|---:|
| K_S from B-produced S | 1.30-1.44 m (median 0.70-0.83 m) | ~50 | 2.7-4.1 % |
| K_S from Higgs-produced S | 12.1-12.2 m (median ~9 m) | ~450 | 67-70 % |
| K_L (both sources) | 725-815 m (B), 6.8 km (H) | ~450 | 99.9-100 % |

The S in the Higgs channel are much more boosted (S energy up to hundreds of GeV), so most of
their K_S decay beyond the front face and are lost; for the B channel the S-vertex
approximation is almost harmless (~2-4 % of K_S lost). In both cases the K_L K_L half of the
neutral KK channel (25 % of the KK channel) is essentially all lost, and its daughters were
previously counted at the S vertex. Sub-channel decomposition at m_S=1.5 (B-source, z-weighted):
K+K- unchanged, K_SK_S 0.0008 -> 0.0008, K_LK_L 0.0006 -> 0.0000.
K_S displacement at these energies is therefore 0.7-1.3 m for the B channel (S-vertex
approximation is a few-% effect) but 9-12 m for the Higgs channel (a large effect).

### Convention for a K_L (or K_S) decaying inside the tracking volume

No special case: the decay point is computed from the sampled flight and the daughters are
transported from that point. If 2.5 m < z_dec < 7.7 m the daughters are inside the S/B decay
volume but still upstream of the SciFi front face, so the kernel transports them through the
field and accepts the pair normally. Only z_dec > frontZ = 7.7 m (BornAfterFront) or a decay
whose tracks leave the |x|<=3.24, |y|<=2.415 envelope is discarded. For the record, the
fraction of K_L decaying before 7.7 m is ~1 % (B source) and <0.2 % (H source), and the
fraction of K_S decaying between 7.7 and 9.4 m is <3 %; neither changes the conclusions.

### Recommended fix

1. In `bc5_lib._K0_charged`/`s_decay_products('KK')` return, together with each daughter, the
   kaon four-momentum and its mass-eigenstate lifetime; in `acceptance_profile` sample one
   decay length per neutral kaon and boost the daughters from `v0 + l*p_hat`.
2. Treat a kaon whose decay point fails the transport kernel as invisible (no pair), i.e.
   call `transport_track` from the true decay point, as done here.
3. Optional (same class, smaller): K+- also has c tau=3.69 m and ~10 % of K+- decay before the
   front face; the current code transports them as stable. This is a second-order correction
   to the charged half of the KK channel and was not applied here.

## F09-1 (kaon parent eta pre-cut)

`samples/minbias_light.txt`: 330,171 rows, eta in [1.50003, 5.49999] (all four PDG codes).
A full-phase-space sample from the same Pythia 8.316 setup exists in
`work/magnus/outputs/nocut2/out/minbias_light.txt.gz` (200k events, no eta cut, contains
321/-321/311/-311); its in-window spectra reproduce the production sample exactly
(eta and pT quantiles agree to 0.001-0.002), so it is a valid truth sample.

Direct estimate (all 1.12M parents, 12 emission directions each, full 2-body K->pi S
kinematics; the two species combined with their br_K_to_S weights):

| m_S | P_eta all | P_eta inside cut | lost fraction of K->piS channel | total-s effect |
|---:|---:|---:|---:|---:|
| 0.1 | 0.2224 | 0.6965 | **7.7 %** | <= 1e-5 |
| 0.2 | 0.2248 | 0.7321 | 4.1 % | <= 1e-5 |
| 0.3 | 0.2258 | 0.7585 | 1.2 % | <= 1e-5 |
| 0.5 / 1.0 | - | - | **0** | 0 |

At m_S >= 0.3541 (K+->pi+ S) / 0.3626 (K0L->pi0 S) and >= 0.2487 (K->SS) the kaon channel is
kinematically closed, so the pre-cut cannot matter there. Kinematically only very soft kaons
can move 2 units of eta: for eta_K=1.4 and m_S=0.1 the reachable eta_S range is
[1.0, 2.0] at p_K=5 GeV but [0.0, 7.2] at p_K<=1 GeV; such soft parents are a small fraction
and their S daughters are soft (mostly below the E>5 GeV birth cut), so the visible-signal loss
is below the 7.7 % quoted.

## F09-2 (K0 flavour state vs K_L)

`run_scan.py:357` maps 311/-311 to which='K0L'; `br_K_to_S(...,'K0L')` uses GAM_KL, i.e. the
full K_L rate is applied to every K0. At production the sample holds strong-interaction
flavour eigenstates; K0 = (K_S+K_L)/sqrt(2), so the K_L component is |K_L|^2 = 1/2 and the
mixing rate is a factor 2 too large. The K0L and K+ contributions are nearly equal in the
flight-suppressed regime (partial widths Gamma(K0L->pi0 S) : Gamma(K+->pi+ S) = 1.05, parent
counts 0.90, so K0L ~48 % of the K channel), hence the K channel drops by ~24 %.
Effect on the total signal = 0.24 x (K/s) <= 3e-5 = 0.003 % at the crossing (K/s = 8e-5-1.2e-4
below 0.36 GeV). Secondary point: the K_S quartic component (BR 7.8e-16, not flight-suppressed)
is not added; it is ~5-10 % of the neutral quartic rate and can be added for consistency.

## F09-3 (K->SS and K->SSpi share one history)

Verified with the production `run_scan.channel_efficiency` on a synthetic table: one parent,
2 scalars per mode, hid=1 for both modes, w=[4.4e-13, 4.4e-13, 6.6e-15, 6.6e-15], p=0.1.
The code returns 1.513e-13 (per-parent expectation) against the correct exclusive-history sum
8.36e-14 + 1.25e-15 = 8.49e-14, i.e. **1.783x**, as claimed (the first row supplies the
branching ratio, so the whole 4-scalar OR is weighted by BR(K->SS)). Effect on the total
signal: the K quartic is at most 2.3e-5 of s (reached only at theta^2 -> 0; at the long-lived
crossing the K quartic is ~1e-19 of the K mixing), so Delta s/s <= 1e-5.

## Recommended fixes (F09)

1. Regenerate the K parent sample without the 1.5<eta<5.5 window (a no-cut sample already
   exists in the magnus outputs, but note the newer v3fix3 minbias file contains **no** 311
   rows at all, so it cannot be used as is for the neutral-kaon channel).
2. Project each 311/-311 parent onto 1/2 K_L + 1/2 K_S: use the K_L lifetimes and BRs for the
   K_L half, add the K_S half (K_S->pi S = 0 by CP; K_S->SS / SSpi0 exist and are not flight
   suppressed), and halve the neutral mixing rate.
3. Give each quartic mode its own hid (histories of one parent add, they do not OR).

Because the K channel is <=0.01 % of the total signal everywhere, all three F09 items are
correctness-only for the present exclusion contour; they would matter only for a benchmark in
which the kaon channel dominates (which is not the case here).
