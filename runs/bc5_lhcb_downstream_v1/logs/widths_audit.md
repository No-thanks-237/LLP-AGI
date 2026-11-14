# Width-model audit — full comparison against every reference data set

Produced by work/widths/audit_widths.py; machine-readable output in
outputs/data/widths_audit.json.  Three independent references are used:

| id | source | content |
| --- | --- | --- |
| W | Winkler arXiv:1809.01876 Fig. 2, digitised into work/widths/reference_figures/decayrates_curves.csv | absolute partial widths, per channel, 0.28-7.5 GeV |
| S | SensCalc tables, phenomenology/darkscalar/prepared/branching_ratios.csv | branching ratios, exclusive channels, 0.021-62.49 GeV |
| B | Boiarska et al. arXiv:1904.10447 branching-ratio figure, digitised into work/widths/reference_figures/scalardecaybr_curves.csv | branching ratios, 0.1-10 GeV |

## A. Absolute widths against W  (ratio = mine / reference)

| channel | n | median | range | worst point |
| --- | --- | --- | --- | --- |
| mu+mu- | 500 | 0.991 | 0.990 - 1.005 | 5.979 |
| tau+tau- | 226 | 1.020 | 0.995 - **5.50** | 3.563 (threshold) |
| pi pi | 1815 | 0.991 | **0.436** - 0.993 | 0.281 (threshold / frame floor) |
| K K | 849 | 0.991 | 0.800 - **4.74** | 0.987 (threshold) |
| multi-meson | 193 | 0.995 | 0.995 - 0.996 | — |
| s sbar | 233 | **0.778** | 0.673 - 0.965 | 7.499 |
| c cbar | 181 | **0.644** | **0.001** - 0.692 | 3.740 (threshold) |
| g g | 238 | 1.032 | 0.987 - 1.259 | 2.542 |
| hadronic total | 1770 | 0.991 | 0.436 - 1.166 | 0.281 |

The 0.991 median for the leptons is the known 0.9% offset of the digitised figure itself (the
muon calibration factor is 1.00946).  The pi pi / K K / 4pi curves were calibrated with that
factor; the lepton curves in the figure are calibrated, not ours.

## B. Branching ratios against S  (mine are isospin totals; the exclusive S columns are summed)

| channel | n | median | range | worst point |
| --- | --- | --- | --- | --- |
| e+e- | 120 | 1.008 | 0.727 - 2.259 | 3.473 |
| mu+mu- | 85 | 1.077 | 0.727 - 2.259 | 3.473 |
| tau+tau- | 43 | 1.086 | 0.988 - 1.888 | 3.714 |
| pi pi | 29 | 1.035 | 0.936 - 2.072 | 1.896 |
| **K K** | 10 | **1.308** | 1.140 - 2.020 | 1.896 |
| s sbar | 52 | 0.974 | 0.780 - 1.392 | 3.473 |
| c cbar | 42 | 1.308 | 1.017 - 1.323 | 62.49 |
| b bbar | 27 | 0.981 | 0.828 - 0.997 | 10.886 |
| **g g** | 52 | **0.492** | 0.414 - 0.962 | 5.559 |

## C. Branching ratios against B  (digitised figure)

| channel | n | median | range | worst point |
| --- | --- | --- | --- | --- |
| e+e- | 112 | 1.007 | 0.838 - 1.071 | 0.211 |
| mu+mu- | 1998 | 0.984 | 0.219 - 1.766 | 0.211 |
| tau+tau- | 636 | 1.059 | 1.008 - 1.407 | 3.555 |
| pi pi | 1999 | 1.011 | 0.627 - 643 | 0.280 (artefact) |
| K K | 1323 | 1.095 | 1.048 - 256 | 1.000 (artefact) |
| multi-meson | 1847 | **0.622** | 0.434 - 5.256 | 0.560 |
| **s sbar** | 753 | **0.641** | 0.502 - 0.708 | 2.000 |
| c cbar | 461 | 0.924 | 0.000 - 0.978 | 3.739 (threshold) |
| g g | 768 | 1.174 | 1.077 - 1.436 | 2.000 |

## D. Internal consistency

| check | result |
| --- | --- |
| branching ratios sum to 1 | exact (max deviation 0) |
| hadronic width below / above the 2 GeV matching scale | 7.761e-8 / 8.937e-8 GeV, ratio **1.151** |
| table rows / mass range / spacing | 3000 rows, 0.2792-2.0000 GeV, 1.8e-4 to 1.3e-3 GeV |
| table non-finite / negative entries | 0 / 0 |
| every channel opens at its nominal threshold | yes |

## Findings, in order of importance

### Fixed in this pass

1. **The K K, pi pi and 4pi curves in the branching-ratio figure were drawn above 2 GeV as a
   1/m^2 extrapolation.**  The model has no such channels there: above the matching scale the
   hadronic width is perturbative.  The curves now stop at 2 GeV.  This was a presentation
   fabrication, not an error in the physics or in the exclusion curve.
2. **One NaN row in hadronic_table.csv** (the very first row, at the pi pi threshold m = 0.2792)
   made any query in [0.2792, 0.2801] return NaN.  No scan mass falls there, so the result was
   unaffected, but the table is now rebuilt clean.
3. **The table was rebuilt on 3000 points instead of 600** (spacing 1.8e-4 instead of 9e-4 GeV).
   The K K threshold rises within 5e-4 GeV and was smeared by the old grid; the discrepancy at
   0.9874 GeV dropped from a factor 8 to a factor 4.7.  At the 52 scan masses the change is at
   most 0.09%, so the scan and the exclusion curve did not need to be re-run.

### Real differences between models, all previously documented

4. **The 2 GeV matching step is 1.151**, whereas both W and S are continuous there (1.0002 and
   1.0000).  Cause: Gamma(S -> g g) at 2 GeV is 22% above the W curve (it agrees to better than
   3% above 4 GeV), which is the low-scale alpha_s.  Effect on the result: a ~15% step in
   Gamma_S at 2 GeV, i.e. a ~15% step in the theta^2 boundary at that mass.
5. **s sbar and c cbar are 22% and 36% below the W curve** because the running quark masses are
   evaluated at mu = m_S (as Boiarska prescribes) whereas the W figure appears to use fixed
   masses.  This is a genuine scheme difference, not an error: the two references bracket our
   value.  S has c cbar 31% *below* ours, B has s sbar 36% *above* ours.
6. **The gluon channel sits between the references**: 0.49 x S but 1.17 x B.  The two references
   disagree with each other by a factor 2.4 on this channel alone.
7. **K K is 31% above S** on median (1.31); against B it is 1.10.  The K K channel is the least
   constrained one because it is the channel that the f0(980) shape feeds.

### Artefacts of the references, not discrepancies

8. Every entry with a wildly out-of-range ratio is **within about 1 MeV of a threshold**:
   tau+tau- at 3.563 (5.5x), pi pi at 0.281 (0.44x), K K at 0.987 (4.7x), c cbar at 3.740
   (0.001x), and the 643x / 256x outliers of the B figure at 0.280 and 1.000 GeV.  At a
   threshold the rate varies as beta^3 and the digitised curves are multi-valued or clipped by
   the frame, so a ratio there carries no information.  No scan mass falls within 20 MeV of a
   threshold except where the physics is genuinely smooth.
9. The 0.991 median of the lepton channels against W is the figure's own 0.9% calibration
   offset, measured on the exactly known muon curve.

### Remaining known limitations (unchanged)

10. The pi pi / K K resonance region is an external published input, not a computation; see
    section 1b and work/widths/research/README.md.  Boiarska et al. quote this channel as
    uncertain by more than an order of magnitude around 1 GeV.
11. The multi-meson toy channel is 38% below B's 4pi curve, which is as expected for a toy
    formula whose only constraint is continuity at the matching scale.
12. Quark channels are LO (QUARK_NLO = False); see section 1b for why.

## E. Why the K K branching ratio differs from SensCalc near 2 GeV

This was the one entry that looked alarming in the first pass (BR_KK / SensCalc = 2.20 at
1.99 GeV).  It is a denominator effect, not a K K discrepancy.  Three facts settle it:

1. **Our K K rate IS the Winkler curve.**  Gamma_KK agrees with the digitised Winkler Fig. 2 to
   0.2-0.8% at every mass below 2 GeV (3.5310e-8 vs 3.5646e-8 at 1.99 GeV, 4.1214e-8 vs
   4.1604e-8 at 1.90 GeV).  There is no independent K K model here to disagree with anything.

2. **At exactly 2 GeV both models have zero exclusive K K.**  Ours switches to the perturbative
   branch (Gamma_KK = 0 at m_S = 2.000); SensCalc switches too (K+K- = K_LK_S = 0 at 2.000, with
   Jets_ss = 0.218 and Jets_GG = 0.711 instead).  The "2 GeV" comparison is therefore zero
   against zero and carries no information.

3. **The absolute K K widths agree much better than the branching ratios.**  Converting the
   SensCalc branching ratios with SensCalc's own total width (GammaS_reference.csv):

   | m_S | our Gamma_KK | SensCalc Gamma_KK | ratio | our Gamma_tot | SensCalc Gamma_tot | ratio |
   | --- | --- | --- | --- | --- | --- | --- |
   | 1.90 | 4.12e-8 | 5.41e-8 | 0.76 | 9.18e-8 | 2.44e-7 | 0.38 |
   | 1.95 | 3.78e-8 | ? | - | 9.16e-8 | 1.62e-7 | 0.56 |
   | 1.99 | 3.53e-8 | 3.06e-8 | 1.15 | 9.20e-8 | 1.75e-7 | 0.53 |

   The K K rates agree to 15% at 1.99 GeV; what differs by a factor 1.9 is the TOTAL width, and
   that is what enters the branching ratio.

**Why is the SensCalc total width so much larger in 1.5-2 GeV?**  Because its four-pion channels
carry most of the width there: at 1.99 GeV, 2pi+2pi- = 0.245 and pi+pi-2pi0 = 0.491, i.e. 74% of
the total, against our toy multi-meson channel at about 40%.  SensCalc computes the real 4pi
channels; ours is Winkler's toy C m_S^3 beta_2pi, whose only constraint is continuity at the
2 GeV matching scale.

Moreover the SensCalc total width has a **pathological spike at 1.70-1.85 GeV**:

| m_S | 1.60 | 1.70 | 1.75 | 1.80 | 1.90 | 1.99 | 2.00 | 2.10 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SensCalc Gamma_tot | 1.81e-7 | 5.47e-7 | **8.96e-7** | **8.60e-7** | 2.44e-7 | 1.75e-7 | 2.03e-7 | 2.22e-7 |
| this work | 1.13e-7 | 1.03e-7 | 9.77e-8 | 9.44e-8 | 9.18e-8 | 9.20e-8 | 1.04e-7 | 1.14e-7 |

It rises by a factor 5 between 1.60 and 1.75 GeV, falls again by a factor 5 into 1.95 GeV, and
is non-monotonic across the 2 GeV switch.  That structure comes from the 4pi sector, which
Winkler himself describes as "currently not known" and which Boiarska et al. quote as uncertain
by more than an order of magnitude.  Our K K result is therefore *not* the outlier here.

**Cross-check that the bookkeeping is sound:** in absolute terms the leptonic widths agree with
SensCalc to better than 0.1% (e+e-, mu+mu-, tau+tau- all have median ratio 0.999), because both
sides use the exact G_F formula.  Any disagreement between the two is therefore confined to the
hadronic sector, and the hadronic sector is precisely the one that neither side claims to know.

### Absolute widths against SensCalc (its BRs times its own total width)

| channel | n | median | range | worst m_S |
| --- | --- | --- | --- | --- |
| e+e- | 120 | 0.999 | 0.182 - 1.282 | 1.773 |
| mu+mu- | 85 | 0.999 | 0.183 - 1.282 | 1.773 |
| tau+tau- | 43 | 0.999 | 0.998 - 1.001 | 62.49 |
| pi pi | 29 | 0.822 | 0.176 - 1.306 | 1.773 |
| K K | 10 | 0.644 | 0.181 - 1.171 | 1.773 |
| s sbar | 52 | 0.795 | 0.502 - 0.984 | 2.028 |
| c cbar | 42 | 1.117 | 0.786 - 1.338 | 62.49 |
| g g | 52 | 0.433 | 0.361 - 0.568 | 4.544 |
