# Decay figures (task-9): three standalone smooth-curve plots

Status: **DONE** -- three figures written, no digitised points and no markers anywhere; every curve is the author Light_scalar_decay program evaluated on 500 log-spaced masses.

Author: decay-figures (teammate), 2026-10-26.  Run: runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream.
Script: work/decay/make_decay_figures.py (sha256 `a67f3d71e37005c97c90dc41f6b713b1e40c32c3496946f7ecde0a3da3af92a8`).

## 1. Deliverables

| figure | content | axes |
|---|---|---|
| `outputs/figures/decay_branching_ratios.png` / `.pdf` | BR(H->XX) vs mH, all 10 channels | log-log, x 0.1-10 GeV, y 1e-9-2 |
| `outputs/figures/decay_partial_widths.png` / `.pdf` | Gamma(H->XX) [GeV] vs mH, all 10 channels | log-log, x 0.1-10 GeV, y 1e-24-1e-9 |
| `outputs/figures/decay_total_width_ctau.png` / `.pdf` | Gamma_tot vs mH for tanb = 5,10,70,200,1000 (left) and c*tau (right) | log-log, x 0.1-10 GeV, y 1e-22-5e-9 GeV / 3.9e-8-2.0e6 m |

Colours and LaTeX labels are identical to make_overlay() in work/decay/validate_decay.py: `gaga`=purple, `ee`=cyan, `mumu`=orange, `tautau`=peru, `pipi`=blue, `KK`=red, `4pi`=green, `cc`=navy, `ss`=darkred, `gg`=darkgreen.
The old outputs/figures/decay_validation.png (1x3 panel with digitised circles) is left untouched; the three new figures replace its left/middle/right panels without any data points.

Two deviations from the old panel ranges, needed so that all ten channels are visible: the BR lower limit is 1e-9 (old 5e-3) because BR(e e) <= 8.8e-5, and the total-width window is 1e-22-5e-9 (old 1e-24-5e-9) because Gamma_tot spans 4.9e-21-9.3e-10 over the five tanb values.

## 2. Generation

python work/decay/make_decay_figures.py -- one batch call to the author executable, 2515 points (500 masses x 5 tanb values + 15 spot/jump-check points), runtime 8.65 s.

| item | value |
|---|---|
| executable | `E:\AAAwork\LLP-AGI\runs\20261025_1325_UTC_2HDM_lightH_lhcb_downstream\work\decay\build\light_scalar_decay.exe` |
| executable sha256 | `5c63124c9b60ead8fd4c68874aa8a5ce916ddada3a2dd9d6ed877ec5ec1c909a` |
| mass grid | 500 log-spaced points, 0.1 ... 10 GeV |
| tanb values | 5, 10, 70, 200, 1000 |
| closest sampled mH to 2.0 GeV | 0.0073 GeV (the mH=2.0 node is never evaluated: the program double counts there) |
| Gamma_tot range on all points | 4.930e-21 ... 9.305e-10 GeV |
| c*tau range on all points | 2.121e-07 ... 4.003e+04 m |

## 3. Sanity checks

| check | result |
|---|---|
| max abs(BR sum - 1) over all 2515 points | **7.46e-11** |
| median abs(BR sum - 1) | 8.24e-12 |
| Gamma_tot / sum(10 channels) min ... max | 0.999999999925 ... 1.000000000074 |
| max rel. deviation of c*tau*Gamma_tot from hbar*c | 2.22e-16 |
| spot values vs logs/decay_validation.md sec. 7: max rel. dev. Gamma_tot | 3.07e-07 |
| same for c*tau | 1.86e-05 |

### Spot values (ours vs the validation log table)

| mH [GeV] | tanb | Gamma_tot [GeV] | rel.diff | c*tau [m] | rel.diff | BR(gaga) | BR(pipi) | BR(gg) | BR(cc) |
|---|---|---|---|---|---|---|---|---|---|
| 0.5 | 10 | 1.272186e-14 | 3.1e-07 | 1.551086e-02 | 8.8e-06 | 0.4814 | 0.4646 | 7.86e-27 | 7.86e-27 |
| 1 | 10 | 1.631479e-13 | 2.7e-08 | 1.209498e-03 | 2.0e-06 | 0.3003 | 0.5481 | 6.13e-28 | 6.13e-28 |
| 2 | 10 | 4.356101e-13 | 9.0e-09 | 4.529899e-04 | 2.1e-07 | 0.8995 | 0.000318 | 0.0325 | 2.3e-28 |
| 5 | 10 | 8.367106e-12 | 1.4e-09 | 2.358366e-05 | 1.4e-05 | 0.7331 | 1.195e-29 | 0.0114 | 0.143 |
| 1 | 5 | 7.610433e-12 | 5.3e-08 | 2.592848e-05 | 1.9e-05 | 0.02528 | 0.7636 | 1.31e-29 | 1.31e-29 |
| 0.5 | 100 | 6.162374e-17 | 8.0e-08 | 3.202126e+00 | 8.2e-06 | 0.9999 | 9.545e-05 | 1.62e-24 | 1.62e-24 |
| 1 | 100 | 4.930578e-16 | 4.6e-08 | 4.002106e-01 | 1.6e-06 | 0.9998 | 0.0001805 | 2.03e-25 | 2.03e-25 |
| 2 | 100 | 3.943782e-15 | 8.6e-08 | 5.003496e-02 | 8.6e-07 | 1 | 3.495e-08 | 3.58e-06 | 2.54e-26 |
| 0.5 | 1000 | 6.162092e-19 | 4.3e-08 | 3.202273e+02 | 8.5e-06 | 1 | 9.545e-09 | 1.62e-22 | 1.62e-22 |
| 1 | 1000 | 4.929737e-18 | 4.0e-08 | 4.002789e+01 | 2.7e-06 | 1 | 1.805e-08 | 2.03e-23 | 2.03e-23 |
| 5 | 1000 | 6.164705e-16 | 4.8e-08 | 3.200915e-01 | 4.8e-06 | 1 | 1.622e-25 | 1.54e-10 | 1.94e-09 |

The tiny-BR entries (gg/cc below threshold) are the 1e-40 program floor divided by Gamma_tot; the validation table prints them with 2 significant digits, so their relative differences are print quantisation, not physics. Gamma_tot reproduces the table to <= 3e-07 and c*tau to <= 2e-05 (limited by the 4-5 significant digits printed there).

## 4. Thresholds and the 2 GeV matching (features visible in the curves)

Channel thresholds hard-coded in the author program (CPevenDecay.cpp), with the first sampled mH at which the width is above the 1e-40 floor (tanb = 10; the 500-point grid resolves a threshold to ~1 % in mass):

| channel | threshold [GeV] | first sampled point [GeV] | in 0.1-10 window |
|---|---|---|---|
| ee | 0.001022 | open from 0.1 | no (below 0.1) |
| mumu | 0.21132 | 0.2131 | yes |
| pipi | 0.26994 | 0.2709 | yes |
| 4pi | 0.53988 | 0.5413 | yes |
| KK | 0.996 | 1.0046 | yes |
| ss | 2 | 2.0073 | yes |
| gg | 2 | 2.0073 | yes |
| tautau | 3.5536 | 3.5572 | yes |
| cc | 3.72969 | 3.7597 | yes |

At the 2 GeV hadron/parton matching (m_H > 2 floors pi pi / K K / 4pi and switches on s sbar / g g) the individual channels jump while the total width stays continuous. Bracketed at mH = 1.999 vs 2.001 GeV, tanb = 10:

| quantity | 1.999 GeV | 2.001 GeV |
|---|---|---|
| Gamma_tot [GeV] | 4.1592e-13 | 4.1621e-13 (ratio 1.0007) |
| BR(pipi) | 0.00279 | 2.403e-28 |
| BR(KK) | 0.02195 | 2.403e-28 |
| BR(4pi) | 0.02583 | 2.403e-28 |
| BR(ss) | 2.404e-28 | 0.01433 |
| BR(gg) | 2.404e-28 | 0.03409 |
| BR(gaga) | 0.9407 | 0.9429 |
| Gamma(pipi) [GeV] | 1.160e-15 | 1.000e-40 |
| Gamma(KK) [GeV] | 9.130e-15 | 1.000e-40 |
| Gamma(4pi) [GeV] | 1.074e-14 | 1.000e-40 |
| Gamma(ss) [GeV] | 1.000e-40 | 5.963e-15 |
| Gamma(gg) [GeV] | 1.000e-40 | 1.419e-14 |

This is the intended matching of the author program (the 4pi coefficient C = 5.1e-9 was tuned so that the hadronic and partonic descriptions meet at 2 GeV), not a plotting artifact; it is why the curves in figures 1-2 have a near-vertical jump there. A thin grey dotted guide line marks mH = 2 GeV in figures 1-2 (a guide line, not a data marker).

Structure at mH ~ 1 GeV (the author form-factor cusp just above the K K threshold, 2 m_K = 0.996 GeV), tanb = 10, 0.9 <= mH <= 1.3 GeV:

| channel | min BR | at mH | max BR | at mH |
|---|---|---|---|---|
| gaga | 0.18 | 0.9772 | 0.7728 | 1.2889 |
| pipi | 0.01966 | 1.2889 | 0.8095 | 0.9772 |
| KK | 3.937e-28 | 0.9772 | 0.2615 | 1.0521 |

Other visible features: gaga BR -> 1 below the pi pi threshold (only gamma gamma and e e are open); abrupt onsets at mH = 0.211 (mu mu), 0.270 (pi pi), 0.540 (4pi), 0.996 (K K), 2.0 (ss, gg), 3.554 (tau tau), 3.730 (c cbar); at large tanb gamma gamma dominates from low mass (figure 3 shows the correspondingly tiny total widths and long c*tau).

## 5. Files

| file | sha256 |
|---|---|
| `outputs/figures/decay_branching_ratios.png` | `15d7e1677effdb9a35d0134340332bb0fc2d1e9aeafb513a169fed9f542e10a1` |
| `outputs/figures/decay_branching_ratios.pdf` | `00cfcdf26b0a01031de43709ea2fdda56e47fa8198d88d93d650a6ad36407bb3` |
| `outputs/figures/decay_partial_widths.png` | `78d2d67c9397d9440bf9f098bb31b6aea366f1f2f8a8c8b5e80f7e32a772b896` |
| `outputs/figures/decay_partial_widths.pdf` | `cff557a1356073df5fd62bce8f86c705703c002635cfdcc032f0e61c07772f7f` |
| `outputs/figures/decay_total_width_ctau.png` | `69df9ffa79c4bad6bccfe668b9251f480cca3a1fb14503e4243d4f6c82985975` |
| `outputs/figures/decay_total_width_ctau.pdf` | `80f4c7c24b847efade83fbbe1045bfde084d8466c69433f051d2c5f609c86697` |

