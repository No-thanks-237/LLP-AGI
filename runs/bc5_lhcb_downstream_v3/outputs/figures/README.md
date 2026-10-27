# outputs/figures — every figure and the script that makes it

All commands are run from the run directory `runs/bc5_lhcb_downstream_v3` (paths below are
relative to it) with a Python 3.11+ interpreter (numpy + matplotlib).
Model tag of the current products: **v3.6-2026-10-25** (Boiarska exclusive B production,
strict-Winkler widths, v3.5 review fixes).  Nothing in this directory is hand-made.

## Exclusion curve (the main result)

| figure | script | needs |
|---|---|---|
| `bc5_lhcb_downstream_exclusion.{png,pdf}` | `work/analysis/make_curves.py` | `outputs/data/signal_scan.json` |
| `bc5_lhcb_downstream_exclusion_refbr.{png,pdf}` | `make_curves.py --suffix _refbr` | `signal_scan_refbr.json` |
| `bc5_lhcb_downstream_exclusion_s3.{png,pdf}` | `make_curves.py --s 3 --suffix _s3` | same scan, 95 % CL |

`signal_scan.json` is produced by `work/analysis/run_scan.py` through
`work/analysis/run_sensitivity.py --adaptive` (~50 min, the expensive step).  Default threshold
s = 2.3 (90 % CL); `--s 3` gives the 95 % CL version.

## Our contour against the digitised reference (BC5.pdf)

| figure | script |
|---|---|
| `bc5_vs_reference.{png,pdf}` | `work/analysis/plot_comparison_reference.py` |
| `bc5_vs_reference_refbr.{png,pdf}` | same, `--suffix _refbr` |
| `bc5_vs_reference_s3.{png,pdf}` | same, `--s 3 --suffix _s3` |

## Decay widths and branching ratios

| figure | script | content |
|---|---|---|
| `bc5_branching_ratios.{png,pdf}` | `work/analysis/plot_branching_ratios.py` | the BC5 decay channels vs m_S (canonical mass window) |
| `bc5_branching_ratios_fullrange.{png,pdf}` | same (written alongside) | full 0.02–62 GeV range |
| `bc5_branching_ratios_with_bb.{png,pdf}` | same (written alongside) | extended to the b bbar threshold |
| `bc5_branching_ratios_senscalc.{png,pdf}` | `work/analysis/plot_branching_ratios_senscalc.py` | our BRs vs the SensCalc BrRatioScalar table |
| `bc5_branching_ratios_senscalc_all.{png,pdf}` | same (written alongside) | all channels, log range |
| `br_validation.png / .pdf` | `work/analysis/plot_br_validation.py` | per-channel ours vs the reference table, absolute widths |
| `widths_validation.png` | `work/widths/validate_widths.py` | our widths vs the digitised Winkler figure and the reference tables |
| `outputs/data/bc5_branching_ratios.csv` | `plot_branching_ratios.py` | the same numbers as the canonical figure |

All of these read `work/widths/scalar_widths.py` — i.e. the production width model — through
`widths()`; none of them recomputes a channel with a different quark-mass convention any more
(that inconsistency was removed in v3.6).

## Production parameter space (parent particles)

| figure | script |
|---|---|
| `production/production_{Bplus,B0,Bs,Kplus,K0L,Higgs}.png / .pdf` | `work/analysis/plot_production.py` |
| `production/production_all_particles.png` | same (2x3 grid) |

2-D histograms of the parent particles in the (theta, p) plane, normalised to the cross sections
of `work/magnus/samples/normalisation.json`.  Reads the gzipped parent samples and, for the
Higgs, keeps only the pre-decay copy (|status| = 62).

## Regenerating everything except the exclusion curve

```bash
python work/analysis/plot_branching_ratios.py
python work/analysis/plot_branching_ratios_senscalc.py
python work/analysis/plot_br_validation.py
python work/analysis/plot_production.py
python work/widths/validate_widths.py
python work/analysis/validate_production.py
python work/analysis/validate_hadronisation.py
python work/analysis/validate_run.py
python work/widths/audit_widths.py
```

Removed in the v3.6 cleanup (they no longer described the code): the Omnes-vs-Winkler figure and
its `work/widths/research/` scripts (the pi pi / K K rates come from the digitised Winkler
curve), `diagnostic_2gev_switch.png`, and the `_quick` / `_test` variants of the exclusion
figures with their data files.
