# outputs/figures — every figure and the script that makes it

All commands are relative to the run directory
runs/20261016_144605_UTC_bc5_lhcb_downstream.

## Exclusion curve (the main result)

| figure | script | command |
| --- | --- | --- |
| bc5_lhcb_downstream_exclusion.png / .pdf | work/analysis/make_curves.py | python work/analysis/make_curves.py |
| bc5_lhcb_downstream_exclusion_s3.png / .pdf | work/analysis/make_curves.py | python work/analysis/make_curves.py --s 3 --suffix _s3 |

make_curves.py reads outputs/data/signal_scan.json, which is produced by
work/analysis/run_scan.py (that is the expensive step, ~20 min).  The default threshold is
s = 2.3 (90% CL); --s 3 gives the 95% CL version with the _s3 suffix.

## Branching ratios

| figure | script | command |
| --- | --- | --- |
| bc5_branching_ratios.png / .pdf | work/analysis/plot_branching_ratios.py | python work/analysis/plot_branching_ratios.py |
| bc5_branching_ratios_fullrange.png / .pdf | same | same (written alongside) |
| bc5_branching_ratios_with_bb.png / .pdf | same | same (written alongside) |
| bc5_branching_ratios_senscalc.png / .pdf | work/analysis/plot_branching_ratios_senscalc.py | python work/analysis/plot_branching_ratios_senscalc.py |
| bc5_branching_ratios_senscalc_all.png / .pdf | same | same (written alongside) |

plot_branching_ratios.py draws our own widths from work/widths/scalar_widths.py; the optional
flags --xmin, --xmax, --xscale, --show-note change the axes and the note box.  The SensCalc
script reads phenomenology/darkscalar/prepared/branching_ratios.csv and computes nothing.

## Widths and validation

| figure | script | command |
| --- | --- | --- |
| widths_validation.png | work/widths/validate_widths.py | python work/widths/validate_widths.py |
| widths_pipi_omnes_vs_winkler.png / .pdf | work/widths/research/plot_omnes_vs_winkler.py | python work/widths/research/plot_omnes_vs_winkler.py |

widths_validation.png compares every computed width with the digitised Winkler curves.
The Omnes figure is a research comparison only; it is not part of the result.

## Digitisation of the reference curves (one-off, already done)

work/widths/reference_figures/digitize.py produced the CSV curves in
work/widths/reference_figures/ from the PDF figures of the Winkler paper.  It does not need to
be re-run unless the digitisation is repeated.
## Comparison with the BC5 reference figure

| figure | script | command |
| --- | --- | --- |
| bc5_vs_reference.png / .pdf | work/analysis/plot_comparison_reference.py | python work/analysis/plot_comparison_reference.py |
| bc5_vs_reference_s3.png / .pdf | same | same with --s3 |

These overlay our contour on the two blue curves of references/model/BC5.pdf, digitised into
work/widths/reference_figures/bc5_reference_curves.csv (vector route, see
bc5_digitisation_report.json).  Each reference curve is the outline of a closed excluded band, so
it is drawn as its separate branches; the CSV stores them as contiguous row blocks and the script
slices them by the ranges recorded in the JSON.  Slicing by row must use
rows[row_start-2 : row_end-1]: row_start/row_end are 1-based CSV line numbers, and an off-by-one
splices the first point of the next branch onto the previous one and draws a fake connector.
