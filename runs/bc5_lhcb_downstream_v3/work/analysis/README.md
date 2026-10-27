# work/analysis — the sensitivity pipeline

## Run the whole thing with one command

    cd runs/bc5_lhcb_downstream_v3
    python work/analysis/run_sensitivity.py

That is the complete sensitivity: it runs the Monte-Carlo scan over the mass grid and then
produces the exclusion curves.  About 20 minutes for the scan, seconds for the curves.

## Options

| command | what it does | time |
| --- | --- | --- |
| python work/analysis/run_sensitivity.py | full scan (52 masses, 1200 samples) and the s = 2.3 curve | ~20 min |
| python work/analysis/run_sensitivity.py --quick | 8 masses, 200 samples, ndecay 5; smoke test, writes to signal_scan_quick.json and the _quick figures and never touches the production outputs | ~2 min |
| python work/analysis/run_sensitivity.py --curves-only | re-plot from the existing outputs/data/signal_scan.json, no scan | seconds |
| python work/analysis/run_sensitivity.py --s 3 | 95% CL curve instead of the default s = 2.3 (90% CL) | seconds |
| python work/analysis/run_sensitivity.py --s3 | both thresholds in one go, the second one with the _s3 suffix | seconds |
| python work/analysis/run_sensitivity.py --nmass 20 --nsamp 600 | custom scan size | scales |
| python work/analysis/run_sensitivity.py --masses 0.5,1,2,5,10 --suffix _test | an explicit mass list, written to its own files | scales |

A preflight check runs first and lists every input it needs, so a missing file is reported
before any CPU is spent.

## Changing the mass points

Two ways:

* **Range and count**: --mmin, --mmax, --nmass build a log-spaced grid. This is what the
  production run uses (0.021 to 62.49 GeV in 52 points).
* **An explicit list**: --masses 0.5,1,2,5,10 runs exactly those masses, in GeV. They are
  sorted automatically.

Things to keep in mind:

1. **Give a --suffix when experimenting.** Any non-empty suffix (or --quick) makes the run write
   signal_scan<suffix>.json and the *_<suffix> figures, so the production outputs are never
   touched. A run with no suffix DOES overwrite them; that is intended for a genuine re-analysis.
2. **The hadronisation library has 52 fixed mass nodes** (log-spaced over the same range). For a
   requested mass the nearest node is used, at most 8% away, and the decay kinematics are those
   of the node. A finer custom grid therefore does not make the decay simulation more accurate,
   it only evaluates the efficiency at different points.
3. Masses must lie inside the model range. Below 0.021 GeV there is no scan range; the pi pi / K K
   table only covers 0.279 to 2.0 GeV, and outside it those channels are correctly zero (below
   their thresholds, and above 2 GeV where the description is perturbative).
4. The scan time is linear in the number of mass points, roughly 20 s per point at the default
   settings (1200 samples, ndecay 8).
5. The curve stage needs at least two masses to draw a curve; a single mass gives one point.

## What the two stages do

**1. Scan — run_scan.py.**  For each mass on a log grid it samples the scalar's parent
production (Higgs decays, B mesons, K mesons), draws the scalar decay with the Pythia
hadronisation library, transports every charged daughter through the LHCb downstream kernel,
and accumulates the event-level efficiency eps(m_S, theta^2).  The output is
outputs/data/signal_scan.json: per mass, per channel, the expected events at 300 fb^-1 as a
function of theta^2, together with the Monte-Carlo errors.

**2. Curves — make_curves.py.**  Turns signal_scan.json into the s = 3-sigma crossing in the
(m_S, theta^2) plane, writes outputs/data/exclusion_boundary.csv and the s-curve CSV, and draws
outputs/figures/bc5_lhcb_downstream_exclusion.png.

## Inputs (all already present in this run)

    work/magnus/samples/bbar_heavy.txt        B mesons
    work/magnus/samples/minbias_light.txt     K mesons
    work/magnus/samples/higgs.txt             SM Higgs
    work/magnus/samples/normalisation.json    n_gen and cross sections for each sample
    work/magnus/hadron_lib/hadron_lib_*.txt.gz   Pythia hadronisation per channel and mass
    work/response/downstream_kernel.dll       the compiled detector kernel
    work/widths/hadronic_table.csv            pi pi / K K widths

The Monte-Carlo normalisation is the one thing that can silently ruin the result: the number of
generated events must match the files, and run_scan.py reads it from normalisation.json
(n_gen B = 373659, K = 200000, H = 200000).  The BC5_NGEN_* environment variables override it if
a different sample is used.

## Other scripts in this directory

| file | role |
| --- | --- |
| run_scan.py | the scan stage |
| make_curves.py | the curve stage |
| run_sensitivity.py | the one-command driver for both |
| bc5_lib.py | kinematics, decay sampling, the transport interface, the acceptance integral |
| bc5_production.py | production cross sections: Higgs, B -> X_s S, K -> pi S, the quartic channels |
| validate_run.py | writes outputs/data/run_validation.json |
| validate_production.py | production cross sections against Boiarska |
| plot_branching_ratios.py | the branching-ratio figure (our own widths) |
| plot_branching_ratios_senscalc.py | the branching-ratio figure from the SensCalc tables |
