# work/widths — decay widths of the BC5 dark scalar

A self-contained program that produces every decay width, branching ratio and lifetime of the
Higgs-portal real scalar S, given its mass.  It does not depend on the rest of the analysis
pipeline: everything it needs is in this directory.

    cd work/widths
    python widths_cli.py --mass 1.0                    # one mass point, human readable
    python widths_cli.py --mass 1.0 --theta2 1e-9      # at a given mixing angle
    python widths_cli.py --mass 2.5 --json             # machine readable
    python widths_cli.py --scan 0.021 62.49 52         # log-spaced mass scan
    python widths_cli.py --check                       # validate against the references

Python 3 with numpy only.

## Model

A real scalar S mixing with the Higgs with angle theta, plus the quartic coupling alpha:

    L = -theta (m_f/v) S fbar f  -  theta (2 m_W^2/v) S W+W-  -  theta (m_Z^2/v) S ZZ
        -  alpha S^2 |H|^2

so for every channel with a Standard-Model Higgs expression

    Gamma(S -> X) = theta^2 * Gamma(h -> X) |_{m_h -> m_S}

and the widths below are all quoted at theta^2 = 1 (except the alpha-induced ones, which are
theta independent and are not part of this module).  Reference: Winkler arXiv:1809.01876 for
the decay rates, Boiarska et al. arXiv:1904.10447 for the perturbative treatment.

## Two regimes

* **m_S < 2 GeV** — the hadronic width is taken from the dispersive analysis of Winkler:
  pi pi and K K from the digitised Fig. 2 (see PROVENANCE below), plus the multi-meson toy
  channel C m_S^3 beta_2pi with C = 5.1e-9 GeV^-2, chosen so that the total width is continuous
  at the 2 GeV matching scale.
* **m_S >= 2 GeV** — the hadronic width is perturbative: s, c, b quarks and gluons, with the
  running quark masses evaluated at mu = m_S and the kinematic threshold of each channel set by
  the lightest hadron containing that quark (K, D, B).  The higher-order QCD correction of
  Boiarska is NOT applied; see QUARK_NLO in scalar_widths.py and section 4 of the run report.

Leptons and photons are present in both regimes.

## Files

| file | role |
| --- | --- |
| **widths_cli.py** | the driver: argument parsing, the report, the provenance table |
| **scalar_widths.py** | the physics: alpha_s and running masses, loop functions, every channel, the assembly |
| **hadronic_table.csv** | the pi pi and K K widths (theta^2 = 1) read off Winkler Fig. 2 |
| **build_hadronic_table.py** | how that table was built, including the self-calibration |
| **hadronic_table_report.json** | the calibration record |
| **validate_widths.py** | validation: against the Winkler curves and the SensCalc tables |
| **reference_comparison.json** | the output of that validation |
| **reference_figures/** | the digitised reference curves and the digitisation code |
| **research/** | Omnes machinery and phase shifts — NOT used in production, see research/README.md |
| **references_thesis/** | the Hoferichter et al. LaTeX source and the extracted equation list |

## PROVENANCE of each channel

| channel | how it is computed | independent? |
| --- | --- | --- |
| S -> e+e-, mu+mu-, tau+tau- | closed form, G_F m_S m_l^2 beta^3 / (4 sqrt2 pi) | yes |
| S -> gamma gamma | SM h -> gamma gamma loop functions with m_h -> m_S | yes |
| S -> g g (m_S >= 2 GeV) | alpha_s^2 m^3 abs(sum F_1/2)^2 / (128 pi^3 v^2) | yes |
| S -> q qbar (m_S >= 2 GeV) | N_c m_S mbar_q^2(m_S) / (8 pi v^2) * beta^3, LO spectator | yes |
| S -> multi-meson | C m_S^3 beta_2pi, C from Winkler | formula, constant from the literature |
| **S -> pi pi, K K (m_S < 2 GeV)** | **digitised from Winkler Fig. 2** | **no — external published input** |

Everything except the pi pi / K K resonance region is computed here from closed-form
expressions.  The pi pi / K K input is not a computation: it is the published dispersive result,
and an independent reproduction is not possible from published inputs alone (the pi pi -> K K
amplitude used by Hoferichter et al. is a private communication and the phase-shift set Winkler
used is not published as a curve).  The same choice is made by the benchmark phenomenology
paper, which quotes this channel as uncertain by more than an order of magnitude around 1 GeV.

## Validation

python widths_cli.py --check replays work/widths/validate_widths.py and writes
reference_comparison.json.  Headline numbers:

| check | result |
| --- | --- |
| S -> mu mu against the exact analytic curve | ratio 1.000, scatter 1e-4 |
| S -> gamma gamma, SM limit at m = 125 GeV | 9.1e-6 GeV reproduced |
| S -> g g against Winkler Fig. 2 | 6% (median), better than 3% above 4 GeV |
| S -> q qbar against Winkler Fig. 2 | strange 3.5% at 2 GeV; charm ~30% low at 5 GeV, within 10% of SensCalc |
| alpha_s(m_Z) | 0.1179 by construction, 3-loop, 3 flavour thresholds |
| heavy-top limit of Gamma(S -> g g) | alpha_s^2 m^3 / (72 pi^3 v^2) reproduced |
| hadronic width continuity at the 2 GeV matching scale | 1.15 (Winkler 1.0002, SensCalc 1.0000) |
| SensCalc cross-reference | leptons < 0.1%; see reference_comparison.json for the rest |
