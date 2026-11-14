# research/ — not part of the production path

These modules are kept for the record.  widths_cli.py and scalar_widths.py do not import
anything from this directory, and the production pi pi / K K widths come from
../hadronic_table.csv (the digitised Winkler figure), not from here.

## Contents

* **phases.py** — the isoscalar S-wave pi pi phase shift and inelasticity in closed form, from
  Garcia-Martin, Kaminski, Pelaez, Ruiz de Elvira, Yndurain, PRD 83 (2011) 074004, Appendix A,
  Eqs. (A1)-(A4) and Table V, CFD column.  Valid for 2 m_pi <= sqrt(s) <= 1.42 GeV.

  Verified: delta and its derivative are continuous at the 0.85 GeV matching point, and
  delta(4 m_K^2) = 226.5 deg is exactly the paper's own d_0.

* **omnes.py** — (i) the two-channel Muskhelishvili-Omnes machinery transcribed from Winkler
  arXiv:1809.01876 and Hoferichter et al. arXiv:1204.6251 (T matrix, unitarity relation,
  form-factor assembly, rate equations).  The two-channel SOLVER IS DISABLED: it raises
  NotImplementedError rather than returning unvalidated numbers.  (ii) omnes_single(), a
  validated single-channel Omnes function for the elastic region.

## What this showed, and why the production path does not use it

Comparing omnes_single() driven by the CFD phase shift with the digitised Winkler curve
(outputs/figures/widths_pipi_omnes_vs_winkler.png) gives agreement at m_S ~ 0.7 GeV, 20-25% low
below 0.5 GeV and up to a factor 3 high near 0.95 GeV.  Two causes, both physical:

1. The phase-shift input disagrees near the f0(980).  CFD gives delta(4 m_K^2) = 226.5 deg;
   Caprini-Colangelo-Leutwyler quote 185 +- 10 deg at the same point and Winkler's own figure
   shows 183 deg at 1 GeV.  This is the known dip-versus-non-dip ambiguity of the S0
   inelasticity, and the CCL curve is not published as a table.

2. The coefficient p = 0.73 of the trace-anomaly form factor cannot be fixed in the
   single-channel approximation: it is defined through Omega'(0) of the full two-channel
   solution, including Omega_12'(0) which is non-zero even though Omega_12(0) = 0.  The p s term
   dominates the rate above ~0.7 GeV.

A faithful independent reproduction is therefore not possible from published inputs alone: the
pi pi -> K K amplitude is a private communication in Hoferichter et al. (their footnote 3), and
the phase-shift set Winkler actually used is not published as a curve.

This is the same conclusion the benchmark phenomenology paper reaches.  Boiarska et al.
arXiv:1904.10447 write that "the reliability of the dispersion relation method is questionable
for m_S >~ 1 GeV because of lack of the experimental data on meson scattering at high energies",
use the Winkler figure as their benchmark instead, and quote the pi pi / K K width as uncertain
by more than an order of magnitude around 1 GeV.

## Reference material

* references_thesis/OMNES_RECIPE.md — all 89 equations of the two-channel Omnes construction in
  Hoferichter et al. (arXiv:1204.6251), extracted verbatim from the LaTeX source, with the
  equation map, the numerical method as stated, and the list of inputs that are needed but not
  printed.
* references_thesis/arXiv-1204.6251v2/ — that paper's full LaTeX source.
* references_thesis/_eps_png/ — previews of its Omega_11, Omega_12, Omega_21, Omega_22 and
  det Omega plots, the natural validation targets for a future two-channel implementation.
