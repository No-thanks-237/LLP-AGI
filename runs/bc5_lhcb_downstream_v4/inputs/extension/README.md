# detector/downstream extension snapshot (this run)

Snapshot of `detector/downstream/` taken at the start of
`runs/20261016_144605_UTC_bc5_lhcb_downstream`, plus the interface notes required by
`docs/downstream_import.md` before `verify_adapter.sh` can be used.

## Contents

| Path | Purpose |
| --- | --- |
| `modules/DownstreamTransport.h` | Header-only deterministic dipole (By = 1 T, 3 < z < 7 m) + SciFi envelope acceptance kernel. No ROOT dependency. |
| `modules/BirthKinematics.h` | Strict birth-point truth cuts, E > 5 GeV and pT > 0.5 GeV. |
| `modules/DownstreamPropagator.{h,cc}` | Delphes module wiring the kernel to `Track` output (route B only). |
| `cards/delphes_card_downstream.tcl` | Production card: full envelope, no UT requirement, no smearing/PID. |
| `cards/delphes_card_downstream_geometry_only.tcl` | Geometry-only test card; must not be used for physics acceptance. |
| `build.sh`, `verify_adapter.sh`, `tests/` | Isolated build, interface tests and packaging. |

## Interface contract used by this run

* Units: GeV for four-momenta, metres for positions **inside the kernel**; HepMC input uses GEV/MM.
* `downstream::Config`: by=1.0, fieldZMin=3.0, fieldZMax=7.0, frontZ=7.7, backZ=9.4,
  halfX=3.24, halfY=2.415, holeRadius=0.09, mode=FullEnvelope.
* `downstream::BirthCuts`: minEnergy=5.0, minPT=0.5 (strict `>`).
* `transport()` judges a single track only. Decay region, S ancestry, opposite-charge pairing,
  collision-level OR, the lifetime integral and the single event-level 0.4 efficiency are the
  responsibility of the response/analysis layer, as `docs/downstream_import.md` requires.
* A track that turns inside the dipole or is tangent at the dipole exit is rejected.
* Vertex mm -> m conversion happens in the calling layer, never inside the kernel.

## Route actually used in this run

Route A (semi-analytic) was taken one step further: the two header files are compiled into a
small ctypes shared library (`work/response/downstream_kernel.cpp`, built by
`work/response/build.ps1`) so that the *same, unmodified* kernel code evaluates responses from
the Python analysis. No detector source in this snapshot was modified. The Delphes route
(route B) was not exercised in this run and remains available for cross-validation.

## Version / provenance

Source snapshot of `detector/downstream` at the start of the run; SHA-256 of every file is
recorded in `inputs/source.sha256`. The upstream repository files themselves were not changed.
