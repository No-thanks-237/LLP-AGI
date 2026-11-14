# T3 — Pythia samples via Magnus (2HDM light-H @ LHCb downstream)

Owner: samples (task-3). Contract: DESIGN.md sections 5-6; interface file names fixed by the Lead.
Run: `20261025_1325_UTC_2HDM_lightH_lhcb_downstream` — everything below lives inside this run.

Service: local Magnus service in container `llp-magnus` (image `sha256:100453c2…3717e6`),
blueprint `llp-job` (registered; pinned tool image
`git.pku.edu.cn/het-agi/collider@sha256:1bf6d8db2cb9348f965c3c6650f1c6af8234c1d8be0b49040ddd55e79a33c1a5`,
Pythia 8.316, ROOT, MG5). The queue is **serial** (one job Running, the rest Pending, see
`deployment/container/check_concurrency.py`); the long job was therefore launched first and the
short hadron-library job immediately behind it, so the queue never idled.

Submission command (from the host, inside the service container):

    docker exec -w /workspace -e PYTHONDONTWRITEBYTECODE=1 llp-magnus \
        python runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream/work/magnus/submit_samples.py

## 1. Jobs

| tag | job id | bundle sha256 | status | wall time | notes |
| --- | --- | --- | --- | --- | --- |
| smoke | f3b8b84cbc78334d | b227b23234c6ddedbfb8a2e37ca594093e2de189a3c6757d524cf05a3b415c81 | Success (script exit 141, diagnostic only) | ~110 s | compile probe + 20 000-event rate probe + 2-node shadlib probe. Measured **540 ev/s**, **sigmaGen = 78.0485 mb**, 45.6 pi/K rows/event. |
| pik | **c588d44d0bd8ac7a** | bfaf8806619f18fbdff3d73f84a3e4f9ffb8939c5def972b03130807216bddd6 | Success (script exit 141, diagnostic only) | **971 s** (launched 13:40:51Z, done 13:57:03Z) | `2hdm_pik`: 500 000 inelastic pp events, seed 20261025 |
| hadronlib | **5b3b8b76029ea0af** | ff7cbd921a476ae1ee9fce3fa03ab94a56f6e9b07013c428a9bb8a5135d848fd | **Success (exit 0)** | **66 s** (launched 13:40:51Z, ran 13:57:03Z → 13:58:09Z) | `2hdm_hadronlib`: 4 channels, 1500 hadronisations/node, seeds 911-914 |

**Exit code 141 note (both minbias jobs).** The job script ends with self-checks; the two commands
`zcat … | awk '…{…; exit}'` and `grep … | head -1 | awk` close the pipe early, so the upstream
`zcat`/`grep` dies with SIGPIPE and `set -o pipefail` turns that into exit 141 **after** all
science output (txt, gzip, sha256-less diagnostics) was written. This is a defect of the
diagnostic tail of *my* script, not of the generator: the gz was gzipped completely and its
integrity + line count + physics content were re-verified on download (see section 6), so the
sample was **not** regenerated. The hadron-library script uses non-exiting consumers and exited 0.
Fixed for future bundles (templates/script_pik.sh keeps the historical form; use `sed -n` instead).

Machine-readable state: `work/magnus/jobs/submitted.json` (job ids, secrets, bundle sha256),
`work/magnus/jobs/bundles.json` (bundle + staged-file sha256), `work/magnus/jobs/*_result.json`,
`work/magnus/jobs/install_report_{pik,hadronlib}.json`, `work/magnus/validation.json`,
`work/magnus/validation_interface.json`, `work/magnus/hadron_lib/manifest.json`.

## 2. Bundles (staged inputs)

`work/magnus/bundles_stage_pik/`:

| file | sha256 | content |
| --- | --- | --- |
| parents_pik.cc | b4982dbf23a261ee11d38232ad6bac23ea664ac38eed3b69abd82b3338d33fe1 | new minbias pi/K generator (adapted from inputs/bc5_reuse/parents.cc) |
| script.sh | 3acb033bbc15c9db40846c9f26f0b376680ab9519237601a51f481cefc3b79a1 | compile + run + gzip + checks |
| config.env | f84ab9a6180cc2cc087290d69e104872981616848721e17d217694fd421a6b09 | NEVENTS=500000, SEED=20261025 |

`work/magnus/bundles_stage_hadronlib/`:

| file | sha256 | content |
| --- | --- | --- |
| shadlib.cc | ebae345ef7cc76a5776b87e37b27ee25114f8bfe4d5951e98c42bf865f94ad28 | **verbatim copy** of inputs/bc5_reuse/shadlib.cc (generator unchanged) |
| mass_list_gg.txt | a87bf430103d2f72a9a2022e9a9a40de1fa5b586d7d71a5d961380831a61d6bc | 21 nodes (verbatim work/common/mass_list_gg.txt) |
| mass_list_ss.txt | a87bf430103d2f72a9a2022e9a9a40de1fa5b586d7d71a5d961380831a61d6bc | 21 nodes (same list as gg, verbatim) |
| mass_list_cc.txt | 01bdacf342d5082abaa6baebead296d48c1813f3767a1b3d78a8b16b9f1a0679 | 13 nodes |
| mass_list_tau.txt | 8e9544373033308b5a12785736f20879baf6e1baea483324bd557187642892c9 | 14 nodes |
| script.sh | e29ff65250f0c93d51f114fe18f185f3ecd41d47e561ae11ed09ab894bbbd21e | compile + 4 channel runs + gzip + node checks |
| config.env | 6a41ce4935fd058910f4c565c3ca4f96a2ea9ba88365a82462f62051a61a7d2a | NSAMPLE=1500 |

Mass-list consistency with `work/common/mass_grid.txt` (67 nodes, 0.021-10 GeV): every list is an
exact subset of the grid; the first node of each list is the first grid point above threshold
(gg/ss 2.00000 after 1.81693 < 2; cc 3.99192 after 3.72870 < 2m_D = 3.729686;
tau 3.72870 after 3.55300 < 2m_tau = 3.5536). Node counts 21/21, 21/21, 13/13, 14/14.

## 3. Generator settings (exact)

### S1 minbias (parents_pik.cc)
Identical to BC5 `parents.cc` mode minbias except one output file and the particle set:

    Beams:eCM = 13000.        Random:setSeed = on   Random:seed = 20261025
    SoftQCD:inelastic = on    SoftQCD:elastic = off
    ParticleDecays:limitTau0 = on   ParticleDecays:tau0Max = 0.0   (keep fragmentation-level hadrons)
    Print:quiet = on          Next:numberShow{Event,Info,Process} = 0

Event loop: status 91/92 only, PDG in {211,-211,321,-321,311,-311}, **no eta cut**, written as
`# event pdg px py pz E x y z t` (x,y,z in mm, t in mm/c as Pythia gives them; 10 columns, %.8g).
K0/K0bar (311/-311) are the fragmentation-level flavour states: the consumer applies the
50 % K_L / 50 % K_S treatment of DESIGN 3.3; no K_S/K_L separation is made at generation.

### S2 hadron library (shadlib.cc, verbatim BC5 generator)
`ProcessLevel:all=off, PartonLevel:all=off, HadronLevel:all=on, ParticleDecays:limitTau0=off,
Beams:eCM=13000, Check:epTolErr=1e12, Random:seed = 911/912/913/914`; for each mass node the
isolated colour-singlet pair (gg / ss / cc / tau+tau-) is injected at rest with the partner along
+-z, hadronised, all decays followed, and every **stable charged** particle is written in the H
rest frame with columns `mass sample pid charge E px py pz` (8 columns). The consumer rotates each
event by a random SO(3) matrix (bc5_lib `HadronLib`), which is required because the injected decay
axis is fixed to z.

## 4. Reused BC5 v3.6 samples (copied verbatim, NOT regenerated)

| file | sha256 | rows | n_gen | source |
| --- | --- | --- | --- | --- |
| samples/bbar_heavy.txt.gz | ac8e1a054b1af0474d80c698b25739afec283edd50c9d7ba97668c343537f54b | 1 509 528 | 800 000 | BC5 v3.6, Magnus job 9d36fa12b8485053; 1.887 B hadrons/event |
| samples/higgs.txt.gz | 7a62b1f01eec0dfa6be7b235a497f812946c1f9b6671e8e8b190ba5cb7e65dc0 | 1 411 171 | 200 000 | BC5 v3.6, Magnus job bc5_higgs4; 11 columns, 200 000 status-62 rows (one per event) |

Provenance (generator settings, pdg/process-code histograms): `samples/provenance_reused.json`.

## 5. Normalisation (`samples/normalisation.json`)

| key | file | n_gen | squared/rows | sigma [pb] | note |
| --- | --- | --- | --- | --- | --- |
| B | bbar_heavy.txt.gz | 800 000 | 1 509 528 | 5.6e8 | 560 ub pp->bbbar at 13 TeV (BC5 v3.6 value) |
| inel | minbias_pik.txt.gz | 500 000 | 22 820 428 | 7.95e10 | 79.5 mb measured inelastic pp (DESIGN 4); Pythia sigmaGen = 78.0485 mb recorded separately (1.8 % difference, not corrected; channel is ~1e-4 of the signal) |
| H | higgs.txt.gz | 200 000 | 1 411 171 | 55.4527 | sum of 13 TeV SM Higgs cross sections |

`s_threshold = 2.302585093`, `luminosities_fb = [25, 300]`, `hadron_lib = ../hadron_lib`,
channels `["gg","ss","cc","tau"]` → files `hadron_lib_{gg,ss,cc,tautau}.txt.gz`
(the decay-module channel `tau` maps to the file channel name `tautau`, exactly as in bc5_lib
`_LIB_MAP`).

## 6. Validation (numbers)

Machine-readable: `work/magnus/validation.json`, `work/magnus/validation_interface.json`.
Scripts: `validate_samples.py`, `check_lib_interface.py`.

### 6.1 S1 minbias pi/K sample (`samples/minbias_pik.txt.gz`, 525 049 850 B)

| quantity | value |
| --- | --- |
| sha256 | cc157eeef958d747da47df6c6e524b26de9eadeeeb73875678ff5cff8a9de41f |
| data rows | 22 820 428 (+2 header lines; gzip integrity OK) |
| n_gen (= max event id + 1) | 500 000 — matches `normalisation.json` |
| columns | 10 in every row, none malformed |
| pdg counts | pi+ 10 106 630, pi- 9 905 875, K+ 743 787, K- 730 920, K0 676 210, K0bar 657 006 |
| charge sum of pi/K | +213 622 on 21.0 M charged pions (0.94 % net charge, as expected for pp) |
| rows / generated event | 45.64 (K0/K0bar counted separately) |
| events with >= 1 pi/K | 490 476 (98.10 %) |
| events with >= 1 pi/K in abs(eta) < 6 | 471 041 (94.21 %) — n_gen = 500 000 was the agreed target, this is the realised | 
| abs(eta) max | 15.69 (no eta cut, as required) |
| multiplicity per event | min 1, max 260, mean 46.53 |
| Pythia sigmaGen | 78.048501 mb, sigmaErr 3.8e-7 mb |
| on-shell mass from (E,p) | 100 % for E < 100 GeV (91.2 % of rows); 63.5 % at 100-200 GeV, 39.9 % at 200-500, 14.7 % above 500 GeV |

The last line is a **text-precision limitation, not a physics defect**: with %.8g, (E-p) of a
multi-TeV pion is smaller than the quantisation step, so sqrt(E^2-p^2) is noise for the 8.8 % of
rows above 100 GeV. Those rows are the extreme-forward parents whose H cannot pass |eta_H| < 5
anyway, and for the moderate-energy rows the consumer should simply use the PDG pi/K masses instead
of sqrt(E^2-p^2) (one-line change; removes the issue entirely). BC5's parents files share the same
format and the same limitation.

### 6.2 S2 hadron library (`work/magnus/hadron_lib/`)

| channel | rows | mass nodes | node match | n_attempt | empty attempts | charge-violating groups | on-shell fraction |
| --- | --- | --- | --- | --- | --- | --- | --- |
| gg | 153 492 | 22 | **exact** | 1500 | 3.48 % | 0 | 1.000 |
| ss | 131 018 | 22 | **exact** | 1500 | 3.57 % | 0 | 1.000 |
| cc | 119 864 | 13 | **exact** | 1500 | 0.046 % | 0 | 1.000 |
| tautau | 55 932 | 14 | **exact** | 1500 | 0.00 % | 0 | 1.000 |

* Mass nodes: identical to `work/common/mass_list_{gg,ss,cc,tau}.txt` (no missing, no extra, no
  borrowed neighbouring node); `sample_index_max = 1499` → `n_attempt = 1500`, so the consumer's
  uniform draw over the 1500 generation attempts (bc5_lib `HadronLib.sample`) is exact.
* Charge conservation per (mass, sample) decay: 0 violations in all four channels.
* All stable daughters are on shell (mass fraction 1.000); no unknown PDG ids.
* Interface test: for every (scan mass, channel) pair where the channel is open — 71 pairs — the
  nearest library node is within the 0.6 % consumer tolerance (bc5_lib NODE_MTOL), so the
  `HadronLib` lookup never raises and never borrows a node from a different mass.
* `numpy.loadtxt` (the exact call in bc5_lib.HadronLib) parses all four files: 8 columns,
  mass range 2.0-10.0 (gg/ss), 3.99192-10.0 (cc), 3.7287-10.0 (tautau).
* Library file sha256: gg e4887686…32be01, ss 41b09bed…5e3a19, cc 7974eab7…50fb88,
  tautau 14b1279b…c40a600 (full values in `hadron_lib/manifest.json`).

### 6.3 Reused BC5 samples

n_gen read back from the files equals `normalisation.json` (B 800 000, H 200 000); the Higgs file
has exactly one status-62 row per generated event (200 000/200 000); 10 columns for B, 11 for H.
