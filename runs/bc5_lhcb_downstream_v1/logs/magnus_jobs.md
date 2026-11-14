# Magnus jobs for this run

Service container: llp-magnus (image sha256:100453c2...3717e6), interface http://localhost:3011/jobs
Project-local blueprint registered for this run: **llp-job**
  source: work/magnus/llp_job_blueprint.py
  sha256: ca4d8ebd85bf396bcb98fabf49895f0c61b42583f108beebeaadd3d7ca2e086f
  image:  git.pku.edu.cn/het-agi/collider@sha256:1bf6d8db2cb9348f965c3c6650f1c6af8234c1d8be0b49040ddd55e79a33c1a5
  purpose: run a shell bundle (script.sh + inputs) inside the pinned collider image and return
           an output tarball file secret; the controlled execution entry point required by
           docs/downstream_import.md.

| Job id | Purpose | Status | Notes |
| --- | --- | --- | --- |
| 376356289cf31f66 | llp-job round-trip smoke test | Success | verified python/numpy/pythia8-config inside the tool image |
| 873591ac7789b36c | first parent generation (minbias 2e6 + bbar + cbar + higgs) | Terminated | killed after 10 min: the bbar mode was too slow and, more importantly, the three modes wrote to the same output filenames. minbias_light.txt and a partial bbar_heavy.txt were salvaged by docker cp before termination. |
| b23bd404fdee756d | second parent generation (minbias 2e5 + bbar 8e5 + cbar + higgs) | Terminated | killed after the bbar rate proved to be ~220 events/s; minbias (200 000 events, sigmaGen = 78.05 mb) and a partial bbar sample (162 530 B mesons from 373 659 events) were salvaged by docker cp before termination. |
| c9d006475b1cebd2 | higgs-only generation (7e5) | Terminated | killed at 213 506 events; the partial file was lost because the container is removed on termination. No science result was taken from it. |
| 3e427489ec85a4fb | higgs-only generation (2e5) | Success | 200 000 events, 199 105 Higgs, sigmaGen = 32.2 pb (Pythia LO); output fetched as bc5_higgs2.tar.gz |

Lesson recorded for future attempts: copy intermediates out of the job container *before*

## Second pass (2026-10-17, closing gaps 1-3)

| Job id | Purpose | Status | Notes |
| --- | --- | --- | --- |
| 7a19951f84e3b87c | first hadronisation-library attempt | Failed | script.sh started with `cd /w`, a leftover from the interactive tests, but the job already runs in /magnus/workspace/job. No output produced. |
| 25a5309c9a4b17fd | Pythia 8.316 hadronisation library | Success | 5 channels x 52 masses x 1500 hadronisations of the isolated colour-singlet S decay; downloaded as bc5_hadronlib.tar.gz and installed in work/magnus/hadron_lib/ |

New sample produced by that job:

    hadron_lib/hadron_lib_{gg,ss,cc,bb,tautau}.txt.gz
       1500 full Pythia hadronisations per mass node, stable charged particles in the
       S rest frame, columns: mass sample pid charge E px py pz
       (gg 35 nodes from 0.30 GeV, ss 27 from 1.06, cc 20 from 3.17, bb 12 from 11.13,
        tautau 19 from 3.71; nodes with no charged particles are absent)

The analysis now samples this library whenever it is present; see logs/run_report.md section 3
for the exact Pythia settings and the colour-flow requirement.
