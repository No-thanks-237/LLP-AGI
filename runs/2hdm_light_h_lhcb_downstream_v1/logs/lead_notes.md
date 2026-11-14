# Lead notes — 2HDM light H at LHCb downstream (run 20261025_1325)

## Authority and scope
- Task: inputs/prompt_2HDM.md. Reference model: inputs/references/FASER_light.tex.
- User instruction 2026-10-25 ~21:50: no FeynRules; the model is handled analytically.
- This run does not modify runs/bc5_lhcb_downstream_v3 or detector/; it reuses their read-only
  artifacts (copied into inputs/bc5_reuse and work/magnus/samples with sha256).

## Team split
- T1 decay (task-1): author Light_scalar_decay compile + wrapper + paper-figure validation.
- T2 production (task-2): analytic 2HDM xi_H, B/K/pi rates, Higgs kappas, h->HH.
- T3 samples (task-3): Magnus Pythia jobs (pi/K minbias + H-decay hadron library), reuse B/higgs.
- T4 response (task-4): per-channel detector response, event OR engine, BC5 regression.
- T5 model (task-5): independent analytic audit (re-implementation + Wolfram symbolic checks).
- Lead: integration, 2D scan, curves, final physics validation, report.

## Key physics expectations (for sanity checks)
- At large tanb the H is gamma-gamma dominated (BR_gaga -> 1), so only the tiny charged channels
  (pi pi, K K, 4pi below 2 GeV; ss/gg above 2; tau tau >3.55; cc >3.73) contribute to the
  downstream signature. Gamma_gamma controls the lifetime and therefore the two edges.
- B -> X_s H: inclusive paper formula, closes at m_b = 4.18 GeV; ~0.2 % at m_H=1, tanb=100.
  K -> pi H is ~20 % of the B rate for m_H < 0.35 GeV (4-quark coherence matters).
  pi -> H e nu is negligible. h -> HH is negligible at large tanb (BR ~1e-8) and at low tanb the
  H is too short-lived (prompt decay); it is kept for completeness.
- s ~ 1e3 L sigma epsilon; the long-lived edge is steep (~1/tanb^8 in the crude long-lifetime
  limit for the hadronic channels), so the MC error on s maps to a small tanb error.
- The BC5 v3.6 regression (same engine, BC5 modules) is the end-to-end acceptance check.

## Open items
- samples ETA and hadron library node coverage.
- response table performance and MC error at the working point.
- decay figure overlay result (author program). Resolve the apparent gg/gaga prefactor convention.
- production validation numbers (xi_H, B/K/pi rates, g_C).
- final single-top treatment wording in the report.

## 22:00 decisions
- B mode: nominal 'foresee', variants 'paper_mb' + 'exclusive'; all three scanned. Primary figure
  uses foresee with a systematic band; paper-Eq and exclusive curves overlaid in a validation figure.
- Decay + production modules accepted after audits (see DESIGN s.10).
- Next: wait for the pi/K minbias sample + hadron library (Magnus), then the response engine
  BC5 regression, then the 67-mass 2-D scan (parallel over masses, 32 logical cores available).

## FINAL (2026-10-25 ~22:45)
- All 5 workstreams complete (decay, production, samples, response, analytic audit).
- 66-mass x 40-tanb x 2-lumi scan done twice (B=3000 then B=6000; endpoints B=20000),
  plus paper_mb and exclusive B-prescription variants. 0 failures.
- Main result: long-lived edge 21.9-357.0 (25 fb^-1) and 29.8-487.5 (300 fb^-1) over
  0.021-4.786 GeV; MC error 3.1-3.9% (<=5%); no reach above m_B-m_K.
- Deliverables: README.md, logs/run_report.md, outputs/data/*, outputs/figures/*,
  outputs/sha256_manifest.txt.
- Caveats: short-lived crossings flagged (EFT-invalid below 2.7 GeV, 43-99% s-level MC);
  B-prescription model uncertainty -8%..0% in tanb for mH<4.2; B->XsH rate is the dominant
  model uncertainty. No detector/theory systematics included (as prompt).

## 2026-10-26 00:20 — photon visibility reverted
- Tried gamma -> e+e- conversion visibility (and a direct-photon variant) for the user; after
  checking SensCalc the user confirmed gamma is NOT a final state, so everything was reverted.
- Code and scan products restored; re-run after revert is bit-identical to the delivered
  charged-only results (max ds/s = 0, 264/264 boundary entries identical).
- Workspace cleaned: temp/ emptied (only the small FORESEE 2HDM_H subset and the author code zip
  moved to references/), run __pycache__/scratch removed, manifest regenerated.
