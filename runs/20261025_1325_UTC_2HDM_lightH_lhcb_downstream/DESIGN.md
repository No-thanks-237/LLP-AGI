# 2HDM light-H @ LHCb Downstream — design and interface contract

Run id: 20261025_1325_UTC_2HDM_lightH_lhcb_downstream
Authoritative task: inputs/prompt_2HDM.md (read first). Reused BC5 reference: inputs/bc5_reuse/.
Everything lives inside this run directory. Do not modify anything outside it except read-only
references and temp/ scratch.

## 0. Conventions

* Units: GeV, metres, c = 1. hbar*c = 1.973269804e-16 GeV m.
* Every produced table records the sha256 of the files it depends on.
* Each stage writes its own validation report under logs/ and a machine-readable JSON.
* Python: the interpreter is D:\py\python.exe (3.14); numpy + matplotlib available.
* C++: g++ 15.2 (MSYS2 ucrt64). Light_scalar_decay needs -DM_PI=3.14159265358979323846.
* Magnus: deployment/container/magnus.ps1 (wraps docker exec llp-magnus magnus ...).
  Blueprint id "llp-job" is registered; submit via a python script inside the container:
    pin = magnus.custody_file(bundle_path, expire_minutes=..., max_downloads=3)
    res = magnus.run_blueprint('llp-job', args={'bundle': sec, 'job_name': name}, ...)
  Template: inputs/bc5_reuse + BC5 run work/magnus/submit_*.py. Job images:
  git.pku.edu.cn/het-agi/collider@sha256:1bf6d8db...  (Pythia 8.316, ROOT, MG5).
  On job submission record job id, bundle sha256, command, log, output sha256 in logs/.

## 1. Physics model (fixed, do not change without telling the Lead)

Type-I 2HDM, light CP-even H, benchmark of Li et al. arXiv:2212.06186 (inputs/references/FASER_light.tex):
  m_h = 125 GeV, v = 246 GeV, m_A = m_Hpm = 600 GeV, lambda v^2 = 0 (m12^2 = mH^2 sinb cosb),
  cos(beta-alpha) = c = 1/tanb, sin(beta-alpha) = s = +sqrt(1-c^2).
  xi_H^V = c, xi_H^f = c - s/tanb  (in code: c - sqrt(1-c^2)/tbeta; algebraically 1/(tbeta^3(1+s))),
  xi_h^V = s, xi_h^f = s + c/tanb.
  lambda_HH+H- = -(2 m_Hpm^2 - m_H^2)/v * c        (paper App. C, lambda v^2 = 0)
  Branching ratios / widths: the AUTHOR'S Light_scalar_decay code, model=1, mA=600, mHC=600,
  lambdav=0, tanb, cosba=1/tanb (see section 2). Never substitute another width model.
  Important: in the author program xi_H^f is computed as cosba - sqrt(1-cosba^2)/tanb (= c - s/tb),
  and gHCC = -(2 mHC^2 - mH^2 + 2 lambdav^2) cosba / v + ...  (lambdav=0: -(2mHC^2-mH^2)c/v).

Parameter plane: 0.021 <= m_H <= 10 GeV, 5 <= tanb <= 1e5, log grids.
Luminosities 25 and 300 fb^-1. Zero background, zero observed events -> 90% CL boundary s_90 = -ln 0.1
= 2.302585093.

## 2. Decay module (work/decay/)  [owner: decay]

Author source: inputs/decay_source/{src,include} (from Light_scalar_decay Combined_Code.zip).
Only allowed adaptations: (i) define M_PI for the Windows compiler; (ii) replace the main() with a
batch driver; (iii) record version/sha256/compile log. No physics edits.

Driver contract (work/decay/build/light_scalar_decay):
  stdin : one line per point  "<mH> <tanb>"
  stdout: one line per point, whitespace separated:
          mH tanb GammaTot Gamma_gaga Gamma_ee Gamma_mumu Gamma_tautau Gamma_pipi Gamma_KK
          Gamma_4pi Gamma_cc Gamma_ss Gamma_gg
          (all in GeV; use %.10e)
Python wrapper work/decay/decay_widths.py:
  widths(mH, tanb) -> dict with keys gaga,ee,mumu,tautau,pipi,KK,4pi,cc,ss,gg  (GeV)
  widths_and_br(mH, tanb) -> (dict, dict) ; br sums to 1.
  Caches batch calls (whole scan grid in one process invocation).
CHANNEL NAMES (fixed, used by response): gaga, ee, mumu, tautau, pipi, KK, 4pi, cc, ss, gg.
Validation (logs/decay_validation.md + work/decay/validation.json):
  1. BR sum = 1 everywhere; all widths >= 0; thresholds: pi pi opens 2m_pi0? note program uses
     mpi=0.13497, K K opens at 2*0.498, cc at 2*m_D=3.729686, tau tau at 2*1.7768, ss/gg at 2.
  2. Reproduce the paper figure inputs/references/mphi_t1Br_tanb102_cba=1tanb.png (tanb=10, BR vs mH,
     0.1..10 GeV) and the left panel of mphi_Type1_cba=1tanb-2.png (total width / ctau). Build an
     overlay PNG outputs/figures/decay_validation.png; digitise at least the gamma-gamma, pi pi,
     K K, 4pi, cc, ss, gg, mu mu curves at ~30 x-positions by colour and report max/median |log ratio|.
  3. Report the 2 GeV matching discontinuity of the program and the absence of b bbar (mH<=10 GeV).
  4. Print a table of (mH, tanb) -> GammaTot, ctau at the sampled validation points.

## 3. Production module (work/production/)  [owner: production]

File: work/production/twohdm_production.py. Pure model physics, no detector, no MC. Own unit tests.

### 3.1 FCNC effective coupling (paper Eq. xi_H_2HDM with App. D functions)
xi_H(i,j) = -(4 sqrt2 G_F)/(16 pi^2) sum_k V*_ki m_k^2 [ g1(x_k,x_H) c - g2(x_k,x_H) s
            - g0(x_k,x_H) (2v/m_W^2) lambda_HH+H- ] V_kj
with (i,j) a down-type pair (d,s),(d,b),(s,b); k=u,c,t. Type-I g0,g1,g2 from App. D:
  g0 = -cot^2 b (3xH^2 -4xH x + x^2 - 2x(2xH-x) ln(xH/x)) / (16 (xH-x)^3)
  g1 = -3/4 + cot^2 b x[5xH^2 -8xH x +3x^2 -2xH(2xH-x) ln(xH/x)] / (4 (xH-x)^3)
  g2 = cot b X1(x,xH) + cot^3 b X2(x,xH)
  X1 = -1/4{ xH/(xH-x) - 6/(x-1)^2 + 3 - xH(3xH-2x)/(xH-x)^2 ln xH
             + [xH(3xH-2x)/(xH-x)^2 + 3(x+1)/(x-1)^3] ln x }
  X2 = x(5xH-3x)/(4(xH-x)^2) - xH x(2xH-x)/(2(xH-x)^3) ln(xH/x)
  x_k=m_k^2/m_W^2, x_H=m_Hpm^2/m_W^2, with m_W=80.385 (record choice), m_u=2.2e-3, m_c=1.28,
  m_t=173.1 (MSbar/PDG choices recorded), CKM Wolfram parameters |V|: use PDG 2022:
    Vud=0.97373 Vus=0.2243 Vub=0.00382, Vcd=0.221 Vcs=0.975 Vcb=0.0408,
    Vtd=0.0086 Vts=0.0415 Vtb=0.999 (complex phases: for xi_H^{bbar s} the combination is
    V*_ts V_tb + V*_cs V_cb + V*_us V_ub; use the standard CKM phase convention with
    gamma=1.144 rad and Vub=|Vub| e^{-i gamma}, Vtd=|Vtd| e^{-i beta} etc. Simpler and sufficient:
    use the unitarity-consistent set from PDG and keep complex V where needed; document.
    If complex phases are omitted the interference between the t and c terms changes by <1% for
    mt >> mc; still, implement with complex V (numpy complex) so the coherent sum is exact.)
Limits to test: for tanb -> infinity xi_H^{ij} -> 0 like 1/tanb; g0,g1,g2 finite.

### 3.2 B -> X_s H (inclusive, paper Eq. after Eq.Lfceff_phi)
Br(B->XsH)/Br(B->Xc e nu) = 12 pi^2 (v^2/m_b^2) (1-mH^2/m_b^2)^2 / f(mc^2/mb^2) |xi_H^{bs}/Vcb|^2
  m_b=4.18, m_c=1.28, Br(B->Xc e nu)=0.104, f(x)=(1-8x+x^2)(1-x^2)-12x^2 ln x.
Nominal total rate = this inclusive formula (NO sum of exclusive channels on top).
Kinematics: the H momentum distribution is taken from Boiarska's exclusive-channel mixture
(work/production/b_channels.py, copied channel list from inputs/bc5_reuse/bc5_production.py
B_CHANNELS but only the b->s channels, i.e. drop the "pi" b->d entry), with the relative channel
weights evaluated at the same mH and renormalised to 1. Physically: the exclusive channels describe
the recoil system X_s that the inclusive rate integrates over; the total is the inclusive rate.
Document this as a kinematic-shape approximation and quantify it by comparing the H lab momentum
distribution and the detector acceptance with (a) all weight on B->K, (b) the channel mixture.
B_s is NOT included (Boiarska/SensCalc exclude it; BC5 v3.6 same).
Return: dict channel -> absolute BR for B+ and B0 (differ only by Gamma_B+/Gamma_B0=0.93),
        total inclusive BR, and channel masses.

### 3.3 K -> pi H (paper Eq. Matrix_4FH, coherent FCNC + four-quark, A_K=0)
M(K->pi H) = G_F^{1/2} 2^{1/4} xi_H^W [ 7 lambda (mK^2+m_pi^2-mH^2)/18 - 7 A_K mK^2/9 ]
             + (xi_H^{ds}/(2v)) m_s (mK^2-m_pi^2)/(m_s-m_d) f_0^{K pi}(q^2)
  lambda = 3.1e-7, A_K=0, f0=0.96, m_s=0.096 (or the program's value - record), m_d=0.0047.
Br = (1/Gamma_K) (2 p*/mK) |M|^2/(16 pi mK), p* = two-body momentum.
K+ -> pi+ H (mK+=0.493677, mpi+=0.13957039, Gamma_K+ = 5.314e-17 GeV)
K0L -> pi0 H (mK0=0.497611, mpi0=0.1349768, Gamma_K0L=1.287e-17 GeV):
  use the same amplitude with the neutral masses; document that the exact SU(2) coefficient for
  K0 is not given in the paper and that we use the K+-normalised isospin relation (factor 1).
K0_S -> pi0 H is CP suppressed -> 0 (as in BC5). PDG 311/-311 rows are 50% K_L, 50% K_S ->
  multiply the K0L rate by 1/2 (flavour projection).
Return absolute BRs.

### 3.4 Semileptonic pi/K -> H e nu (paper Eq. after n_h,n_l)
Br(X->H e nu) = sqrt2 G_F mX^4 |xi_H^W|^2 /(96 pi^2 m_mu^2 (1-m_mu^2/mX^2)^2)
                * Br(X->mu nu) * f(mH^2/mX^2) * (1-2n_h/(33-2n_l))^2
  n_h=n_l=3 -> (7/9)^2. Br(pi->mu nu)=0.999877, Br(K->mu nu)=0.6356.
Return absolute BRs for pi+-, K+-.

### 3.5 Higgs sector: h couplings, h -> HH, and 2HDM Higgs production cross sections
g_C = s/(2v)[ (mH^2 - 3 lambda v^2 - mh^2)(2 cot(2b) s c - c^2 + s^2) + (lambda v^2 - mH^2) ]
    with lambda v^2 = 0  (App. C, Eq. ghHH)
g_hHH^(3) = -2 g_C,  L = -1/2 g_hHH^(3) h H^2
Gamma(h->HH) = |g_hHH^(3)|^2/(32 pi mh) sqrt(1-4mH^2/mh^2)
Gamma_h^other = Gamma_h^SM * sum_i BR_i^SM * kappa_i^2   (kappa from the same H couplings:
   kappa_f = xi_h^f, kappa_V = xi_h^V, kappa_g from the top loop, kappa_gamma from the full loop)
   with Gamma_h^SM = 4.07e-3 GeV, mh = 125 GeV. BR_i^SM = {bb 0.5824, WW 0.2137, gg 0.0818,
   tautau 0.0627, cc 0.0288, ZZ 0.0262, gamma gamma 0.00227, Zgamma 0.00153, ss 0.0003,
   mumu 0.00022} (HXSWG; normalised to 1). Document the approximation.
BR(h->HH) = Gamma(h->HH)/(Gamma_h^other + Gamma(h->HH))  (NOT fixed to 0.01).
SM Higgs production cross sections at 13 TeV [pb] (same as BC5, inputs/bc5_reuse):
   ggF 48.58, VBF 3.782, WH 1.369, ZH 0.883, bbH 0.2384, ttH 0.5065, tHj 0.0740, tHW 0.0152, tHb 0.0046
2HDM scaling (kappa framework, per process):
   ggF: |kappa_g|^2 with kappa_g = [xi_h^t A12(t) + xi_h^b A12(b)] / [A12(t)+A12(b)] (top loop
        dominated; A12 = paper App. A, tau = mh^2/(4 m^2))
   VBF, WH, ZH: (xi_h^V)^2
   ttH, tHj, tHb: (xi_h^f)^2   (tHj = t-channel, tHb = s-channel: Higgs couples to the top)
   bbH: (xi_h^b)^2
   tHW: |kappa_t A_t + kappa_V A_V|^2 / |A_t + A_V|^2 -- an interference-aware form; the two
        amplitudes are NOT known analytically here, so use the conservative bracket
        sigma = sigma_SM * [kappa_t^2 + kappa_V^2 + 2 kappa_t kappa_V cos(phi)]/4 with phi=0 as
        the baseline and evaluate the maximal spread over phi in [0,pi] as an uncertainty. The
        channel is 2.7e-4 of the total; document this and, if a parton-level decomposition becomes
        available, replace it. (Do NOT apply one global kappa to all processes.)
Return dict: sigma per process (pb), sigma_total, BR_hHH, Gamma_h_other, g_C, g3, kappas.

Validation (logs/production_validation.md + work/production/validation.json):
  a. xi_H^{bs} at tanb=100, mH=1 reproduces an independent hand calculation (print all terms).
  b. Br(B->XsH) at mH->0 for tanb=100: compare with Boiarska's inclusive/exclusive relation
     (inclusive ~ 2x exclusive sum) and with the FASER-paper B-physics picture; flag the tbeta
     range where Br>0.1 (model saturated/outside validity) explicitly.
  c. K->pi H at mH->0, tanb=100: compare the FCNC-only value with Boiarska's K->pi S master
     formula rescaled by (xi_H^{ds}/xi_Boiarska_ds)^2; report the four-quark/FCNC ratio.
  d. Br(pi->He nu), Br(K->He nu) at mH->0.
  e. g_C at tanb=100, mH=1 vs the paper's approximation -mh^2 c^2/(4v); BR(h->HH) vs 4700 (g/v)^2.
  f. All kappas -> 1 as tanb -> inf; sigma_total -> 55.45 pb.

## 4. Detector response and event engine (work/response/, work/analysis/)  [owner: response]

Reuse inputs/bc5_reuse/bc5_lib.py verbatim where possible (kinematics, nbody, tau/kaon treatment,
HadronLib, transport, p_window). The only new layer is that the response must be per DECAY CHANNEL
(the 2HDM BRs depend on tanb, so the old "sample a channel with BR weights" is not reusable):
   A_ch(mH; i, z_j): for H candidate i, probability that its decay into channel ch at decay point
   z_j (z grid 1.0..2.5 m, 13 nodes) yields >=1 opposite-charge daughter pair passing
   BirthKinematics + DownstreamTransport (full envelope, By=1T 3..7m, SciFi 7.5..9.4? -> keep
   exactly the BC5 config: frontZ 7.7, backZ 9.4, halfX 3.24, halfY 2.415, hole 0.09).
   The p4/vertex of candidate i is the same as in BC5: production point + parent flight.
Then for each (mH, tanb):
   p_i(tanb) = sum_ch BR_ch(mH,tanb) * p_window(mH, pz_i, z0_i, A_ch(i), ctau(mH,tanb))
   collision OR: o_event = 1 - prod_{parents in event} (1 - BR_prod(parent type; mH,tanb) * p_i)
   eps_sample = 0.4/N_gen * sum_events o_event ;  s = 1e3 L * sum_samples sigma_sample * eps_sample
The 0.4 reconstruction efficiency is applied exactly once per collision (inside eps).
Channels and their response: gaga -> no charged daughters (A=0); ee,mumu -> two-body e/mu
(E>5 GeV cut makes them 0 below 10 GeV; still computed); tautau -> tau decays; pipi -> 1/3 pi0pi0
invisible + 2/3 pi+pi-; KK -> CP-even split 1/2 K+K- + 1/2 K0K0bar (K_S/K_L flight); 4pi -> BC5
isospin composition (4/9 4 charged, 4/9 pi+pi-pi0pi0, 1/9 invisible); cc,ss,gg -> Pythia hadron
library on the exact mass grid (new library, section 5).
Parent samples and normalisation (each sample contributes sigma_c eps_c):
   B  : bbar_heavy.txt.gz from BC5 v3 (800k events, sigma_bb=560 ub=5.6e8 pb), B+/B0 only
        (drop Bs/Bc rows), weight = Br(B->XsH) * channel-mixture sampling.
   inel: NEW Pythia SoftQCD:inelastic sample (section 5) with pi+-, K+-, K0/K0bar rows,
        sigma_inel = 79.5 mb (Pythia sigmaGen recorded; normalisation uses 79.5 mb measured).
        K+ and K0L contribute Br(K->pi H) + Br(K->H e nu) (one H per parent, mutually exclusive
        histories -> sum the BRs); pi+- contribute Br(pi->H e nu).
   h  : higgs.txt.gz from BC5 v3 (200k events, status 62, process code col 2), reweighted per
        process to the 2HDM sigma_i mixture; weight = sigma_i^2HDM/(sigma_h^2HDM * frac_i);
        o_event = BR(h->HH) * (1-(1-p1)(1-p2)) with the two H from h->HH.
   s = 1e3 L [ sigma_bb*eps_B + sigma_inel*eps_inel + sigma_h*eps_h ]  (L in fb^-1, sigma in pb)
     = 1e3 * L * sum ... ; document the 1e3 factor (1 fb^-1 = 1e3 pb^-1).
MC uncertainty: delta_s/s from the event-level variance; requirement <=5% near the boundary;
print per-mass convergence.
Cross-checks: (i) kernel unit tests from inputs/bc5_reuse/test_kernel.py; (ii) p_window vs analytic
A=1; (iii) event-OR algebra toy test; (iv) reproduce the BC5 v3.6 25 fb^-1 curve (same machinery,
BC5 model) to <5% -- this is the strongest end-to-end regression and MUST be done.

## 5. Pythia samples (work/magnus/)  [owner: samples]

Generated fresh in THIS run (record job ids, bundle sha256, sigmaGen, file sha256):
  S1. minbias: Pythia 8.316, pp 13 TeV, SoftQCD:inelastic=on, 2e6 events (adjust for statistics;
      target >=5e5 events containing at least one pi/K in |eta|<6), dump rows for pi+-, K+-,
      K0/K0bar (311/-311) at status 91/92 with x,y,z,t [mm] (same 10-column format as BC5
      parents.cc). No eta pre-cut. Use seeds recorded.
  S2. hadron library: Pythia 8.316 isolated colour-singlet H decay at rest for channels
      gg, ss, cc, tautau on the exact mass grid (work/analysis/grids/mass_grid.txt, see section 6),
      per-channel mass lists (gg/ss only m>=2.0, cc m>3.73, tau m>3.553), 1500 hadronisations
      per node, same generator as inputs/bc5_reuse/shadlib.cc (it already accepts a mass file),
      output 8-column format (mass sample pid charge E px py pz) in the H rest frame.
Reuse (copy into work/magnus/samples/ with sha256 + provenance, do not regenerate):
  R1. bbar_heavy.txt.gz (BC5 v3.6, job 9d36fa12b8485053, 800k events, sigmaGen recorded)
  R2. higgs.txt.gz (BC5 v3.6, job bc5_higgs4, 200k events, status column)
Validation: verify row counts, pdg content, event multiplicity per event, sigmaGen vs normalisation;
record everything in logs/samples.md.

## 6. Mass grid (work/analysis/grids/mass_grid.txt)  [Lead fixes it before S2]

48 log-spaced points from 0.021 to 10 GeV + physical nodes:
0.1396 (pi->H e nu closes), 0.2118 (2 m_mu), 0.2700 (2 m_pi), 0.3540 (K->pi H closes),
0.4937 (K->H e nu closes), 0.5399 (4 pi), 0.9960 (2 m_K), 1.0000, 1.1, 1.2, 2.0, 2.2,
3.553 (2 m_tau), 3.729 (2 m_D, c cbar), 4.18 (m_b, inclusive B->XsH closes), 4.784 (B->K H closes),
5.279 (m_B), 6.0, 8.0.
Merge: a log point within 0.6% (in ln m) of a node is dropped. The hadron library is generated on
exactly this list (per-channel subsets). The scan uses exactly this list. Do not edit it after S2
starts; if a node must change, regenerate the library.

## 7. Scan and outputs (work/analysis/)  [Lead + response]

launch_scan.py (grid = work/analysis/grids/scan_mass_grid.txt) --tanb-grid <log grid 5..1e5 ~40 pts>
 -> outputs/data/signal_scan.json with, for every (mH,tanb) and L:
    s, delta_s, s/B components (B, inel, h), epsilon, BR_prod used, Gamma_H, ctau, BR_hHH
make_curves.py -> contour of s(mH,tanb)=2.302585 for L=25,300; both crossings per mass;
 -> outputs/data/exclusion_boundary.csv (+ MC errors), outputs/data/exclusion_s_curves.csv
 -> outputs/figures/exclusion_2hdm_lhcb_downstream.{png,pdf} with MC error bands.
Also output width/BR/lifetime tables (outputs/data/twohdm_widths_br.csv at reference tanb).
Validation: boundary MC error <=5%; lifetime integral convergence; compare with BC5 regression;
physics review of the shape (B-channel closure at 4.18, K-channel kink at 0.354, pi at 0.1396,
KK/4pi/pipi thresholds, h->HH tail).

## 8. Task board

Tasks are on the shared board; claim before working, keep write scopes disjoint:
 decay -> work/decay/**, logs/decay*
 production -> work/production/**, logs/production*
 samples -> work/magnus/**, logs/samples*
 response -> work/response/**, work/analysis/**, logs/response*
 Lead -> work/analysis/** (grids/, engine, scan, curves), outputs/**, logs/**
Report to the Lead with: what changed, validation numbers, open questions. Never edit another
owner's files; put interfaces in work/analysis/grids/ or tell the Lead.

## 9. Lead decisions log (2026-10-25, binding)

* **No FeynRules / UFO / MG5-for-model.** User decision: "我们不用feynrules，解析处理模型".
  The 2HDM is implemented analytically (production module) and with the author's C++ decay
  program. The prompt_2HDM sentence about FeynRules/UFO is superseded. The Pythia8 HiggsSM
  sample (BC5 v3.6, 200k events) is reused for the h kinematics with per-process reweighting to
  the 2HDM cross sections (shape approximation, documented; BC5 precedent). No MG5 event sample.
* **tH/tHW/tHb channels**: they are 0.17 % of the SM h cross section and are absent from the
  Pythia8 HiggsSM sample. They are included in the *reported* 2HDM sigma_h and in the kappa
  table, but the signal estimator normalises the h-channel MC to the sum of the sampled process
  groups only (ggF, VBF, WH, ZH, bbH, ttH). The missing tH contribution is bounded separately
  (<0.2 % of an already small h->HH component) and reported. Do not fabricate kinematics for them.
* **tHW interference**: no parton-level 2HDM tool; implement a documented kappa-interference
  bracket sigma = sigma_SM*(kt^2 + kV^2 + 2 kt kV)/4 and report that the whole channel is
  2.7e-4 of sigma_h so the choice changes s by <1e-3 relative. Never use one global kappa.
* **B normalisation**: the 2HDM B channel uses sigma_bb * (B+-/B0 rows per generated event) with
  a *derived* yield correction: truth total B+-/B0/Bs = 1.918 per bb event (Pythia 8.316 truth,
  BC5 v3.6), measured Bs fraction in the sample f_Bs, so truth B+-/B0 = 1.918*(1-f_Bs) and
  corr = truth/(N_BpmB0/N_gen) (about 1.018). Apply it once to the B-channel epsilon. Keep BC5's
  own 1.0172 only inside the BC5 regression adapter.
* **K0L semileptonic**: not defined in the paper formula (which normalises to Br(X->mu nu) for
  pi+-/K+-); set to 0. K0L contributes 0.5*BR(K0L->pi0 H) (flavour projection applied by the
  response layer, NOT in the production module).
* **Hadron library channels**: gg, ss, cc, tautau only (m_H <= 10 GeV < 2 m_b, so no b bbar;
  the author program has no b bbar at all).
* **Sample names**: work/magnus/samples/{bbar_heavy,higgs,minbias_pik}.txt.gz and
  work/magnus/hadron_lib/hadron_lib_{gg,ss,cc,tautau}.txt.gz.
* **Mass grid**: work/analysis/grids/mass_grid.txt, 67 nodes, fixed before the hadron library job.
* **No FeynRules audit**: the "model" teammate is repurposed to an independent analytic audit.

## 10. Lead decisions — 2026-10-25 22:00 (after T1/T2/T5 audits)

* **Decay module ACCEPTED.** work/decay/build/light_scalar_decay.exe (sha256 5c63124c…)
  reproduces the authors' own FORESEE 2HDM_H decay tables (references/tools/FORESEE/Models/2HDM_H/model)
  on 12261 (mH,tanb) nodes x 10 channels at median |ln ratio| 3.2e-6, max 4.7e-5 (= the tables'
  printing precision); ctau agrees after the FORESEE hbar*c = 2.0e-16 convention (our 1.97327e-16,
  prompt value). Paper-figure overlay agrees at the 0.6 % level. The gamma-gamma "1/4" and gg
  "9/4" apparent factors are conventions, not errors; the author's gg is 9/16 of the standard LO
  (a real defect shared with the paper figure) but has zero impact on this curve because the
  hadronic/gg channels are 1/tanb^3 suppressed and gamma-gamma is verified to 0.4 %.
* **Production module ACCEPTED** (MODEL_TAG 2HDM-lightH-production-1.0.0, sha256 4cc687ad…).
  Independent audit agrees to 4.4e-16 in the paper_mb mode.
* **B -> X_s H prescription (FINAL).** Nominal = mode 'foresee': exact 2HDM xi_H^{bs} with the
  meson phase space (1 - mH^2/mB^2)^2, as in the authors' own FORESEE analysis (their numerical
  model is what produced the paper's reach figure). The physical mesonic endpoint mB - mK =
  4.78566 GeV is imposed (no strangeness recoil channel above it). Systematic variants, all run:
  'paper_mb' (literal paper Eq. 2.2 with m_b, closes at 4.18) and 'exclusive' (Boiarska channel
  sum with the 2HDM xi, closes per channel). At mH = 1 the foresee/paper_mb/exclusive rates are
  3.846/3.678/2.086 e-3 (x1.05/x0.57); at mH = 3 they are 1.897/0.973/1.647 e-3 (x1.95/x1.69).
  The three curves are reported; the primary figure uses 'foresee' with a systematic band.
* **BR(B0->XsH) = 0.92735 BR(B+->XsH)** (Gamma_B+/Gamma_B0), applied in the response layer.
* **Validity flags**: Br(B->XsH) > 0.1 for tanb < ~19.6 (mH=1), >1 for tanb < ~6.2; the scan
  records them per point and the report marks that region as outside the effective description.
* **Probability algebra independently validated** (p_window analytic 1.1e-16; OR toy MC
  0.4-1.2 sigma). The 13-node A(z) grid is the dominant response approximation (<=0.4 % relative).
* **No FeynRules/UFO/MG5.** The T5 workstream is an analytic audit; its FeynRules artefacts deleted.
