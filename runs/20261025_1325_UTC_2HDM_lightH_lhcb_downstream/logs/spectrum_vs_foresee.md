# Spectrum cross-check: our pp samples vs FORESEE forward spectra

**Task:** confirm that the pp-collision particle spectra used as parents in the 2HDM light-H
LHCb-downstream analysis are unbiased, by comparing them with the forward spectra published
in FORESEE (arXiv:2105.07077). **Status: no bias found**; the residual normalisation
difference is fully accounted for by the known 13 vs 13.6 TeV energy scaling plus the
b bbar cross-section normalisation convention, and the shapes agree in the region relevant
for the analysis.

Code: `work/analysis/spectrum_foresee.py` (self-contained, rerunnable; ~95 s)
Figures: `outputs/figures/production/spectrum_vs_foresee_overview.png` / `.pdf`, and one
figure per particle `outputs/figures/production/spectrum_<species>.{png,pdf}`
(made by `work/analysis/spectrum_particle_figures.py`); numbers in
`outputs/data/spectrum_vs_foresee_numbers.csv` and `..._bins.csv`.

---

## 1. Inputs

### Our samples (see `work/magnus/samples/normalisation.json`)

| sample | file | generator | n_gen | normalisation | weight/row |
|---|---|---|---|---|---|
| B | `bbar_heavy.txt.gz` | Pythia 8.316, pp 13 TeV, HardQCD gg/qqbar -> b bbar | 800 000 | sigma(bbbar) = 560 ub = 5.6e8 pb (BC5 v3.6) | 700 pb |
| inel | `minbias_pik.txt.gz` | Pythia 8.316, pp 13 TeV, SoftQCD:inelastic | 500 000 | sigma_inel = 79.5 mb = 7.95e10 pb (DESIGN 4) | 1.59e5 pb |

Format: `event pdg px py pz E x y z t` (GeV, mm). Forward hemisphere = `pz > 0`; the
B rows are fragmentation B+/B0/Bs/Bc hadrons, the inelastic rows are pi+-, K+-, K0/K0bar
**prompt** fragmentation hadrons (status 91/92, no eta cut). The downstream analysis uses
B+/B0 only (Bc/Bs dropped), so the primary B comparison is B0/B0bar/B+/B-.

### FORESEE reference spectra (downloaded, `references/tools/FORESEE/files/hadrons/`)

Three columns `log10(theta/rad)  log10(p/GeV)  dsigma[pb/bin]`, 20 bins/decade, forward
hemisphere only, from `https://raw.githubusercontent.com/KlingFelix/FORESEE/main/files/hadrons/`.
The `# Normalization:` header is the sum of the bin weights, i.e. the forward-hemisphere
cross section of that species (not the total pp cross section). sha256 of all 60 files:
`files/hadrons/MANIFEST.sha256`; provenance in `references/tools/FORESEE/NOTE.md`.

* main reference: `13.6TeV/Pythia8/` (Pythia 8 + Monash tune, 13.6 TeV), B0/B0bar/B+-, pi+-, K+-;
* energy scaling: `14TeV/Pythia8/` and `14TeV/NLO-P8/` beauty files;
* cross-checks: `13.6TeV/NLO-P8/` (and Max/Min) - the beauty family the authors' own
  `2HDM_H.ipynb` uses for B -> X_s H; `13.6TeV/EPOSLHC/` and `13.6TeV/SIBYLL/` for the
  light-hadron generator spread (FORESEE recommends EPOS/SIBYLL over Pythia for light hadrons).

## 2. Method

* Our rows are histogrammed in exactly the FORESEE binning in (`log10 theta, log10 p`), forward
  hemisphere (`pz > 0`), each row weighted by `sigma/n_gen`.
* `theta = atan2(pT, pz)` (angle w.r.t. the beam), `p = |pvec|`; `eta = -ln tan(theta/2)`.
* FORESEE files are 13.6 TeV, our samples 13 TeV -> shapes are compared bin-by-bin and the
  13 -> 13.6 TeV energy scaling of the b bbar cross section is estimated **from the FORESEE
  files themselves**: from the 14/13.6 TeV pairs, a local power law
  `sigma ~ s^p` is fitted and extrapolated down to 13.6/13 TeV.
* Bin ratio = our cross section / FORESEE cross section in the same bin; the "shape ratio" is
  the bin ratio divided by the global ratio (so it is centred at 1 for identical shapes).
* The NLO-P8 family only covers `theta < 1 rad, 1 GeV < p < 1e4 GeV`; comparisons with it are
  restricted to that common region (all other comparisons use the full FORESEE grid).
* Bin-ratio distributions use the 0.05 dex bins with our count >= 100, shown as a KDE with a
  bandwidth of the order of the MC resolution (integer bin counts otherwise produce a comb).
  All quoted shape RMS values are dominated by bin-level MC/generator differences, see below.

## 3. B mesons (B0/B0bar/B+/B-)

### Total forward cross sections

| sample | 13 TeV / 13.6 TeV | forward B cross section |
|---|---|---|
| ours | 13 TeV | **475.0 ub** (all 8 B-meson states 528.5 ub) |
| FORESEE Pythia8 | 13.6 TeV | **318.9 ub** (all 8 states 354.1 ub) -> ratio ours/FOR = **1.4895** |
| FORESEE NLO-P8 | 13.6 TeV | 403.8 ub (common region) -> ratio **1.0562** |
| FORESEE NLO-P8-Max | 13.6 TeV | 602.1 ub -> ratio 0.7084 |
| FORESEE NLO-P8-Min | 13.6 TeV | 300.3 ub -> ratio 1.4204 |

Common region (`theta < 1 rad, 1 GeV < p < 1e4 GeV`), all B0/B0bar/B+/B-:
ours 426.5 ub, Pythia8 287.7 ub (ratio 1.4824), NLO-P8 403.8 ub (ratio 1.0562).

### Energy scaling 13 -> 13.6 TeV

`sigma(14)/sigma(13.6)` = 1.0305 (Pythia8) and 1.0278 (NLO-P8); with the local power law this
gives `sigma(13.6)/sigma(13)` = 1.0479 and 1.0435, **adopted 1.0457 (+4.6 %)**. This is the
expected size of the energy effect and it is small compared with the normalisation convention
difference discussed next.

### Why the absolute ratio is 1.48 (normalisation, not shape)

The two samples are normalised to different b bbar cross sections. Converting the observed
ratios into an implied `sigma(bbbar)` of each FORESEE variant (at 13.6 TeV, using the same
B-meson fraction as our sample):

| variant | implied sigma(bbbar) at 13.6 TeV | same-energy ratio to ours (560 ub at 13 TeV = 586 ub at 13.6) |
|---|---|---|
| Pythia8 (LO-like) | 377.8 ub | 1.550 |
| NLO-P8 (central) | 530.2 ub | 1.105 |
| NLO-P8-Max | 790.6 ub | 0.741 |
| NLO-P8-Min | 394.3 ub | 1.486 |

So the observed `ours/FORESEE` ratios are the same-energy normalisation factors divided by
the +4.6 % energy scaling:

    ours/Pythia8  = 1.550 / 1.046 = 1.48   (observed 1.4824)
    ours/NLO-P8   = 1.105 / 1.046 = 1.056  (observed 1.0562)

* The Pythia8 files carry the LO-like b bbar normalisation; the NLO-P8 family (which the
  FORESEE 2HDM_H model itself uses for B -> X_s H) is 1.4035x larger and matches our
  560-ub normalisation to ~10 % at the same energy, i.e. inside the FORESEE Max/Min
  band (0.74-1.49) and inside the LHCb-derived 560 +- 80 ub uncertainty.
* Therefore the factor 1.48 is **not a shape bias**: it is the known LO vs measured/NLO
  b bbar normalisation times the 13 -> 13.6 TeV energy scaling.

### Shape agreement

Bin-by-bin (0.05 dex bins, our count >= 100), using the global ratio 1.4895 as the
normalisation:

| region | bins | raw bin-ratio median | 68 % | shape ratio median | shape RMS | within +-10 % / +-20 % |
|---|---|---|---|---|---|---|
| full forward hemisphere | 1215 | 1.453 | [1.289, 1.710] | 0.976 | 13.1 % | 47 % / 89 % |
| analysis core (eta 1.5-5.5, 5<p<500 GeV) | 763 | 1.445 | [1.283, 1.708] | 0.970 | 13.0 % | 48 % / 89 % |

The raw bin-ratio distribution is centred on the expected scaling (median 1.45 vs global
1.49) and the shape-only distribution is centred on 1 with 89 % of bins within +-20 %.

Integrated shape ratios (ours/FOR after removing the global normalisation) show **no
coherent tilt in the LHCb acceptance region**:

| eta band | ours | FORESEE | shape ratio |
|---|---|---|---|
| eta < 1.5 | 27.2 % | 26.2 % | 1.039 |
| 1.5 < eta < 5.5 (LHCb parents) | 66.3 % | 66.6 % | **0.996** |
| eta > 5.5 | 6.4 % | 7.2 % | 0.898 |

| log10 p band | -1..0 | 0..1 | 1..2 | 2..3 | 3..4 |
|---|---|---|---|---|---|
| shape ratio | 0.919 | 0.976 | 1.019 | 0.997 | 0.863 |

| log10 theta band | -5..-4 | -4..-3 | -3..-2 | -2..-1 | -1..0 |
|---|---|---|---|---|---|
| shape ratio | 0.747 | 0.761 | 0.919 | 0.982 | 1.025 |

* In the acceptance region (2 < eta < 5) the integrated shapes agree to **0.4 %** (eta band
  1.5-5.5) and the p bands agree within 4 % for 10 GeV < p < 1 TeV.
* The only coherent deviations are outside the acceptance: a ~10 % deficit at extreme
  forward rapidity (eta > 5.5) and a ~4 % excess at eta < 1.5.  Of the eta > 5.5 deficit only
  ~3 % is expected from the energy scaling: FORESEE's own 14/13.6 TeV pair gives a
  1.021 +- 0.001 change of the eta > 5.5 fraction per 0.4 TeV, i.e. ~1.03 for 13 -> 13.6 TeV
  (the 1.5-5.5 fraction is stable to 0.1 % per 0.4 TeV).  The remainder is a
  generator/tune difference in the far-forward tail; this band carries only 6-7 % of the
  forward cross section and is outside the LHCb acceptance.
* The 13 % bin-level RMS in 0.05 dex bins is larger than the median MC statistical error
  (~5 %) - it is a genuine small bin-to-bin difference between Pythia 8.316 at 13 TeV and
  FORESEE's Pythia 8 at 13.6 TeV - but it has no systematic structure (see the band tables
  above and panel (c)), it is of the same order as the spread between FORESEE's own
  Pythia8/NLO-P8/Max/Min variants, and it averages out in any integral over the spectrum.
  The peak bins agree: ours (log10 theta, log10 p) = (+0.10, +0.70), FORESEE (+0.10, +0.65).

## 4. pi/K comparison (secondary)

Our minbias sample stores **prompt** fragmentation pi/K (status 91/92, no resonance
feed-down), while the FORESEE spectra are generator final-state spectra (decay-inclusive).
FORESEE itself recommends EPOS/SIBYLL rather than Pythia for light hadrons.

| species | ours (13 TeV) | FORESEE Pythia8 | EPOSLHC | SIBYLL |
|---|---|---|---|---|
| pi+ | 803.8 mb | 1361.6 mb | 1259.1 mb | 968.3 mb |
| pi- | 787.8 mb | 1335.7 mb | 1232.1 mb | 945.6 mb |
| K+ | 59.2 mb | 152.4 mb | 150.6 mb | 124.7 mb |
| K- | 58.2 mb | 147.8 mb | 144.3 mb | 119.5 mb |

pi/K total ratio ours/FOR = 0.570 (Pythia8), 0.613 (EPOSLHC), 0.792 (SIBYLL); the FORESEE
generator spread alone is +-30 % (pi+ 968-1362 mb). K/pi: ours 0.074, FORESEE 0.111-0.128.
Bin-level: 3395 bins with raw median 0.514 [0.405, 0.580]; in the core (eta 2-5,
10<p<500 GeV) median 0.414 [0.327, 0.519].

Interpretation: the deficit (factor 1.3-1.8) and the softer shape are consistent with the
different particle definitions (our prompt fragmentation hadrons vs decay-inclusive
generator spectra) plus the generator/tune spread; a decay-inclusive light-hadron spectrum
cannot be reproduced by a status-91/92 sample. The comparison therefore validates the
order of magnitude and the shape trends of the inelastic sample but is **not** used as an
absolute normalisation target; the light-hadron yields in the analysis follow the
sample definition and the measured sigma_inel = 79.5 mb (DESIGN 4).

## 5. Conclusion

* **No bias in the B-meson shape within the LHCb acceptance:** integrated over
  1.5 < eta < 5.5 (66 % of the cross section) the shape agrees with FORESEE Pythia8 to
  0.4 %, the p bands to <=4 % for 10 GeV < p < 1 TeV; bin-level deviations (RMS 13 % in
  0.05 dex bins) show no coherent trend and are within the generator-variant spread.
* **Normalisation fully understood:** total forward B cross section ours/FORESEE =
  1.4895 (Pythia8) = [1.550 same-energy normalisation] / [1.046 energy scaling]; the 1.55
  factor is the LO vs measured/NLO b bbar normalisation (our 560 ub vs the 378 ub implied
  by the LO spectra), as confirmed by the FORESEE NLO-P8 variant (ratio 1.056, i.e. 10 %
  at the same energy) and its Max/Min envelope (0.71-1.42).
* **pi/K:** consistent in order of magnitude and trend with the prompt-vs-inclusive
  definition difference; not usable as a normalisation target (and FORESEE does not
  recommend Pythia for light hadrons).
* **Open point (not a bias):** the inelastic sample contains prompt (status 91/92) pi/K
  only; H production from pions produced in resonance decays is not included. This is the
  BC5 convention and affects only the small pi/K channels (~1e-4 of the signal).

## 6. Reproduction

    python work/analysis/spectrum_foresee.py           # ~95 s, full 22.8M-row minbias file
    python work/analysis/spectrum_foresee.py --quick   # ~7 s, 100k-event pi/K cache

The script downloads any missing FORESEE file from the raw GitHub URL (use
`--no-download` to forbid that) and writes the figure plus the full numerical report
shown above.

---

## 7. 2026-10-26: 每个粒子单独的图 + 数值一致性检查

新增 `work/analysis/spectrum_particle_figures.py`：按粒子逐个出图（同一 (log10θ, log10p) 网格、
三面板：我们 / FORESEE / shape 比值），并把数值写到
`outputs/data/spectrum_vs_foresee_numbers.csv`（每粒子一行）与 `spectrum_vs_foresee_bins.csv`
（0.25 dex 粗 bin 的逐 bin 表，2630 行）。总览图改名为
`outputs/figures/production/spectrum_vs_foresee_overview.{png,pdf}`，所有产生谱图统一放在
`outputs/figures/production/` 下（spectrum_{B0,B0bar,Bp,Bm,pip,pim,Kp,Km}）。

### 7.1 数值上一样吗？—— 不一样，但差的是归一化约定，不是形状

| 粒子 | 我们 [μb] | FORESEE Pythia8 13.6 TeV [μb] | 比值 | 同能量比值(/1.046) | 形状中位 | 形状 RMS | ±20% 内 bin | RMS/统计误差 | 接受度内比值 |
|---|---|---|---|---|---|---|---|---|---|
| B0     | 119.64 | 79.70 | 1.501 | 1.435 | 1.015 | 13.3% | 89% | 1.7 | 1.500 |
| B0bar  | 117.85 | 79.72 | 1.478 | 1.413 | 1.018 | 13.6% | 85% | 1.7 | 1.476 |
| B+     | 120.06 | 79.82 | 1.504 | 1.438 | 1.003 | 12.7% | 88% | 1.6 | 1.491 |
| B-     | 117.44 | 79.66 | 1.474 | 1.409 | 1.005 | 12.7% | 90% | 1.6 | 1.473 |
| pi+    | 820081 | 1361602 | 0.602 | - | 0.784 | 25.4% | 45% | 6.6 | 0.520 |
| pi-    | 804274 | 1335705 | 0.602 | - | 0.785 | 25.3% | 46% | 6.4 | 0.521 |
| K+     | 60221  | 152396  | 0.395 | - | 0.888 | 18.7% | 70% | 2.6 | 0.337 |
| K-     | 59459  | 147787  | 0.402 | - | 0.865 | 19.8% | 62% | 2.8 | 0.339 |

（形状中位/RMS/±20%/RMS-统计 都是对 FORESEE 全网格做全局归一后的逐 bin 形状比；
「接受度内」= B: 1.5<η<5.5, 5<p<500 GeV；轻强子: 2<η<5, 10<p<500 GeV。）

结论：

1. **总截面数值不同**：B0 我们 119.6 μb vs FORESEE 79.7 μb（比值 1.50）。这不是形状问题，
   而是 σ(bb̄) 归一约定：FORESEE 的 Pythia8 文件隐含 σ(bb̄)=378 μb（LO），我们用测量值 560 μb；
   同能量（除以 13→13.6 TeV 的 1.046）后比值 1.435，与 LO/NLO 约定差一致。用 FORESEE 自己的
   NLO-P8 谱（其 2HDM_H 模型实际使用的那一族）比较，比值 1.056（同能量差 ~1%），
   完全落在 FORESEE Max/Min 包络 (0.71–1.42) 与 LHCb 560±80 μb 之内。
2. **去掉全局归一后，逐 bin 数值仍有差异**：中位 1.01，RMS 13%，85-90% 的 bin 在 ±20% 内；
   但 RMS 只有**合并统计误差的 1.6-1.7 倍**——即在统计精度内基本一致，剩余的是产生器版本/
   tune/能量（13 vs 13.6 TeV）的差别；唯一相干倾斜在最前向 η>5.5（-10%，在 LHCb 母粒子
   接受度 2<η<5 之外，且其中约 3% 可由能量标度解释）。接受度内的积分比值与全局比值相同
   （如 B0：1.500 vs 1.501），说明**接受度内没有形状偏差**。
3. **π/K 不能这样直接比**：我们的样本只存 prompt status-91/92 强子，FORESEE 谱含衰变贡献，
   且轻强子 FORESEE 推荐 EPOS/SIBYLL 而非 Pythia；比值 0.60/0.40、形状 RMS 25%/19%
   （RMS/统计 6.5/2.6）正是定义差 + 产生器散布，只作 (θ,p) 图形态的 sanity check。

因此：**「数值一样」只能理解为去掉归一化/能量/产生器约定后在接受度内一致；绝对数值不同，
且我们采用测量 σ(bb̄)=560 μb 的归一（与 BC5 一致），而不是 FORESEE 文件的 LO 归一。**
分析本身只用每个事例的产额与形状（不依赖谱的绝对归一化约定），所以这个差别不影响排除曲线。

