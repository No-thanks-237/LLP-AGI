# 如何运行 / 如何把这一版发给别人

本目录自带输入、可独立重跑（约 0.6 GB，主要是一个 500k 事例的 Pythia pi/K 样本）。
衰变与产生代码分别在 work/decay/ 与 work/production/，分析在 work/analysis/。

## 1. 环境

| 项 | 要求 |
|---|---|
| OS | Windows x64（work/response/downstream_kernel.dll 为预编译 Windows 二进制） |
| Python | 3.11+（本版在 3.14.3 上跑过）；numpy + matplotlib |
| C++ | 仅重建作者衰变程序时需要 g++（MSYS2）；work/decay/build/light_scalar_decay.exe 已含 |
| Magnus/Docker | 仅重新生成 Pythia 样本时需要；现有样本/强子化库已在 work/magnus/ 内 |

## 2. 最快自检（不重跑扫描，约 1 分钟）

```bash
python work/analysis/final_analysis.py     # 从 outputs/data/signal_scan.json 重算边界/表
python work/analysis/final_figure.py       # 重画主图（闭合轮廓）与诊断图
python work/analysis/signal_map.py         # 显著度图
python work/analysis/make_manifest.py      # 重生成 sha256 清单
```

## 3. 完整重跑（首次约 35-45 分钟，8 workers；命中响应表缓存后约 5-8 分钟）

```bash
python work/analysis/launch_scan.py --workers 8 --nsamp '{"B":17293,"inel":2500,"h":4620}' --b-mode foresee
python work/analysis/aggregate_scan.py --indir work/analysis/scan_2hdm --out outputs/data/signal_scan.json --b-mode foresee
python work/analysis/final_analysis.py
python work/analysis/final_figure.py
python work/analysis/signal_map.py
```

- 母样本抽样按 |pz| 分层（尾部普查 + Horvitz-Thompson 权重）。**--nsamp 现在真正决定样本量**：
  B_STRAT_ALLOC / H_STRAT_ALLOC 只是相对分配形状，_scaled_alloc 按 nsamp['B']/nsamp['h']
  等比缩放（默认 17293 / 4620 = 当前网格实际使用的样本量）。
- 逐质量 JSON 的 skip 判定（launch_scan.done_ok）核对 b_mode、nsamp、nparam、ndecay 和响应代码
  指纹；响应表缓存（_cache_matches）核对代码哈希、DLL、detector_config、birth_cuts、zgrid 与
  母样本 (size,mtime)，任何一项变化都会强制重建。
- 母粒子（B/K/pi）的衰变位置在引擎里**解析积分掉**（Rao-Blackwell：p = ∫A W dz，
  W = C[e^{-z/lam_P}-e^{-alpha z}]，lam_P = (pz/m)ctau 带符号；见
  logs/physics_code_review.md 第 5 节）。响应表 v3 记录每 history 的 lam_P。
- 旧 cache/*.npz（v2 及更早）会在第一次重跑时自动重建；v3 全网格首次约 30-35 分钟（8 workers），
  之后变体命中缓存约 5-8 分钟/模式。
- 扫描范围内**没有短寿交点**（低 tanb 信号远高于阈值），主图的下边界因此是扫描下界 tanb=5。
- 端点质量（4.55151, 4.784 GeV）已包含在扫描网格内，使用同一 nsamp。

## 4. 系统性变体（B->X_s H 处方）

```bash
for MODE in paper_mb exclusive; do
  python work/analysis/launch_scan.py --workers 8 --nsamp '{"B":17293,"inel":2500,"h":4620}' --b-mode $MODE --outdir work/analysis/scan_2hdm_$MODE
  python work/analysis/aggregate_scan.py --indir work/analysis/scan_2hdm_$MODE --out work/analysis/scan_2hdm_$MODE/signal_scan.json --b-mode $MODE
done
```
然后重跑 final_analysis.py（会自动读取变体并画 exclusion_systematics_bmode）。

## 5. 证据 vs 运行必需

| 目录 | 内容 | 可删? |
|---|---|---|
| work/magnus/samples/ | Pythia 母样本（B、higgs、minbias pi/K） | 不能 |
| work/magnus/hadron_lib/ | H 衰变强子化库（精确网格）+ manifest.json（缓存校验用） | 不能 |
| work/analysis/grids/ | 质量网格与 tanb 批量输入 | 不能（网格变了要重生成 library） |
| work/analysis/cache/ | 响应表缓存 | 可删（重跑时重建） |
| work/review_check/、logs/ | 审计/验证/历史记录 | 建议保留 |
| outputs/ | 最终图与数据 | 不能（signal_scan.json 是重画输入） |

## 6. 数据来源

* 衰变：作者 Light_scalar_decay（Combined_Code.zip，见 work/decay/PATCHES.md），与作者 FORESEE
  表逐点一致性见 logs/decay_validation.md。
* 产生：解析公式（见 work/production/ 与 logs/production_validation.md、logs/production_formula_audit.md）。
* 母样本：B/higgs 复用 BC5 v3.6（作业号 9d36fa12b8485053 / bc5_higgs4）；pi/K 与强子化库为本 run
  Magnus 作业 c588d44d0bd8ac7a / 5b3b8b76029ea0af 生成（logs/samples.md）。
* 参考谱：FORESEE 公开前向谱（GitHub files/hadrons，13.6 TeV Pythia8/NLO-P8/EPOS/SIBYLL）
  已下载到 references/tools/FORESEE/files/hadrons/（见 NOTE.md 与 MANIFEST.sha256）。
