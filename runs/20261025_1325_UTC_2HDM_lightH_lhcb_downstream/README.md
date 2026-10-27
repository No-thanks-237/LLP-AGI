# 2HDM light-H @ LHCb Downstream — run 20261025_1325_UTC_2HDM_lightH_lhcb_downstream

Type-I 2HDM 轻 CP-even 标量 H（cos(beta-alpha)=1/tanb, m_A=m_Hpm=600 GeV, lambda v^2=0）在 LHCb
Downstream 的 90% CL 预期排除曲线，按 inputs/prompt_2HDM.md 执行。零背景、零观测；
L=25, 300 fb^-1；pp, sqrt(s)=13 TeV；扫描 0.021<=mH<=10 GeV, 5<=tanb<=1e5。

**无 FeynRules/UFO/MG5**：模型解析实现（产生用解析 2HDM 公式，衰变用论文作者的
Light_scalar_decay C++ 程序）。

## 一眼看懂

| 项 | 值 |
|---|---|
| 模型 | Type-I 2HDM，CP-even 轻 H：m_h=125、m_A=m_H±=600 GeV、λv²=0、cos(β−α)=1/tanβ |
| 探测 | LHCb Downstream：1<z_dec<2.5 m、2<η<5、反向双径迹 E>5 GeV / pT>0.5 GeV、ε_rec=0.4 |
| 统计 | 零背景零观测 s₉₀=−ln 0.1=2.3026（90% CL）；L = 25 / 300 fb⁻¹ |
| 扫描 | 0.021 ≤ m_H ≤ 10 GeV（66 点）× 5 ≤ tanβ ≤ 10⁵（40 点），三种 B→X_sH 处方 |
| **结果** | 长寿边（tanβ 上限）25 fb⁻¹：22.2–358.0；300 fb⁻¹：30.3–488.9；**扫描范围内无短寿交点** |
| 精度 | 长寿边统计误差 3.6% 中位、最大 4.3%/4.8%；边界位置误差中位 0.5% |

**关键文档**：[RUNNING.md](RUNNING.md)（怎么跑）· [DESIGN.md](DESIGN.md)（设计）·
[logs/physics_code_review.md](logs/physics_code_review.md)（全流程审查 + 第三方独立审计逐项处置）·
[logs/production_formula_audit.md](logs/production_formula_audit.md)（产生公式逐条审计）·
[logs/spectrum_vs_foresee.md](logs/spectrum_vs_foresee.md)（产生谱 vs FORESEE）·
[logs/run_report.md](logs/run_report.md)（物理报告）·
[outputs/sha256_manifest.txt](outputs/sha256_manifest.txt)（sha256 清单）

## 2026-10-26 审查轮（全流程物理/代码审查 + MC 误差修复 + 新图）

1. **产生谱 vs FORESEE**：母样本（B 强子）的 (theta,p) 前向谱与 FORESEE 公开谱逐 bin 对比，
   LHCb 母粒子接受度区间 1.5<eta<5.5 内的**积分**形状一致到 **0.4%**（逐 bin RMS 13.1%，
   89% 的 bin 在 ±20% 内——积分一致不等于全形状精度）；整体归一差 1.49 已分解为
   LO-vs-测量/NLO 的 sigma(bbbar) 约定差（用 FORESEE NLO-P8 比对同能量差 10%）。
   -> outputs/figures/production/（总览 + 每个粒子单独一张），logs/spectrum_vs_foresee.md
2. **产生道公式逐条对照**（论文 arXiv:2212.06186）：所有公式恒等式一致到 3.1e-16；
   唯一 HIGH 级发现是 B->X_s H 名义处方用 m_B 相空间（FORESEE 数值）而非论文印刷公式的
   m_b（'paper_mb' 变体给出字面结果）。-> logs/production_formula_audit.md
3. **MC 误差修复（第二版，结论修正）**：短寿交点 50-100% MC 误差的真正根因不是能量尾部，
   而是**抽样母粒子（B）的衰变顶点**：短寿 H 只有在 B 恰好飞到 1-2.5 m 衰变窗内时才可能被探测，
   P ~ exp(-1/lam_B) ~ 1e-3…1e-5，旧实现抽样的 17293 个顶点里几乎从不落入窗口（估计量无偏但
   单次实现远低于期望，且报告误差本身低估约 1.5 倍）。修复：把母粒子衰变位置**解析积分掉**
   （Rao-Blackwell：p = ∫A W dz，W = C[e^{-z/lam_P}-e^{-alpha z}]，带符号 lam_P；v3 响应表），
   并用暴力/重要性抽样验证到 0.1-0.8%（主导事例）。后果：**短寿交点消失**（低 tanb 信号比
   阈值大 1e3-1e6 倍），排除区从扫描下界 tanb=5 一直延伸到长寿边；长寿边不变（±1%）。
4. **图**：主图改为每个亮度一条闭合轮廓（300 fb^-1 实线、25 fb^-1 虚线，上下限相连，
   无竖线/散点）；衰变验证拆成三张独立平滑曲线图；产生谱对比图新增。
5. 详见 [logs/physics_code_review.md](logs/physics_code_review.md)。

## 2026-10-27 第三方独立审计（logs/project_audit_20261026.md）处置

独立审计提了 7 条（F3 已撤回）；核对后 F1/F2/F4/F5/F6/F7 全部成立，已逐条修复并**重跑三种处方**：

1. **F1（高）数值精度**：样本按 float32 载入 + boost 用 beta=p/E，高 gamma 下轻 H 会偏离质量壳
   （mH=0.021 时 26% 候选非类时）。修复：用已知母粒子质量由三动量重建 E（_boost_on_shell），
   并把末态重新置壳。修复后最大 |m^2-m_H^2|/m_H^2 = 8.8e-6、非类时 0 个。
2. **F2（高）中性 B 端点**：B0 的 recoil 应为 K0（4.782039 GeV）而非 K+（4.785663）。
   修复后 mH=4.784 处 B0 贡献严格为零（原先把已关闭的 B0->K0H 计入信号），端点区分 B+/B0；
   exclusive 变体也改为按 meson 直接求总率（不再用 "B+ 基 × 寿命比"）。
3. **F4**：新增 outputs/data/exclusion_crossings.csv（每交点一行），boundary CSV 增加 all_crossings 列。
4. **F5**：done_ok 与响应表缓存现在核对完整配置/依赖指纹（nparam/ndecay/代码/DLL/探测器/样本），
   代码指纹纳入 b_channels 等，--nsamp 真正控制分层样本量，聚合脚本校验依赖指纹与质量网格。
5. **F6**：两个 sha256 清单都能 100% 校验通过（排除自身与 pyc）。
6. **F7**：可去奇点改用精确极限（不再扰动物理参数），文档极限式补 alpha，OR 的"严格无偏"措辞收窄。
   详见 [logs/physics_code_review.md](logs/physics_code_review.md) 第 8 节。

## 主要产物

| 路径 | 内容 |
|---|---|
| outputs/figures/exclusion_2hdm_lhcb_downstream.{png,pdf} | 主图：闭合排除轮廓（300 fb^-1 实线，25 fb^-1 虚线；上=长寿边，下=扫描下界 tanb=5，因扫描范围内无短寿交点） |
| outputs/figures/production/ | pp 母粒子产生谱与 FORESEE 对比：总览 spectrum_vs_foresee_overview + 每个粒子单独 spectrum_{B0,B0bar,Bp,Bm,pip,pim,Kp,Km}（9 张 x png/pdf） |
| outputs/figures/decay_branching_ratios.{png,pdf} | BR(H->XX) vs mH，tanb=10，10 通道（平滑曲线） |
| outputs/figures/decay_partial_widths.{png,pdf} | 分宽度 vs mH，tanb=10 |
| outputs/figures/decay_total_width_ctau.{png,pdf} | Gamma_tot（tanb=5..1000）与 c*tau |
| outputs/figures/signal_significance_map.{png,pdf} | log10(s/s90) 二维显著度图 + 等值线 |
| outputs/figures/exclusion_crossings_mc.{png,pdf} | 长寿边交点的 tanb/s MC 误差诊断图（修正后短寿交点不再存在） |
| outputs/figures/exclusion_systematics_bmode.{png,pdf} | B->X_s H 三种处方（foresee/paper_mb/exclusive）的长寿边对比；变体的多交点结构见 exclusion_crossings.csv（20-30% MC，未达标） |
| outputs/figures/decay_validation.png | 旧版衰变验证叠加图（历史产物，本次已被三张新图取代） |
| outputs/figures/comparison_paper_reach.png | 与论文 FASER/FASER2 reach 的对比（context） |
| outputs/data/exclusion_boundary.csv | 长寿/短寿交点、MC 误差、EFT 有效性标记、b_mode、all_crossings（完整交点串） |
| outputs/data/exclusion_crossings.csv | **每个交点一行**（三套处方共 340 个），含 tanb/误差/rel_ds/局部显著性/resolved 标记 |
| outputs/data/exclusion_s_curves.csv | s(mH,tanb) 网格（25/300 fb^-1 及误差） |
| outputs/data/signal_scan.json | 完整逐点记录（Gamma、ctau、BR、eps、sigma、s、误差、指纹） |
| outputs/data/twohdm_widths_br.csv | tanb=10.7/105/1082 的宽度/BR/ctau 表（完整表: decay_widths_grid.csv） |
| outputs/sha256_manifest.txt | 全部输出与代码的 sha256 清单（work/analysis/make_manifest.py 生成） |
| logs/physics_code_review.md | **本次全流程审查报告**（谱对比、公式审计、MC 误差、结果） |
| logs/production_formula_audit.md | 产生道公式逐条审计（含 7 条发现 F1-F7） |
| logs/spectrum_vs_foresee.md | 产生谱 vs FORESEE 的完整数字 |
| logs/run_report.md | 详细物理报告（模型、产生、验证、结果、系统误差、限制） |

## 结果摘要（长寿边 = 物理排除上限，tanb）

| mH [GeV] | 25 fb^-1 | 300 fb^-1 |
|---|---|---|
| 0.021 | 22.2 | 30.3 |
| 0.195 | 38.8 | 52.9 |
| 0.290 (pi pi 打开) | 178 | 243 |
| 0.494 | 232 | 316 |
| 0.827 | 291 | 397 |
| 1.000 | **358** | **489** |
| 2.072 | 299 | 414 |
| 3.501 | 276 | 415 |
| 4.180 | 304 | 469 |
| 4.784 (B+→K+H 端点) | **190** | **353** |

mH>4.786 GeV 无排除能力（B->X_sH 运动学关闭；h->HH 可忽略）——如实报告。
mH=4.784 GeV 只有 B+→K+H 还开着（B0→K0H 在 4.7820 已关闭），边界因此比上一版低 18%/12%。
长寿边 MC 误差 3.6% 中位 / 4.3% 和 4.8% 最大（<=5%），边界位置误差中位 0.5%。
低 tanb 侧：修正后的估计量给出 s(tanb=5…10) 比阈值大 1e3-1e6 倍（例：mH=4.18, tanb=5:
s=2.5e4；mH=1, tanb=5: s=3.0e5），因此**扫描范围内不存在短寿交点**，排除区从 tanb=5
（扫描下界）延续到长寿边；低 tanb 侧没有被本次分析闭合，且 tanb 远低于长寿边时模型本身
落在 Br(B->X_sH)≳0.1 的 EFT 失效区（CSV 的 validity 标记）。低 tanb 各点 s 的残差 MC 误差
30-45%，建议保守乘以 ~1.5。

## 复现

依赖：Python 3.11+（本 run 用 3.14）+ numpy + matplotlib；Windows x64 预编译探测器内核
（work/response/downstream_kernel.dll）；C++ g++（仅重建衰变程序时需要）。
输入基本都在本 run 内：母样本 `work/magnus/samples/{bbar_heavy,higgs}.txt.gz` 与强子化库随
git 仓库提供；只有 500k 事例的 `minbias_pik.txt.gz`（525 MB）因 GitHub 单文件限制放在
Release `samples-2026-10` 里（`gh release download samples-2026-10 -R No-thanks-237/LLP-AGI
-p minbias_pik.txt.gz -D work/magnus/samples`）。完整网格首次重跑约 35-45 分钟（8 workers）。

    python work/analysis/launch_scan.py --workers 8 --force --nsamp '{"B":17293,"inel":2500,"h":4620}' --b-mode foresee
    python work/analysis/aggregate_scan.py --indir work/analysis/scan_2hdm --out outputs/data/signal_scan.json --b-mode foresee
    python work/analysis/final_analysis.py
    python work/analysis/final_figure.py
    python work/analysis/signal_map.py
    python work/analysis/make_manifest.py

变体：把 --b-mode foresee 换成 paper_mb / exclusive，--outdir work/analysis/scan_2hdm_<mode>。
注意：响应表 v3 记录每个 history 的母粒子飞行长度 lam_P（顶点解析积分），v2 及更早的缓存会
自动重建；完整网格（66 质量点）首次重跑约 30-35 分钟（8 workers），之后命中缓存约 5-8 分钟/变体。

## 验证（详见 logs/）

* 衰变: 与作者自己的 FORESEE 2HDM_H 表逐点比 12261 节点 x 10 通道，|log10| 中位 1.4e-6；
  论文图数字化叠加中位 0.1-0.7%；三张新图 BR 和 = 1 到 7.5e-11。
* 产生: 与独立重写的论文公式一致到 3.1e-16；与 FORESEE 解析式 bracket 比 0.7%。
* 产生谱: 与 FORESEE 前向谱在 LHCb 母粒子接受度内的**积分**形状一致 0.4%（逐 bin RMS 13.1%；
  整体归一见上）。
* 事例 OR/归一化: 与朴素 MC 玩具 0.4-1.2 sigma；s 分量重组精确到 2e-16。
* **BC5 v3.6 端到端回归**: 25 fb^-1 长寿边界 mine/published 中位 1.0144（中位偏差 2.2%，
  最大 7.6%，在 BC5 自身 MC 散布内）；本次审计修复前后逐位相同，说明未改动 BC5 共享路径。
* MC 采样: 分层抽样估计量无偏（HT 权重与 n_tot/n_smp 归一严格配套）；母粒子顶点解析积分与
  暴力 z0 平均/重要性抽样一致（主导事例 0.1-0.8%）；prompt 极限与旧 p_window 逐位相同。

## 已知限制（量化见 logs/physics_code_review.md 与 run_report.md 第 6 节）

1. B->X_sH 处方: 名义 'foresee'（m_B 相空间，与论文 FORESEE 数值一致；2026-10-26 用户确认
   主结果采用此名义处方，主图标题已注明 "FORESEE (m_B) prescription"）vs 论文印刷公式
   'paper_mb'（m_b=4.18）。长寿边比 paper_mb/foresee ≈ 1.00（mH≤4.18；mH>4.18 印刷公式无交点），
   exclusive/foresee ≈ 0.92–0.95（端点 4.784 处 0.11，因两处方在运动学末端行为不同）。
   三条曲线都给出（exclusion_systematics_bmode），不合并成误差带。
2. 低 tanb 侧未被闭合：修正后 s(tanb=5) 已远高于阈值，短寿交点落在扫描范围之外（tanb<5）；
   低 tanb 各点 s 的残差 MC 误差 30-45%（高 pz 母粒子统计），且该区域 Br(B->X_sH)≳0.1
   属 EFT 失效区。若要把它变成可引用的物理边界，需要专门的高 pz b 样本（10-30M 事例）。
3. 宽度模型按作者程序原样（无 b bbar、无带电轻子 gg 圈、alpha_em=1/137、作者 gg 9/16 缺陷，
   对本曲线无影响）；κ_γγ 未含 H± 圈（对 h 宽度 0.08%）。
4. h 通道样本缺 tH (0.17%)；单顶干涉用 kappa 括号。
5. 探测器/统计沿用 BC5（13 点 A(z) 分段线性、0.4 效率一次、无 PID/UT），不计探测器/理论系统误差。

## 目录

inputs/ (任务书、探测器快照、论文、作者源码、复用代码、sha256) · work/{decay,production,response,
analysis,magnus,review_check} · logs/ · outputs/{data,figures} · 本 README 与 DESIGN.md
