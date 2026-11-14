# BC5 暗标量 @ LHCb Downstream（v4.0）

> **版本**：4.0.0（2026-11-14）　**状态**：当前代码基线　**物理基线**：v3.6-2026-10-25（未改动）
>
> 本目录是 `bc5_lhcb_downstream_v3` 的结构重构版：把全部代码按物理用途整理成 `bc5/` 包
> （model / production / decay / response / analysis），入口脚本放进 `scripts/`，数据放进
> `work/`。**物理、数值、样本、随机种子都没有改**，重跑结果与 v3.6 逐位一致（见第 6 节）。

这是 `inputs/prompt.md` 的实现：Higgs 混合单态实标量 `S`（PBC benchmark **BC5**，
`BR(h→SS)=0.01`）在 LHCb Downstream 的 `(m_S, θ²)` 90% CL 排除曲线，`L = 25 / 300 fb⁻¹`，
零背景、零观测。

## 1. 结果

排除区在每个质量上都是 θ² 上的一个**带**：θ² 高于长寿命边时 S 在 `z_dec > 2.5 m` 之外衰变，
低于短寿边时在 `z_dec < 1 m` 之内衰变。只有一个交点落在扫描范围（θ² ≤ 1e-3）内时，
带的上侧**保持开放**（`exclusion_boundary.csv` 的 `band_open_above=True`），不伪造曲线。

代表性边界（s = 2.3，90% CL）：

| m_S [GeV] | 25 fb⁻¹ | 300 fb⁻¹ | 主导道 |
|---|---|---|---|
| 0.021 | > 2.84×10⁻⁴（开放） | > 8.18×10⁻⁵（开放） | h → SS |
| 0.302 | 3.23×10⁻⁸ … 8.79×10⁻⁴ | > 3.86×10⁻⁹（开放） | h → SS |
| 0.484 | 5.30×10⁻⁹ … 9.97×10⁻⁵ | 4.65×10⁻¹⁰ … 1.69×10⁻⁴ | h → SS |
| 1.059 | 2.83×10⁻¹⁰ … 4.47×10⁻⁶ | 2.36×10⁻¹¹ … 7.54×10⁻⁶ | h → SS |
| 1.983 | 4.51×10⁻¹⁰ … 5.97×10⁻⁶ | 3.76×10⁻¹¹ … 1.14×10⁻⁵ | h → SS |
| 3.175 | 7.00×10⁻¹¹ … 1.17×10⁻⁶ | 5.81×10⁻¹² … 2.07×10⁻⁶ | h → SS |
| 5.082 | 1.95×10⁻¹² … 2.27×10⁻⁸ | 1.62×10⁻¹³ … 4.02×10⁻⁸ | h → SS（B→X_sS 刚关闭） |

* 67 个质量点；25 fb⁻¹ 与 300 fb⁻¹ 各有 47 个点是双边、其余为开放带。
* 两个运动学台阶是物理的：m_S ≈ 0.28 GeV（ππ 打开）、m_S ≈ 5.13 GeV = m_B − m_π（B→X_sS 关闭，
  只剩 h→SS）；f₀(980) 在 0.98 GeV 附近给出可见的窄结构。
* 长寿边 MC 误差：25 fb⁻¹ 中位 4.8%、300 fb⁻¹ 中位 4.7%，最大 5.0%，0/67 超过 5%；
  边界位置误差中位 3.4%，0/134 超过 5%。短寿边仍受重尾限制（已知限制，见第 8 节）。
* 早期 s = 3（95% CL）的结果原样保留（文件名后缀 `_s3`），与 s = 2.3 来自同一次信号扫描。

## 2. 模型处方（与 v3.6 相同）

| 环节 | 处方 |
|---|---|
| S 宽度 | **strict Winkler**（arXiv:1809.01876）：色散 ππ/KK、C m³β₂π 多介子、固定质量微扰夸克道、LO gg、α_em(0) 的 γγ |
| B→X_sS | **Boiarska 逐道**（arXiv:1904.10447），B±/B⁰ 十道 + b→d，B_s 无混合道 |
| K→πS、四点道 | Boiarska；α 由 `BR(h→SS)=0.01` 反解 |
| h→SS | Pythia 8 13 TeV 混合产生 × 逐过程重加权 |
| 探测器 | 1 < z_dec < 2.5 m、2 < η_S < 5、对向双径迹 E > 5 GeV / pT > 0.5 GeV、B_y = 1 T、SciFi 3.24×2.415 m、ε_rec = 0.4 |
| 统计 | s₉₀ = 2.3（= −ln 0.1），事例级 OR，母粒子分层 + HT 权重 |

## 3. 代码结构

代码只有一处可导入入口：`bc5/` 包（每个子包对应一个物理阶段）。

```
bc5_lhcb_downstream_v4/
├─ bc5/                        Python 包（唯一代码源）
│  ├─ paths.py                 运行目录/数据位置（所有路径都从这里取）
│  ├─ model/                   S 的衰变：宽度、BR、ctau
│  │  ├─ widths.py             Winkler 模型（含数字化 ππ/KK 表）
│  │  ├─ reference.py          validation-only：从 SensCalc 表反推的对照模型
│  │  ├─ validation.py         宽度验证（图 + JSON）
│  │  └─ audit.py              逐通道审计
│  ├─ production/              S 的产生
│  │  └─ rates.py              B→X_sS、K→πS、四点道、h→SS、σ(h)
│  ├─ decay/                   S 的衰变运动学
│  │  ├─ kinematics.py         boost、精确 n 体相空间、喷注碎裂
│  │  ├─ hadronlib.py          Pythia 强子化库随机访问
│  │  └─ channels.py           道组成、带电末态、库路由
│  ├─ response/                探测器响应
│  │  ├─ kernel.py             downstream_kernel.dll 的 ctypes 封装
│  │  ├─ lifetime.py           A(z)、精确的衰变窗积分、选择
│  │  └─ native/               内核 C++ 源、预编译 DLL、build.ps1
│  └─ analysis/                统计与出图
│     ├─ scan.py               MC 扫描引擎（s(m_S, θ²)）
│     ├─ curves.py             排除边界提取与主图
│     ├─ sensitivity.py        一键驱动（scan + curves + 对比）
│     └─ plots/                5 个作图模块
├─ scripts/                    薄 CLI 入口（只做 sys.path + main()）
├─ tests/                      pytest 单元/回归测试（22 项，约 8 s）
├─ inputs/                     冻结输入：任务书、论文源、模型表、探测器头文件快照、sha256
├─ work/
│  ├─ magnus/                  Pythia 母样本、强子化库、作业脚本与 bundle
│  ├─ model/                   hadronic_table.csv + 数字化参考曲线
│  └─ checks/                  验证脚本产出的 JSON（reference_comparison 等）
├─ outputs/data/               signal_scan.json、边界 CSV、宽度表
├─ outputs/figures/            主图、参考对比、BR/宽度、产生相空间
└─ logs/                       验收记录、物理报告、v3 审查历史（见 logs/README.md）
```

### 数据流

```
inputs/  ──►  bc5.model.widths        Γ_S(m_S), BR, ctau
              │
              ├──► bc5.production.rates   BR(B→X_sS), BR(K→πS), quartic, α(h→SS)
              │
Pythia 样本 ──► bc5.analysis.scan    每个质量：抽母粒子 → S 衰变 → 探测器 → ε(m_S, θ²)
强子化库   ──►   │                        │
                  │                        └──► bc5.decay.channels + response.lifetime
                  ▼
            outputs/data/signal_scan.json
                  │
                  ▼
            bc5.analysis.curves  ──►  outputs/data/exclusion_boundary.csv + 主图
```

## 4. 快速开始

环境：Windows x64（预编译 `downstream_kernel.dll`）、Python ≥3.11、numpy + matplotlib
（`pandas` 只有 `scripts/plot_production.py` 需要；`pytest` 只有测试需要）。

```powershell
cd runs\bc5_lhcb_downstream_v4

# 1) 不重跑 MC，只从已提交的扫描重画主图与参考对比（约 5 s）
python scripts/run_sensitivity.py --curves-only

# 2) 单元/回归测试（约 8 s）
python -m pytest tests -q

# 3) 完整重跑（scan + 自适应 + 曲线 + 参考对比）
python scripts/run_sensitivity.py --adaptive --target-rel-err 0.05 \
    --adaptive-rounds 6 --nsamp-max 25000 --force
```

更多命令（参考模型变体、全部验证、重建内核、重新生成样本）见 [RUNNING.md](RUNNING.md)。

## 5. 产物

| 路径 | 内容 |
|---|---|
| `outputs/data/signal_scan.json` | 67 个质量 × 240 个 θ² 的逐点记录（含输入 sha256 指纹、nsamp、模型标签） |
| `outputs/data/exclusion_boundary.csv` | 每个质量的寿命边、MC 误差、`band_open_above`、道分解 |
| `outputs/data/exclusion_s_curves.csv` | s(m_S, θ²) 网格（25/300 fb⁻¹） |
| `outputs/figures/bc5_lhcb_downstream_exclusion.{png,pdf}` | 主图：25/300 fb⁻¹ 排除带 |
| `outputs/figures/bc5_vs_reference.{png,pdf}` | 与 BC5 参考曲线（数字化的 PDF）对比 |
| `outputs/figures/bc5_branching_ratios*.{png,pdf}` | BR vs m_S（模型、SensCalc 参考、验证） |
| `outputs/figures/production/` | 各母粒子在 (θ, p) 平面的产生谱 |
| `outputs/data/scalar_widths.csv`、`widths_audit.json` | 宽度表与逐通道审计 |
| `work/checks/*.json` | 各验证脚本的机器可读结果 |
| `outputs/figures/README.md` | 每张图的生成命令 |

## 6. 验证：与 v3.6 逐位一致

重构后用同一套输入重跑，产物对比（v4 vs v3）：

| 检查 | 结果 |
|---|---|
| 宽度 `Gamma(m_S)`（15 个质量 × 11 道） | 逐位相同 |
| 衰变产物 `s_decay_products`（102 组，含强子化库） | 逐位相同 |
| `acceptance_profile`、`p_window` | 逐位相同（A=1 解析积分到 1e-16） |
| 扫描引擎（3 个质量 × 240 个 θ²，固定 nsamp） | 7,920 个值，0 个不同 |
| **完整自适应扫描**（67 个质量 × 240 个 θ²，名义 + reference-BR 两套） | **各 176,880 个值，0 个不同**；nsamp 映射逐质量相同 |
| 自适应轮次 | base 67 质量 → 67 → 24 → 6，全部交点 ≤ 5% 后停止（与 v3 一致） |
| `exclusion_boundary.csv`、`exclusion_s_curves.csv`、`exclusion_report.json` | 逐字节相同 |
| 主图、参考对比图、BR 图、产生谱图（PNG） | 逐字节相同 |
| `reference_comparison.json`、`scalar_widths.csv`、`widths_validation.png` | 逐字节相同 |
| 22 项 pytest（相空间、寿命积分、宽度、产生率、曲线回归、输入指纹） | 全部通过 |

> PDF 只差内部时间戳元数据；PNG/CSV/JSON 的内容逐字节一致。
> `run_validation.json` 与 v3 的差异仅在于 v3 内嵌的 `hadronisation_validation` 是旧副本
> （v3 的独立文件本身与 v4 一致）；物理数字不变。

## 7. 与 v3 的关系

* v3 的 269 个一次性审查脚本（`work/review_check/`）与 93 MB 审查中间产物
  （`work/physics_review_20261025/`）**留在 v3**，没有进 v4；v4 只保留可复跑的自检
  （`tests/`、`scripts/validate_*`、`scripts/audit_widths.py`）与关键报告（`logs/`）。
* 审查结论、v3.6 验收记录和逐项修改理由仍在 v3：
  [`../bc5_lhcb_downstream_v3/logs/v36_acceptance.md`](../bc5_lhcb_downstream_v3/logs/v36_acceptance.md)、
  [`../bc5_lhcb_downstream_v3/logs/run_report.md`](../bc5_lhcb_downstream_v3/logs/run_report.md)。
* 文件级迁移表见 [CHANGELOG.md](CHANGELOG.md)。

## 8. 已知限制

* 固定夸克质量、无 b bbar、2 GeV 强子/部分子衔接、α_em = 1/137 等作者/参考程序的近似原样保留。
* 无真实背景、无探测器/理论系统误差：这是"零背景假设下的预期灵敏度"。
* γ 不作为可见末态（S→γγ 只进总宽度）。
* 短寿边 MC 误差仍大（25 fb⁻¹ 中位 25%，300 fb⁻¹ 中位 35%）：这是重尾统计问题，需要
  `>>10⁴` 母事例才能收敛，超出当前设置的 `--nsamp-max`；长寿边（物理上限）不受影响。
* 详细量化见 [`logs/run_report.md`](logs/run_report.md) 与 v3 的
  `logs/v36_acceptance.md` 第 5 节。

## 9. 文档索引

| 文档 | 内容 |
|---|---|
| [RUNNING.md](RUNNING.md) | 环境、全部命令、数据来源、指纹与合并规则 |
| [CHANGELOG.md](CHANGELOG.md) | v4.0 相对 v3.6 的逐文件迁移与验证记录 |
| [bc5/README.md](bc5/README.md) | 包结构、各阶段接口、扩展方式 |
| [inputs/README.md](inputs/README.md) | 冻结输入清单与 sha256 校验 |
| [work/README.md](work/README.md) | 样本、模型表、检查产物的位置 |
| [outputs/figures/README.md](outputs/figures/README.md) | 每张图对应哪个脚本 |
| [logs/README.md](logs/README.md) | v3 审查历史在 v4 中的位置说明 |
| [`inputs/prompt.md`](inputs/prompt.md) | 原始任务书（研究框架、验收标准） |
