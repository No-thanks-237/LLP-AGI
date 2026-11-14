# v4.0 变更记录（2026-11-14）

本版是 `bc5_lhcb_downstream_v3`（v3.6 验收基线）的**纯结构重构**：只移动、重命名、拆分代码与
文档。**没有修改任何物理公式、数值常数、样本、随机种子或扫描设置**。所有物理产物与 v3.6
逐位一致（第 6 节）。

## 1. 目录迁移

| v3 路径 | v4 路径 |
|---|---|
| `work/widths/scalar_widths.py` | `bc5/model/widths.py` |
| `work/widths/reference_br_model.py` | `bc5/model/reference.py` |
| `work/widths/validate_widths.py` | `bc5/model/validation.py` + `scripts/validate_widths.py` |
| `work/widths/audit_widths.py` | `bc5/model/audit.py` + `scripts/audit_widths.py` |
| `work/widths/widths_cli.py` | `scripts/widths_cli.py` |
| `work/widths/build_hadronic_table.py` | `scripts/build_hadronic_table.py` |
| `work/widths/hadronic_table*.csv/json` | `work/model/` |
| `work/widths/reference_figures/` | `work/model/reference/` |
| `work/widths/references_thesis/` | `work/model/reference/thesis/` |
| `work/analysis/bc5_production.py` | `bc5/production/rates.py` |
| `work/analysis/bc5_lib.py` | `bc5/decay/{kinematics,hadronlib,channels}.py` + `bc5/response/lifetime.py` |
| `work/analysis/run_scan.py` | `bc5/analysis/scan.py`（CLI 仍在，另有 `scripts/run_scan.py`） |
| `work/analysis/make_curves.py` | `bc5/analysis/curves.py` + `scripts/make_curves.py` |
| `work/analysis/run_sensitivity.py` | `bc5/analysis/sensitivity.py` + `scripts/run_sensitivity.py` |
| `work/analysis/plot_*.py` | `bc5/analysis/plots/*.py` + 5 个 `scripts/plot_*.py` 薄包装 |
| `work/analysis/validate_*.py` | `scripts/validate_*.py` |
| `work/response/downstream_ctypes.py` | `bc5/response/kernel.py` |
| `work/response/downstream_kernel.{cpp,dll}`、`build.ps1` | `bc5/response/native/` |
| `work/response/test_kernel.py` | `tests/kernel_suite.py`（+ `tests/test_kernel.py` pytest 包装） |
| `work/widths/README.md` | `bc5/model/README.md` |
| `work/analysis/README.md` | `bc5/analysis/README.md`（重写） |
| `work/review_check/`（269 个一次性脚本） | **留在 v3**（原始审查证据） |
| `work/physics_review_20261024/`、`physics_review_20261025/`（93 MB） | **留在 v3** |
| `inputs/README.md`（误放的项目根 README） | 重写为真正的输入清单 |

## 2. 代码改动（语义等价）

* **bc5_lib.py 拆分**：运动学/相空间 → `decay/kinematics.py`；强子化库 → `decay/hadronlib.py`；
  道组成与带电末态 → `decay/channels.py`；A(z)、寿命积分、选择 → `response/lifetime.py`。
  逐行搬移，没有改写算法；102 组衰变抽样与 v3 逐位相同。
* **导入方式**：删除所有 `sys.path.insert` 拼接与顶层重名模块（`scalar_widths`、
  `downstream_ctypes`、`run_scan` 等），改为 `bc5.*` 包内导入；路径统一走 `bc5/paths.py`。
* **探测器内核**：`kernel.py` 在 `native/` 找 DLL；`build.ps1` 用 `inputs/extension/modules`
  快照编译（不再引用 `E:/AAAwork/LLP-AGI/detector/...`），内核测试移到 `tests/`。
* **模型表自动装载**：`bc5.model.widths` 导入时从 `work/model/hadronic_table.csv` 装载
  数字化 Winkler 表；`scripts/widths_cli.py --hadronic chpt` 仍可显式切回 ChPT。
* **扫描指纹**：`_FINGERPRINT_TARGETS` 改为 8 个新源文件；`MODEL_TAG` 从
  `v3.6-2026-10-25` 改为 `v4.0-2026-11-14`（物理未变，仅代码指纹/布局变了）。
* **脚本入口**：`scripts/*.py` 都是薄包装（加 run 根到 `sys.path` 后调用包内 `main()`），
  不再包含物理代码。
* **新增测试**：`tests/` 22 项 pytest（相空间、寿命积分、宽度、产生率、曲线回归、输入指纹、
  内核套件），约 8 s 跑完。

## 3. 删除/清理

* 删除重复的 `runs/bc5_lhcb_downstream_v3.zip`（168 MB，未跟踪，与 v3 目录相同）。
* 不复制 `__pycache__`、`.dsh-ops-*` 锁文件。
* 不复制 v3 的一次性审查脚本与 93 MB 中间产物（见上表）。
* 不复制 v3 的收敛性扫描变体（`signal_scan_{base,kink*,nd*,np*,ns*,seed*}.json`，
  约 4 MB）：它们对应 v3 的 `work/review_check/` 研究，原文与脚本都留在 v3。

## 4. 文档

* 新增 `README.md`（结构 + 结果 + 验证）、`RUNNING.md`（全部命令）、本文件、
  `VERSION`、`run.json`、`bc5/README.md`、`work/README.md`、`scripts/README.md`、
  `tests/README.md`、`logs/README.md`。
* 重写 `inputs/README.md`（原文件是 2026-10-16 的项目根 README，与 inputs 无关）。
* 更新 `outputs/figures/README.md`、`bc5/model/README.md`、`bc5/analysis/README.md` 的命令与路径。
* v3 的 `logs/` 原样复制，作为历史证据；其中出现的 `work/analysis/...` 路径指 v3 布局。

## 5. 版本与命名

* 所有 run 目录统一为 `<程序>_v<版本>`：`bc5_lhcb_downstream_v1..v4`、
  `2hdm_light_h_lhcb_downstream_v1`、`deployment_smoke_v1`。
* v4 的内部版本为 **4.0.0**，模型标签 `v4.0-2026-11-14`；物理基线仍标注为 v3.6。

## 6. 回归验证（v4 vs v3.6）

同一套输入、同一种子重跑后：

| 检查 | 结果 |
|---|---|
| 宽度（15 质量 × 11 道） | 逐位相同 |
| `s_decay_products`（102 组，含库） | 逐位相同 |
| `acceptance_profile` / `p_window` | 逐位相同 |
| `exclusion_boundary.csv` / `exclusion_s_curves.csv` / `exclusion_report.json` | 逐字节相同 |
| 主图、参考对比、BR、产生谱（PNG） | 逐字节相同 |
| `reference_comparison.json` / `scalar_widths.csv` / `widths_validation.png` | 逐字节相同 |
| 22 项 pytest | 通过 |
| 完整扫描 `signal_scan.json` / `signal_scan_refbr.json` | 由 v4 代码重新生成；与 v3 逐值相同（各 176,880 个值，0 个不同；见 `logs/v4_regression.md`） |

已知的**非物理**差异：

* PDF 文件的内部创建时间元数据不同（内容相同）。
* `outputs/data/run_validation.json`：v3 内嵌的 `hadronisation_library` 是旧副本，
  v4 重新生成了当前版本；v3 的独立文件 `work/analysis/hadronisation_validation.json`
  与 v4 逐字节相同。
