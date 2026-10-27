# 如何运行 / 如何把这一版发给别人

这个目录是**自带输入、可独立运行**的：把整个 `bc5_lhcb_downstream_v3/`（约 240 MB，含 Pythia
母样本与强子化库；这些数据不在 git 仓库内）拷到任何位置（Windows x64）就能重跑分析、重画所有图。
目录内所有路径都由脚本自身的 `__file__` 推出，不依赖仓库根目录（参考数据已拷进
`inputs/phenomenology/`）。

## 1. 环境要求

| 项 | 要求 |
|---|---|
| OS | Windows x64（`work/response/downstream_kernel.dll` 是预编译的 Windows 二进制） |
| Python | 3.11+（本版在 **3.14.3** 上跑过）；只用标准库 + numpy + matplotlib |
| 可选 | `pandas` —— 只有 `work/analysis/plot_production.py` 需要 |
| 磁盘 | ~240 MB（含 1.5 M 行的 B 母样本、Higgs 母样本、67 点衰变库） |

```powershell
python -m pip install numpy matplotlib pandas
```

换机器/换 Python 版本时要重新编译探测器核（`work/response/build.ps1`，需要 g++/MSVC）；
源码在 `work/response/downstream_kernel.cpp` + `detector/` 的头文件快照。

## 2. 最快的一次自检（不重跑昂贵的扫描，约 1 分钟）

```powershell
cd bc5_lhcb_downstream_v3\work\analysis
python run_sensitivity.py --curves-only
```

它会从 `outputs/data/signal_scan.json` 重画排除曲线与参考对比，并打印 MC 验收
（应看到 `0/67 above 5%`）。

## 3. 完整重跑（约 50 分钟）

```powershell
cd bc5_lhcb_downstream_v3\work\analysis
python run_sensitivity.py --adaptive --target-rel-err 0.05 --adaptive-rounds 6 --nsamp-max 25000 --force
# 验证扫描（喂参考分支比表），另 ~50 分钟
$env:BC5_WIDTH_MODEL="reference_br"; python run_sensitivity.py --adaptive --target-rel-err 0.05 --adaptive-rounds 6 --nsamp-max 25000 --force --suffix _refbr
```

## 4. 重画所有图（排除曲线之外）

见 `outputs/figures/README.md` 最后一段的命令清单（8 条，合计约 2 分钟）。

## 5. 其它入口

| 想做什么 | 命令 |
|---|---|
| 单质量点扫描 | `python run_scan.py --only-masses 1.0,5.0 --nsamp 400` |
| 宽度表/分支比表 | `python work/widths/widths_cli.py --help` |
| 逐道产生率核对 Boiarska tab:BR | `python work/analysis/validate_production.py` |
| 宽度自洽 / 2 GeV 接口 | `python work/widths/audit_widths.py` |
| 强子化库 vs 解析模型 | `python work/analysis/validate_hadronisation.py` |
| 探测器响应/寿命积分一致性 | `python work/analysis/validate_run.py` |
| 数字化 Winkler 曲线（需要原始 PDF） | `set BC5_REFS=<放着 papers 的目录>` 然后 `python work/widths/reference_figures/digitize.py` |

## 6. 目录里哪些是"证据"，不是运行必需

| 目录 | 内容 | 可以删吗 |
|---|---|---|
| `work/magnus/` | 母样本、衰变库、Magnus 作业 bundle 与摘要 | **不能删**（样本和库是运行输入） |
| `work/review_check/` | 三轮审查的复核脚本与 JSON（我写的） | 可删，只影响复核 |
| `work/physics_review_2026102{4,5}/` | GPT 审查报告的数据产物（参考接受度表、参考谱等） | 可删，只影响复核 |
| `logs/` | 审查报告、核对记录、验收记录、运行日志 | 建议保留（说明"为什么这么做"） |
| `outputs/` | 最终图与数据（`signal_scan*.json` 是重画曲线的输入） | **不能删** |

删掉两个 review 目录后约 **185 MB**。

## 7. 部分重跑与指纹

`run_scan.py --only-masses` 会把新结果**合并**进已有的 `signal_scan.json`，合并前检查
`model_tag`、σ(bb̄)、B 产额因子**以及 12 个文件的 sha256**（4 个生产源文件 + 3 个母样本 +
5 个衰变库）。任何一项变了都会拒绝合并并提示改用 `--force-mix`（或全网格重跑）。
这是有意的：一个扫描不能混用两套代码/两套样本的结果。

## 8. 数据来源与可追溯性

* `outputs/data/signal_scan.json` 里记录了 4 个生产源文件、3 个母样本、5 个衰变库文件的
  **sha256 指纹**；`--only-masses` 部分重跑时指纹不符会拒绝合并。
* 每个母样本/衰变库都记录在 `work/magnus/samples/normalisation.json` 与 `logs/` 里，
  并对应一个 Magnus 作业 id（见 `logs/v36_acceptance.md` 第 2 节）。
* 参考曲线来自 `references/model/BC5.pdf` 的数字化结果，在
  `work/widths/reference_figures/` 下（CSV + 数字化报告）。