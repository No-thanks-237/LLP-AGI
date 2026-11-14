# 如何运行 BC5 v4

这个目录是**自带输入、可独立运行**的：把整个 `bc5_lhcb_downstream_v4/`（约 170 MB，含 Pythia
母样本与强子化库）拷到任何位置（Windows x64）都能重跑分析和重画所有图。目录内所有路径都由
`bc5/paths.py` 从包自身的 `__file__` 推出，不依赖仓库根目录。

## 1. 环境

| 项 | 要求 |
|---|---|
| OS | Windows x64（`bc5/response/native/downstream_kernel.dll` 是预编译二进制） |
| Python | 3.11+（本版在 **3.14.3** 上跑过）；标准库 + numpy + matplotlib |
| 可选 | `pandas`（只有 `scripts/plot_production.py` 需要）；`pytest`（只有 `tests/` 需要） |
| 磁盘 | ~170 MB（含 116 MB 母样本与 29 MB 强子化库） |

```powershell
python -m pip install numpy matplotlib pandas pytest
```

换机器/换 Python 版本时重建探测器核：

```powershell
pwsh -File bc5/response/native/build.ps1
# 需要 g++（MSYS2 ucrt64）或 -Compiler <g++.exe>；头文件用 inputs/extension/modules 快照
```

## 2. 最快的自检（不重跑 MC，约 10 s）

```powershell
cd runs\bc5_lhcb_downstream_v4

python -m pytest tests -q                      # 22 项：相空间、寿命积分、宽度、产生率、曲线回归、指纹
python scripts/run_sensitivity.py --curves-only # 从已提交的 signal_scan.json 重画主图
```

`--curves-only` 会打印长寿边的 MC 误差，应看到 `0/67 above 5%`。

## 3. 完整重跑

### 3.1 名义模型（Winkler 宽度，约 50–90 分钟）

```powershell
cd runs\bc5_lhcb_downstream_v4
python scripts/run_sensitivity.py --adaptive --target-rel-err 0.05 `
    --adaptive-rounds 6 --nsamp-max 25000 --force
```

它依次执行：预检 → 67 个质量的基扫描（`--nsamp 1200 --ndecay 8`）→ 自适应加密
（只重跑交点 MC 误差 > 5% 的质量）→ 曲线（s = 2.3）→ 与数字化参考曲线对比。

### 3.2 参考宽度模型（validation，另 ~50–90 分钟）

```powershell
$env:BC5_WIDTH_MODEL = "reference_br"
python scripts/run_sensitivity.py --adaptive --target-rel-err 0.05 `
    --adaptive-rounds 6 --nsamp-max 25000 --force --suffix _refbr
Remove-Item Env:BC5_WIDTH_MODEL
```

### 3.3 95% CL 曲线

```powershell
python scripts/make_curves.py --s 3 --suffix _s3
python scripts/plot_comparison_reference.py --s 3 --suffix _s3
```

## 4. 重画图 / 重跑验证（不需要重跑 MC）

```powershell
python scripts/plot_branching_ratios.py
python scripts/plot_branching_ratios_senscalc.py
python scripts/plot_br_validation.py
python scripts/plot_production.py
python scripts/validate_widths.py        # -> work/checks/reference_comparison.json + widths_validation.png
python scripts/audit_widths.py           # -> outputs/data/widths_audit.json
python scripts/validate_production.py    # -> work/checks/production_validation.json
python scripts/validate_hadronisation.py # -> work/checks/hadronisation_validation.json
python scripts/validate_run.py           # -> outputs/data/run_validation.json
```

`validate_run.py` 会读取 `work/checks/` 里前四个脚本的 JSON，所以先跑它们。

## 5. 其它入口

| 想做什么 | 命令 |
|---|---|
| 单质量点扫描 | `python scripts/run_scan.py --masses 1.0,5.0 --nsamp 400 --ndecay 5` |
| 宽度/BR 表（单点或扫描） | `python scripts/widths_cli.py --mass 1.0 --provenance` |
| 只用 ChPT 的 ππ/KK（对比） | `python scripts/widths_cli.py --mass 1.0 --hadronic chpt` |
| 逐道产生率核对 Boiarska tab:BR | `python scripts/validate_production.py` |
| 宽度自洽 / 2 GeV 接口 | `python scripts/audit_widths.py` |
| 强子化库 vs 解析模型 | `python scripts/validate_hadronisation.py` |
| 探测器响应/寿命积分一致性 | `python scripts/validate_run.py` |
| 探测器内核自检 | `python tests/kernel_suite.py` |
| 重新生成 hadronic_table.csv | `python scripts/build_hadronic_table.py` |

`run_scan.py` 的 `--masses` 必须是强子化库的质量节点（生产网格已经与库对齐）。给一个库外
质量时程序会**拒绝**而不是借用邻近节点。库节点就是生产网格本身；`work/checks/hadronisation_validation.json`
的 `library_node` 字段给出了验证时实际使用的节点。

## 6. 部分重跑与指纹

`outputs/data/signal_scan.json` 记录 8 个源码文件（`bc5/model/widths.py`、
`bc5/production/rates.py`、`bc5/analysis/scan.py`、`bc5/decay/*`、`bc5/response/*`）、
3 个母样本、5 个强子化库文件的 **sha256**。`run_scan.py --only-masses` 会检查：

* `model_tag`、σ(bb̄)、B 产额因子；
* 上述全部文件指纹。

任何一项变了都会拒绝合并，并提示全网格重跑（或 `--force-mix`，危险）。

逐值比较两套扫描（例如 v4 与 v3.6 的回归）：

```powershell
python scripts/compare_scans.py outputs/data/signal_scan.json `
    ../bc5_lhcb_downstream_v3/outputs/data/signal_scan.json
```

退出码 0 表示物理负载逐值相同（只允许 `model_tag`/指纹等元数据不同）。

## 7. 数据来源

| 输入 | 位置 | 说明 |
|---|---|---|
| `bbar_heavy.txt.gz` | `work/magnus/samples/` | 1,509,528 行 b 强子；Magnus 作业 `9d36fa12b8485053` |
| `minbias_light.txt.gz` | 同上 | K 母样本；作业 `bc5_minbias4` |
| `higgs.txt.gz` | 同上 | 200,000 个 Higgs；作业 `bc5_higgs4` |
| `hadron_lib/*.txt.gz` | `work/magnus/hadron_lib/` | 67 个扫描质量点；作业 `bc5_hadronlib3` |
| `hadronic_table.csv` | `work/model/` | 数字化 Winkler Fig. 3（ππ/KK） |
| 数字化参考曲线 | `work/model/reference/` | BC5.pdf、Winkler Fig. 2、Boiarska BR 图 |
| 探测器头文件快照 | `inputs/extension/modules/` | 与 `detector/downstream/modules` 逐字节相同 |

样本与强子化库随 git 提供（v4 与 v3 的内容相同）。若从零重建，按 `work/magnus/` 里的
`submit_*.py` 与 `logs/magnus_jobs.md` 记录的作业 ID 重新提交（需要 Docker + Pythia 8.316）。

## 8. 目录里哪些是"证据"，不是运行必需

| 目录 | 内容 | 可以删吗 |
|---|---|---|
| `work/magnus/` | 母样本、强子化库、作业 bundle 与摘要 | **不能删**（运行输入） |
| `work/model/` | ππ/KK 表 + 数字化参考曲线 | **不能删**（模型/验证输入） |
| `work/checks/` | 验证脚本的 JSON 结果 | 可删，重跑验证即重建 |
| `tests/` | 单元/回归测试 | 建议保留 |
| `logs/` | v3 验收/审查历史（见 `logs/README.md`） | 建议保留（说明"为什么这么做"） |
| `outputs/` | 最终图与数据 | **不能删**（`signal_scan.json` 是重画曲线的输入） |

## 9. 排错

| 现象 | 原因/处理 |
|---|---|
| `could not find downstream_kernel.dll` | 没编译或不在 `bc5/response/native/`；跑 `build.ps1` |
| `hadron library has no node ...` | 用了库外质量；换生产网格节点，或重建库 |
| `refusing to merge into signal_scan.json` | 源码/样本/模型标签变了；全网格重跑或 `--force-mix` |
| `ModuleNotFoundError: bc5` | 从 run 根目录运行，或用 `scripts/` 下的入口（会自己加路径） |
| 图里中文/符号乱码 | 图内文字都是 LaTeX/英文；检查 matplotlib 版本 |
