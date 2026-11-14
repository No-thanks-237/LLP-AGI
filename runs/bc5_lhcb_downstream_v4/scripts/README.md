# scripts — 薄 CLI 入口

每个脚本只做两件事：把 run 根目录加入 `sys.path`，调用 `bc5.*` 的 `main()`。物理/统计代码
都在包里。所有命令都在 run 根目录下运行。

| 脚本 | 调用 | 主要产物 | 耗时 |
|---|---|---|---|
| `run_sensitivity.py` | `bc5.analysis.sensitivity` | 扫描 + 主图 + 参考对比 | 秒～小时 |
| `run_scan.py` | `bc5.analysis.scan` | `outputs/data/signal_scan*.json` | ~20 min+ |
| `make_curves.py` | `bc5.analysis.curves` | 边界 CSV + 主图 | 秒 |
| `plot_comparison_reference.py` | `plots.reference` | `bc5_vs_reference.*` | 秒 |
| `plot_branching_ratios.py` | `plots.branching` | `bc5_branching_ratios*.csv/png` | 秒 |
| `plot_branching_ratios_senscalc.py` | `plots.senscalc` | `bc5_branching_ratios_senscalc.*` | 秒 |
| `plot_br_validation.py` | `plots.br_validation` | `br_validation.*` | 秒 |
| `plot_production.py` | `plots.production` | `outputs/figures/production/` | ~30 s（读大样本） |
| `widths_cli.py` | `bc5.model.widths` | 终端/JSON | 毫秒 |
| `validate_widths.py` | `bc5.model.validation` | `work/checks/reference_comparison.json` + 图 | ~20 s |
| `audit_widths.py` | `bc5.model.audit` | `outputs/data/widths_audit.json` | ~20 s |
| `validate_production.py` | `bc5.production.rates` | `work/checks/production_validation.json` | 秒 |
| `validate_hadronisation.py` | `bc5.decay` | `work/checks/hadronisation_validation.json` | ~30 s |
| `validate_run.py` | `bc5.analysis.scan` + response | `outputs/data/run_validation.json` | 秒 |
| `build_hadronic_table.py` | `bc5.model.widths` | `work/model/hadronic_table.csv` | 秒 |
| `check_inputs.py` | — | 校验 `inputs/source.sha256` | 毫秒 |
| `compare_scans.py` | — | 逐值比较两个 `signal_scan.json`（回归检查） | 秒 |

`validate_run.py` 依赖前三个验证 JSON，先跑它们。完整命令见 `../RUNNING.md`。
