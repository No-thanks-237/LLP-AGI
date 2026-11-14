# outputs/figures — 每张图与它的生成脚本

所有命令都在 run 根目录 `runs/bc5_lhcb_downstream_v4` 下执行，Python 3.11+（numpy +
matplotlib；`plot_production.py` 另需 pandas）。当前产物的模型标签：**v4.0-2026-11-14**，
物理与 v3.6 相同。**没有一张图是手工做的**。

## 主结果：排除曲线

| 图 | 脚本 | 需要 |
|---|---|---|
| `bc5_lhcb_downstream_exclusion.{png,pdf}` | `scripts/make_curves.py` | `outputs/data/signal_scan.json` |
| `bc5_lhcb_downstream_exclusion_refbr.*` | `make_curves.py --suffix _refbr` | `signal_scan_refbr.json` |
| `bc5_lhcb_downstream_exclusion_s3.*` | `make_curves.py --s 3 --suffix _s3` | 同一扫描，95% CL |

`signal_scan.json` 由 `scripts/run_sensitivity.py --adaptive ...` 生成（约 50–90 分钟，
最贵的一步）。默认阈值 s = 2.3（90% CL）；`--s 3` 给 95% CL。

## 与参考曲线对比（数字化 BC5.pdf）

| 图 | 脚本 |
|---|---|
| `bc5_vs_reference.{png,pdf}` | `scripts/plot_comparison_reference.py` |
| `bc5_vs_reference_refbr.*` | 同上，`--suffix _refbr` |
| `bc5_vs_reference_s3.*` | 同上，`--s 3 --suffix _s3` |

## 宽度与分支比

| 图 | 脚本 | 内容 |
|---|---|---|
| `bc5_branching_ratios.*` | `scripts/plot_branching_ratios.py` | 模型 BR vs m_S |
| `bc5_branching_ratios_fullrange.*` | 同上 | 0.02–62 GeV 全区间 |
| `bc5_branching_ratios_with_bb.*` | 同上 | 含 b bbar 阈值 |
| `bc5_branching_ratios_senscalc.*` | `scripts/plot_branching_ratios_senscalc.py` | SensCalc BrRatioScalar 表 |
| `bc5_branching_ratios_senscalc_all.*` | 同上 | 所有道，对数坐标 |
| `br_validation.*` | `scripts/plot_br_validation.py` | 我们 vs 参考表，逐道绝对宽度 |
| `widths_validation.png` | `scripts/validate_widths.py` | 宽度 vs 数字化 Winkler/参考表 |

对应的数据：`outputs/data/bc5_branching_ratios.csv`、`bc5_branching_ratios_senscalc.csv`。

## 产生相空间（母粒子）

| 图 | 脚本 |
|---|---|
| `production/production_{Bplus,B0,Bs,Kplus,K0L,Higgs}.{png,pdf}` | `scripts/plot_production.py` |
| `production/production_all_particles.png` | 同上（2×3 拼图） |

(`theta, p`) 平面上的二维直方图，用 `work/magnus/samples/normalisation.json` 的截面归一。
Higgs 只取衰变前的副本（|status| = 62）。

## 一键重画（不重跑 MC）

```powershell
python scripts/run_sensitivity.py --curves-only
python scripts/plot_branching_ratios.py
python scripts/plot_branching_ratios_senscalc.py
python scripts/plot_br_validation.py
python scripts/plot_production.py
python scripts/validate_widths.py
python scripts/audit_widths.py
```

## 说明

* PDF 与 PNG 内容相同；PDF 内部有生成时间戳，所以不同次生成的 PDF 不会逐字节相同，PNG 会。
* v3 的 `_quick` / `_test` 变体与 Omnes-vs-Winkler 图已在 v3.6 清理中删除，v4 不保留。
