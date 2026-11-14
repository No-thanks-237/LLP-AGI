# v4 重构回归记录（2026-11-14）

目标：证明 `bc5_lhcb_downstream_v4` 相对 `bc5_lhcb_downstream_v3`（v3.6）**只改了结构**，
物理与数值逐位不变。

## 1. 模型/衰变/响应（函数级）

| 检查 | 方法 | 结果 |
|---|---|---|
| 宽度 `widths(m_S)` | 15 个质量 × 11 个道，v3 `scalar_widths` vs v4 `bc5.model.widths` | 0 个不同 |
| 衰变产物 `s_decay_products` | 5 个库节点 × 11 个道 × 2 个种子 = 102 组，含 Pythia 强子化库 | 0 组不同 |
| `acceptance_profile` | 6 组（质量/种子），同一 A(z) 网格 | 0 组不同 |
| `p_window` | 27 组（m_S, pz, z0），A 为线性 | 0 组不同 |
| `s_decay_products` 库外质量 | 应与 v3 一样拒绝（不借用邻近节点） | 行为一致 |

## 2. 扫描引擎（固定 nsamp）

v3 与 v4 各自跑 `--masses 1.0,2.0,5.082255 --nsamp 400 --ndecay 8`（同一随机种子、
同一批样本/强子化库），比较 `signal_per_fb_per_channel` 的每个字段：

    compared values: 7920   differences: 0   max relative difference: 0.0

参考宽度模型（`BC5_WIDTH_MODEL=reference_br`）同样跑 3 个质量：N_S 计数逐质量相同
（224 / 187 / 105），说明 `WIDTH_OVERRIDE` 路径也没有改。

## 3. 完整生产扫描（自适应）

用 v3 的 RUNNING 命令重跑（`--adaptive --target-rel-err 0.05 --adaptive-rounds 6
--nsamp-max 25000 --force`）。自适应轮次与 v3 完全一致：

| 轮次 | 重跑质量数 | nsamp 总量 |
|---|---|---|
| base | 67 | 80,400 |
| 1 | 67 | 159,700 |
| 2 | 24 | 59,500 |
| 3 | 6 | 15,200 |
| 之后 | 0（全部交点 ≤ 5%） | — |

参考宽度模型（`--suffix _refbr`，同一命令）也完整重跑：

| 轮次 | 重跑质量数 | nsamp 总量 |
|---|---|---|
| base | 67 | 80,400 |
| 1 | 67 | 162,600 |
| 2 | 25 | 59,400 |
| 3 | 2 | 5,200 |
| 之后 | 0 | — |

最终 `signal_scan.json`：

    masses: 67  compared values: 176880  differing: 0  max relative difference: 0.000e+00
    nsamp_by_mass: 与 v3 完全相同（min 1400 / median 2300 / max 4200）
    sigma_bb_pb, b_yield_corr: 相同
    model_tag: v4.0-2026-11-14（v3 是 v3.6-2026-10-25；仅标签/指纹不同）

参考模型扫描（`signal_scan_refbr.json`）同样是 67 个质量、240 个 θ²：

    compared values: 176880  differing: 0  max relative difference: 0.000e+00
    nsamp_by_mass: 与 v3 完全相同（min 1500 / median 2300 / max 4000）

## 4. 曲线与图（逐字节）

| 文件 | 结果 |
|---|---|
| `outputs/data/exclusion_boundary.csv` | 相同 |
| `outputs/data/exclusion_s_curves.csv` | 相同 |
| `outputs/data/exclusion_report.json` | 相同 |
| `outputs/figures/bc5_lhcb_downstream_exclusion.png` | 相同 |
| `outputs/figures/bc5_vs_reference.png` | 相同 |
| `work/checks/reference_comparison.json`（在 v3 是 `work/widths/`） | 相同 |
| `outputs/data/scalar_widths.csv` | 相同 |
| `outputs/figures/widths_validation.png` | 相同 |
| `work/checks/hadronisation_validation.json` | 相同 |
| `outputs/figures/production/*.png`、BR 图 PNG | 相同 |
| `exclusion_boundary_refbr.csv`、`exclusion_s_curves_refbr.csv`、`exclusion_report_refbr.json` | 相同 |
| `bc5_lhcb_downstream_exclusion_refbr.png`、`bc5_vs_reference_refbr.png` | 相同 |
| `exclusion_boundary_s3.csv`、`exclusion_s_curves_s3.csv`、`bc5_lhcb_downstream_exclusion_s3.png`、`bc5_vs_reference_s3.png` | 相同 |

PDF 只有内部时间戳元数据不同（内容相同）。

## 5. 唯一的内容差异

`outputs/data/run_validation.json`：v3 内嵌的 `hadronisation_library` 是旧副本，
v4 重新生成的是当前版本（与 v3 的独立文件
`work/analysis/hadronisation_validation.json` 逐字节相同）。这是 v3 的文件不同步，
不是物理变化。

## 6. 测试

    python -m pytest tests -q     # 22 passed

包括：相空间不变量、A=1 的解析寿命积分、BR 和/阈值、Boiarska 产生率、曲线回归
（从已提交 `signal_scan.json` 重算并与 `exclusion_boundary.csv` 逐字节比较）、
探测器内核 33 项自检。

## 7. 结论

v4.0 是 v3.6 的等价实现：代码按用途重排、导入改为包内、路径集中到 `bc5.paths`，
所有物理产物逐位/逐字节一致。以后修改 v4 时，可以用本文件的命令与 v3.6 对比。
