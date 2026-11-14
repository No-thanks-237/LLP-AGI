# work/model — 模型数据表与数字化参考曲线

代码在 `bc5/model/`，这里只放它读取的数据。

| 文件 | 内容 |
|---|---|
| `hadronic_table.csv` | ππ、KK 宽度（theta² = 1，GeV），3000 行，质量 0.2792–2.0 GeV |
| `hadronic_table_600.csv` | 同一张表的 600 点版本（审计/交叉检查用） |
| `hadronic_table_report.json` | 数字化与 μ 子道自校准记录 |
| `reference/decayrates_curves.csv` | Winkler Fig. 2 的数字化分宽度曲线 |
| `reference/scalardecaybr_curves.csv` | Boiarska BR 图的数字化曲线 |
| `reference/tauscalar_curves.csv` | 寿命图的数字化曲线 |
| `reference/bc5_reference_curves.csv` | BC5.pdf 中参考排除轮廓的数字化 |
| `reference/bc5_digitisation_report.json` | BC5.pdf 的分段/行号记录（画图时按段切开） |
| `reference/overlay_*.png`、`bc5_digitisation_check.png` | 数字化校验叠加图 |
| `reference/digitization_report.md`、`OMNES_RECIPE.md` | 数字化方法与 Omnes 配方 |
| `reference/thesis/` | Hoferichter 等相位移动论文源（历史参考） |
| `reference/digitize.py` | 数字化脚本的**残缺快照**（仓库首次导入时就不完整）；以已提交的 CSV 为准 |

## hadronic_table.csv 的来源

ππ/KK 在 f₀(980) 附近需要 Winkler 的双通道 Muskhelishvili-Omnes 形状因子；相位移动表不在
仓库里，所以直接用 Winkler 自己 Fig. 2 的**发表曲线**（数字化），并用同图的 S→μμ 曲线
（有精确解析式）做自校准。生成脚本：`scripts/build_hadronic_table.py`。

格式：`mass_GeV,gamma_pipi_GeV,gamma_KK_GeV`，单位 GeV，theta² = 1。

## 约定

* 这些是**派生**数据（从论文 PDF 提取），不是原始输入；原始 PDF 在仓库 `references/model/`。
* 换数字化结果时，必须更新 `hadronic_table_report.json` 并重跑
  `scripts/validate_widths.py` 与 `python -m pytest tests -q`。
