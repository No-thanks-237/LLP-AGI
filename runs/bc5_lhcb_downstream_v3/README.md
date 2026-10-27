# BC5 暗标量 @ LHCb Downstream（v3.6 验收基线）

[`prompt.md`](../../prompt.md) 的实现：Higgs 混合单态实标量 `S`（PBC benchmark **BC5**，
`BR(h→SS)=0.01`）在 LHCb Downstream 的 `(m_S, θ²)` 90% CL 排除曲线，`L = 25 / 300 fb⁻¹`，
零背景、零观测。

本目录是**模型 + 程序 + 产物**的验收快照（2026-10-25，v3.6）；验收记录见
[`logs/v36_acceptance.md`](logs/v36_acceptance.md)，物理报告见 [`logs/run_report.md`](logs/run_report.md)。

## 结果

排除区在每个质量上都是 θ² 上的一个**带**：θ² 高于长寿命边时 S 在 `z_dec > 2.5 m` 之外衰变，
低于短寿命边时在 `z_dec < 1 m` 之内衰变。只有一个交点落在扫描范围（θ² ≤ 1e-3）内时，
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

* 交点结构：300 fb⁻¹ 有 34 个质量是双边、18 个开放（全部 m_S < 0.35 GeV）；25 fb⁻¹ 为 35 / 17。
* 两个运动学台阶是物理的：m_S ≈ 0.28 GeV（ππ 打开）、m_S ≈ 5.13 GeV = m_B − m_π（B→X_sS 关闭，
  只剩 h→SS）；f₀(980) 在 0.98 GeV 附近给出可见的窄结构。
* 早期 s = 3（95% CL）的结果原样保留（文件名后缀 `_s3`），与 s = 2.3 来自同一次信号扫描。

## 模型处方（v3.6）

| 环节 | 处方 |
|---|---|
| S 宽度 | **strict Winkler**（arXiv:1809.01876）：色散 ππ/KK、C m³β₂π 多介子、固定质量微扰夸克道、LO gg、α_em(0) 的 γγ |
| B→X_sS | **Boiarska 逐道**（arXiv:1904.10447），B±/B⁰ 十道 + b→d，B_s 无混合道 |
| K→πS、四点道 | Boiarska；α 由 `BR(h→SS)=0.01` 反解 |
| h→SS | Pythia 8 13 TeV 混合产生 × 逐过程重加权 |
| 探测器 | 1 < z_dec < 2.5 m、2 < η_S < 5、对向双径迹 E > 5 GeV / pT > 0.5 GeV、B_y = 1 T、SciFi 3.24×2.415 m、ε_rec = 0.4 |
| 统计 | s₉₀ = 2.3（= −ln 0.1），事例级 OR，母粒子分层 + HT 权重 |

## 目录

| 路径 | 内容 |
|---|---|
| `inputs/` | 带 sha256 的输入快照（任务书、模型数据、探测器快照、复用代码） |
| `work/widths/` | 宽度/BR/ctau（Winkler + Boiarska）、数字化表与审计 |
| `work/production/` | 产生道（B→X_sS、K→πS、四点道、h→SS）与验证 |
| `work/response/` | Downstream 传播内核封装、响应表 |
| `work/analysis/` | 扫描、曲线、拟合与出图（`run_sensitivity.py` 是入口） |
| `work/review_check/` | 自检脚本（相空间、寿命积分、转动不变性、参考对比） |
| `outputs/data/` | `exclusion_boundary.csv`、`exclusion_s_curves.csv`、`signal_scan.json` 等 |
| `outputs/figures/` | 排除曲线（s=2.3 与 s=3）、分支比、与参考对比、产生谱 |

## 复现

完整步骤见 [`RUNNING.md`](RUNNING.md)。最短自检（从已提交的聚合结果重画，约 1 分钟）：

```powershell
cd work\analysis
python run_sensitivity.py --curves-only
```

完整重跑（约 50 分钟，含自适应采样到 MC 误差 < 5%）：

```powershell
python run_sensitivity.py --adaptive --target-rel-err 0.05 --adaptive-rounds 6 --nsamp-max 25000 --force
```

母样本与强子化库不在仓库内（GB 级），按 [`logs/magnus_jobs.md`](logs/magnus_jobs.md) 里记录的
Magnus 作业 ID 重新生成；`signal_scan.json` 记录了全部输入的 sha256 指纹。

## 验证

* 与参考实现逐点对比：分支比十道 0.3–6%；B⁰/B± = 0.9275；合计 3.356。
* 程序自检（[`work/review_check/`](work/review_check/)）：相空间 ⟨E_i/M⟩、⟨s23⟩、阈值成功率、
  转动不变性、寿命积分 vs 解析、中性 K 共享顶点、强子化库质量节点覆盖。
* MC：边界相对误差中位 ≈ 4.8%，最大 < 5%；`rel_s_mc_err_at_long_edge` 见 CSV。

## 已知限制

* 固定夸克质量、无 b bbar、2 GeV 强子/部分子衔接、α_em = 1/137 等作者/参考程序的近似原样保留。
* 无真实背景、无探测器/理论系统误差：这是"零背景假设下的预期灵敏度"。
* γ 不作为可见末态（S→γγ 只进总宽度）。
