# bc5/analysis — 扫描、曲线与出图

统计层：把模型/产生/衰变/响应拼成 s(m_S, θ²)，再从扫描提取排除边界并出图。

## 一条命令跑完

```powershell
cd runs\bc5_lhcb_downstream_v4
python scripts/run_sensitivity.py                 # 完整 scan + s = 2.3 曲线 + 参考对比
python scripts/run_sensitivity.py --quick          # 8 质量 × 200 样本，冒烟测试（写 *_quick）
python scripts/run_sensitivity.py --curves-only    # 只从 signal_scan.json 重画（秒级）
```

`--s 3`（95% CL）与 `--suffix` 的行为与 v3 相同；`--adaptive` 的说明见下。

## 两个阶段

**1. scan — [`scan.py`](scan.py)。** 对每个质量抽样母粒子产生（Higgs、B、K）、用强子化库
抽 S 衰变、把带电末态送进 LHCb Downstream 内核，汇总成事例级效率 ε(m_S, θ²)。
输出 `outputs/data/signal_scan.json`：每个质量、每个道在 300 fb⁻¹ 下的期望事例数
（theta² 扫描网格）与 MC 误差；同时记录输入 sha256、nsamp、模型标签。

**2. curves — [`curves.py`](curves.py)。** 把 s 曲线与阈值 s_cut 相交，写
`outputs/data/exclusion_boundary.csv`、`exclusion_s_curves.csv`、`exclusion_report.json`
并画 `outputs/figures/bc5_lhcb_downstream_exclusion.{png,pdf}`。只有长寿命边在扫描窗口内时，
`band_open_above=True`，不伪造闭合曲线。

**3. 驱动 — [`sensitivity.py`](sensitivity.py)。** 预检输入、跑 scan、可选自适应加密、
跑 curves、跑参考对比，最后打印交点处的 MC 误差。

## 扫描参数

| 参数 | 含义 |
|---|---|
| `--nmass 52` | 对数质量点数（0.021–62.49 GeV） |
| `--masses 0.5,1,2` | 显式质量列表（会覆盖 nmass；必须是强子化库节点） |
| `--nsamp 1200` | 每个质量、每个道抽多少母粒子 |
| `--ndecay 8` | 每个 S 抽样多少次衰变来估 A(z) |
| `--adaptive` | 先跑基扫描，再只对交点 MC 误差 > `--target-rel-err` 的质量加密 |
| `--adaptive-rounds 6` | 最多几轮 |
| `--nsamp-max 25000` | 单质量样本上限 |
| `--suffix _refbr` | 所有输出加后缀（不会覆盖生产产物） |

* 阈值质量节点（2m_μ、2m_π、K 关道、2m_K、1 GeV、2 GeV、2m_τ、2m_D、2m_b）会自动并入
  对数网格（`merge_mass_grid`），所以 52 个对数点会变成 67 个质量点。
* 扫描时间与质量点数近似线性；低质量点最贵。
* 曲线阶段至少需要两个质量点。

## 输入

```
work/magnus/samples/{bbar_heavy,minbias_light,higgs}.txt(.gz)  母样本
work/magnus/samples/normalisation.json                          n_gen / σ
work/magnus/hadron_lib/hadron_lib_*.txt.gz                     Pythia 强子化库
work/model/hadronic_table.csv                                  ππ/KK 宽度
bc5/response/native/downstream_kernel.dll                      探测器内核
```

归一化只从 `normalisation.json` 读（B = 800000、K = 200000、H = 200000）；`BC5_NGEN_*`
环境变量可覆盖。

## 出图模块（[`plots/`](plots/)）

| 模块 | 产物 | 入口 |
|---|---|---|
| `plots/branching.py` | 模型 BR vs m_S（含全区间/含 bb 版本） | `scripts/plot_branching_ratios.py` |
| `plots/senscalc.py` | SensCalc 表的 BR 对比 | `scripts/plot_branching_ratios_senscalc.py` |
| `plots/br_validation.py` | 我们 vs 参考表的逐道验证图 | `scripts/plot_br_validation.py` |
| `plots/production.py` | 各母粒子 (θ, p) 产生谱 | `scripts/plot_production.py` |
| `plots/reference.py` | 我们的轮廓 vs 数字化 BC5 参考 | `scripts/plot_comparison_reference.py` |

每个模块暴露 `main()`，可以直接 `python -m bc5.analysis.plots.<name>`。

## 与 v3 的差异

只改了路径与导入：函数签名、算法、随机流、输出格式都没变。v3 的
`work/analysis/README.md` 里对 `run_scan.py` / `make_curves.py` 的描述仍适用，只是文件名
换成了 `scan.py` / `curves.py`，入口在 `scripts/`。
