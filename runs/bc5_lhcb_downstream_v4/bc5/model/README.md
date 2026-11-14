# bc5/model — S 的衰变宽度、分支比与寿命

BC5 Higgs-portal 实标量 S 的完整衰变描述：给定 m_S 返回每个道的宽度、BR 和 c*tau。
本子包不依赖探测器或 MC，只依赖 numpy。

```powershell
cd runs\bc5_lhcb_downstream_v4
python scripts/widths_cli.py --mass 1.0                 # 单点，人类可读
python scripts/widths_cli.py --mass 1.0 --theta2 1e-9   # 给定混合角
python scripts/widths_cli.py --mass 2.5 --json          # 机器可读
python scripts/widths_cli.py --scan 0.021 62.49 52      # 对数质量扫描
python scripts/widths_cli.py --check                    # 对参考曲线做验证
```

## 模型

S 以角度 theta 与 Higgs 混合，另有四次耦合 alpha：

```
L = -theta (m_f/v) S fbar f  -  theta (2 m_W^2/v) S W+W-  -  theta (m_Z^2/v) S ZZ
    -  alpha S^2 |H|^2
```

因此对每个有 SM Higgs 表达式的道，`Gamma(S→X) = theta^2 * Gamma(h→X)|_{m_h→m_S}`；
本模块给出的宽度都在 theta^2 = 1。alpha 诱导的道不属于这里（在 `bc5.production`）。

## 两个区间

* **m_S < 2 GeV**：强子宽度用 Winkler 的色散结果——ππ、KK 来自数字化的 Fig. 2
  （`work/model/hadronic_table.csv`），多介子道用 `C m_S^3 beta_2pi`，C = 5.1e-9 GeV^-2；
  在 2 GeV 与微扰区衔接。
* **m_S ≥ 2 GeV**：微扰夸克/胶子，跑动夸克质量在 μ = m_S，阈值由最轻的含该夸克的强子
  （K、D、B）给出。Boiarska 的高阶 QCD 修正**没有**加（见 `scalar_widths.py` 的
  `QUARK_NLO` 与 `logs/run_report.md`）。

轻子和光子在两个区间都存在。

## 模块

| 文件 | 作用 |
|---|---|
| [`widths.py`](widths.py) | 物理实现：α_s、跑动质量、圈函数、每个道、汇总；导入时自动装载 `work/model/hadronic_table.csv` |
| [`reference.py`](reference.py) | **仅验证**：从 SensCalc 的 BR 表反推宽度模型，用于"复现参考"检查 |
| [`validation.py`](validation.py) | 宽度验证（图 + JSON），入口 `scripts/validate_widths.py` |
| [`audit.py`](audit.py) | 逐通道审计（Winkler/SensCalc/Boiarska + 内部自洽），入口 `scripts/audit_widths.py` |

数据与生成脚本：

| 位置 | 内容 |
|---|---|
| `work/model/hadronic_table.csv` | ππ/KK 宽度（theta^2 = 1），数字化 Winkler Fig. 3 并做 μ 子道自校准 |
| `work/model/hadronic_table_report.json` | 上述表的校准记录 |
| `work/model/hadronic_table_600.csv` | 600 点版本（审计用） |
| `work/model/reference/` | 数字化参考曲线与数字化报告 |
| `scripts/build_hadronic_table.py` | 从 `work/model/reference/decayrates_curves.csv` 重建表 |

## 各道的来源

| 道 | 计算方式 | 独立？ |
|---|---|---|
| e+e-、μ+μ-、τ+τ- | 闭式：G_F m_S m_l^2 β^3 / (4√2 π) | 是 |
| γγ | SM h→γγ 圈函数，m_h → m_S | 是 |
| gg（m_S ≥ 2 GeV） | α_s^2 m^3 |ΣF_1/2|^2 / (128 π^3 v^2) | 是 |
| q qbar（m_S ≥ 2 GeV） | N_c m_S mbar_q^2(m_S) β^3 / (8 π v^2)，LO spectator | 是 |
| 多介子 | C m_S^3 β_2π，C 来自 Winkler | 公式；常数取文献 |
| **ππ、KK（m_S < 2 GeV）** | **数字化 Winkler Fig. 2** | **否——外部发表输入** |

## 验证

```powershell
python scripts/validate_widths.py   # -> work/checks/reference_comparison.json + outputs/figures/widths_validation.png
python scripts/audit_widths.py      # -> outputs/data/widths_audit.json
```

要点（详见上面的 JSON）：

| 检查 | 结果 |
|---|---|
| S→μμ 对解析曲线 | 比值 1.000，散布 1e-4 |
| S→γγ 在 m = 125 GeV 的 SM 极限 | 复现 9.1e-6 GeV |
| S→gg 对 Winkler Fig. 2 | 中位 6%，4 GeV 以上好于 3% |
| S→q qbar 对 Winkler Fig. 2 | 2 GeV 的 s 道 3.5%；5 GeV 的 c 道低约 30%（与 SensCalc 在 10% 内） |
| α_s(m_Z) | 0.1179（3-loop，3 个味阈值） |
| 重顶极限 | α_s^2 m^3 / (72 π^3 v^2) 复现 |
| 2 GeV 强子宽度衔接 | 0.981（Winkler 1.0002；SensCalc 1.0000） |
| SensCalc 交叉参考 | 轻子 < 0.1%；其余见 JSON |
