# bc5 — 包结构

BC5 暗标量分析的唯一代码源。按物理阶段分成 5 个子包；**依赖只能向下**，没有循环导入。

```
bc5/
├─ paths.py            运行目录与数据位置（所有路径的唯一入口）
├─ model/              衰变：宽度、BR、ctau
│  ├─ widths.py        物理实现（导入时装载 work/model/hadronic_table.csv）
│  ├─ reference.py     validation-only 的 SensCalc 反推模型
│  ├─ validation.py    宽度验证
│  └─ audit.py         逐通道审计
├─ production/         产生：B→X_sS、K→πS、四点道、h→SS
│  └─ rates.py
├─ decay/              衰变运动学与道组成
│  ├─ kinematics.py    boost、精确 n 体相空间、喷注碎裂
│  ├─ hadronlib.py     Pythia 强子化库
│  └─ channels.py      道列表、带电末态、库路由、channel_branching
├─ response/           探测器响应
│  ├─ kernel.py        downstream_kernel.dll 的 ctypes 封装
│  ├─ lifetime.py      A(z)、精确衰变窗积分、选择
│  └─ native/          C++ 源、预编译 DLL、build.ps1
└─ analysis/           统计与出图
   ├─ scan.py          扫描引擎 + CLI
   ├─ curves.py        排除边界 + 主图
   ├─ sensitivity.py   一键驱动
   └─ plots/           5 个作图模块
```

## 依赖方向

```
paths  ◄── 所有子包
model  ◄── production
model  ◄── decay
decay  ◄── response
model/production/decay/response  ◄── analysis
```

任何子包都**不**导入 `analysis`；`response` 是唯一导入编译内核的地方；`analysis` 是唯一
读写 `outputs/` 的地方（`model.validation` 与 `model.audit` 也会写图/JSON，属于显式的
验证入口）。

## 公共 API

| 子包 | 常用名字 |
|---|---|
| `bc5.paths` | `ROOT`、`SAMPLES`、`HADRON_LIB`、`MODEL_DATA`、`REFERENCE_DATA`、`OUT_DATA`、`OUT_FIG`、`LOGS` |
| `bc5.model` | `widths.widths(ms)`、`widths.gamma_total`、`widths.ctau_m`、`widths.load_hadronic_table`、`reference.widths` |
| `bc5.production` | `br_B_to_S`、`br_K_to_S`、`quartic_channels`、`quartic_kaon_channels`、`sample_dipion_mass`、`alpha_from_bh`、`SIGMA_H`、`SIGMA_H_TOTAL` |
| `bc5.decay` | `two_body`、`nbody`、`boost`、`s_decay_products`、`channel_branching`、`load_hadron_lib` |
| `bc5.response` | `acceptance_profile`、`p_window`、`lifetime_probability`、`kernel`、`RECO_EFF`、`ZGRID` |
| `bc5.analysis` | `scan.main`、`curves.main`、`sensitivity.main`、`plots.*.main` |

示例：

```python
import random

from bc5.model import widths
from bc5.production import rates
from bc5.decay import nbody
from bc5.response import acceptance_profile, p_window

print(widths.gamma_total(1.0), rates.br_B_to_S(1e-6)['K'])
print(nbody(1.0, [0.13957039] * 3, random.Random(1)))
```

## 约定

* **单位**：GeV、米、c = 1；样本顶点是 mm，进入响应层前统一乘 1e-3。
* **宽度**一律在 theta² = 1；`ctau_m(ms, theta2)` 再乘 1/theta²。
* **随机数**：所有抽样函数显式接收 `rng`（`random.Random`）；模块内不使用全局随机源，
  只有 `scan.py` 的母粒子子采样用固定种子的 numpy Generator，以保持可复现。
* **路径**：不写死 `E:/...`；一律通过 `bc5.paths`。样本/强子化库可用
  `BC5_SAMPLES`、`BC5_HADRON_LIB` 覆盖。
* **表**：`bc5.model.widths` 导入时自动装载 `work/model/hadronic_table.csv`；要 ChPT 近似时
  显式 `widths.HADRONIC_TABLE = None`（CLI 的 `--hadronic chpt`）。
* **入口**：`scripts/*.py` 只做 `sys.path` + `main()`，不含物理。

## 扩展方式

* **加一个衰变道**：`decay/channels.py` 的 `CHANNELS`、`s_decay_products`、`model/widths.py`
  的 `widths()` 三处同步；再在 `tests/test_widths.py` 加 BR 和/阈值断言。
* **加一个产生道**：`production/rates.py` 加矩阵元/BR，`scan.sample_meson_decay_S` 或
  `sample_quartic_decay_S` 里接上；`tests/test_production.py` 加对照。
* **加一张图**：`analysis/plots/<name>.py` 暴露 `main()`，在 `scripts/plot_<name>.py` 加
  5 行包装，并更新 `outputs/figures/README.md`。
* **加一个验证**：脚本放 `scripts/`，结果写 `work/checks/*.json` 或 `outputs/data/*`，
  在 `RUNNING.md` 第 4 节和 `logs/` 里登记。
