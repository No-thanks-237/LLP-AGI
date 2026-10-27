# work/analysis — 扫描、响应、统计、作图（本 run 的"分析层"）

对应 BC5 v3.6 的 work/analysis 位置。衰变和产生是**独立的**代码目录：
  ../decay/      作者 Light_scalar_decay 编译、批处理驱动、python 封装、验证（T1）
  ../production/ 解析 2HDM 产生模块（xi_H、B/K/pi、Higgs kappa、h->HH）与验证（T2）
  ../response/   LHCb-downstream 探测器内核 + 逐通道响应库 twohdm_lib.py（T4）
  ../magnus/     Pythia 样本、H 衰变强子化库、Magnus 作业记录
  ../review_check/ 独立解析审计脚本（T5）

## 文件

| 文件 | 作用 |
|---|---|
| `models.py` | 2HDM 与 BC5 两套模型适配（产生权重、宽度、sigma、逐过程重加权） |
| `response_table.py` | 逐质量响应表：母样本读取/事件子采样、H 候选者生成、逐通道 A_ch(z)、npz 缓存 |
| `event_engine.py` | 概率代数（p_i -> 历史 OR -> 母粒子 OR -> 碰撞 OR）、eps/s、MC 误差、交点二分；CLI 逐质量运行 |
| `launch_scan.py` | 并行逐质量扫描驱动（默认 8 workers） |
| `aggregate_scan.py` | 合并逐质量 JSON -> `outputs/data/signal_scan.json` |
| `final_analysis.py` | 提取长寿/短寿交点、MC 误差、EFT 有效性；写边界 CSV、s 网格、宽度/BR 表 |
| `final_figure.py` | 主排除图（只画曲线，BC5 风格）+ 交点/MC 诊断图 |
| `signal_map.py` | log10(s/s90) 二维显著度图 |
| `bc5_regression.py` | BC5 v3.6 端到端回归（同一引擎 + BC5 模块/样本） |
| `test_engine.py` | 14 项引擎单元测试 |
| `collect_validation.py` | 汇总响应/引擎验证 JSON |
| `bc5ref/` | BC5 v3.6 的只读参考副本（回归用，sha256 记录在 validation JSON） |
| `grids/` | 质量网格与逐通道强子化库节点表、tanb=40 点批量衰变输入 |
| `cache/` | 逐质量响应表 npz（版本化，含 DLL/代码 sha256）；可删除后重建 |
| `scan_2hdm*/` | 逐质量扫描 JSON 与日志（foresee / paper_mb / exclusive） |

## 复现

```bash
python work/analysis/launch_scan.py --workers 8 --force --nsamp '{"B":17293,"inel":2500,"h":4620}' --b-mode foresee
python work/analysis/aggregate_scan.py --indir work/analysis/scan_2hdm --out outputs/data/signal_scan.json --b-mode foresee
python work/analysis/final_analysis.py     # 边界表 + exclusion_crossings.csv + 系统误差图
python work/analysis/final_figure.py       # 主图与交点/MC 诊断图
python work/analysis/signal_map.py         # 显著度图
python work/analysis/make_manifest.py      # sha256 清单
```

变体把 `--b-mode` 换成 `paper_mb` / `exclusive` 并指定 `--outdir work/analysis/scan_2hdm_<mode>`，
聚合时对应改 `--indir/--out/--b-mode`。

`--nsamp` 的 B/h 值决定实际分层样本量（`models._scaled_alloc` 按比例缩放
`B_STRAT_ALLOC`/`H_STRAT_ALLOC`）；`cache/` 与 `scan_2hdm*/` 是自动重建的中间产物（不随仓库分发）。
