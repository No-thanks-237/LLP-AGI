# BC5 暗标量（第二版，历史）

> **状态：历史记录。** 当前 BC5 结果以 [`../bc5_lhcb_downstream_v3/`](../bc5_lhcb_downstream_v3/)
> 的 v3.6 验收基线为准。本目录保留用于对照 v1 → v2 → v3 的演化，其中的数字不作为最终结果。

| | |
|---|---|
| 前身 | [`../20261016_144605_UTC_bc5_lhcb_downstream/`](../20261016_144605_UTC_bc5_lhcb_downstream/)（v1） |
| 后继 | [`../bc5_lhcb_downstream_v3/`](../bc5_lhcb_downstream_v3/)（v3.6 验收基线） |

本版相对于 v1 的修改清单见 [`CHANGES_vs_v1.md`](CHANGES_vs_v1.md)。
v3 相对本版的最终修正（strict Winkler 宽度、Boiarska 逐道产生率、归一化修复、四体/强子化
补齐等）见 [`../bc5_lhcb_downstream_v3/logs/v3_fixes.md`](../bc5_lhcb_downstream_v3/logs/v3_fixes.md)
与 [`../bc5_lhcb_downstream_v3/logs/v36_acceptance.md`](../bc5_lhcb_downstream_v3/logs/v36_acceptance.md)。

目录结构与其他 run 相同：`inputs/`（输入快照）、`work/`（widths/production/response/analysis）、
`logs/`（本版报告与审查）、`outputs/`（`data/` + `figures/`）。数据（母样本/强子化库/缓存）
不在此目录中，按 `logs/` 里记录的 Magnus 作业 ID 重建。

完整重跑方式与 v3 相同（`work/analysis/run_sensitivity.py`），可参考
[`../bc5_lhcb_downstream_v3/RUNNING.md`](../bc5_lhcb_downstream_v3/RUNNING.md)。
