# BC5 暗标量（第一版，历史）

> **状态：历史记录。** 这是 BC5 分析的第一次完整运行（2026-10-16），后续已被
> [`../bc5_lhcb_downstream_v2/`](../bc5_lhcb_downstream_v2/) 与
> [`../bc5_lhcb_downstream_v3/`](../bc5_lhcb_downstream_v3/)（v3.6 验收基线）取代。
> 本目录只作为过程记录，其中的曲线与数字不要引用。

第一次运行完成了从 Pythia 母样本、宽度/分支比、探测器响应到排除曲线的端到端链路，
并在随后几轮审查中暴露并修复了若干问题（归一化、真实强子化与 τ 衰变、四点产生道、
宽度模型的处方选择等）。这些问题与修复的演化见：

* [`../bc5_lhcb_downstream_v2/CHANGES_vs_v1.md`](../bc5_lhcb_downstream_v2/CHANGES_vs_v1.md)
* [`../bc5_lhcb_downstream_v3/logs/v3_fixes.md`](../bc5_lhcb_downstream_v3/logs/v3_fixes.md)
* [`../bc5_lhcb_downstream_v3/logs/v36_acceptance.md`](../bc5_lhcb_downstream_v3/logs/v36_acceptance.md)

目录结构：`inputs/`、`work/`（widths/production/response/analysis）、`logs/`、`outputs/`。
数据（母样本/强子化库）不在仓库内。最终结果请看 v3。
