# logs — 报告与历史

## v4 的日志

| 文件 | 内容 |
|---|---|
| `v4_regression.md` | v4 完整重跑与 v3.6 的逐点对比（重构验收） |

## 从 v3 复制的历史（原样保留）

v3.6 的验收记录、物理报告、审查报告与运行日志整体复制到这里，作为"当时为什么这么做"的
证据。**这些文件里的 `work/analysis/...`、`work/widths/...` 路径指 v3 的布局**，在 v4 中
对应的位置见 `../CHANGELOG.md` 第 1 节。

重点文档：

| 文件 | 内容 |
|---|---|
| `v36_acceptance.md` | v3.6 验收记录（模型选择、自检、MC 误差、限制） |
| `run_report.md` | BC5 物理报告（最完整） |
| `v3_task.md` | v3 的任务书（大致目标，历史） |
| `v3_plan.md` / `v3_fixes.md` | v3 的设计与修正清单 |
| `physics_review_20261024.md` / `physics_review_v34_20261025.md` | 两轮物理审查 |
| `review_response_v3.md` / `review_response_v35.md` | 对审查的逐项回复 |
| `mc_convergence_v3.md` | MC 收敛研究 |
| `magnus_jobs.md` | 样本生成的 Magnus 作业 ID |
| `production_*.log`、`scan_*.log` | 生产运行日志 |

v3 的一次性审查脚本与 93 MB 中间产物**没有**复制到 v4，仍在
`../../bc5_lhcb_downstream_v3/work/review_check/` 与 `.../physics_review_20261025/`。
