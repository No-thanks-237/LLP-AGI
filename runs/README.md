# 运行索引

每次运行一个目录，互不覆盖。目录内统一结构：`inputs/`（带 sha256 的输入快照）、
`work/`（各阶段程序：生产/衰变/响应/分析）、`logs/`（报告、审查、验证记录）、
`outputs/`（`data/` + `figures/`）、`README.md`（结果与复现）、`RUNNING.md`（怎么跑）。

| 目录 | 状态 | 内容 |
|---|---|---|
| [`20261025_1325_UTC_2HDM_lightH_lhcb_downstream/`](20261025_1325_UTC_2HDM_lightH_lhcb_downstream/) | **当前主结果**（2026-10-27） | Type-I 2HDM 轻 CP-even H 的 (m_H, tanβ) 90% CL 排除曲线；解析产生 + 作者衰变程序 + BC5 探测器内核 |
| [`bc5_lhcb_downstream_v3/`](bc5_lhcb_downstream_v3/) | **BC5 验收基线**（v3.6，2026-10-25） | Higgs 混合单态暗标量的 (m_S, θ²) 排除曲线 + 完整自检记录 |
| [`bc5_lhcb_downstream_v2/`](bc5_lhcb_downstream_v2/) | 历史（已被 v3 取代） | BC5 的第二版，保留作对照；差异见 `CHANGES_vs_v1.md` |
| [`20261016_144605_UTC_bc5_lhcb_downstream/`](20261016_144605_UTC_bc5_lhcb_downstream/) | 历史（已被 v3 取代） | BC5 第一版完整运行，含早期修正记录 |
| `20261016_051119_UTC_deployment-smoke/` | 环境自检 | ColliderAgent/Magnus 容器部署的 smoke test 记录 |

## 阅读顺序建议

1. 根 [`README.md`](../README.md)：项目总览、两次分析的结论表、复现入口。
2. 当前主结果：[`20261025_1325_UTC_2HDM_lightH_lhcb_downstream/README.md`](20261025_1325_UTC_2HDM_lightH_lhcb_downstream/README.md)，
   审查细节见同目录 `logs/physics_code_review.md`。
3. BC5 基线：[`bc5_lhcb_downstream_v3/README.md`](bc5_lhcb_downstream_v3/README.md)，
   验收记录见同目录 `logs/v36_acceptance.md`。

## 数据

母样本、强子化库、响应表缓存、逐质量扫描 JSON 不在仓库内（体积 GB 级或可自动重建），
各 run 的 `logs/samples.md` / `logs/magnus_jobs.md` 记录了生成它们的 Magnus 作业 ID。
