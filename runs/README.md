# 运行索引

每个 run 一个目录，互不覆盖。命名规则：**`<程序>_v<版本>`**（2026-11-14 统一）。
目录内统一结构：`inputs/`（带 sha256 的输入快照）、`work/`（数据/作业/缓存）、
`logs/`（报告、审查、验证记录）、`outputs/`（`data/` + `figures/`）、
`README.md`（结果与复现）、`RUNNING.md`（怎么跑）。代码位置因版本而异：

* **BC5 v4**：`bc5/` 包 + `scripts/` 入口（按用途分包：model/production/decay/response/analysis）。
* **BC5 v1–v3、2HDM v1**：`work/<阶段>/` 脚本（widths/decay/production/response/analysis/magnus）。

## 目录

| 目录 | 版本 | 日期 | 状态 | 内容 |
|---|---|---|---|---|
| [`bc5_lhcb_downstream_v4/`](bc5_lhcb_downstream_v4/) | **4.0.0** | 2026-11-14 | **当前代码基线** | BC5 的结构重构版：`bc5/` 包 + `scripts/`；物理与 v3.6 逐位相同 |
| [`bc5_lhcb_downstream_v3/`](bc5_lhcb_downstream_v3/) | 3.6 | 2026-10-20 → 10-25 | **v3.6 验收基线**（物理参考） | BC5 的完整历史证据：269 个审查脚本、93 MB 审查产物、验收记录 |
| [`bc5_lhcb_downstream_v2/`](bc5_lhcb_downstream_v2/) | 2 | 2026-10-18 | 历史 | BC5 第二版，差异见 `CHANGES_vs_v1.md` |
| [`bc5_lhcb_downstream_v1/`](bc5_lhcb_downstream_v1/) | 1 | 2026-10-16 | 历史 | BC5 第一次完整运行（原名 `20261016_144605_UTC_bc5_lhcb_downstream`） |
| [`2hdm_light_h_lhcb_downstream_v1/`](2hdm_light_h_lhcb_downstream_v1/) | 1 | 2026-10-25 → 10-27 | **当前 2HDM 结果** | Type-I 2HDM 轻 CP-even H 的 (m_H, tanβ) 90% CL 排除曲线 |
| [`deployment_smoke_v1/`](deployment_smoke_v1/) | 1 | 2026-10-16 | 环境自检 | ColliderAgent/Magnus 容器部署 smoke test（原名 `20261016_051119_UTC_deployment-smoke`） |

## 阅读顺序

1. 根 [`README.md`](../README.md)：项目总览、两次分析的结论表、复现入口。
2. 当前 BC5 代码基线：[`bc5_lhcb_downstream_v4/README.md`](bc5_lhcb_downstream_v4/README.md)，
   结构见同目录 `bc5/README.md`，逐文件迁移见 `CHANGELOG.md`。
3. 当前 2HDM 结果：[`2hdm_light_h_lhcb_downstream_v1/README.md`](2hdm_light_h_lhcb_downstream_v1/README.md)，
   审查细节见同目录 `logs/physics_code_review.md`。
4. BC5 物理与审查历史：[`bc5_lhcb_downstream_v3/README.md`](bc5_lhcb_downstream_v3/README.md)，
   验收记录 `logs/v36_acceptance.md`，物理报告 `logs/run_report.md`。

## 历史路径说明

2026-11-14 统一命名时改了三个目录：

| 原名 | 现名 |
|---|---|
| `20261016_051119_UTC_deployment-smoke` | `deployment_smoke_v1` |
| `20261016_144605_UTC_bc5_lhcb_downstream` | `bc5_lhcb_downstream_v1` |
| `20261025_1325_UTC_2HDM_lightH_lhcb_downstream` | `2hdm_light_h_lhcb_downstream_v1` |

`bc5_lhcb_downstream_v2/v3` 的名字本来就符合规则，未改。**历史日志/审查脚本里出现的旧目录名
按上表对应**；这些文件作为当时的证据原样保留，没有批量改写。

## 数据

* Pythia 母样本（`work/magnus/samples/*.txt.gz`）与强子化库（`work/magnus/hadron_lib/`）
  **随仓库提供**；旧版 run（v1/v2）里字节相同的重复副本不重复入库（见 `.gitignore`）。
* BC5 v4 与 v3 的样本/强子化库内容相同（git 按内容去重，不会重复存储）。
* 2HDM 的 500k 事例 pi/K minbias 样本（`minbias_pik.txt.gz`，525 MB）超过 GitHub 单文件
  限制，以 **Release 资产** `samples-2026-10` 发布，取用方法见根 README 的"数据与第三方材料"。
* 响应表缓存与逐质量扫描 JSON 不入库，可自动重建；BC5 的缓存现在是 `work/cache/`
  （v3 及更早是 `work/analysis/cache/`），2HDM 仍是 `work/analysis/cache/`。
* `logs/samples.md` / `logs/magnus_jobs.md` 记录了生成样本的 Magnus 作业 ID。
