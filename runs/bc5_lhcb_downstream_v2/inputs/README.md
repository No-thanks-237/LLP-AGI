# LLP-AGI

研究 Higgs 混合单态 dark scalar（BC5）及 Type-I 2HDM 的 CP-even light scalar，在 LHCb Downstream 探测器上的参数空间灵敏度。

2026-10-16 已重新部署官方 ColliderAgent 与 Magnus：新服务容器为 `llp-magnus`，界面为 <http://localhost:3011/jobs>。当前入口、版本、验证范围及恢复方法见 [容器部署说明](docs/container_deployment.md)，环境锁为 `configs/container_deployment.lock.json`。旧部署配置与旧锁只作历史材料，不据此恢复旧环境。

## 文件管理

| 目录               | 用途                                          |
| ---------------- | ------------------------------------------- |
| `configs/`       | 研究参数与环境配置；旧配置须核对后才能使用                       |
| `deployment/`    | 部署源码与资产；重新部署前核对版本和依赖                        |
| `detector/`      | Downstream 探测器卡、传播与选择代码                     |
| `docs/`          | 今后确认的设计和使用说明；旧项目说明清理后不保留                    |
| `phenomenology/` | 按模型组织的物理输入、数据与处理程序；历史派生数据不直接用于新计算           |
| `references/`    | 论文 PDF 等原始参考资料                              |
| `runs/`          | 每次运行建立独立子目录，集中保存输入、代码快照、日志、作业 ID 和结果；不覆盖旧运行 |
| `temp/`          | 临时脚本和可重建文件，用完清理，不存唯一副本                      |

2026-10-16 已完成工作区清理：根目录只保留本 README 和上述目录，原 `scripts/` 已改为 `temp/`。论文 PDF、物理数据、配置和源码保留，旧项目说明已删除。此次新部署源码和适配器集中于 `deployment/sources/`、`deployment/container/`，部署测试与临时记录分别进入 `runs/`、`temp/`。

后续每次尝试使用独立的 `runs/<时间戳_任务名>/`，按阶段保存输入、代码快照、日志和结果；不同尝试不共用可写目录，不覆盖旧产物。断线后先查询或收集已有作业，不自动重复提交。模型数据按 `phenomenology/<模型名>/` 管理，并区分原始输入与派生数据；保留的历史派生数据和结果不作为已验收的新计算输入。

