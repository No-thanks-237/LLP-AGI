# inputs — 本 run 的冻结输入（只读）

这里的每个文件都是计算开始前冻结的输入。运行中的代码**只读**这个目录，不写它；所有派生产物
（扫描、响应表、图）都在 `outputs/`、`work/` 下。`source.sha256` 记录全部输入的 sha256。

## 清单

| 路径                                                       | 内容                                                                        | 用途                                                                                                           |
| -------------------------------------------------------- | ------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| `prompt.md`                                              | 通用任务书（研究框架、探测器、验收标准）                                                      | 物理设定与验收依据                                                                                                    |
| `downstream_import.md`                                   | LHCb Downstream 几何/选择导入说明                                                 | 探测器实现说明                                                                                                      |
| `container_deployment.lock.json`                         | 2026-10-16 的 ColliderAgent/Magnus 环境锁                                     | 环境可追溯（其中 `validation_run` 仍写当时的 `runs/20261016_051119_UTC_deployment-smoke`，即现在的 `runs/deployment_smoke_v1`） |
| `extension/modules/`                                     | `DownstreamTransport.h`、`BirthKinematics.h`、`DownstreamPropagator.{h,cc}` | 探测器内核头文件快照（与仓库 `detector/downstream/modules/` 逐字节相同）                                                         |
| `extension/cards/`                                       | Delphes 探测器卡（完整 / 仅几何）                                                    | 探测器响应卡                                                                                                       |
| `extension/tests/`                                       | 边界/传输/出生运动学的 C++ 与 Python 测试                                              | 内核验证（历史）                                                                                                     |
| `extension/build.sh`、`verify_adapter.sh`                 | 扩展构建/校验脚本                                                                 | 探测器接入                                                                                                        |
| `phenomenology/darkscalar/prepared/branching_ratios.csv` | SensCalc `BrRatioScalar` 转换表（8115 个质量）                                    | reference-BR 验证模型与 BR 对比图                                                                                    |
| `phenomenology/darkscalar/prepared/lifetime.csv`         | SensCalc 寿命表                                                              | 宽度/寿命交叉验证                                                                                                    |
| `phenomenology/darkscalar/GammaS_reference.csv`          | SensCalc 总宽度                                                              | 宽度审计                                                                                                         |
| `phenomenology/darkscalar/BrRatioScalar_channels.csv`    | SensCalc 分道表                                                              | 通道对照                                                                                                         |
| `physics/winkler1809_scalar.tex`                         | Winkler, arXiv:1809.01876 源码                                              | 宽度/色散 ππ/KK 的第一性依据                                                                                           |
| `physics/boiarska1904_main.tex`                          | Boiarska et al., arXiv:1904.10447 源码                                      | 产生道 / 形状因子                                                                                                   |
| `physics/boiarska1904_formfactors.json`                  | 论文附录 C 的形状因子参数                                                            | 产生率实现                                                                                                        |
| `physics/boiarska1904_production_reference.md`           | 产生道逐条对照笔记                                                                 | 实现核对                                                                                                         |
| `source.sha256`                                          | 以上文件的 sha256                                                              | 完整性校验                                                                                                        |

> 数字化的参考曲线（Winkler Fig. 2、Boiarska BR 图、BC5.pdf）是从论文 PDF 提取的**派生**
> 数据，放在 `work/model/reference/`，不在 inputs 里；生成脚本与报告也在那里。

## 校验

```powershell
cd runs\bc5_lhcb_downstream_v4
python scripts/check_inputs.py
# 期望：checked 27 file(s), 0 problem(s)
```

在 Git Bash 里 `sha256sum -c inputs/source.sha256` 也可以（清单用正斜杠路径）。

## 约定

* 这个目录**只读**。需要新输入时新建一个 run 版本，不要把运行期产物写进来。
* 修改任何输入都必须同时更新 `source.sha256`，并在 `CHANGELOG.md` 记录原因。
* `inputs/README.md` 与 `source.sha256` 本身不参与哈希清单。
