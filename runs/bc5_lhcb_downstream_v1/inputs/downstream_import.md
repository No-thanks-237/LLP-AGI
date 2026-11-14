# Downstream 导入流程

适用于当前 `prompt.md` 与新容器部署。由 ColliderAgent 经 Magnus 完成接入；本文仅规定流程，不代表已经编译、注册或验证通过。

## 1. 当前状态与接入方式

`detector/downstream/` 保存了专用传播内核、Delphes 模块、card 和测试源码。新部署的 8 个官方蓝图**尚未接入这些扩展**，旧测试锁已清除，旧运行结论不能作为新环境的验收。

当前半解析计算优先直接调用 C++ 传播内核计算末态响应；Delphes 模块用于 HepMC/ROOT 接口及交叉验证。两者复用同一内核，避免另写一套磁场传播。

| 文件 | 作用 |
| --- | --- |
| `modules/DownstreamTransport.h` | $B_y$ 磁场传播与 SciFi 几何接受度，不依赖 ROOT |
| `modules/BirthKinematics.h` | 出生点真值能量和横动量选择 |
| `modules/DownstreamPropagator.{h,cc}` | 将上述内核接入 Delphes，生成 `Track` |
| `cards/delphes_card_downstream.tcl` | 正式轨迹选择与几何配置 |
| `cards/delphes_card_downstream_geometry_only.tcl` | 仅几何测试，关闭出生动量阈值，不能用于正式接受度 |
| `build.sh`、`verify_adapter.sh`、`tests/` | 隔离编译、接口测试及打包 |

**仅复制 card 不够。** 官方 Delphes 未注册 `DownstreamPropagator`；使用该 card 前，必须编译新模块并生成 ROOT 字典，或在半解析驱动中直接调用内核。

## 2. 固定配置与职责

| 项目 | 本轮定义 |
| --- | --- |
| $S$ 衰变区域 | $1<z_{\rm dec}<2.5$ m，$2<\eta_S<5$ |
| 磁场 | $B_y=1$ T，$3<z<7$ m；其余区域为零 |
| SciFi 包络 | $7.7\le z\le9.4$ m，$|x|\le3.24$ m、$|y|\le2.415$ m，束流孔外侧 $r\ge0.09$ m；边界接受 |
| 带电后代选择 | 出生点真值 $E>5$ GeV、$p_T>0.5$ GeV |
| 轨迹接受模式 | `full_envelope`，不要求穿过 UT |
| 事件选择 | 同一 $S$ 至少一对相反电荷后代通过；同一碰撞事件只计一次 |
| 额外重建效率 | 事件级 $0.4$，只乘一次 |

传播内核只判断单条轨迹。衰变区域、$S$ 祖先匹配、异号配对、碰撞级 OR、寿命积分及重建效率由响应／分析层处理，不能靠 card 自动完成。

SciFi 外边缘和束流孔边界沿用 `DownstreamTransport.h` 的包含边界规则；完整 SciFi 段均须满足包络与避孔条件。衰变区和出生动量 cut 仍为严格不等号。

输入 HepMC 使用 `GEV MM`。`DownstreamPropagator` 将顶点 mm 转为内核的 m，再将输出位置转回 mm；直接调用内核时四动量用 GeV、位置和 $ct$ 用 m。`Track.PT` 是传播后的值，不能用它重做出生点 cut；`D0/DZ` 等未定义字段不是重建顶点。

## 3. 隔离准备与编译

### 第一步：建立运行快照

每次新尝试建立独立目录，例如：

```text
runs/<run-id>/
  inputs/extension/        # detector/downstream 的源码、card、测试快照
  inputs/source.sha256     # 输入文件哈希
  work/build/              # 本次独占编译目录
  logs/                    # 命令、版本、job ID、编译与测试日志
  outputs/                 # 编译包、校验和、验收结果
```

复制前检查 `build.sh` 与 `verify_adapter.sh` 的依赖。当前 `verify_adapter.sh` 最后会复制 `/extension/README.md`，但清理后原文件不存在：**须在本 run 的 `inputs/extension/README.md` 写入包用途、版本与接口说明，或复制本文作为该文件，再封存输入哈希。**不必为此修改原始扩展或静默恢复旧说明。

### 第二步：在工具作业内编译

固定使用当前 `collider` 镜像：

```text
git.pku.edu.cn/het-agi/collider@sha256:1bf6d8db2cb9348f965c3c6650f1c6af8234c1d8be0b49040ddd55e79a33c1a5
```

接入驱动须将本 run 的源码快照提供为 `/extension`，将独占编译工作目录提供为 `/build`。这些路径属于临时作业容器，不是服务容器中的共享编译区。容器内执行：

```bash
source /opt/root/bin/thisroot.sh
export BUILD_JOBS=4
bash /extension/build.sh
sha256sum /build/downstream-delphes.tar.gz
```

`build.sh` 会复制镜像内 `/opt/MG5_aMC_v3_7_0/Delphes` 到 `/build/Delphes`，加入模块和 card，更新该副本的 `ModulesLinkDef.h`，然后编译共享库及 HepMC2/3 reader。`make clean` 仅作用于该隔离副本，不修改镜像内的原 Delphes。

完成后将 `/build/downstream-delphes.tar.gz`、SHA-256、编译／测试日志收集到本 run。原镜像、官方源码及其他 run 不改写；编译包只在本轮内复用。

## 4. 导入与调用

### 路线 A：半解析响应，优先采用

薄驱动包含：

```cpp
#include "DownstreamTransport.h"
#include "BirthKinematics.h"
```

对每个候选 $S$ 衰变顶点生成后代四动量，将每条带电后代写为 `downstream::State`。显式设置与正式 card 相同的 `Config`、`BirthCuts`，先调用 `passesBirthKinematics(...)`，再调用 `transport(state, config).accepted()`。

驱动按祖先关系配对，得到给定顶点的响应 $A_a(\ell)$；随后按 `prompt.md` 积分：

$$
p_S=\sum_a\mathrm{BR}_a\int_0^\infty
\frac{d\ell}{\lambda_S}e^{-\ell/\lambda_S}A_a(\ell),
\qquad \lambda_S=\frac{|\mathbf p_S|}{m_S}c\tau_S.
$$

最后在共同产生状态与互斥衰变历史下求事件 OR，并乘一次 $0.4$。不要在此之前随机抽一次 $S$ 飞行距离、筛掉未衰变事件，否则会重复施加衰变概率。完整 HepMC/ROOT 文件不是该路线的逐点生产中间量；保留紧凑事件及响应缓存，并记录其配置与输入哈希。

**当前仍需编写并接入响应／积分薄驱动，已有头文件本身不等于完整半解析入口。**

### 路线 B：Delphes 接口与交叉验证

将本 run 编译包解压到新的作业目录 `package/`，验证 SHA-256，设置：

```bash
source /opt/root/bin/thisroot.sh
export LD_LIBRARY_PATH="$PWD/package/lib:${LD_LIBRARY_PATH:-}"
export ROOT_INCLUDE_PATH="$PWD/package/source:/opt/MG5_aMC_v3_7_0/Delphes:/opt/MG5_aMC_v3_7_0/Delphes/external:${ROOT_INCLUDE_PATH:-}"
./package/bin/DelphesHepMC2 ./package/cards/delphes_card_downstream.tcl output.root input.hepmc
```

HepMC3 使用包内 `DelphesHepMC3`。必须加载包内匹配的 `libDelphes.so`、字典和 card，不能用镜像原版 reader 加新 card。`ROOT_INCLUDE_PATH` 保留基础头文件与新增模块头文件，以支持移位后的字典加载。

输入事件必须已经包含正确的位移顶点及后代；模块不生成 $S$ 衰变。ROOT 中 `Track.Particle` 关联原始真值粒子，分析层据完整祖先链识别其所属 $S$。

## 5. Magnus 接入与验收

当前官方 `madgraph-launch` 没有专门传入 Downstream 编译包或加载该模块的接口，不能仅设置 `detector=Delphes` 就认为扩展生效。后续须落实受控的编译、合成测试及响应执行入口：输入使用本 run 的文件包，固定镜像，输出返回可核验结果和日志；客户端断线后按 job ID 收集，不重复提交。

`detector/downstream/magnus_boundary_blueprint.py` 和 `tests/run_boundary_job.py` 可作为合成测试参考，但它们未在新服务中注册。后者的默认流程主要测试几何，不覆盖正式动量选择，正式 card 必须另行验收。

验收至少包含：

1. `build.sh` 中的独立传播／出生动量测试与 `verify_adapter.sh` 的模块测试通过。
2. 用合成 HepMC2/3 检查 GeV/mm 转换、真值关联、正负电荷弯曲方向、零场极限及 ROOT 字典加载。
3. 检查 $z=3,7$ m 的磁场分段、SciFi 边缘／束流孔、反向或转弯轨迹，以及严格 $E>5$、$p_T>0.5$ 的阈值行为。
4. 分析层另测衰变区域边界、同一 $S$ 异号配对、双 $S$ 只计一次及一次 $0.4$；半解析层另测闭式寿命积分、长短寿命极限和 MC 收敛。

作业 `Success`、返回结果中的 `success`、实际输出检查须分别通过。保存新 run 的验收报告；软件测试通过不等于两种标量模型的物理灵敏度已经验证。

当前剩余工作是：补齐输入快照的打包说明文件，落实 Magnus 编译／响应入口，并在新环境完成上述验收。本文没有执行这些步骤。
