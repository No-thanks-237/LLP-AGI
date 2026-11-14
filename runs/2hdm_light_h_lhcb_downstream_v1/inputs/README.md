# inputs — 本 run 的输入快照（带 sha256，只读）

这里的每个文件都是计算开始前冻结的输入；`source.sha256` 记录全部输入的 sha256，可用
`cd inputs && sha256sum -c source.sha256` 校验。运行中的代码**不写**这个目录。

| 路径 | 内容 |
|---|---|
| `prompt.md` | 通用任务书（研究框架、探测器、验收标准） |
| `prompt_2HDM.md` | 本次 2HDM 运行的任务书（模型、产生道、扫描范围、验收指标） |
| `detector/` | LHCb Downstream 探测器卡与传播内核快照（`cards/`、`modules/`、`tests/`、Magnus 边界作业蓝图） |
| `decay_source/` | 论文作者 `Light_scalar_decay` C++ 源码（`src/`、`include/`、`makefile`、`manual.pdf`），本 run 仅加批处理驱动补丁 |
| `bc5_reuse/` | 从 BC5 v3.6 复用的代码与数据：`bc5_lib.py`、`bc5_production.py`、`downstream_ctypes.py` + 内核 DLL、强子化表、样本归一化、Magnus 作业蓝图、Pythia 母样本生成代码（`parents.cc`、`shadlib.cc`） |
| `references/` | 模型论文 LaTeX 源（`FASER_light.tex`，arXiv:2212.06186）与插图 |
| `container_deployment.lock.json` | ColliderAgent/Magnus 容器环境锁 |
| `source.sha256` | 以上全部输入的 sha256 清单 |

产生/衰变/响应的具体处方与验证见上级目录的 `README.md`、`DESIGN.md` 与 `logs/`。
