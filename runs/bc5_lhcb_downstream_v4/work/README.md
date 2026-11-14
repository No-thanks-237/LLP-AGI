# work — 运行数据（不含代码）

代码全部在 `bc5/`；这个目录只放数据与作业产物。

| 路径 | 内容 | 是否入库 |
|---|---|---|
| `magnus/samples/` | Pythia 母样本（B、K、Higgs）与 `normalisation.json` | 是（gz） |
| `magnus/hadron_lib/` | 67 个扫描质量点的 Pythia 强子化库 | 是（gz） |
| `magnus/bundles*/`、`magnus/jobs/` | Magnus 作业 bundle、脚本、状态、摘要 | 是（小文件） |
| `model/` | `hadronic_table.csv`（ππ/KK）与数字化参考曲线 | 是 |
| `checks/` | 验证脚本的 JSON 结果（reference/production/hadronisation） | 是（小） |

## 约定

* 所有样本/库文件在 `bc5.analysis.scan` 的指纹里（sha256），部分重跑会校验。
* 换样本时用 `BC5_SAMPLES` / `BC5_HADRON_LIB` 指向别的目录，不要覆盖这里的生产输入。
* 缓存目录可以随时删除；下一次扫描会重建。
* 这里**不放代码**。脚本放 `scripts/`，库代码放 `bc5/`。

## 样本来源（Magnus 作业）

| 文件 | 作业 | 说明 |
|---|---|---|
| `samples/bbar_heavy.txt.gz` | `9d36fa12b8485053` | 1,509,528 行 b 强子 |
| `samples/minbias_light.txt.gz` | `bc5_minbias4` | K 母样本 |
| `samples/higgs.txt.gz` | `bc5_higgs4` | 200,000 Higgs（含 status 列） |
| `hadron_lib/*` | `bc5_hadronlib3` | 67 个扫描质量点 |

作业脚本在 `magnus/`（`submit_*.py`、`get_*.py`、`install_parents4.py`、`llp_job_blueprint.py`）；
更早的提交记录见 `logs/magnus_jobs.md`。
