# BC5 LHCb-downstream 验收记录（v3.6-2026-10-25）

这一版是【模型 + 程序 + 产物】的验收快照。下面每一条都能用仓库里的脚本复现。

## 1. 模型（与参考工具对齐的最终选择）

| 环节 | 处方 | 依据 |
|---|---|---|
| S 的衰变宽度 | **strict Winkler**（arXiv:1809.01876）：色散 ππ/KK（数字化 Fig. 3）、C m³β_2π 多介子、固定质量微扰夸克道、LO gg、α_em(0) 的 γγ | 用户指定"完全按 Winkler" |
| B→X_sS | **Boiarska 逐道**（arXiv:1904.10447 tab:BR，ξ_sb = 3.6e-4），B±/B⁰ 十道 + b→d；**B_s 无混合道**（原文明确排除；SensCalc 同上） | prompt.md §3.1；SensCalc `SelectedBrChoice="Exclusive"` |
| K→πS | Boiarska（K± 1.61e-3、K⁰_L 6.96e-3 @m→0，表值 1.7e-3 / 7e-3） | 同上 |
| 四点道 | Boiarska tab:b-production-quartic（B_s→SS、B→X'SS、K→SS…），α 由 BR(h→SS)=0.01 反解 | 同上 |
| h→SS | Pythia8 13 TeV 混合产生 × SIGMA_H 逐过程重加权；b_h = 0.01 | — |
| 探测器 | 1 < z_dec < 2.5 m、2 < η_S < 5、对向双径迹 E>5 GeV / pT>0.5 GeV、By=1 T（3–7 m）、SciFi 3.24×2.415 m、ε_rec = 0.4 | 与 SensCalc 的接受度卡逐项一致 |
| 统计 | s = 2.3（90% CL），L = 25 / 300 fb⁻¹ | prompt.md |

## 2. 输入（全部可追溯到 Magnus 作业）

| 输入 | 作业 | 状态 |
|---|---|---|
| `samples/bbar_heavy.txt.gz`（1,509,528 行，1.8855 B/event） | `9d36fa12b8485053` | 已核 |
| `samples/minbias_light.txt.gz`（1,122,896 行，无 η 预筛） | `bc5_minbias4` | v3.5 |
| `samples/higgs.txt.gz`（1,411,171 行，含 status 列，200,000 条 \|status\|=62） | `bc5_higgs4` | v3.5 |
| `hadron_lib/*`（67 个扫描质量点） | `bc5_hadronlib3` | v3.4 |
| `widths/hadronic_table.csv`（数字化 Winkler Fig. 3） | — | 已核 |

`signal_scan.json` 里记录了全部源文件/样本/库的 sha256 指纹；部分重跑时指纹不符会拒绝合并。

## 3. 程序自检（`work/review_check/` 下的脚本）

| 检验 | 结果 |
|---|---|
| 无质量四体 ⟨E_i/M⟩ | 0.250 / 0.252 / 0.251 / 0.247 |
| 有质量三体 ⟨s23⟩（3π, 1 GeV） | 0.3538 ± 0.0012（精确 0.3528） |
| 四体阈值附近生成成功率 | 100 %（m = 4m_π(1+1e-6) … ） |
| 转动不变 ⟨cos²θ_z⟩ | 0.331–0.335（应 1/3；修前 0.52–0.99） |
| 寿命积分 vs 解析（A=1，含 z0 在窗内） | 全部 1.000000000（50 位对照） |
| 中性 K 共享顶点 | 4000 次 KK 抽样，2976/2976 对同源后代共享同一出生点 |
| 逐道产生率 vs Boiarska tab:BR | 十道 0.3–6 %；B⁰/B± = 0.9275；合计 3.356 |
| 强子化库质量节点覆盖 | 67 个扫描点全部精确命中（0 个需要借用） |
| 宽度自洽 | BR 和 = 1.000000；2 GeV 接口比 0.981；阈值全部正确 |

## 4. 结果

| 量 | 值 |
|---|---|
| 25 fb⁻¹ 长寿边 MC 误差 | 中位 **4.8 %**，最大 5.0 %，**0/67 超 5 %** |
| 300 fb⁻¹ 长寿边 MC 误差 | 中位 4.7 %，最大 5.0 %，**0/67 超 5 %** |
| 边界位置误差（长寿边） | 中位 3.4 %，最大 5.0 %，**0/134 超 5 %** |
| 生产曲线 / 参考（中位，\|ln\| 中位） | **1.061 / 0.402** |
| 验证曲线（参考 BR/宽度）/ 参考 | **0.955 / 0.217** |
| 输入指纹 | 重扫后与包内源文件一致（`--only-masses` 合并检查通过） |
| 短寿边 | 47/47 点 > 5 %（已知重尾，需 5e4–1e5 母事例，超出上限） |

## 5. 已知限制（全部有量化）

1. **宽度模型**：固定夸克质量（Winkler）vs 参考的跑动质量，40–60 GeV 处 Γ_S 差 1.32–1.51 倍 → 高 mass 段 0.68–0.76 的系统偏移；
2. **参考产物自身不闭合**：本地 SensCalc 曲线只有原 PDF 的 0.55–0.62；其预计算标量谱与自己母谱差 17–24 %；
3. B_s 只有四点道（Boiarska 排除 B_s 混合）；B_s 混合的物理贡献被牺牲，量级 ≤10 % 的 B 道；
4. K 道在长寿边占 ~1e-4；D、B_c 道未实现（D/B ~ 5e-9）；
5. 强子化库无次级顶点、无显式 QCD shower（LO 条件末态近似）；重味注入 p=E, m_q=0；
6. 短寿边 MC 未收敛；62 GeV 端点比较受参考封口与 m_h=125 vs 125.25 影响；
7. σ(bb̄) = 560 μb 的实验外推 ±15 %；Higgs 过程混合已按截面重加权，tH（0.17 %）未生成。

## 6. 产物清单

* 图：`outputs/figures/`（说明见该目录的 README；排除曲线 + 参考对比 + 宽度/分支比 + 生产相空间）
* 数据：`outputs/data/`（signal_scan.json、signal_scan_refbr.json、exclusion_boundary{,_refbr,_s3}.csv、exclusion_s_curves*.csv、bc5_branching_ratios*.csv、scalar_widths.csv、widths_audit.json）
* 日志：`logs/`（v3_plan、v3_model_summary、production_audit_v3、widths_audit_v3、mc_convergence_v3、physics_review_20261024 + 核对、physics_review_v34_20261025 + 核对、high_mass_diagnosis_v35、v36_boiarska_revert、本文件）
* 审查/验证脚本：`work/review_check/`（每个结论都有对应脚本与 JSON）
## 7. 可移植性（2026-10-25 实测）

目录**自带全部运行输入，可独立运行**（见 `RUNNING.md`）：

* 参考表已拷进 `inputs/phenomenology/darkscalar/`，脚本不再引用跑目录之外的路径；
  `digitize.py` 的原始 PDF 路径可用 `BC5_REFS` 覆盖；
* 实测：整目录拷到另一路径后，`run_sensitivity.py --curves-only`、6 张生产相空间图、全部
  验证脚本（validate_widths / validate_production / validate_hadronisation / validate_run /
  audit_widths）全部跑通；重扫后 `run_scan.py --only-masses` 的指纹检查也通过（"provenance ok"）；
* 依赖：Windows x64（预编译的 `downstream_kernel.dll`）、Python ≥3.11、numpy + matplotlib
  （pandas 仅 `plot_production.py`）；体积 237 MB（删掉两个审查产物目录后约 185 MB）。

## 8. 最后一次重扫（清理死代码之后）

为让扫描记录的源文件指纹与交付包一致，清理死代码后用最终代码把两套扫描各重跑一遍：
数字与清理前**完全一致**（MC 4.8%/5.0%、0/67；生产 1.061/0.402；验证 0.955/0.217），
即清理确实是纯删除、无物理影响。
