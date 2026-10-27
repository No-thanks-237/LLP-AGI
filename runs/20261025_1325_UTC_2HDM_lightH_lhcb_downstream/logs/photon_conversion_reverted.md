# γ→e⁺e⁻ 转换方案：实现、对比与回退（2026-10-26 00:20）

## 决定

用户查证 SensCalc 的约定：**γ 不作为末态**。因此最终结果维持"γγ 只计入总宽度、在带电轨迹
选择下响应为零"（prompt_2HDM 原始处方）。本文件记录期间做过的转换方案与回退，避免以后再试。

## 曾实现的方案（已回退）

* 在 `work/response/twohdm_lib.py` 中增加 γ 在探测器材料中的转换模型：
  沿光子直线路径的材料按指数分布抽样转换点（有效 X/X0 = 0.10 over 6.7 m 参考路径，
  μ=(7/9)X0⁻¹），e⁺e⁻ 共享光子能量（Bethe–Heitler 权重 1−4/3x(1−x)、共线近似 m_ee=0），
  从真实转换点经 B_y 场传播并过 SciFi 包络；任一光子转换且 e⁺e⁻ 对通过选择即该衰变可见。
  另有"直接光子（ECAL 型）"和"不可见"两种可切换模式。
* 代码、缓存配置和扫描脚本均已回退，工作区无残留（grep PHOTON_MODE / CONVERSION_ 为空）。
* 参考数据：`references/tools/FORESEE/Models/2HDM_H`（作者衰变表/解析式），
  作者代码 zip 在 `references/model/Light_scalar_decay_Combined_Code.zip`。

## 对比数字（仅记录，未进入最终结果）

mH = 1 GeV，长寿边 tanβ（25 fb⁻¹ / 300 fb⁻¹）：
| 方案 | 25 fb⁻¹ | 300 fb⁻¹ |
|---|---|---|
| γ 不可见（最终采用，prompt/SensCalc） | 357.0 | 487.5 |
| γ→e⁺e⁻ 转换（X/X0=0.10） | 5083 | 9461 |
| γ 直接可见（ECAL 型，测试用） | 11018 | — |

γγ 可见会使 $BR_\gamma\gamma\to1$ 的通道从"寿命控制"变成"信号主体"，长寿命边上移 1–1.5 个
数量级；这正是不能用错可见性假设的原因。

## 回退后的验证

清空扫描产物、用回退后的代码完整重跑 foresee/paper_mb/exclusive 后，与回退前快照逐点比较：
2640 个点 s(mH,tanb) 最大相对差 **0**，264 个边界条目 0 处不一致（logs/restructure_reproducibility.md
的同一方法）。最终交付数字与 figure 与回退前完全一致。
