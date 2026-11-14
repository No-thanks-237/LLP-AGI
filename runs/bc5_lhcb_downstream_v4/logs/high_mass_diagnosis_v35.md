# 为什么高质量段和参考曲线有差别 —— 全面筛查（2026-10-25）

对象：$m_S \gtrsim 5$ GeV 的长寿命边，该段信号 100% 来自 $h \to SS$。
我们的曲线在那里是参考曲线的 0.45–0.65（25 fb⁻¹）。本文回答三件事：**SensCalc 用什么产生处方、
我们的 h 产生是否有问题、探测器是否有问题**，并把残差逐项分解。

## 1. SensCalc 默认用"独占道"，不是包容公式

本地 SensCalc 的 `codes/EventCalc/Production` 里（`work/review_check/senscalc/codes_EventCalc_ProductionProbabilities.nb.txt`）：

    (*Exclusive/inclusive production of dark scalar in decays B->Xs+S*)
    SelectedBrChoice="Exclusive";
    PmotherToLLP[mLLP_,θ2_,BrhToSS_,"B-mixing"] = θ2*(fbtob0[Facility]+fbtobplus[Facility])
        * intProd[If[SelectedBrChoice=="Exclusive","BtoXsStotal.dat","BtoXsSinclusive.dat"],
                  mSM["Bcharged"], mSM["PiCharged"]];

即 **默认 `BtoXsStotal.dat` = Boiarska 的逐道（独占）产生率**，包容公式 `BtoXsSinclusive.dat`
只是可切换的备选。四点道同样是 Boiarska 的表（`BplustoXsSStotal.dat`、`BstoSS.dat`）。

推断：v3.4 把 B→X_sS 换成 Winkler 包容公式（当时基于"SensCalc 好像也用包容"的假设）**方向反了**；
这也解释了为什么换之前 v3.3/v3.4 的低质量段中位比值是 1.03，换之后变成 0.84。

**已按此更正（v3.6-2026-10-25）**：B→X_sS 回到 Boiarska 逐道 tab:BR（ξ_sb = 3.6e-4，十个 b→s 道 + b→d），
并删掉 B_s 混合道——Boiarska 明确只讨论 D⁰、D⁺、B⁻、B⁰（main.tex L682 明确排除 B_s），
SensCalc 的 "B-mixing" 也只有 `fbtob0+fbtobplus`。B_s 只经四点道 B_s→SS 贡献。

## 2. 我们的 h 产生没有问题（与 FORESEE 对比）

FORESEE（github.com/KlingFelix/FORESEE）自带 Pythia8 的 14 TeV 强子谱，含 Higgs（PDG 25）
`files/hadrons/14TeV/Pythia8/Pythia8_14TeV_25.txt`。用完全相同的判据（$h\to SS$ 各向同性二体、
$2<\eta_S<5$、$m_S/p_z$ 权重，不含探测器核）比较长寿命系数

$$K=\left\langle \mathbf 1_{2<\eta_S<5}\,\frac{m_S}{p_z}\right\rangle$$

| $m_S$ (GeV) | 我们 | FORESEE Pythia8 | 比值 |
|---:|---:|---:|---:|
| 5.08 | 0.002503 | 0.002545 | **0.983 ± 0.009** |
| 13.02 | 0.006623 | 0.006769 | **0.978 ± 0.008** |
| 33.37 | 0.022678 | 0.023044 | **0.984 ± 0.008** |
| 53.42 | 0.063773 | 0.065449 | **0.974 ± 0.006** |

只看 ggF 子样本时比值 0.997–1.004。**我们的 Higgs 产生谱是标准 Pythia8，与 FORESEE 一致到 2–3%。**

但参考(SensCalc)的表格化谱给出的 K 比我们大 20–25%（见第 3 节），即**参考的谱比 Pythia8 软**，
这是他们输入的性质，不是我们产生的错误。

## 3. 探测器没有问题

| 检验 | 结果 |
|---|---|
| 相同产生谱下，我们的响应 vs 参考接受率表（诊断报告 §3） | 1.041 / 1.013 / 0.990 / 0.955 (5/13/33/53 GeV) |
| 把磁铁从 3–7 m 移到论文的 3.5–7.5 m | 信号的 0.1–0.7% |
| 62.49 GeV 的 bb 条件末态库：我们 vs 参考 | 1.003（vs 2023 公开版 0.979） |

## 4. 残差分解（25 fb⁻¹ 长寿边，$\theta^2_{ours}/\theta^2_{ref}$）

$R \approx (\Gamma_{\rm ref}/\Gamma_{\rm ours})\,(K_{\rm ours}/K_{\rm ref})\,(\mathrm{resp}_{\rm ours}/\mathrm{resp}_{\rm ref})\,(\varepsilon_{\rm ref}/\varepsilon_{\rm ours})$

| $m_S$ | 实测 | $\Gamma$ 比 | $K$ 比 | resp 比 | $\varepsilon_{S/B}=0.87$ | 预计 | 残差 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 5.08 | 0.653 | 0.951 | 0.830 | 1.041 | 0.87 | 0.715 | 0.91 |
| 13.02 | 0.600 | 0.926 | 0.810 | 1.013 | 0.87 | 0.661 | 0.91 |
| 33.37 | 0.491 | 0.758 | 0.796 | 0.990 | 0.87 | 0.520 | 0.95 |
| 53.42 | 0.446 | 0.684 | 0.761 | 0.955 | 0.87 | 0.433 | 1.03 |
| 62.49 | 0.325 | 0.665 | 0.796 | 0.951 | 0.87 | 0.438 | 0.74 |

不带 $\varepsilon_{S/B}$ 时残差是 0.79–0.90；带上它残差是 0.91–1.03（62.5 GeV 端点另有参考封口
与 $m_h=125$ vs 125.25 的运动学差异，见诊断报告 §5）。

结论：**5–53 GeV 的差别由三件已知的事完整解释——宽度模型、参考的 h 谱、可能的 $\varepsilon_{S/B}$——
没有剩下任何可以归给"探测器 bug"或"h 产生 bug"的部分。**

## 5. 三件因素的性质与可改性

1. **宽度模型（最大项，×0.68–0.95）**：我们用 Winkler 的固定夸克质量（$m_b=4.18$），参考用
   Boiarska 的跑动质量（$m_b(m_S)\approx3.0$ GeV @40 GeV）→ $\Gamma\propto m_q^2$，40–60 GeV 处
   我们的总宽度是他们的 1.32–1.47 倍，寿命更短、信号更小。这是**用户选定的 strict-Winkler 模型**，
   要消掉它就得放弃 Winkler 的固定质量。
2. **参考的 h 谱（×0.76–0.83）**：我们的 Pythia8 谱与 FORESEE 一致，参考的表格谱更软。改它意味着
   把参考的谱当输入，属于模型选择。
3. **$\varepsilon_{S/B}=0.87$（×0.87）**：LHCb 论文 Eq. (1) 里与 $\varepsilon_{\rm rec}=0.4$ 并列的
   信号保留效率（NN 背景压低后），我们只乘了 0.4。**但本地 SensCalc 代码里只有 0.4，没有 0.87**
   （审查第三份报告 §1），所以若原图出自 SensCalc，则该因子未被使用。这一项有歧义。
4. **参考产物自身不闭合（0.55–0.62）**：本地 SensCalc 保存的曲线本身只有原 PDF 的 0.55–0.62，
   这是"参考"本身的 40% 级不确定度，无法归因。
---

# 附：SensCalc 审计结果（2026-10-25，`work/review_check/senscalc/report.md`）

| 项 | 参考工具 SensCalc | 我们（v3.6） |
|---|---|---|
| B→X_sS 处方 | `SelectedBrChoice="Exclusive"`（`3. Scalar sensitivity.nb:506`），读 `BtoXsStotal.dat`，BR/θ²(0) = **3.3096**，通道表止于 5.13 GeV | Boiarska 逐道，**3.356**（+1.4%），止于 m_B−m_π = 5.14 GeV ✓ |
| 4.7 GeV 开关 | 两个选项都**没有** | 无（逐道自然截止）✓ |
| K→πS | 只在 SPS 设施里定义，**LHCb 曲线里没有** | 我们有（长寿边占 ~1e-4，可忽略） |
| σ(pp→h) | **55 pb** 单一数字 | 55.45 pb（+0.8%）✓ |
| BR(h→SS) | 2×10⁻²α²β，BR = 0.01 | 同 ✓ |
| 效率 | **只有 ε_rec = 0.4**，全仓库没有 0.87，`NevInt` 里乘一次 | 同 ✓ |
| 接受度 | A(z) **微分**进入 z 积分；z 1.0–2.5 m、2<η<5 环柱、平面 9.4 m、偶极 3–7 m 共 4 T·m；γ/K_L 不可见 | 完全同一张卡（frontZ 7.7、backZ 9.4、场 3–7 m、By=1 T×4 m = 4 T·m）✓ |
| 亮度/判据 | N_pot = 1.75e15 = 25 fb⁻¹×70 mb，N_ev = 2.3 | L = 25 fb⁻¹，s = 2.3 ✓ |
| 总宽度 | 9.09e-6 / 1.43e-4 / 8.82e-4 / 1.39e-3 / 1.61e-3 GeV @5.08/13.0/33.4/53.4/62.5 | 9.58e-6 / 1.55e-4 / 1.16e-3 / 2.04e-3 / 2.42e-3 → **我们大 1.05/1.09/1.32/1.46/1.51 倍** |
| 母 h 谱 | `spectra/SM particles/DoubleDistr_LHC_h.txt`（无生成器说明） | Pythia8 13 TeV（与 FORESEE 14 TeV 一致 2–3%）|

结论补充：

* **探测器卡完全一致**（这是这轮审计最有价值的一条：我们的 z/η 窗口、磁铁 4 T·m、μ 可见而 γ/K_L 不可见，与参考逐项对上），所以"探测器有问题"可以彻底排除。
* **总宽度**：参考的 BR 表也自称来自 Winkler（`BrRatioScalar.mx` 默认描述 1809.01876），但高 mass 处比我们小 32–51%。差异只能来自夸克质量约定（我们的 b 用固定 m_b = 4.18 GeV；他们的等效值约 3.5 GeV）。这是**我们这一侧最大的一项高 mass 系统误差**，但 Winkler 原文只规定了 m_s、m_c，b 的处理没有唯一答案。
* **母 h 谱不是驱动项**：按接受度加权的 K 与我们的差 2–6%（尽管 dN/dη 形状差很多）。之前诊断里 0.76–0.83 的谱因子来自 SensCalc 的**预计算标量谱**（`DoubleDistr_Scalar_h-quartic_LHC.m`），它与他们自己的母谱不一致——即参考产物内部的又一处不闭合。
* SensCalc 的本地副本是 **2025 重写版**，不是产生 BC5 的 2023 版（commit 6664a92）；其自带 json 曲线只有 PDF 的 0.55–0.62。六个产生表与 2023 版字节一致，BR/ctau 表被重算过。
