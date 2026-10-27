# Type-I 2HDM 轻 CP-even 标量 H 在 LHCb Downstream 的排除曲线 — 运行报告

run-id: 20261025_1325_UTC_2HDM_lightH_lhcb_downstream
日期: 2026-10-25 (Asia/Shanghai)
任务书: inputs/prompt_2HDM.md
主图: outputs/figures/exclusion_2hdm_lhcb_downstream.{png,pdf}
边界表: outputs/data/exclusion_boundary.csv

--------------------------------------------------------------------------------------------------
## 0. 摘要

在 Type-I 2HDM 的 light-H 基准（cos(beta-alpha)=1/tanb, m_A=m_Hpm=600 GeV, lambda v^2=0）
下计算了 LHCb Downstream 探测器对长寿命 CP-even H 的 90% CL 预期排除曲线
(pp, 13 TeV, L=25/300 fb^-1, 零背景零观测, s90=-ln 0.1=2.302585)。H 的宽度/分支比来自
论文作者原版 C++ 程序 Light_scalar_decay（仅做 M_PI/批处理驱动适配），产生截面用解析 2HDM
公式（App. D 的 g0,g1,g2 + 复 CKM；B->X_sH、K->piH 相干振幅、pi/K->He nu、h->HH、9 个 Higgs
过程的逐过程 kappa 标度）。探测器响应、寿命积分和事例级 OR 继承并复用了 BC5 v3.6 的
LHCb-downstream 传播内核。

> **2026-10-26 审查轮更新**：本报告第 0 节的结果已按审查轮的最终扫描更新；
> 详细审查（产生谱 vs FORESEE、产生道公式逐条审计、MC 误差诊断与修复、结果变化）
> 见 **logs/physics_code_review.md**，公式审计见 logs/production_formula_audit.md，
> 产生谱对比见 logs/spectrum_vs_foresee.md。主图改为闭合轮廓（300 fb^-1 实线、
> 25 fb^-1 虚线，上下限相连），衰变验证拆成三张独立平滑曲线图。

主要结果（物理上限 = 长寿边）：
  * 25 fb^-1 : tanb 排除上限 22.0 (mH=0.021 GeV) -> 358.0 (mH=1.0) -> 232.0 (mH=4.784)；
  * 300 fb^-1: 30.0 -> 488.9 -> 403.1；
  * 灵敏度质量范围 0.021 <= mH <= 4.786 GeV；mH > m_B-m_K 后 B->X_sH 运动学闭合，
    h->HH 贡献可忽略，无排除能力（如实报告，无交点）。
  * 长寿边 MC 误差: 中位 3.6%，最大 4.0%（两亮度），满足 <=5% 要求；
    边界位置误差中位 0.5%。
  * 低 tanb 侧（2026-10-26 修正后）：母粒子衰变位置解析积分去掉顶点抽样后，
    s(tanb=5…10) 比阈值大 1e3-1e6 倍（例 mH=4.18, tanb=5: s=2.5e4；mH=1, tanb=5: s=3.0e5），
    因此**扫描范围内（tanb>=5）不存在短寿交点**，排除区从 tanb=5 一直延伸到长寿边。
    低 tanb 侧没有被本分析闭合；tanb 远低于长寿边时 Br(B->X_sH)≳0.1，属 EFT 失效区
    （CSV 中 validity 标记）。低 tanb 各点 s 的残差 MC 误差 30-45%（高 pz 母粒子统计），
    建议保守乘 ~1.5。

--------------------------------------------------------------------------------------------------
## 1. 模型与固定参数

Type-I 2HDM，只有 Phi_2 耦合费米子；CP 守恒、软破缺 Z2 势能（论文 Eq. higgspotential）。
  m_h=125 GeV, v=246 GeV, m_A=m_Hpm=600 GeV, lambda v^2=0 (m12^2=mH^2 sinb cosb),
  c=cos(beta-alpha)=1/tanb, s=+sqrt(1-c^2).
  xi_H^V=c, xi_H^f=c-s/tanb (=1/(tanb^3(1+s))), xi_h^V=s, xi_h^f=s+c/tanb.
  lambda_HH+H- = -(2 m_Hpm^2-m_H^2) c / v.
扫描: 0.021 <= mH <= 10 GeV（66 点，见 work/analysis/grids/scan_mass_grid.txt；2.0 GeV 节点被移除，
  因为作者程序在 mH=2.0 恰好同时计算强子与部分子道，总宽度偏大 4.9%(tanb=10)/43%(tanb=5)/85%(tanb=1)），
  5 <= tanb <= 1e5（40 点对数网格）。

衰变: 作者程序 model=1, mH, mA=600, mHC=600, lambdav=0, tanb, cosba=1/tanb；
  通道: gaga, ee, mumu, tautau, pipi, KK, 4pi, cc, ss, gg（程序在 2 GeV 处做强子/部分子衔接；
  0 < mH <= 10 GeV 内无 b bbar，与程序实现一致）。Gamma_H=sum Gamma_a, BR_a=Gamma_a/Gamma_H,
  ctau=hbar c/Gamma_H, hbar c=1.973269804e-16 GeV m。

--------------------------------------------------------------------------------------------------
## 2. 产生

### 2.1 FCNC 有效耦合（论文 Eq. xi_H_2HDM + App. D）
xi_H^{ij} = -(4 sqrt2 G_F)/(16 pi^2) sum_k V*_ki m_k^2 [g1 c - g2 s - g0 (2v/m_W^2) lambda_HH+H-] V_kj
g0,g1,g2 为 App. D 的 Type-I 函数（X1,X2 形式）。CKM 用 (|Vus|,|Vcb|,|Vub|,gamma) 的标准
幺正参数化（复相位保留；与"|Vtd|,|Vts| 混合测量值"不可同时幺正，差 13% 的 K->piH，
采用幺正方案并在 validation 中记录两种）。

### 2.2 B -> X_s H
名义处方 mode='foresee'：保留论文 Eq. 的比值结构
  Br(B->XsH)/Br(B->Xc e nu)=12 pi^2 (v^2/mb^2)(1-mH^2/mB^2)^2/f(mc^2/mb^2) |xi_H^{bs}/Vcb|^2
（论文数值分析/FORESEE 用 m_B 作为相位空间端点；物理介子端点 mB-mK=4.7857 GeV 以上置零）。
运动学：把包容率按 Boiarska 的 b->s 排他道相对权重分配（K, K*, K0*, K1, K2*；去掉 b->d），
总率仍用上面的包容式（不与排他率相加；见 prompt "包容与排他道不重复相加"）。
系统变体：'paper_mb'（论文公式字面 m_b，4.18 GeV 关闭）与 'exclusive'（Boiarska 排他道求和，
各道运动学端点）——两者都单独扫描并给出曲线（outputs/figures/exclusion_systematics_bmode.png）。
BR(B0->XsH)=0.92735 BR(B+->XsH)（Gamma_B+/Gamma_B0=1.519/1.638）。
Br>0.1 (Br>1) 的区域标为 EFT 失效：mH=1 时 tanb<19.6 (<6.21)。

### 2.3 K -> pi H（FCNC 与四夸克振幅相干）
M = G_F^{1/2}2^{1/4} xi_H^W [7 lambda (mK^2+mpi^2-mH^2)/18 - 7 A_K mK^2/9]
    + (xi_H^{ds}/(2v)) m_s (mK^2-mpi^2)/(m_s-m_d) f0^{Kpi}(q^2)
Br=(1/Gamma_K)(2p*/mK)|M|^2/(16 pi mK)。K+->pi+H 与 K0L->pi0+H；K0S 取 0；PDG 311 的
1/2 K_L 投影由响应层乘。数值：FCNC^2 项占主导，|4q|^2/|FCNC|^2=1.66e-24/1.395e-21=0.12%，
相干干涉 -6.8%（mH->0, tanb=100）；与 Boiarska 主公式重标度一致到 0.3%。

### 2.4 pi/K -> H e nu
Br(X->H e nu)=sqrt2 G_F mX^4 |xi_H^W|^2/(96 pi^2 m_mu^2 (1-m_mu^2/mX^2)^2)
             * Br(X->mu nu) f(mH^2/mX^2) (7/9)^2。
运动学用 3 体相空间（电子质量 m_e，不是 m_mu；m_mu 只通过 Br(X->mu nu) 归一化进入公式）。

### 2.5 Higgs 产生与 h -> HH
SM 13 TeV 截面 [pb]: ggF 48.58, VBF 3.782, WH 1.369, ZH 0.883, bbH 0.2384, ttH 0.5065,
  tHj 0.0740, tHW 0.0152, tHb 0.0046（与 BC5 相同 HXSWG 值）。
逐过程 kappa 标度（不做全局统一缩放）：ggF: kappa_g^2（top/bottom 圈，Type-I 下=xi_h^f）；
  VBF/WH/ZH: (xi_h^V)^2；ttH/tHj/tHb: (xi_h^f)^2；bbH: (xi_h^f)^2；
  tHW: |kt A_t + kV A_V|^2 的干涉括号，基线 (kt+kV)^2/4，并按 phi in [0,pi] 给出范围。
  tH 三个过程只占 sigma_h 的 0.17%，Pythia 样本没有它们；报告里保留其截面，
  信号估计只用样本覆盖的 6 个过程（s_h 因此低 0.17%，远小于 MC）。
g_C = s/(2v)[(mH^2-mh^2)(2 cot2b s c - c^2+s^2) - mH^2]（lambda v^2=0），g_hHH^(3)=-2g_C,
Gamma(h->HH)=|g3|^2 beta/(32 pi mh), Gamma_h^other=Gamma_h^SM sum_i BR_i^SM kappa_i^2,
Gamma_h^SM=4.07e-3 GeV, BR(h->HH)=Gamma_HH/(Gamma_other+Gamma_HH)。数值:
  BR(h->HH)=3.035% (mH=1,tanb=5), 1.918e-3 (mH=1,tanb=10), 1.545e-8 (mH=1,tanb=100)。
h 的运动学复用 BC5 v3.6 的 Pythia8 HiggsSM 样本（status=62），按过程重加权到 2HDM 的 sigma_i 混合
（单顶 0.17% 未生成，记录为近似）。

--------------------------------------------------------------------------------------------------
## 3. 探测器响应与统计

固定选择与 BC5 v3.6 完全一致（detector/downstream 完整包络传播内核，work/response 下 verbatim 副本）：
  1 < z_dec < 2.5 m, 2 < eta_H < 5,
  B_y=1 T (3<z<7 m)，其余为 0,
  SciFi 7.7 <= z <= 9.4 m, |x|<=3.24 m, |y|<=2.415 m, sqrt(x^2+y^2)>=0.09 m,
  后代对：来自同一 H 的异号带电径迹，出生点 E>5 GeV, pT>0.5 GeV,
  事件级重建效率 0.4 只乘一次；不要求 UT/PID/质量窗。
响应按 10 个衰变通道分别计算 A_ch(i,z)：逐通道比逐 tanb 采样更精确，因为 2HDM 的分支比随 tanb 变化。
  p_i(tanb)=sum_ch BR_ch(mH,tanb) * int dz (m/(ctau pz)) exp(-(z-z0)m/(ctau pz)) A_ch(i,z)
寿命积分用分段线性 A(z) 的精确指数积分（BC5 v3.6 F10 版本；与解析 A=1 相符到 1e-12）。
事例级概率代数: o_hist=w_hist*(1-prod(1-p_i)); 同一母粒子的互斥历史相加；同一碰撞的母粒子 OR；
  eps_c = yield_corr_c * 0.4/N_gen * sum_events o; s = 1e3 L sum_c sigma_c eps_c
  (1 fb^-1 = 1e3 pb^-1)。B 组 σ=560 μb，B 产额修正按样本 Bs 份额运行时导出 1.018；
  inel 组 σ=79.5 mb（Pythia sigmaGen 78.05 记录但不使用，1.8% 差异记录）；
  h 组 σ_h=sum sigma_i，逐过程权重 = BR(h->HH) sigma_i/(sigma_h frac_i)。

母样本：
  B:      BC5 v3.6 bbar_heavy.txt.gz（800k 事例，1.5M B 行；作业 9d36fa12b8485053）
  inel:   本 run 新生成 minbias_pik.txt.gz（Pythia 8.316 SoftQCD:inelastic，500k 事例，
          22.8M 行 pi+-/K+-/K0；Magnus job c588d44d0bd8ac7a；σGen=78.05 mb）
  h:      BC5 v3.6 higgs.txt.gz（200k 事例，status 行；作业 bc5_higgs4）
  H 衰变强子化库: 本 run 新生成 hadron_lib_{gg,ss,cc,tautau}（按精确扫描质量网格，1500/节点；
          Magnus job 5b3b8b76029ea0af）。
扫描: 66 质量 x 40 tanb x 两个亮度；8 进程并行，两次 pass（B=3000 初扫 + B=6000 终扫，
  端点质量 4.55/4.78 GeV 用 B=20000），总共约 10 分钟墙钟。

--------------------------------------------------------------------------------------------------
## 4. 验证（按 prompt 要求的顺序）

### 4.1 作者衰变程序 vs 论文 light-H 衰变图（第一优先）
* 与作者自己的 FORESEE 2HDM_H 表（references/tools/FORESEE/Models/2HDM_H/model，同一 C++ 程序产生）逐点比：
  10 个通道 x 12261 (mH,tanb) 节点：|log10(ours/author)| 中位 1.4e-6，最大 2.7e-5，100% 在 0.1% 内；
  ctau 用作者的 hbar c=2.0e-16 约定后中位 5e-6，物理区最大 1.3e-4。结论: 我们的可执行文件
  等价于作者程序（仅 M_PI 与批处理驱动适配）。
* 论文 Fig. mphi_t1Br_tanb102_cba=1tanb.png（tanb=10）和总宽度/ctau 图颜色数字化后叠加：
  各通道中位 |log10| 0.001-0.007（小于线宽）；大值只出现在阈值/共振的竖直段（像素量化）。
  图: outputs/figures/decay_validation.png。
* 复核了 gg/gaga 预因子的"9/4"/"1/4"疑问：是 xi^g=(3/2)sum xi A 与 App. A 的 1/2 约定，
  代码=论文。作者 gg 是标准 LO 的 9/16（与论文图共享的真实缺陷），但对本曲线无影响
  （强子/gg 道 1/tanb^3 压低，tanb>=30 时 BR<1e-6；总宽度由 gamma-gamma 主导）。
  记录的近似: 作者未计带电轻子 gamma-gamma 圈（-0.27% at tanb=10）且 alpha_em=1/137。

### 4.2 产生模型
* 独立实现（work/review_check）与生产模块: xi_H^bs(1,100)=-9.037848e-6+1.756021e-7i，
  两套代码相对差 4.4e-16；App. D g0,g1,g2 经 Wolfram 15 复算，最大相对差 2e-12；
  bracket(tanb=5/100)=-0.360919/-0.018013，与 Lead 独立手算一致。
* 与作者 FORESEE 的 B->XsH 解析式比: bracket 比值 0.99294（tanb 5..1e5 常数）；
  绝对率比 0.919 (mH->0)、0.962 (1)、1.19 (2)、1.88 (3)、24 (4)——差异完全由
  (1-m^2/mB^2)^2 vs (1-m^2/mb^2)^2 解释。
* K->piH: FCNC-only 与 Boiarska 主公式重标度一致到 0.28%；4-夸克项占率 0.12%，相干干涉 -6.8%。
* h->HH: g_C 与论文近似 -mh^2 c^2/(4v) 关系在 mH->0 时比 0.9995；4700(g/v)^2 关系复现 <0.4%。
* 全部 12 个单元测试通过（work/production/tests/）。

### 4.3 事例归一化与 OR（第二优先）
* 探测器内核: 32/33 通过（唯一失败 a2 是 BC5 已文档化的规格冲突：轴向径迹落在束流孔）。
* 逐通道响应与 BC5 bc5_lib 在重叠通道上逐事例一致（bit-identical）。
* p_window 与解析 A=1 积分: 1.1e-16；与高精度 quadrature 3.7e-9。
* 事例/历史/母粒子 OR: 与朴素逐事例 MC 玩具一致（0.4-1.2 sigma）；0.4 只乘一次、1e3 pb->fb 已核。
* 独立审计: 交叉点位置用 ln s-ln tanb 插值复现 CSV 到 1.6e-6；MC 误差指标精确复现；
  s 的分量重组 s=1e3 L sum sigma_g eps_g 精确到 2e-16；b-mode 变体比与产生模块比一致到 1e-7。
* **BC5 v3.6 端到端回归**（同一引擎+BC5 模型/样本/宽度/产生，25 fb^-1，16 质量点）:
  长寿边 mine/published 中位 1.017，|dev| 中位 2.5%，最大 9.0%，13/16 在 5% 内；
  曲线比中位 0.984。相对偏差在 BC5 自身种子散布（5.5%）内；±9% 的离群点已用三种子重跑
  证明是引擎自身 MC 散布（1.090/0.929/1.016）。报告: logs/response_validation.md。

--------------------------------------------------------------------------------------------------
## 5. 主要结果表（长寿边，物理上限；单位 tanb）

mH[GeV]   25 fb^-1 (tanb +- err)   300 fb^-1          边界 MC(rel)
0.021     22.25 +- 0.10             30.35 +- 0.14       3.6%
0.0526    27.80 +- 0.13             37.92 +- 0.18       3.7%
0.101     33.24 +- 0.15             45.35 +- 0.20       3.5%
0.195     38.80 +- 0.18             52.94 +- 0.24       3.6%
0.2895    177.8 +- 0.83             242.6 +- 1.13       3.7%
0.4937    231.8 +- 1.04             316.3 +- 1.42       3.6%
0.827     290.7 +- 1.34             396.8 +- 1.83       3.7%
1.000     358.0 +- 1.58             488.9 +- 2.15       3.5%
1.594     307.9 +- 1.34             422.6 +- 1.82       3.4%
2.072     298.5 +- 1.32             414.3 +- 1.79       3.3%
3.071     298.6 +- 1.59             432.5 +- 2.01       3.3%
3.992     304.4 +- 2.00             462.9 +- 2.37       3.4%
4.551     276.7 +- 2.06             448.6 +- 2.50       3.4%
4.784     190.1 +- 2.84             352.9 +- 3.34       4.3%   (只剩 B+->K+H；B0->K0H 在 4.7820 关闭)

完整表（交点、EFT 有效性标记、亮度、MC 误差）: outputs/data/exclusion_boundary.csv；
逐交点完整列表（含变体多交点）: outputs/data/exclusion_crossings.csv。
信号网格: outputs/data/exclusion_s_curves.csv；s(mH,tanb) 全记录: outputs/data/signal_scan.json。
宽度/BR/寿命表: outputs/data/twohdm_widths_br.csv（tanb=10.7/105/1082 三列；完整 67x40 表:
outputs/data/decay_widths_grid.csv）。

物理形状（重要）:
* 长寿边不是"衰变长度边"。在长寿边处 ctau*beta*gamma 已远大于窗口（mH=1, tanb=357:
  ctau=5.1 m, beta*gamma 中位 10.7 -> Lambda~55 m，约 34 倍窗口）；寿命最优点在 tanb~63，
  深处于已排除区。真正的驱动是 s ~ tanb^-8：产生率 ~tanb^-2，可见分支比 ~tanb^-4
  （H 逐渐变成 gamma-gamma 为主），1/Lambda ~tanb^-2。实测对数斜率 -7.97。
* 阈值结构清晰: pi pi 打开 (0.27 GeV) 使上限从 ~39 跳到 ~179；K K (0.996)、tau (3.55)、
  c (3.73) 阈值处有小折点；mB-mK=4.786 处 B 道关闭，曲线终止。
* 低质量端 (mH<0.2) 只有 ee 可见道 -> 上限只有 22-39（模型有效区 19.6<tanb<上限 很窄）。

--------------------------------------------------------------------------------------------------
## 6. 系统误差与已知限制

1. B->X_s H 处方（最大模型不确定）: 以长寿边 tanb 比值计（variants/foresee，108-114 个交点）：
   paper_mb/foresee 中位 0.999（mH<4.18 时 0.98-1.00），但在 mH>4.2 后随 m_b 端点迅速降到 0.48-0.6；
   exclusive/foresee 中位 0.921（0.85-0.99），高质量端 0.53-0.8。即 mH<4.2 时 B 处方不确定在
   tanb 上约 -8%..0%，4.2-4.8 区间为 -50%..-15%（且 paper_mb 无曲线）。三种处方都给出了完整扫描，
   见 outputs/figures/exclusion_systematics_bmode.png。
2. 宽度模型: 作者程序原样（固定夸克质量、2 GeV 强子/部分子衔接、无 b bbar、无带电轻子 gg 圈、
   alpha_em=1/137）。这些都是 prompt 要求保留的作者实现；未做修改。
3. 作者程序 mH=2.0 GeV 双重计数 -> 该节点从扫描中移除（66 点）。
4. 产生归一化: sigma(bb)=560 ub（BC5 由 LHCb 数据外推，+-15%）；sigma_inel=79.5 mb 而
   Pythia sigmaGen=78.05 mb（1.8%）；B 产额修正 1.018（Bs 份额导出）。
5. h 通道: 样本覆盖 6/9 过程（tHj/tHW/tHb 0.17% 未生成），s_h 低 0.17%；单顶干涉用
   kappa 括号（全通道权重 0.027%，对 s 影响 <1e-3）。
6. 探测器: 完全沿用 BC5 的传播内核和选择；13 点 A(z) 分段线性是主要响应近似（<=0.4% relative）；
   0.4 重建效率只乘一次；未加 PID/质量窗/UT。
7. MC: 长寿边统计误差 3.6% 中位 / 3.9% 最大（<=5% 达标）；这只是统计方差，重尾情形解析
   方差低估约 1.5 倍，且不含理论/探测器系统误差。低 tanb=5、B 可产生区域的报告误差中位 6.0%、
   最大 42.1%（高 pz 母粒子统计），但该区域 s 远高于阈值，不影响排除结论。
8. 模型有效性: Br(B->XsH)>0.1 (>1) 的区域在 CSV 中用 validity 列标记；修正后低 tanb 侧没有
   物理交点（s 远高于阈值），主图下边界只是扫描下界 tanb=5。
9. 未做: 背景模拟（b=0 的研究假设）、探测器系统误差、理论系统误差、h->HH 的 NLO 修正。
10. γ 可见性：按 SensCalc/用户约定，γ 不作为末态；H→γγ 只进总宽度、响应为零。曾实现
    γ→e+e- 材料转换方案并做了单点对比（mH=1 长寿边会从 357/487 移到 5083/9461），
    确认约定后已整体回退，见 logs/photon_conversion_reverted.md。
11. 与论文 FASER/FASER2 比较（context，非同一探测器）: 我们的 25/300 fb^-1 长寿边上限
    在 mH=0.3/1/3 GeV 分别为 179/357/298 与 245/487/432，高于论文 digitised FASER
    (133/444/561) 和 1-2 个数量级低于 FASER2 (1330/4324/1e4)，阈值折点结构一致
    (outputs/figures/comparison_paper_reach.png)。

--------------------------------------------------------------------------------------------------
## 6bis. 工作区整理与重跑一致性（2026-10-25 23:20）

按用户要求将 work 目录整理为与 BC5 v3.6 同构：衰变/产生各占独立目录
(work/decay, work/production)，分析层 work/analysis（含 grids/ 与 cache/），另有
work/response、work/magnus、work/review_check（审计）。整理细节见 logs/restructure.md。
清空响应表缓存与逐质量结果后按新结构完整重跑 foresee/paper_mb/exclusive（0 失败），
与整理前快照逐点比较：**s(mH,tanb) 最大相对差 0（2640 点），264 个边界条目 0 处不一致**
（logs/restructure_reproducibility.md）。主排除图改为只画边界曲线（无阈值竖线、无区域阴影、
无 MC 带），MC 误差/EFT 标记保留在 CSV、诊断图与本报告。

## 7. 文件与可追溯性

inputs/        prompt、探测器快照、论文 tex/图、作者衰变源码、BC5 复用代码、sha256 清单
work/decay/    作者程序编译+批处理驱动+python 封装+验证 — 独立的衰变代码目录
work/analysis/grids/   mass grid、聚合、曲线、最终图、MC 误差表脚本
work/decay/    作者程序编译+批处理驱动+python 封装+验证（validation.json, PATCHES.md）
work/production/ 解析产生模块、b 道表、单元测试、验证（validation.json, cross_check_scan.md）
work/response/ 探测器内核 verbatim、逐通道 A_ch 库
work/analysis/ 引擎、模型适配、逐质量扫描 JSON（scan_2hdm/, scan_2hdm_paper_mb/, scan_2hdm_exclusive/）、
               BC5 回归、验证
work/magnus/  Magnus 作业（bundle/sha256/job id）、样本与强子化库、normalisation.json
work/review_check/ 独立解析审计、论文 reach 数字化叠加
logs/         decay_validation.md, production_validation.md, response_validation.md,
              model_audit.md, final_review.md, samples.md, lead_notes.md, run_report.md(本文件)
outputs/      data/ + figures/

关键 Magnus 作业:
  c588d44d0bd8ac7a  minbias_pik (500k events, SoftQCD:inelastic, sigmaGen 78.05 mb) Success
  5b3b8b76029ea0af  H decay hadron library on the exact 66-point (per-channel) grid Success
  (复用) 9d36fa12b8485053 bbar_heavy BC5 v3.6; bc5_higgs4 higgs BC5 v3.6

复现命令（在本 run 目录）:
  python work/analysis/launch_scan.py --workers 8 --force --nsamp '{"B":17293,"inel":2500,"h":4620}' --b-mode foresee
  python work/analysis/aggregate_scan.py --indir work/analysis/scan_2hdm --out outputs/data/signal_scan.json --b-mode foresee
  python work/analysis/final_analysis.py && python work/analysis/final_figure.py && python work/analysis/signal_map.py
  (变体: --b-mode paper_mb/exclusive --outdir work/analysis/scan_2hdm_<mode>)
