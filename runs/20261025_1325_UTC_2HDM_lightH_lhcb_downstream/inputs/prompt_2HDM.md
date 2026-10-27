# 1. 目标

计算 Type-I 2HDM 的长寿命 CP-even light scalar $H$ 在 LHCb Downstream 的预期排除曲线：

$$
pp\ (\sqrt{s}=13\ \mathrm{TeV}),\qquad
\mathcal L=25,\ 300\ \mathrm{fb}^{-1},\qquad
(m_H,\tan\beta)\text{ 平面，零背景，90\% CL}.
$$

模型采用 [Li 等，Light Scalars at FASER，2212.06186v1](https://arxiv.org/abs/2212.06186v1) 的 light-$H$ 基准，依据第 2、4.1、4.4 节及附录 A、C、D；本地源码为 `references/model/arXiv-2212.06186v1/FASER_light.tex`。产生过程在论文介子道之外加入模型预测的 $h\to HH$。pp 模拟、末态演化和探测器实现参照 `prompt.md`、`runs/bc5_lhcb_downstream_v3/` 与 `detector/downstream/`。

# 2. 模型与参数

两个标量双重态 $\Phi_{1,2}\sim(\mathbf1,\mathbf2,1/2)$，只有 $\Phi_2$ 耦合所有 SM 费米子。采用论文的 CP 守恒、软破缺 $Z_2$ 势能：

$$
\begin{aligned}
V={}&m_{11}^2\Phi_1^\dagger\Phi_1+m_{22}^2\Phi_2^\dagger\Phi_2
-m_{12}^2(\Phi_1^\dagger\Phi_2+\mathrm{h.c.})\\
&+\frac{\lambda_1}{2}(\Phi_1^\dagger\Phi_1)^2
+\frac{\lambda_2}{2}(\Phi_2^\dagger\Phi_2)^2
+\lambda_3(\Phi_1^\dagger\Phi_1)(\Phi_2^\dagger\Phi_2)\\
&+\lambda_4|\Phi_1^\dagger\Phi_2|^2
+\frac{\lambda_5}{2}\big[(\Phi_1^\dagger\Phi_2)^2+\mathrm{h.c.}\big].
\end{aligned}
$$

标量拉氏量为 $\sum_i(D_\mu\Phi_i)^\dagger D^\mu\Phi_i-V$ 加 Type-I Yukawa 项。令 $\Phi_i^0=(v_i+\rho_i+i\eta_i)/\sqrt2$，$h=-\sin\alpha\,\rho_1+\cos\alpha\,\rho_2$、$H=\cos\alpha\,\rho_1+\sin\alpha\,\rho_2$。$h$ 为 125 GeV 的 SM-like Higgs，$H$ 为待探测的中性自共轭 CP-even 标量；另有 CP-odd $A$ 与 $H^\pm$。

扫描参数和固定基准为

$$
0.021\le m_H\le10\ \mathrm{GeV},\qquad 5\le t_\beta\equiv\tan\beta\le10^5,
$$

$$
m_h=125\ \mathrm{GeV},\quad v=246\ \mathrm{GeV},\quad
m_A=m_{H^\pm}=600\ \mathrm{GeV},\quad
c\equiv\cos(\beta-\alpha)=t_\beta^{-1},\quad s\equiv\sin(\beta-\alpha)=+\sqrt{1-c^2},
$$

$$
\lambda v^2\equiv m_H^2-\frac{m_{12}^2}{\sin\beta\cos\beta}=0,
\qquad m_{12}^2=m_H^2\sin\beta\cos\beta.
$$

这里的 $\lambda v^2$ 与势能中的 $\lambda_{1,\ldots,5}$ 区分；后者由物理质量、混合角和真空条件确定。耦合修正为

$$
\xi_H^V=c,\qquad
\xi_H^f=c-\frac{s}{t_\beta}=\frac{1}{t_\beta^3(1+s)},\qquad
\xi_h^V=s,\qquad\xi_h^f=s+\frac{c}{t_\beta}.
$$

$f$ 为所有 SM 夸克和带电轻子；$\xi_H^f$ 的最后一种写法避免大 $t_\beta$ 下的相消误差。带电 Higgs 环所需的三线性系数按论文约定为

$$
\mathcal L\supset\lambda_{HH^+H^-}HH^+H^-,\qquad
\lambda_{HH^+H^-}=-\frac{2m_{H^\pm}^2-m_H^2}{v}\,c.
$$

# 3. 衰变：直接接入 Light_scalar_decay

直接使用论文作者的 [Light_scalar_decay](https://github.com/shiggs90/Light_scalar_decay) 原版 C++ 程序计算 $H$ 的部分宽度、总宽度与分支比，不重写衰变模型，也不使用 BC5 或 SensCalc 的衰变表。按其 `manual.pdf` 初始化 `BSModel`：

```cpp
BSModel model;
model.set_model(1);       // Type-I；不要照搬示例中的其他 model 编号
model.set_mH(mH);
model.set_mA(600.0);
model.set_mHC(600.0);
model.set_lambdav(0.0);
model.set_tanb(tanb);
model.set_cosba(1.0 / tanb);
```

读取 `DecayWidth_HTotal()`、各 `DecayWidth_H...()` 与 `BranchRatio_H...()`，保留原程序的强子形状因子、约 2 GeV 的强子/夸克匹配以及 $W$、费米子和 $H^\pm$ 环贡献。统一使用

$$
\Gamma_H=\sum_a\Gamma_a,\qquad
\mathrm{BR}_a=\frac{\Gamma_a}{\Gamma_H},\qquad
c\tau_H=\frac{\hbar c}{\Gamma_H},\qquad
\hbar c=1.973269804\times10^{-16}\ \mathrm{GeV\,m}.
$$

$H\to\gamma\gamma$ 必须计入总宽度，但在本轮带电轨迹选择下响应为零。宽度和分支比均依赖 $(m_H,t_\beta)$，不能整体套用 BC5 的耦合缩放。原版程序的 CP-even 总宽度不含 $b\bar b$，因此本轮质量上限取 10 GeV，严格保留其原有通道与物理实现。记录程序版本、源码/压缩包校验值及必要的编译适配。

# 4. 产生与碰撞模拟

介子道沿用论文第 2.2、4.4 节的数值分析范围，包含电荷共轭：

$$
pp\to B,K,\pi+X,\qquad
B\to X_sH,\quad K\to\pi H,\quad K^\pm,\pi^\pm\to He^\pm\nu.
$$

介子跃迁使用论文的完整 2HDM 有效耦合，不能用 $\xi_H^f$ 缩放 BC5 的 FCNC 系数：

$$
\xi_H^{ij}=-\frac{4\sqrt2G_F}{16\pi^2}
\sum_{k=u,c,t}V_{ki}^*m_k^2
\left[g_1c-g_2s-g_0\frac{2v}{m_W^2}\lambda_{HH^+H^-}\right]V_{kj}.
$$

$g_{0,1,2}(m_k^2/m_W^2,m_{H^\pm}^2/m_W^2)$ 按附录 D；介子率和形状因子按第 2.2 节，$K\to\pi H$ 保留 FCNC 与四夸克振幅的相干贡献。包容 $B\to X_sH$ 与对应排他道不重复相加。保留母介子飞行、祖先关系和中性 kaon 演化。

另外纳入 `prompt.md` 中的全部 Higgs 产生类别及 $h\to HH$，产生率按上述 $\xi_h^{f,V}$ 修正；单顶伴随 Higgs 保留干涉，不能统一乘一个缩放因子。论文附录 C 的三线性系数在本基准下为

$$
g_C=\frac{s}{2v}\left[(m_H^2-m_h^2)
\left(2\cot2\beta\,sc-c^2+s^2\right)-m_H^2\right].
$$

将其换算为 $g_{hHH}^{(3)}=-2g_C$，统一顶角归一化：

$$
\mathcal L\supset-\frac12g_{hHH}^{(3)}hH^2,\qquad
\Gamma(h\to HH)=\frac{|g_{hHH}^{(3)}|^2}{32\pi m_h}
\sqrt{1-\frac{4m_H^2}{m_h^2}},\qquad
\mathrm{BR}(h\to HH)=\frac{\Gamma(h\to HH)}{\Gamma_h^{\rm other}+\Gamma(h\to HH)}.
$$

$\Gamma_h^{\rm other}$ 按同一模型的 $h$ 耦合计算，SM 参考宽度为 $4.07\times10^{-3}$ GeV；不固定 $\mathrm{BR}(h\to HH)=0.01$。

pp 生成与归一化参照 BC5：介子母样本采用 Pythia8 `SoftQCD:inelastic=on`，Higgs 硬过程采用 MG5+Pythia8，记录 PDF、尺度、tune 与截面来源；统一 GeV、m 单位。按本模型生成并验证 FeynRules/UFO，沿用本地 ColliderAgent/Magnus 流程。可复用模型无关的生成和响应代码，核验母样本是否覆盖所需 $\pi/K/B$、事件编号与总率，避免不同样本重复计数。$H$ 各通道的衰变运动学、强子化和次级衰变沿用 BC5 方法，以本模型的分支比加权。

# 5. Downstream 接受率与统计

完全沿用 BC5 的带电轨迹选择及 `detector/downstream/` 完整包络传播内核：

$$
1<z_{\rm dec}<2.5\ \mathrm m,\qquad 2<\eta_H<5,\qquad
E_{\rm track}>5\ \mathrm{GeV},\quad p_{T,\rm track}>0.5\ \mathrm{GeV}.
$$

至少一对来自同一 $H$ 的相反电荷后代在各自出生点通过动量 cuts，并通过以下磁场与 SciFi 几何选择：

$$
B_y=1\ \mathrm T\ (3<z<7\ \mathrm m),\quad B_y=0\text{（其余区域）},\qquad
z_{\rm SciFi}=7.7\text{--}9.4\ \mathrm m,
$$

$$
|x|\le3.24\ \mathrm m,\quad |y|\le2.415\ \mathrm m,\quad
\sqrt{x^2+y^2}\ge0.09\ \mathrm m.
$$

不要求 UT，不增加 PID、质量窗或顶点拟合 cuts；次级粒子与全中性末态的处理沿用 BC5。额外事件重建效率为 0.4，只乘一次。

对每个 $H$ 从实际出生点起做半解析寿命积分：

$$
p_k=\sum_a\mathrm{BR}_a\int_0^\infty\frac{d\ell}{\Lambda_k}e^{-\ell/\Lambda_k}A_{ka}(\ell),
\qquad\Lambda_k=\frac{|\mathbf p_k|}{m_H}c\tau_H.
$$

同一碰撞中至少一个 $H$ 通过即计一次。先对完整、互斥的产生与衰变历史 $\mathcal H$ 求事件 OR，再平均：

$$
\epsilon_c=0.4\left\langle1-\prod_{k\in\mathcal H}(1-p_k)\right\rangle_{c,\mathcal H},
\qquad
s=10^3\mathcal L\sum_c\sigma_c\epsilon_c,
\quad[\sigma_c]=\mathrm{pb},\quad[\mathcal L]=\mathrm{fb}^{-1}.
$$

产生分支概率包含在历史权重中，只计一次；$A_{ka}$ 不含额外 0.4 效率。零背景、零观测，取

$$
s_{90}=-\ln0.1=2.302585\ldots.
$$

# 6. 扫描与输出

对 $(m_H,\tan\beta)$ 使用对数网格，在衰变阈值、共振和边界附近加密。复用同一质量下的运动学与逐通道响应，逐参数点调用衰变程序并重算产生权重、寿命和事件 OR。对两种亮度寻找全部 $s=s_{90}$ 交点，保留所有边界；无交点时如实报告。边界附近要求 $\delta s_{\rm MC}/s\le5\%$，并检查寿命积分收敛。

输出双对数排除图 PNG/PDF、边界及 MC 误差 CSV、信号网格和宽度/分支比/寿命表；先核验作者程序输出与论文 light-$H$ 衰变图，再核验事件归一化和 OR。使用独立 `runs/<run-id>/` 保存输入、程序版本、种子、作业 ID、日志与结果，不覆盖 BC5 产物；必要输入或论文约定不明确时向用户确认。
