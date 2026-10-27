# 1. 目标

计算 BC5 Higgs 门户暗标量在 LHCb Downstream 的预期排除曲线：

$$
pp\ (\sqrt{s}=13\ \mathrm{TeV}),\qquad
\mathcal L=25,\ 300\ \mathrm{fb}^{-1},\qquad
(m_S,\theta^2)\text{ 平面，零背景，95\% CL}.
$$

---

# 2. BC5 暗标量模型

## 2.1 拉氏量

依据 [Boiarska 等，1904.10447](https://arxiv.org/abs/1904.10447)，在小混合角下取

$$
\begin{aligned}
\mathcal L={}&\mathcal L_{\rm SM}
+\frac12(\partial_\mu S)^2-\frac12m_S^2S^2
-\theta\sum_f\frac{m_f}{v}S\bar f f\\
&+\frac{2\theta m_W^2}{v}SW_\mu^+W^{-\mu}
+\frac{\theta m_Z^2}{v}SZ_\mu Z^\mu
+\frac{\alpha}{2}hS^2+\frac{\alpha}{4v}h^2S^2 .
\end{aligned}
$$

$S$ 为实的 CP-even 单态，$h$ 为 SM-like Higgs；四点项按 $\alpha_2 S^2H^\dagger H$ 展开，约定 $\alpha=2\alpha_2v$。环诱导耦合与介子有效跃迁按文献一致匹配。

## 2.2 参数

自由参数：

$$
0.021\le m_S\le62.49\ \mathrm{GeV},\qquad
10^{-20}\le\theta^2\le10^{-3}.
$$

固定参数：$v=246\ \mathrm{GeV}$、$m_h=125\ \mathrm{GeV}$、$\Gamma_h^{\rm SM}=4.07\times10^{-3}\ \mathrm{GeV}$、$\theta>0$、$\alpha>0$，且每个质量点满足

$$
b_h\equiv\mathrm{BR}(h\to SS)=0.01,\qquad
\Gamma(h\to SS)=\frac{\alpha^2}{32\pi m_h}
\sqrt{1-\frac{4m_S^2}{m_h^2}}
=\frac{b_h}{1-b_h}\Gamma_h^{\rm SM}.
$$

所有开放衰变道共同定义

$$
\Gamma_S=\sum_a\Gamma_a,\qquad
\mathrm{BR}_a=\frac{\Gamma_a}{\Gamma_S},\qquad
c\tau_S=\frac{\hbar c}{\Gamma_S}.
$$

依据 [Winkler，1809.01876](https://arxiv.org/abs/1809.01876) 及 Boiarska 的描述重新计算衰变宽度；低能强子道使用相应形状因子与匹配，不以自由夸克宽度替代。SensCalc 的分支比表及其派生数据只作交叉参考，不作为生产输入；寿命与分支比须来自同一套宽度。必要物理输入缺失时明确报告。

## 2.3 粒子

按论文第 1 节，规范基底中的实单态 $S_0\sim(\mathbf1,\mathbf1,0)$ 与 Higgs 双重态 $H\sim(\mathbf1,\mathbf2,1/2)$ 的中性 CP-even 分量混合；括号依次表示 $SU(3)_c$、$SU(2)_L$ 表象及超荷 $Y$。

| 质量本征态 | 自旋 | 电荷 | 性质 |
| --- | --- | --- | --- |
| $S$ | $0$ | $0$ | CP-even、自共轭，主要来自 $S_0$ |
| $h$ | $0$ | $0$ | CP-even，125 GeV 的 SM-like Higgs |

拉氏量中的 $f$ 为 SM 夸克和带电轻子，$W^\pm$、$Z$ 为 SM 弱规范玻色子；不引入额外暗扇区粒子或不可见暗衰变。

---

# 3. 碰撞模拟

## 3.1 物理过程

介子产生与稀有衰变：

$$
pp\to M+X,\qquad
M\to X_M S,\quad M\to X_M SS,\quad M\to SS.
$$

按 Boiarska 文献逐道纳入运动学允许的 $K$、$D$、$B$ 等介子过程及电荷共轭；保留母粒子飞行、祖先关系与必要的中性介子混合。排他道与对应包容率不重复相加。不考虑直接产生 $S$。

Higgs 产生与衰变：

$$
pp\to h\ (\mathrm{ggF}),\ hjj\ (\mathrm{VBF}),\ Wh,\ Zh,
\ t\bar th,\ b\bar bh,\ thj,\ th\bar b,\ tWh,
$$

$$
h\to SS.
$$

上述 $Zh$ 仅含非 $gg$ 分量，$thj$、$th\bar b$ 分别为 $t$、$s$ 道，包含电荷共轭。$b\bar bh$ 取 4FS，$tWh$ 取 5FS；LO+Pythia 运动学配权威总截面，明确标注形状近似。

标量衰变及后续演化：

$$
S\to\ell^+\ell^-,\ \pi\pi,\ K\bar K,\ 4\pi,\ q\bar q,\ gg,\ \gamma\gamma,\ldots
\quad\longrightarrow\quad\text{可追踪末态}.
$$

逐道覆盖开放通道；四 pion 采用均匀四体相空间近似，喷注强子化和次级衰变由 MG5+Pythia 链安排。中性 kaon 对按 CP-even 映射为等比例 $K_SK_S$、$K_LK_L$。

## 3.3 碰撞设置

- $pp$，$\sqrt{s}=13\ \mathrm{TeV}$；动量单位 GeV、传播长度单位 m。
- ColliderAgent 子代理经本地 Magnus 执行；重新建立并验证 FeynRules 模型与 UFO。
- 强子母样本统一由 Pythia8 `SoftQCD:inelastic=on` 生成；MG5 重味样本仅作对照。Higgs 硬产生使用 MG5+Pythia8。
- PDF、尺度、tune 和总率来源由 ColliderAgent 选择、核验并记录；不同来源不得重复计数。

## 3.4 背景样本

$$
b=0.
$$

这是研究假设，本轮不模拟背景。

## 3.5 信号样本

按质量分批生成母样本和衰变运动学，复用于该质量下的混合角扫描；不逐混合角保存完整事件文件。模型、宽度、样本和响应均重新计算，不复用首轮产物，也不采用 SensCalc 母谱。

样本量自适应增加，使排除曲线附近预测信号数的 MC 相对统计误差不超过 $5\%$。每次尝试使用独立 `runs/<run-id>/`，记录输入、种子、软件版本、作业 ID、日志与结果；样本分块处理，断线只收集、不自动重交。

---

# 4. 数值分析

## 4.1 观测量：通过选择的碰撞事件数

同一碰撞事件中，至少一个 $S$ 满足条件即计一次；两个 $S$ 都通过也只计一次。

### 步骤 1：计算逐标量响应

衰变顶点和方向满足

$$
1<z_{\rm dec}<2.5\ \mathrm m,\qquad 2<\eta_S<5.
$$

至少一对来自同一 $S$ 的相反电荷后代满足出生点真值

$$
E>5\ \mathrm{GeV},\qquad p_T>0.5\ \mathrm{GeV},
$$

并通过固定 Downstream 传播与 SciFi 几何选择：

$$
B_y=1\ \mathrm T\quad(3<z<7\ \mathrm m),\qquad B_y=0\text{（其余区域）},
$$

$$
z_{\rm SciFi}=7.7\text{--}9.4\ \mathrm m,\quad
|x|\le3.24\ \mathrm m,\quad |y|\le2.415\ \mathrm m,\quad
\sqrt{x^2+y^2}\ge0.09\ \mathrm m.
$$

复用 `detector/downstream/` 的完整包络传播内核，不要求穿过 UT，不增加 PID、质量窗或顶点拟合 cuts；带电次级粒子的飞行中衰变沿用忽略近似。全中性末态本身不算可见；中性母粒子的后续衰变按统一后代规则处理。

对给定运动学和衰变道，MC 求响应 $A_{ka}(\ell)$，半解析积分标量飞行距离：

$$
p_k=\sum_a\mathrm{BR}_a\int_0^\infty
\frac{d\ell}{\lambda_k}e^{-\ell/\lambda_k}A_{ka}(\ell),
\qquad \lambda_k=\frac{|\mathbf p_k|}{m_S}c\tau_S.
$$

$\ell$ 从 $S$ 的实际出生点起算；$A$ 包含衰变区域、末态运动学与轨迹选择，不含额外重建效率。

### 步骤 2：得到任意混合角的事件接受率

在小混合角基准下使用 $\Gamma_a\propto\theta^2$ 与相应产生率的耦合依赖；固定 $b_h$ 时 $\alpha$ 随 $m_S$ 调整，不随 $\theta$ 扫描。

$$
\epsilon_c(m_S,\theta^2)
=0.4\left\langle1-\prod_{k\in\mathcal S(\mathcal H)}(1-p_k)\right\rangle_{c,\mathcal H}.
$$

$\mathcal H$ 是共同产生状态及完整互斥衰变历史，包含产生分支概率；先在共同条件下求事件 OR，再平均，不能把同母互斥道独立相加或把平均单标量效率相乘。$0.4$ 为额外事件重建效率，只乘一次。

## 4.2 统计分析：零背景泊松计数

对去重后的碰撞事件类别 $c$，预期信号数为

$$
s(m_S,\theta^2;\mathcal L)
=10^3\mathcal L\sum_c\sigma_c\epsilon_c,
\qquad [\sigma_c]=\mathrm{pb},\quad [\mathcal L]=\mathrm{fb}^{-1}.
$$

### 95% CL 排除

取背景预期和观测事件数均为零：

$$
P(n=0\mid s,b=0)=e^{-s},\qquad
s_{95}=-\ln0.05\simeq3.
$$

本轮以 $s=3$ 为排除边界，分别计算 $25$、$300\ \mathrm{fb}^{-1}$。只输出中央结果与 MC 统计误差；以完整碰撞事件为统计单位，计入产生和响应抽样误差，另核验寿命积分收敛；不加入理论或重建系统误差。

## 4.3 绘图

**步骤 1：提取边界。** 每个质量点扫描 $\theta^2$，寻找全部 $s=3$ 交点；保留长、短寿命两侧边界，无交点时不补造曲线。验证曲线附近 $\delta s_{\rm MC}/s\le5\%$。

**步骤 2：绘制曲线。**

- 横轴：$m_S\ [\mathrm{GeV}]$，$[0.021,62.49]$，对数坐标。
- 纵轴：$\theta^2$，$[10^{-20},10^{-3}]$，对数坐标。
- 两条亮度曲线分别标注 $25$、$300\ \mathrm{fb}^{-1}$，显示 MC 误差带。
- 输出 PNG、PDF、边界与误差 CSV；核验宽度、归一化、事件 OR 和寿命极限，并记录近似及未解决的物理缺口。
