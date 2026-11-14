# Scalar production from meson decays — Boiarska et al. (arXiv:1904.10447v4)

**Source file (read-only, not modified):** `references/model/arXiv-1904.10447v4/main.tex` (1482 lines)

**Paper:** I. Boiarska, K. Bondarenko, A. Boyarsky, V. Gorkavenko, M. Ovchynnikov, A. Sokolenko, *"Phenomenology of GeV-scale scalar portal"*, arXiv:1904.10447v4 (JHEP).

**Purpose:** implementation-ready transcription of every formula and number needed to compute the branching ratios of meson decays producing a real scalar $S$ that mixes with the Higgs boson.

**Reading conventions used below**

* Equation tags are the paper's own `\label{}` names. "$Lnnnn$" = line number in `main.tex`.
* Every formula is transcribed verbatim in LaTeX. Where the printed source contains an inconsistency or a typo, it is quoted as printed and flagged in the **Gaps / ambiguities** section.
* Numbers are given with their units and uncertainties exactly as printed. Nothing has been recomputed or filled in from outside the paper.

> **Task/paper mismatch (read first).** The request asked for a table containing $\chi_B,\chi_K,\chi_D$. **This paper contains no such table and no per-meson abundances $\chi_B,\chi_K,\chi_D$.** The only abundance table is Table `tab:meson-amounts` and it contains only $\chi_{c\bar c}$ and $\chi_{b\bar b}$ (Sec. A.3). The kaon-to-$B$ multiplicity is only discussed qualitatively (L261: *"the number of kaons is at least $10^3$ times larger than the number of $B$ mesons"*). Likewise, the paper never prints numerical values of $v$, $m_h$ or $G_F$, and it never lists meson total widths $\Gamma_X$ (it cites PDG 2018 for them). See Sec. E.

---

## 0. Conventions and definitions

### 0.1 The mixing angle

$\theta$ is **the mixing angle itself**, not a separately defined $\sin\theta$. The field shift is written (L75)

$$
h \to h + \theta S .
$$

All couplings of $S$ to SM particles are linear in $\theta$, hence every production rate/BR $\propto\theta^2$. The paper uses **no** symbol $s_\theta$, and it **never prints** the relation of $\theta$ to the fundamental parameters ($\alpha_1,v,m_h,m_S$). Treat the printed $\theta$ as $\sin\theta$ at leading order in the small mixing (they coincide to $\mathcal{O}(\theta^2)$).

### 0.2 Lagrangian

**Eq. `eq:L1`** (L69–73), the scalar-portal Lagrangian:

$$
\mathcal{L} = \mathcal{L}_{SM}  + \frac{1}{2} \partial_{\mu} S \partial^{\mu} S +
(\alpha_1 S + \alpha_{2} S^2) (H^{\dagger} H) -\frac{m_{S}^{2}}{2}S^2 .
$$

**Eq. `g01`** (L76–82), the $S$–SM couplings after EWSB and the shift $h\to h+\theta S$:

$$
\mathcal{L}_{SM}^{S} = -\theta\frac{m_f}{v} S \bar{f}f +
2\theta\frac{ m_W^2}{v} S W^+ W^- +
\theta\frac{m_Z^2}{v} S Z^2 + \alpha\left(\frac{1}{4v}S^2 h^2 + \frac{1}{2}S^2 h\right),
$$

with the **definition of the quartic coupling**

$$
\alpha \equiv 2\alpha_{2} v .
$$

The quartic interaction is isolated as (Eq. `eq:quartic-coupling-1`, L347–350, and identically Eq. `eq:quartic-coupling`, L1352–1355)

$$
\mathcal{L}_{\text{quartic}} = \frac{\alpha}{2} S^2 h \quad\Big(=\tfrac{\alpha}{2}hS^2\Big).
$$

### 0.3 Momentum transfer $q^2$ (all form factors)

* **Two-body decay** $X\to X'S$: $q^2 \equiv (p_X-p_{X'})^2$ and the form factors are evaluated at $q^2 = m_S^2$. This is explicit in the BR formula `eq:productionBRmesons`, which calls $|M_{XX'}(m_S^2)|^2$; the matrix-element definition `eq:hadronic-matrix-element-1` gives it the argument $(p_X-p_{X'})^2$.
* **Three-body quartic decay** $X\to X'SS$: $q^2$ is the **invariant mass squared of the $SS$ pair**, integrated from $4m_S^2$ to $(m_X-m_{X'})^2$ (Eq. `eq:br-quartic-3-body`).
* Only $q^2$ (not the sign of $q^\mu$) enters the form factors. The paper is internally inconsistent about the sign/order of $q$: Eq. `eq:trick` defines $q_\mu \equiv p_{X'\mu}-p_{X\mu}$ (L1105), the pseudoscalar subsection says "$q = p_X - p_P$" (L1118), and `eq:hadronic-matrix-element-1` uses $(p_X-p_{X'})^2$. **Use $q^2$ positive (timelike)**; the pole factors are written $1-q^2/m^2$.

### 0.4 The hadronic matrix element $M_{XX'}$

Defined in two (slightly different) ways in the paper:

**Eq. `eq:hadronic-matrix-element-1`** (main text, L725–728):

$$
M_{XX'}((p_{X}-p_{X'})^{2}) \equiv \langle X'(p_{X'})|\bar{Q}_{i}(1+\gamma_{5})Q_{j}|X(p_{X})\rangle
$$

**Eq. `eq:meson-decay-scalar-matrix-element`** (App., L1082–1085), split by parity:

$$
M_{XX'}^{P=P'} =\langle X'(p_{X'})|\bar{Q}_{j}Q_{i}|X(p_{X})\rangle, \quad
M_{XX'}^{P \neq P'}= \langle X'(p_{X'})|\bar{Q}_{j}\gamma_{5}Q_{i}|X(p_{X})\rangle .
$$

The quark-bilinear ordering ($\bar Q_i(\ldots)Q_j$ vs $\bar Q_j(\ldots)Q_i$) differs between the two equations. Only $|M_{XX'}|^2$ enters the BRs, so an overall sign is irrelevant; keep the ordering in mind if amplitudes are summed.

**Which one applies** (L1109): $X$ is taken pseudoscalar, therefore

* $X'=$ pseudoscalar $P$ or pseudovector $A$ $\Rightarrow$ parity-even transition $\Rightarrow M^{P=P'}$ (scalar density; gives $f_0^{XP}$ / $V_0^{XA}$);
* $X'=$ scalar $\tilde S$, vector $V$ or tensor $T$ $\Rightarrow$ parity-odd transition $\Rightarrow M^{P\neq P'}$ (pseudoscalar density; gives $f_+^{X\tilde S}$, $A_0^{XV}$, $A_0^{XT}$).

The relations to the weak (charged-current) form factors are derived in the paper from Eq. `eq:meson-decay-scalar-pseudoscalar-matrix-element`, `eq:meson-decay-vector-matrix-element` and the "trick" Eq. `eq:trick` (verbatim, L1092–1108):

$$
M_{XX'}^{P = P'} \equiv \frac{1}{m_{Q_{i}} - M_{Q_{j}}}\langle X'(p_{X'})|\bar{Q}_{i}\slashed{p}_{Q}Q_{j}|X(p_{X})\rangle, \qquad
M_{XX'}^{P \neq P'} \equiv \frac{1}{m_{Q_{i}} + m_{Q_{j}}}\langle X'(p_{X'})|\bar{Q}_{i}\gamma_{5}\slashed{p}_{Q}Q_{j}|X(p_{X})\rangle,
$$

with $p_{Q}^{\mu} \equiv p^{\mu}_{Q_{i}} - p^{\mu}_{Q_{j}}$, then

$$
M_{XX'}^{P = P'} = -\frac{1}{m_{Q_{i}} - m_{Q_{j}}}q^{\mu}M_{\mu}, \qquad q_{\mu} \equiv p_{X'\mu} - p_{X\mu},
$$
$$
M_{XX'}^{P \neq P'} = -\frac{1}{m_{Q_{i}} + m_{Q_{j}}}q^{\mu}M_{XX\mu}^{5}.
$$

(The derivation of the pseudoscalar case follows **Bobeth:2001sq**; the paper generalises it to arbitrary final-state mesons — the generalised derivation is reproduced here, the original pseudoscalar-only one is not.)

### 0.5 Three-momentum of the scalar

**Unlabelled equation L695–697:**

$$
|\bm{p}_{S}| = \frac{\sqrt{(m_{X}^{2}-(m_{S}+m_{X'})^{2})(m_{X}^{2}-(m_{S}-m_{X'})^{2})}}{2m_{X}},
$$

where $\bm p_S$ is the $S$ momentum in the rest frame of $X$. By momentum conservation $|\bm p_{X'}|=|\bm p_S|$.

---

## A. Mixing with the Higgs boson (Sec. 2.1, L107–334)

### A.1 Production channels considered

Three proton-proton CM energies are used throughout (L115): $\sqrt{s_{pp}}\approx16$ GeV (DUNE), $\sqrt{s_{pp}}\approx28$ GeV (SHiP), $\sqrt{s_{pp}}=13$ TeV (LHC). The meson-production channels (this reference) dominate at SHiP/LHC energies (L233).

### A.2 Reference numbers for the direct channels (needed for the comparison, not for meson BRs)

**Eq. `eq:Pbrem`** (L130):

$$
P_{\text{brem}} = \theta^{2} g_{SNN}^{2} \mathcal{P}_{\text{brem}}(m_{S},s_{pp}),
$$

with $g_{SNN}\approx 1.2\cdot10^{-3}$ (L626) and $\mathcal{P}_{\text{brem}}$ the probability at $\theta=g_{SNN}=1$ (ranges $10^{-2}$ for DUNE/SHiP to $10^{-1}$ for LHC).

**Eq. `eq:partonic-weight`** (L841):

$$
W_{XX}(s,m_{S}) \equiv \int\limits_{m_S^2/s}^1 \frac{dx}{x} f_{X}(m_{S},x) f_{X}\left(m_{S}, \frac{m_S^2}{sx} \right).
$$

**Eq. `eq:DISsigma`** (L153–157):

$$
\sigma_{\text{DIS},G} \sim \frac{\theta^2 \alpha_{s}^2(m_S) m_S^2}{s_{pp} v^2} W_{GG}(m_S,s_{pp}),
\qquad
\sigma_{\text{DIS},q} \sim \frac{\theta^2 m_q^2}{s_{pp} v^2} W_{q\bar{q}}(m_S,s_{pp}).
$$

**Table `tab:WDIS`** — dimensionless gluon-gluon partonic weight $W_{GG}$ (Eq. `eq:partonic-weight` above; counts parton pairs in two protons that can fuse to an $S$), per proton-proton collision at the stated CM energy:

| $W_{GG}$ | DUNE ($\sqrt{s_{pp}}=16$ GeV) | SHiP ($\sqrt{s_{pp}}=28$ GeV) | LHC ($\sqrt{s_{pp}}=13$ TeV) |
|---|---|---|---|
| $W_{GG}(m_S=2\text{ GeV})$ | $1.2\cdot 10^{3}$ | $1.4\cdot 10^{3}$ | $6.2\cdot 10^{10}$ |
| $W_{GG}(m_S=5\text{ GeV})$ | $1.5\cdot 10^{-1}$ | $7.9$ | $3.9 \cdot 10^{10}$ |

### A.3 Table of meson amounts — `tab:meson-amounts` (L239–253)

**Normalisation (stated in the caption):** $\chi_{q\bar{q}} = \sigma_{q\bar{q}}/\sigma_{pp}$, i.e. the production fraction of a $q\bar{q}$ pair in $pp$ collisions (dimensionless, per $pp$ collision) at the given CM energy.

| Experiment | DUNE ($\sqrt{s_{pp}}=16$ GeV) | SHiP ($\sqrt{s_{pp}}=28$ GeV) | LHC ($\sqrt{s_{pp}}=13$ TeV) |
|---|---|---|---|
| $\chi_{c\bar{c}}$ | $1.0\cdot 10^{-4}$ | $3.9\cdot 10^{-3}$ | $2.9\cdot 10^{-2}$ |
| $\chi_{b\bar{b}}$ | $1.0\cdot 10^{-10}$ | $2.7\cdot 10^{-7}$ | $8.6 \cdot 10^{-3}$ |

Sources: DUNE/SHiP values from **Lourenco:2006vw**, **CERN-SHiP-NOTE-2015-009**; LHC values computed with FONLL (**Cacciari:1998it**) for the $B/D$ cross sections and PDG (**Tanabashi:2018oca**) for the total $pp$ cross section.

**Raw-table ambiguity / gap.** The requested per-meson fractions $\chi_B,\chi_D,\chi_K$ do **not** exist in this paper. Only the $q\bar q$-pair fractions above are tabulated. $\chi_h\sim10^{-9}$ (LHC Higgs production fraction) is quoted in passing (L371, citing **Heinmeyer:2013tqa**); $\chi_{b\bar b}\sim10^{-2}$ is quoted as a round number there.

### A.4 Inclusive production estimate (quark-level)

**Eq. `eq:prodBRinclusiveB`** (L201–204):

$$
\text{BR}_{\text{incl}}(X_{Q_{i}} \to X_{Q_{j}}S) \simeq \frac{\Gamma(Q_{i}\to Q_{j}S )}{\Gamma(Q_{i} \to Q' e\bar{\nu}_{e})}\times\text{BR}_{\text{incl}}(X\to X_{Q'}e\bar{\nu}_{e}),
$$

following **Chivukula:1988gp**, **Curtin:2018mvb**. Here $\Gamma(Q_i\to Q'e\bar\nu_e)$ is the free-quark semileptonic width (Fermi theory) and $\text{BR}_{\text{incl}}(X\to X_{Q'}e\bar\nu_e)$ is the **experimentally measured** inclusive semileptonic BR (so QCD corrections largely cancel).

**Unlabelled width** (L687–689):

$$
\Gamma_{Q_{i}\to Q_{j} +S} = \frac{|\mathcal{M}_{Q_{i}\to Q_{j}S}|^{2}}{8\pi m_{b}}\frac{|\bm{p}_{S}|}{m_{b}} \approx |\xi_{bs}|^{2}\frac{m_{b}^{3}\left(1-\frac{m_{S}^{2}}{m_{b}^{2}} \right)^{2}}{32 \pi v^{2}}\theta^{2}.
$$

(As printed the prefactor uses $m_b$; this equation is written for the $b$-quark case. Only $|\xi_{bs}|$ appears.)

**Unlabelled semileptonic width** (L690–693):

$$
\Gamma_{Q_{i}\to Q_{k} +e+\bar{\nu}_{e}} = \frac{1}{(2\pi)^{3}}\int \limits_{m_{Q_{k}}^{2}}^{m_{Q_{i}}^{2}}ds_{Q_{k}e} \int \limits_{s_{e\nu,\text{min}}}^{s_{e\nu,\text{max}}} ds_{e\nu} \frac{|\mathcal{M}_{Q_{i}\to Q_{k}+e+\bar{\nu}_{e}}|^{2}}{32m_{Q_{i}}^{3}} \approx
\frac{G_{F}^{2}|V_{Q_{i}Q_{k}}|^{2}m_{Q_{i}}^{5}}{192\pi^{3}}\times f(m_{Q_{k}}/m_{Q_{i}}),
$$

with integration limits (L699–701)

$$
s_{e\nu,\text{min}} = 0, \quad s_{e\nu,\text{max}}= m_{Q_{i}}^{2}+m_{Q_{i}}^{2}-s_{Q_{k}e}-\frac{m_{Q_{i}}^{2}m_{Q_{k}}^{2}}{s_{Q_{k}e}},
$$

and the phase-space factor (L703–705)

$$
f(m_{Q_{k}}/m_{Q_{i}}) = \bigg(1-8\frac{m_{Q_{k}}^{2}}{m_{Q_{i}}^{2}}-24\frac{m_{Q_{k}}^{4}}{m_{Q_{i}}^{4}}\ln\left(\frac{m_{Q_{k}}}{m_{Q_{i}}}\right)+8\frac{m_{Q_{k}}^{6}}{m_{Q_{i}}^{6}}-\frac{m_{Q_{k}}^{8}}{m_{Q_{i}}^{8}} \bigg) \approx 1/2 .
$$

**Inclusive production probabilities, Eqs. `eq:PD` and `eq:PB`** (L207–215):

$$
P_{D} \sim 2\chi_{c\bar{c}} \times \text{BR}(c\to Su) \sim 6\cdot 10^{-11} \ \theta^{2} \chi_{c\bar{c}} \left( 1-\frac{m_{S}^{2}}{m_{c}^{2}}\right)^{2},
$$

$$
P_{B} \sim 2\chi_{b\bar{b}} \times \text{BR}(b\to Ss) \sim 13 \ \theta^{2}\chi_{b\bar{b}}\left( 1-\frac{m_{S}^{2}}{m_{b}^{2}}\right)^{2}.
$$

Interpretation: $P_{D}$, $P_{B}$ are **scalar-production probabilities per $pp$ collision** from the decays of the produced $D$/$B$ hadrons. The factor $2\chi_{q\bar q}$ is the number of heavy hadrons per $pp$ collision (each $q\bar q$ pair gives two), and $\text{BR}(b\to Ss)$ is the inclusive partonic branching per heavy quark.

**Source of the $10^{-11}$ vs $\mathcal{O}(1)$ hierarchy** (L215): mostly $(m_b/m_t)^4\sim10^{-7}$ and $|V_{ub}|^2/|V_{ts}|^2\sim10^{-2}$.

**Cross-check values quoted in the text:**

* $D$ mesons: the leptonic channel $D\to S e\nu$ with $\text{BR}(D\to S e\nu)\sim 5\cdot 10^{-9}\theta^{2}$ is *more important* than `eq:PD` (L216).
* **Eq. `eq:Pbrem2B`** (L219–221):

$$
\frac{P_{\text{brem}}}{P_{B}} \sim \frac{g_{SNN}^2}{\text{BR}(b\to Ss)} \frac{\mathcal{P}_{\text{brem}}}{\chi_{b\bar{b}}} \sim 10^{-7} \frac{\mathcal{P}_{\text{brem}}}{\chi_{b\bar{b}}}.
$$

* **Eq. `eq:pDIS-to-B`** (L223–227):

$$
\frac{P_{\text{DIS}}}{P_{B}} \sim 10^{-6} \frac{1}{s_{pp} \sigma_{pp}} \left(\frac{m_S}{1\text{ GeV}}\right)^2
\frac{W_{GG}(m_S,s_{pp})}{\chi_{b\bar{b}}}.
$$

* **Eq. `eq:total-exclusive`** (L318–321), sum over all $K$ final states at $m_S\ll m_B$:

$$
\text{BR}(B \to SX_{s}) \equiv \sum_{X_{s}}\text{BR}(B\to SX_{s}) \approx 3.3^{+0.8}_{-0.7}\ \theta^{2}.
$$

* **Eq. `eq:exclusive-to-inclusive`** (L323–326):

$$
\text{BR}(B \to SX_{s})/ \text{BR}_{\text{incl}}(B \to SX_{s}) \approx 0.5 .
$$

### A.5 Table of two-body BRs via mixing — `tab:BR` (L263–285)

Column 2 is $\text{BR}(m_S=0)/\theta^2$ (i.e. the BR **normalised by $\theta^2$**, evaluated at $m_S=0$, using Eq. `eq:productionBRmesons`). Values are for $B^\pm$; for $B^0$ multiply by $0.93$ (see A.6). Column 3 is the *effective closing mass*: the $m_S$ at which the BR has dropped by a factor 10. Column 4 points to the form-factor appendix.

| Process | $\text{BR}(m_S=0)/\theta^2$ | Closing mass [GeV] | Appendix |
|---|---|---|---|
| $K^\pm\to S\pi^{\pm}$ | $1.7 \cdot 10^{-3}$ | $0.350$ | `app:pseudoscalar` |
| $K^0_L\to S\pi^0$ | $7 \cdot 10^{-3}$ | $0.360$ | `app:pseudoscalar` |
| $B^{\pm} \to S K^{\pm}_{1}(1270)$ | $(9.1^{+3.6}_{-4.0})\cdot 10^{-1}$ | 3.82 | `app:pseudovector` |
| $B^{\pm} \to S K^{*,\pm}_{0}(700)$ | $7.6\cdot 10^{-1}$ | 4.27 | `app:scalar` |
| $B^{\pm} \to S K^{*,\pm}(892)$ | $(4.7^{+0.9}_{-0.8})\cdot 10^{-1}$ | 4.29 | `app:vector` |
| $B^{\pm} \to S K^{\pm}$ | $(4.3^{+1.1}_{-1.0})\cdot 10^{-1}$ | 4.79 | `app:pseudoscalar` |
| $B^{\pm} \to S K^{*,\pm}_{2}(1430)$ | $3.0\cdot 10^{-1}$ | 3.85 | `app:tensor` |
| $B^{\pm} \to S K^{*,\pm}(1410)$ | $(2.1^{+0.6}_{-1.1})\cdot 10^{-1}$ | 3.57 | `app:vector` |
| $B^{\pm} \to S K^{*,\pm}(1680)$ | $(1.3^{+0.5}_{-0.4})\cdot 10^{-1}$ | 3.26 | `app:vector` |
| $B^{\pm} \to S K^{*,\pm}_{0}(1430)$ | $8.1\cdot 10^{-2}$ | 3.82 | `app:scalar` |
| $B^{\pm} \to S K^{\pm}_{1}(1400)$ | $(1.6^{+0.6}_{-1.1})\cdot 10^{-2}$ | 2.28 | `app:pseudovector` |
| $B^{\pm} \to S \pi^{\pm}$ | $(1.3^{+0.3}_{-0.3})\cdot 10^{-2}$ | 5.14 | `app:pseudoscalar` |

> **Raw-source note:** in `main.tex` the row `K^0_S\to S\pi^0` is **commented out** (L271): `% $K^0_S\to S\pi^0$ & $1.2 \cdot 10^{-5}$ & $0.360$ & \ref{app:pseudoscalar}`. It is not part of the published table because the $K_S^0\to\pi S$ rate is assumed to vanish (see B.4).

Caption additions: uncertainties come from meson transition form-factor uncertainties; the $m_S\to0$ total exclusive BR sum is Eq. `eq:total-exclusive`.

### A.6 $B^\pm$ vs $B^0$ normalisation

Quoted repeatedly (L283, L313, L418): the form factors are the same for $B^+$ and $B^0$, so

$$
\text{BR}(B^0\to X'^0 S) = \frac{\Gamma_{B^{+}}}{\Gamma_{B^{0}}}\,\text{BR}(B^+\to X'^+ S) \approx 0.93\,\text{BR}(B^+\to X'^+ S).
$$

(The text at L313 mislabels the first factor as $B^0$; the intended statement is the ratio of $B^+\to B^0$ BRs.) This ratio is taken from PDG lifetimes (**Tanabashi:2018oca**).

### A.7 Production via the quartic coupling — summary table `tab:b-production-quartic` (L394–420)

**Normalisation:** column 2 is $\text{BR}(m_S=0)$ **at $\alpha=1$ GeV** (no $\theta$ normalisation; these channels are $\theta$-independent and scale as $\alpha^2$). Values are for $B^\pm$; multiply by $0.93$ for $B^0$. Column 3 is the effective closing mass (BR drops by factor 10). Formulas: `eq:bstossbranching`, `eq:br-quartic-3-body` (Sec. D).

| Process | $\text{BR}(m_S=0)$ at $\alpha=1$ GeV | Closing mass [GeV] | Appendix |
|---|---|---|---|
| $K^0_L\to SS$ | $4.4\cdot 10^{-13}$ | $0.252$ | `sec:quartic-coupling` |
| $K^0_L\to SS\pi^0$ | $6.6\cdot 10^{-15}$ | $0.140$ | `app:pseudoscalar` |
| $K^\pm\to SS\pi^{\pm}$ | $1.4\cdot 10^{-15}$ | $0.136$ | `app:pseudoscalar` |
| $K^{0}_{S}\to SS$ | $7.8\cdot 10^{-16}$ | $0.252$ | `sec:quartic-coupling` |
| $K^0_S\to SS\pi^0$ | $1.2\cdot 10^{-17}$ | $0.140$ | `app:pseudoscalar` |
| $B_{s}\to SS$ | $4.0\cdot 10^{-10}$ | $2.670$ | `sec:quartic-coupling` |
| $B^{\pm} \to SS K^{\pm}$ | $1.4\cdot 10^{-10}$ | $1.998$ | `app:pseudoscalar` |
| $B^{\pm} \to SS K^{*,\pm}_{0}(700)$ | $1.2\cdot 10^{-10}$ | $1.621$ | `app:scalar` |
| $B^{\pm} \to SS K^{\pm}_{1}(1270)$ | $(1.2^{+0.5}_{-0.5})\cdot 10^{-10}$ | $1.478$ | `app:pseudovector` |
| $B^{\pm} \to SS K^{*,\pm}(892)$ | $9.1\cdot 10^{-11}$ | $1.701$ | `app:vector` |
| $B^{\pm} \to SS K^{*,\pm}_{0}(1430)$ | $3.5\cdot 10^{-11}$ | $1.621$ | `app:scalar` |
| $B^{\pm} \to SS K^{*,\pm}(1410)$ | $(1.9^{+0.6}_{-0.5})\cdot 10^{-11}$ | $1.358$ | `app:vector` |
| $B^{\pm} \to SS K^{*,\pm}_{2}(1430)$ | $2.5\cdot 10^{-11}$ | $1.499$ | `app:tensor` |
| $B_{0}\to SS$ | $1.1\cdot 10^{-11}$ | $2.624$ | `sec:quartic-coupling` |
| $B^{\pm} \to SS K^{*,\pm}(1680)$ | $(9.9^{+0.4}_{-0.3})\cdot 10^{-12}$ | $1.240$ | `app:vector` |
| $B^{\pm} \to SS \pi^{\pm}$ | $(4.7^{+1.2}_{-1.1})\cdot 10^{-12}$ | $2.149$ | `app:pseudoscalar` |
| $B^{\pm} \to SS K^{\pm}_{1}(1400)$ | $(7.3^{+0.3}_{-0.3})\cdot 10^{-13}$ | $1.545$ | `app:pseudovector` |

Related order-of-magnitude comparisons for $m_S\ll m_X$ (L376–387):

$$
\frac{\text{BR}(X\to X'SS)}{\text{BR}(X\to X'S)} \sim \frac{\alpha^2 m_X^2}{\theta^2 m_h^4} \sim 10^{-9} \frac{1}{\theta^2} \left(\frac{\alpha}{1\text{ GeV}}\right)^{2}\left(\frac{m_X}{1\text{ GeV}}\right)^{2},
$$

$$
\frac{\text{BR}(X\to SS)}{\text{BR}(X\to X'S)} \sim \frac{\alpha^2 f_X^2}{\theta^2 m_h^4} \sim 10^{-9} \frac{1}{\theta^2} \left(\frac{\alpha}{1\text{ GeV}}\right)^{2}\left(\frac{f_X}{1\text{ GeV}}\right)^{2},
$$

and the Higgs-mediated production (Eq. `eq:BR-h-to-SS`, `eq:Pss-to-PB`, L354–370):

$$
\text{BR}(h\to S S) = \frac{\alpha^2|\bm{p}_S|}{16\pi m_h^2 \Gamma_h} \approx 2.0\cdot 10^{-2}\left(\frac{\alpha}{1\text{ GeV}}\right)^{2}\sqrt{1 - \frac{4m_{S}^{2}}{m_{h}^{2}}},
\qquad \Gamma_h\approx4\text{ MeV},
$$

$$
\text{BR}(h \to \text{invis.}) < 0.19 \ (2\sigma,\ \text{CMS}) \ \Rightarrow\ \alpha< 3.1\text{ GeV for } m_S<m_h/2,
$$

$$
\frac{P_{h\to SS}}{P_B} \sim 10^{-3} \frac{1}{\theta^2} \left(\frac{\alpha}{1\text{ GeV}}\right)^{2} \frac{\chi_h}{\chi_{b\bar{b}}}
    \sim 10^{-10} \frac{1}{\theta^2} \left(\frac{\alpha}{1\text{ GeV}}\right)^{2}.
$$

---

## B. Flavour-changing effective Lagrangian and meson-decay BRs (L639–767)

### B.1 Effective Lagrangian and $\xi_{ij}$

**Eq. `eq:g05`** (L650–653):

$$
\mathcal{L}_{eff}^{Sqq}= \theta\ \frac{S}{v}\sum_{i,j}\xi_{ij}m_{Q_j}\bar{Q}_{i}P_{R}Q_{j}+h.c.,
$$

where $Q_i,Q_j$ are both up-type or both down-type quarks and $P_R\equiv(1+\gamma_5)/2$. **Note** the mass multiplying the operator is $m_{Q_j}$ (the *final*, lighter quark in this indexing).

**Definition of $\xi_{ij}$** (unlabelled, L655–657):

$$
\xi_{ij} = \frac{3G_{F}\sqrt{2}}{16\pi^{2}}\sum_{k}V_{ki}^{*}m_{k}^{2}V_{kj},
$$

where the sum runs over the quarks $Q_k$ of the **opposite** type to $Q_i,Q_j$ (i.e. the virtual $t$ for down-type transitions, the virtual $b$ for up-type), $V$ is the CKM matrix, and $G_F$ is the Fermi constant. One power of the quark mass comes from the $h\bar qq$ coupling, the second from the helicity flip. Transitions among down-type quarks (through virtual $t$) dominate.

**Table `tab:xi`** — numerical $\xi_{ij}$ (dimensionless):

| $\xi_{ij}$ | $\xi_{ds}$ | $\xi_{uc}$ | $\xi_{db}$ | $\xi_{sb}$ |
|---|---|---|---|---|
| Value | $3.3 \cdot 10^{-6}$ | $1.4 \cdot 10^{-9}$ | $7.9 \cdot 10^{-5}$ | $3.6 \cdot 10^{-4}$ |

No uncertainties are given for these. Related flavour-changing amplitudes were computed in **Willey:1982ti, Willey:1986mj, Grzadkowski:1983yp, Leutwyler:1989xj, Haber:1987ua, Chivukula:1988gp**.

### B.2 Matrix element (general 2-body transition)

**Eq. `eq:matrixelementhhS`** (L719–723):

$$
\mathcal{M} (X_{Q_{i}} \to S X'_{Q_{j}}) =
\frac{\theta}{2}\frac{m_{Q_{i}}}{v}\times  \xi_{ij}\times M_{XX'}(m_{S}^{2}).
$$

Here $m_{Q_i}$ is the **initial (heavier) quark mass** and $M_{XX'}(m_S^2)$ is the hadronic matrix element of Sec. 0.4 evaluated at $q^2=m_S^2$.

### B.3 Two-body branching ratio

**Eq. `eq:productionBRmesons`** (L730–733) — the master formula for every $X_{Q_i}\to X'_{Q_j}S$ channel (all resonances in Table `tab:BR`), from **Tanabashi:2018oca**:

$$
\text{BR}(X_{Q_{i}}\to X'_{Q_{j}}+S) = \frac{1}{\Gamma_{X}}\theta^{2}\frac{|\xi_{ij}|^{2}m_{Q_{i}}^{2}|M_{XX'}(m_{S}^{2})|^{2}}{32\pi v^{2}}\frac{|\bm{p}_{S}|}{m_{X}^{2}},
$$

where $\Gamma_X$ is the **total** decay width of $X$ ("We use the lifetimes of mesons from [PDG 2018]" — L734; the numerical lifetimes are *not* reproduced in this paper). $|\bm p_S|$ is Eq. L695–697 (Sec. 0.5).

**The paper provides no resonance-specific width formula.** Every channel of Table `tab:BR` (and Table `tab:b-production-quartic`) uses this single expression together with the corresponding matrix element $M_{XX'}$ of Sec. C. The decay width follows trivially as $\Gamma=\text{BR}\cdot\Gamma_X$.

### B.4 The $K\to\pi S$ case (L736–740)

The only possible 2-body kaon decay is

$$
K \to \pi +S .
$$

Three kaon states exist: $K^\pm$, $K^0_L$, $K^0_S$. All are governed by the **same loop factor $\xi_{sd}$**, but the BRs differ because (i) the total widths differ, and (ii) $K^0_S$ is approximately the $CP$-even eigenstate, so $K^0_S\to\pi S$ is proportional to the CKM $CP$-violating phase and is strongly suppressed (cite **Leutwyler:1989xj** — the suppression argument is *not* reproduced in this paper). **The paper assumes $\text{BR}(K^0_S\to\pi S)=0$.** Only $K^0_L$ and $K^\pm$ appear in Table `tab:BR`.

### B.5 Three-body leptonic decay $X\to S e\nu$

**Unlabelled formula, L745–751** (from **Dawson:1989kr, Cheng:1989ib**):

$$
\text{BR} (X \rightarrow S e \nu ) =
  \dfrac{\sqrt{2} G_F m_X^4}{96 \pi^2 m_{\mu}^2 (1 - m_{\mu}^2/m_h^2)^2}
  \times \text{BR} (X \rightarrow \mu \nu) \left(\frac{7}{9}\right)^2 f\left( \frac{m_S^2}{m_X^2}\right),
$$

$$
f(x) = (1 - 8x + x^2) (1 - x^2) - 12 x^2 \ln x .
$$

**Normalisation / ambiguity (flagged):** as printed, this equation contains **no explicit $\theta^2$**, whereas Table `tab:BR-3-body` divides by $\theta^2$ and the text (L216) writes "$\text{BR}(D\to S e\nu)\sim 5\cdot 10^{-9}\theta^{2}$". So the physical scalar BR is intended to carry $\theta^2$; treat the printed formula as the $\theta$-independent factor and multiply by $\theta^2$ (or, equivalently, read the table entries as $\text{BR}/\theta^2$ per unit $f$). The printed formula is also specific to the electron channel and uses $m_\mu$ from the reference $\text{BR}(X\to\mu\nu)$ normalisation; $m_\mu$ has no printed numerical value.

**Table `tab:BR-3-body`** (L753–766). Column header as printed: $\text{BR}(h\to Se\nu)/f(x)\theta^2$ (the "$h$" is a typo for the meson $X$):

| Meson | $\text{BR}(X \rightarrow S e \nu )/ f\left( x\right) \theta^2$ |
|---|---|
| $D \to S e \nu$ | $5.2 \cdot 10^{-9}$ |
| $K \rightarrow S e \nu$ | $4.1 \cdot 10^{-8}$ |
| $B \rightarrow S e \nu$ | $< 7.4 \cdot 10^{-10}$ |

Caption note: only an experimental upper bound exists on $\text{BR}(B\to\mu\nu)$, hence the bound on $B\to Se\nu$. Text (L767): for $D$ this channel enhances production by $\mathcal{O}(100)$, but $D$ production remains sub-dominant.

---

## C. Form factors for the FCNC meson decays (App., L1079–1348)

General setup: Eq. `eq:meson-decay-scalar-matrix-element`, `eq:meson-decay-vector-matrix-element` in Sec. 0.4. For every transition below, $M_{XX'}$ means the scalar object entering `eq:productionBRmesons` as $|M_{XX'}(m_S^2)|^2$; $q^2=(p_X-p_{X'})^2=m_S^2$ for the 2-body decay.

### C.1 Pseudoscalar final state $X'=P$ — `app:pseudoscalar`

Weak matrix element (L1114–1118, from **Ebert:1997mg**):

$$
M^{\mu}_{XP} = \langle P(p_{P})|\bar{Q}_{i}\gamma^{\mu}Q_{j}|X(p_{X})\rangle =
\left[(p_{X}+p_{P})^{\mu} - \frac{m_{X}^{2} - m_{P}^{2}}{q^{2}}q^{\mu}\right]f^{XP}_{1}(q^{2})+\frac{m_{X}^{2} - m_{P}^{2}}{q^{2}}q^{\mu}f^{XP}_{0}(q^{2}),
\qquad q = p_{X} - p_{P}.
$$

Contraction (L1121–1126):

$$
q_{\mu}M^{\mu}_{XP} = (m_{X}^{2} - m_{P}^{2})f^{XP}_{0}(q^{2})
\ \Longrightarrow\
M_{XP} = \frac{m_{X}^{2} - m_{P}^{2}}{m_{Q_{j}} - m_{Q_{i}}}f^{XP}_{0}(q^{2}).
$$

**Eq. `eq:pseudoscalar-form-factor`** (L1129–1131, from **Ball:2004ye**), single-pole form:

$$
f^{XP}_{0}(q^{2})=\frac{F_{0}^{XP}}{1-q^{2}/(m^{X}_{\text{fit}})^{2}} .
$$

**Table `tab:pseudoscalar-form-factor-parameters`** (L1134–1147; sources **Ball:2004ye, Marciano:1996wy**):

| $X, P$ | $B^{+/0}, K^{+/0}$ | $B^{+/0},\pi^{+/0}$ | $K,\pi$ |
|---|---|---|---|
| $m_{\text{fit}}^{X}$ [GeV] | $6.16$ | $6.16$ | $\infty$ |
| $F_{0}^{XP}$ | $0.33\pm 0.04$ | $0.258\pm 0.031$ | $0.96$ |

For $K\to\pi$ the pole is removed ($m_{\text{fit}}=\infty$), i.e. $f_0^{K\pi}=F_0^{K\pi}=0.96$ = constant.

### C.2 Scalar final state $X'=\tilde S$ — `app:scalar`

Weak matrix element (L1152–1156, from **Sun:2010nv**; the paper used $f_+(q^2)=-f_-(q^2)$ in Eq. (6) of Sun:2010nv):

$$
M^{\mu}_{X\tilde{S}} = -i\left[(p_{X}+p_{\tilde{S}})^{\mu} -q^{\mu}\right]f^{X\tilde{S}}_{+}(q^{2}).
$$

Matrix element (L1157–1159):

$$
M_{X\tilde{S}} =  i\frac{m_{X}^{2} - m_{\tilde{S}}^{2} - q^2}{m_{Q_{j}} + m_{Q_{i}}}f^{X\tilde{S}}_{+}(q^{2}).
$$

**Eq. `eq:scalar-form-factor`** (L1162–1167) — pole-like fit to **Sun:2010nv**:

$$
f^{BK_{0}^{*}}_{+}(q^{2}) = \frac{F^{BK_{0}^{*}}_{0}}{1-a\frac{q^{2}}{m_{B}^{2}}+b\left(\frac{q^{2}}{m_{B}^{2}}\right)^{2}},
\qquad m_{B} = 5.3\text{ GeV (mass of the } B^+ \text{ meson).}
$$

**Table `tab:scalar-form-factor-parameters`** (L1169–1182; source **Sun:2010nv**). The caption writes $\tilde S = K_0^{*0}(700),\,K_0^*(1430)$ but the values are the $B\to K_0^{*}$ (charged) parameters used for $B^\pm$:

| $\tilde{S}$ | $F_{0}^{B\tilde{S}}$ | $a$ | $b$ |
|---|---|---|---|
| $K_{0}^{*}(700)$ | $0.46$ | $1.6$ | $1.35$ |
| $K_{0}^{*}(1430)$ | $0.17$ | $4.4$ | $6.4$ |

Modelling assumption (L1160): $K_0^*(700)$ is taken as a di-quark state and $K_0^*(1430)$ as its excited state; there are no observed $B\to K_0^*(700)X$ decays, so the form-factor uncertainty is large (see **Issadykov:2015iba**, **Daldrop:2012sr**, **Cheng:2013fba**). The results for $K_0^*(1430)$ agree with $B\to K_0^*(1430)\eta'$ data.

### C.3 Vector final state $X'=V$ — `app:vector`

Weak matrix elements (L1189–1201, from **Ebert:1997mg, Ball:2004rg**):

$$
\langle V(p_{V})| \bar{Q}_{i}\gamma^{\mu}\gamma_{5}Q_{j} |X(p_{X})\rangle = (m_{X}+m_{V})\epsilon^{\mu *}(p_{V})A_{1}(q^{2}) -
(\epsilon^{*}(p_{V})\cdot q)(p_{X} + p_{V})^{\mu}\frac{A_{2}(q^{2})}{m_{X}+m_{V}}-2m_{V}\frac{\epsilon^{*}(p_{V})\cdot q}{q^{2}}q^{\mu}(A_{3}(q^{2})-A_{0}(q^{2})),
$$

$$
\langle V(p_{V})| \bar{Q}_{i}\gamma^{\mu}Q_{j} |X(p_{X})\rangle= \frac{2V(q^{2})}{m_{X}+m_{V}}i\epsilon^{\mu\nu\rho\sigma}\epsilon^{*}_{\nu}(p_{V})p_{X, \rho}p_{V, \sigma},
$$

$$
A_{3}(q^{2}) = \frac{m_{X}+m_{V}}{2m_{V}}A_{1}(q^{2})-\frac{m_{X} - m_{V}}{2m_{V}}A_{2}(q^{2}).
$$

The vector part vanishes upon contraction with $q_\mu$; the axial part gives (L1203–1209)

$$
M_{XV} = \langle V(p_{V})| \bar{Q}_{i}\gamma_{5}Q_{j} |X(p_{X})\rangle = -\frac{(\epsilon^{*}(p_{V}) \cdot p_{X})}{m_{Q_{i}}+m_{Q_{j}}}2m_{V}A^{XV}_{0}(q^{2})
\ \xrightarrow[\text{X rest frame}]{} \
M_{XV} = -\frac{2m_{X}|\bm{p}_{V}|}{m_{Q_{i}}+m_{Q_{j}}}A_{0}(q^{2}),
$$

using $\epsilon_{\mu}^{L,*}(p_{V}) = \left(\frac{|\bm{p}_{V}|}{m_{V}},\frac{\bm{p}_{V}}{|\bm{p}_{V}|}\frac{E_{V}}{m_{V}}\right)$ (only the longitudinal polarisation contributes).

**$B\to K^*(892)$**, Eq. `eq:vector-form-factor-1` (L1211–1215, from **Ball:2004rg**):

$$
A^{BK^{*}(892)}_{0}(q^{2}) = \frac{r_{1}}{1-q^{2}/m_{R}^{2}}+\frac{r_{2}}{1-q^{2}/(m_{\text{fit}}^{A_{0}})^{2}} .
$$

**Table `tab:vector-form-factor-parameters1`** (L1229–1240):

| $V$ | $r_{1}$ | $r_{2}$ | $m_{R}$ [GeV] | $m_{\text{fit}}$ [GeV] | $A^{BV}_{0}(0)$ |
|---|---|---|---|---|---|
| $K^{*}(892)$ | $1.364$ | $-0.99$ | $m_{B}^{+}$ ($=5.279$ per the paper's $B^\pm$ mass) | $\sqrt{36.8}\ (=6.0663)$ | $0.374^{+0.033}_{-0.033}$ |

**$B\to K^*(1410),\,K^*(1680)$**, Eq. `eq:vector-form-factor-2` (L1218–1226, from **Lu:2011jm, Hatanaka:2009gb**):

$$
A^{BV}_{0}(q^{2}) = \left( 1-\frac{2m_{V}^{2}}{m_B^2 + m_V^2 - q^2}\right)\xi_{||}(q^{2})+\frac{m_{V}}{m_{B}}\xi_{\perp}(q^{2}),
\qquad
\xi_{\perp/||}(q^{2}) = \frac{\xi_{\perp/||}(0)}{1-q^{2}/m_{B}^{2}} .
$$

**Table `tab:vector-form-factor-parameters2`** (L1242–1255):

| $V$ | $\xi_{\perp}(0)$ | $\xi_{||}(0)$ | $A^{BV}_{0}(0)$ |
|---|---|---|---|
| $K^{*}(1410)$ | $0.28^{+0.04}_{-0.04}$ | $0.22^{+0.03}_{-0.03}$ | $0.3^{+0.036}_{-0.036}$ |
| $K^{*}(1680)$ | $0.24^{+0.05}_{-0.05}$ | $0.18^{+0.03}_{-0.03}$ | $0.22^{+0.04}_{-0.04}$ |

($m_B$ in this pair of equations is not redefined; the natural reading is the same $m_B=5.3$ GeV used in C.2. There is **no explicit $m_R$ or $m_{\text{fit}}$** for $K^*(1410)/K^*(1680)$ — the pole is at $q^2=m_B^2$.)

### C.4 Pseudovector (axial-vector) final state $X'=A$ — `app:pseudovector`

Matrix elements (L1261–1268, from **Bashiry:2009wq, Hatanaka:2008gu**): same as the vector case with $A_i\leftrightarrow V_i$ interchanged,

$$
\langle A(p_{A})| \bar{Q}_{i}\gamma^{\mu}Q_{j} |X(p_{X})\rangle = (m_{X}+m_{A})\epsilon^{\mu *}(p_{A})V_{1}(q^{2}) -
(\epsilon^{*}(p_{A})\cdot q)(p_{X} + p_{A})^{\mu}\frac{V_{2}(q^{2})}{m_{X}+m_{A}}-2m_{A}\frac{\epsilon^{*}(p_{A})\cdot q}{q^{2}}q^{\mu}(V_{3}(q^{2})-V_{0}(q^{2})),
$$

$$
\langle A(p_{A})| \bar{Q}_{i}\gamma^{\mu}\gamma_{5}Q_{j} |X(p_{X})\rangle= \frac{2A(q^{2})}{m_{X}+m_{A}}i\epsilon^{\mu\nu\rho\sigma}\epsilon^{*}_{\nu}(p_{A})p_{X, \rho}p_{A, \sigma},
$$

$$
V_{3}(q^{2}) = \frac{m_{X}+m_{A}}{2m_{A}}V_{1}(q^{2})-\frac{m_{X} - m_{A}}{2m_{A}}V_{2}(q^{2}) \quad\text{(same relation as \ref{eq:a3-a0})},
$$

$$
M_{XA} = \frac{2m_{X}|\bm{p}_{A}|}{m_{Q_{j}}-m_{Q_{i}}}V^{XA}_{0}(q^{2}).
$$

**$K_1(1270)$–$K_1(1400)$ mixing** (L1274–1284), physical states are mixtures of $K_{1A},K_{1B}$:

$$
\begin{pmatrix} |K_{1}(1270) \rangle \\ |K_{1}(1400) \rangle \end{pmatrix} =  \begin{pmatrix} \sin(\theta_{K_{1}}) & \cos(\theta_{K_{1}}) \\ \cos(\theta_{K_{1}}) & -\sin(\theta_{K_{1}})\end{pmatrix}
\begin{pmatrix} |K_{1A} \rangle \\ |K_{1B}\rangle \end{pmatrix},
$$

$$
V^{BK_{1}(1270)}_{0}(q^{2}) =\frac{1}{m_{K_{1}(1270)}}\left[\sin(\theta_{K_{1}})m_{K_{1A}}V_{0}^{A}(q^{2})+\cos(\theta_{K_{1}})m_{K_{1B}}V_{0}^{B}(q^{2})\right],
$$

$$
V^{BK_{1}(1400)}_{0}(q^{2}) =\frac{1}{m_{K_{1}(1400)}}\left[\cos(\theta_{K_{1}})m_{K_{1A}}V_{0}^{A}(q^{2})-\sin(\theta_{K_{1}})m_{K_{1B}}V_{0}^{B}(q^{2})\right],
$$

**Eq. `eq:pseudovector-form-factors-unphys`** (L1286–1288):

$$
V_{0}^{A/B}(q^{2}) = \frac{F^{A/B}_{0}}{1-a_{A/B}\frac{q^{2}}{m_{B}^{2}}+b_{A/B}\left(\frac{q^{2}}{m_{B}^{2}}\right)^{2}} .
$$

**Table `tab:pseudovector-form-factors-parameters1`** (L1292–1302; source **Bashiry:2009wq**):

| $F_{0}^{A}$ | $F_{0}^{B}$ | $a_{A}$ | $a_{B}$ | $b_{A}$ | $b_{B}$ |
|---|---|---|---|---|---|
| $0.22^{+0.04}_{-0.04}$ | $-0.45^{+0.12}_{-0.08}$ | $2.4$ | $1.34$ | $1.78$ | $0.69$ |

**Table `tab:pseudovector-form-factors-parameters2`** (L1304–1314):

| $\theta_{K_1}$ | $m_{K_{1A}}$ | $m_{K_{1B}}$ | $V_{0}^{BK_{1}(1270)}(0)$ | $V_{0}^{BK_{1}(1400)}(0)$ |
|---|---|---|---|---|
| $-34^{\circ}\pm 13^{\circ}$ | $1.31$ GeV | $1.34$ GeV | $-0.52^{+0.13}_{-0.09}$ | $-0.07^{+0.033}_{-0.012}$ |

$m_{K_1(1270)}$ and $m_{K_1(1400)}$ themselves are **not tabulated** here (they appear only as normalisation denominators); the paper's meson list (L677–680) does not include them.

### C.5 Tensor final state $X'=T$ — `app:tensor`

Matrix element (L1322–1325, from **Li:2010ra, Cheng:2010yd**):

$$
\langle T(p_{T})| \bar{Q}_{i}\gamma^{\mu}\gamma_{5}Q_{j} |X(p_{X})\rangle = (m_{X}+m_{T})\epsilon^{\mu *,s}_{T}(p_{T})A_{1}(q^{2}) -
(\epsilon^{*,s}_{T}(p_{T})\cdot q)(p_{X} + p_{T})^{\mu}\frac{A_{2}(q^{2})}{m_{X}+m_{T}}-2m_{T}\frac{\epsilon^{*,s}_{T}(p_{T})\cdot q}{q^{2}}q^{\mu}(A_{3}(q^{2})-A_{0}(q^{2}))
$$

with (L1326–1337)

$$
\epsilon_{T\mu}^{s}(p_{T}) \equiv  \frac{1}{m_{X}}\epsilon_{\mu\nu}^{s}(p_{T})p_{X}^{\nu},
$$

$$
\epsilon^{\pm 2}_{T\mu} = 0, \quad \epsilon_{T\mu}^{\pm 1} = \frac{1}{m_{h}\sqrt{2}}(\epsilon^{0}\cdot p_{X})\epsilon^{\pm 1}_{\mu}, \quad \epsilon^{0}_{T\mu} = \sqrt{\frac{2}{3}}\frac{\epsilon^{0}\cdot p_{X}}{m_{X}}\epsilon_{\mu}^{0},
$$

$$
\epsilon^{\pm 1}_{\mu} = \frac{1}{\sqrt{2}}(0, \mp 1, i, 0), \quad \epsilon_{\mu}^{0} = \frac{1}{m_{T}}(|\bm{p}_{T}|,0,0,E_{T}),
$$

and $\epsilon_{\mu\nu}^{s}$ the polarisation tensor ($p_\mu\epsilon^{\mu\nu,s}=0$, symmetric, traceless). Only $s=0$ contributes to $q_\mu M^{\mu,s}$:

$$
M_{XT}  = -\frac{q_{\mu}M^{\mu, 0}_{XT}}{m_{Q_{i}}+m_{Q_{j}}}= -\frac{1}{m_{Q_{i}}+m_{Q_{j}}}\sqrt{\frac{2}{3}}\frac{m_{X}|\bm{p}_{T}|^{2}}{m_{T}}2 A^{XT}_{0}(q^{2}).
$$

Form-factor parametrisation (L1342–1346, from **Li:2010ra, Cheng:2010yd**):

$$
A^{XT}_{0}(q^{2}) = \frac{F^{XT}_{0}}{\left( 1-\frac{q^{2}}{m_{X}^{2}}\right)\left( 1-a_{T}\frac{q^{2}}{m_{X}^{2}}+b_{T}\left(\frac{q^{2}}{m_{X}^{2}}\right)^{2}\right)} .
$$

For $B\to K_{2}^{*}(1430)$ (**Cheng:2010yd**): $F^{BK_{2}^{*}}_{0} = 0.23$, $a_{T} = 1.23$, $b_{T} = 0.76$ (no uncertainties given). Note the first pole is at $q^2=m_X^2$ and the second uses $m_X$ (the $B$ mass), **not** $m_B=5.3$ GeV — as printed.

### C.6 Two-body width/BR for all of the above

There is **one** formula, Eq. `eq:productionBRmesons` (Sec. B.3), applied per transition with the $M_{XX'}$ of C.1–C.5 and $q^2=m_S^2$. Explicitly:

$$
\Gamma(X_{Q_i}\to X'_{Q_j}S) = \theta^2\,\frac{|\xi_{ij}|^{2}m_{Q_{i}}^{2}|M_{XX'}(m_S^2)|^{2}}{32\pi v^{2}}\frac{|\bm{p}_{S}|}{m_{X}^{2}},
\qquad
\text{BR}=\Gamma/\Gamma_X .
$$

The paper gives **no** separate closed-form width per resonance (no sum over polarisations is written explicitly for $V$, $A$, $T$ beyond the $M_{XX'}$ above).

### C.7 Meson masses and total decay widths used

Masses **explicitly printed** by the paper (L676–681, format "hadron (quark content, mass in MeV)"):

| Hadron | Quark content | Mass |
|---|---|---|
| $K^-$ | $s\bar u$ | 494 MeV |
| $K^0_{S,L}$ | $s\bar d$ | 498 MeV |
| $D^0$ | $c\bar u$ | 1865 MeV |
| $D^+$ | $c\bar d$ | 1870 MeV |
| $D_s$ | $c\bar s$ | 1968 MeV |
| $J/\psi$ | $c\bar c$ | 3097 MeV |
| $B^-$ | $b\bar u$ | 5279 MeV |
| $B^0$ | $b\bar d$ | 5280 MeV |
| $B_s$ | $b\bar s$ | 5367 MeV |
| $B_c$ | $b\bar c$ | 6276 MeV |
| $\Upsilon$ | $b\bar b$ | 9460 MeV |

Other printed masses: $m_B=5.3$ GeV (L1167, "$m_B=5.3$ GeV is the mass of the $B^{+}$ meson", used in C.2/C.3/C.4), $m_{K_{1A}}=1.31$ GeV and $m_{K_{1B}}=1.34$ GeV (C.4), $\Gamma_h\approx4$ MeV (L358), $\overline m_c(2\text{ GeV})=1.23$ GeV and $\overline m_s(2\text{ GeV})=0.0924$ GeV (L1469, scalar decays).

**Total widths / lifetimes $\Gamma_X$ for $B^+,B^0,B_s,K^+,K_L,K_S,D^+,D^0,D_s$: NOT given in this paper.** The paper states only "We use the lifetimes of mesons from [PDG 2018]" (**Tanabashi:2018oca**, L734) and quotes the derived ratio $\Gamma_{B^+}/\Gamma_{B^0}\approx0.93$ (A.6). The pion and kaon masses per se appear only through the hadron list above and through the kinematic closures; the only additional printed numbers are $m_\pi$ implicitly (via the $K\to\pi S$ closing masses 0.350/0.360 GeV) and the scalar-decay thresholds (e.g. $2m_\pi$). **Take $\Gamma_X$ and any precision meson masses from PDG 2018 (Tanabashi:2018oca); they are not reproducible from this source.**

---

## D. Production from mesons through the quartic coupling (App., L1349–1408)

The interaction (Eq. `eq:quartic-coupling`):

$$
\mathcal{L}_{\text{quartic}} = \frac{\alpha}{2}hS^{2}
$$

generates (Eq. `eq:quartic-coupling-processes`):

$$
X_{Q_i} \to X'_{Q_j}SS, \qquad X\to SS .
$$

### D.1 Three-body $X_{Q_i}\to X'_{Q_j}SS$

Matrix element (L1363–1367):

$$
\mathcal{M}(X_{Q_{i}} \to X'_{Q_{j}}SS) \approx \frac{\alpha}{m_{h}^{2}}\frac{m_{Q_i}}{2v}\xi_{ij}M_{XX'}(q^{2}),
$$

where $q^{2}$ is the invariant mass of the scalar pair and $M_{XX'}(q^{2})$ is the same hadronic matrix element as in the mixing case (Eq. `eq:hadronic-matrix-element-1`).

Decay width, **Eq. `eq:br-quartic-3-body`** (L1401–1404):

$$
\Gamma_{X_{Q_{i}} \to X'_{Q_{j}}SS} = \frac{|\xi_{ij}|^{2}m_{Q_{i}}^{2}\alpha^{2}}{512\pi^{3}m_{X}^{3}v^{2}m_{h}^{4}}\int \limits_{4m_{S}^{2}}^{(m_{X}-m_{X'})^{2}} |M_{XX'}(q^{2})|^{2}\sqrt{(E_{2}^{*})^{2}-m_{S}^{2}}\sqrt{(E_{3}^{*})^{2}-m_{X'}^{2}}dq^{2},
$$

with (Eq. L1406–1408)

$$
E_{2}^{*} =\frac{\sqrt{q^{2}}}{2}, \quad E_{3}^{*} = \frac{m_{X}^{2}-q^{2}-m_{X'}^{2}}{2\sqrt{q^{2}}} .
$$

$\alpha$-dependence: $\Gamma\propto\alpha^{2}$ (and $\propto m_h^{-4}$, $v^{-2}$).

### D.2 Two-body $X_{Q_iQ_j}\to SS$

Decay constant definition (L1370–1377):

$$
\langle 0|\bar{Q}_{i}\gamma_{\mu}\gamma_{5}Q_{j} |X(p)\rangle \equiv if_{X}p_{\mu}
\ \Longrightarrow\
\langle 0|\bar{Q}_{i}\gamma_{5}Q_{j} |X(p)\rangle \equiv -\frac{if_{X}m_{X}^{2}}{m_{Q_{i}} - m_{Q_{j}}} .
$$

Matrix element (L1378–1381):

$$
\mathcal{M}(X_{Q_{i}Q_{j}} \to SS)= \frac{m_{Q_{i}}\xi_{ij}}{2v m_{h}^{2}} \langle 0|\bar{Q}_{i}\gamma_{5}Q_{j} |X(p)\rangle \approx i\frac{\alpha f_{X}m_{X}^{2}}{2v m_{h}^{2}}\xi_{ij},
$$

Decay width, **Eq. `eq:bstossbranching`** (L1395–1399):

$$
\Gamma (X_{Q_{i}Q_{j}}\to SS) = \frac{m_{X}^3}{v^2} \frac{|\xi_{ij}|^2 f_{X}^2 \alpha^2}{128\pi m_h^4}\sqrt{1-\frac{4m_{S}^{2}}{m_{X}^{2}}} .
$$

$\alpha$-dependence: $\Gamma\propto\alpha^{2}$.

### D.3 Meson decay constants $f_X$

**Table `tab:meson-decay-constants`** (L1383–1393; source **Chang:2018aut** and references therein):

| Meson $X$ | $B_0$ | $B_s$ | $K_0$ |
|---|---|---|---|
| $f_{X}$ [GeV] | $0.19$ | $0.23$ | $0.16$ |

No uncertainties are quoted. No decay constants are given for $D$, $K^\pm$, etc.

**Printed inconsistency (flagged):** in the matrix element at L1380 the prefactor $\frac{m_{Q_i}\xi_{ij}}{2vm_h^2}$ has **no factor $\alpha$** and no $\frac{m_{Q_i}}{m_{Q_i}-m_{Q_j}}$ (the latter is $\approx1$ for $b\to d,s$), yet the same line's "$\approx$" result and Eq. `eq:bstossbranching` both contain one power of $\alpha$. Since $\Gamma\propto|\mathcal{M}|^2\propto\alpha^2$ as printed, the intended prefactor must carry $\alpha$ (as in the 3-body case, L1365). Read the matrix element as $\mathcal{M}\approx i\,\alpha f_X m_X^2\,\xi_{ij}/(2vm_h^2)$, up to $\mathcal{O}(m_{Q_j}/m_{Q_i})$ corrections.

---

## E. Gaps, ambiguities, and things NOT reproduced in this paper

1. **No $\chi_B,\chi_K,\chi_D$ table.** Only $\chi_{c\bar c},\chi_{b\bar b}$ (Sec. A.3). Per-meson abundances must come from elsewhere; the paper measures heavy-hadron yields via $\chi_{q\bar q}$ and the factor 2 for two hadrons per pair.
2. **No numerical $v$, $m_h$, $G_F$, $m_\mu$, $m_\pi$, $m_c$, $m_b$.** These are symbolic in the formulas. Only $\Gamma_h\approx4$ MeV, $\overline m_c=1.23$ GeV, $\overline m_s=0.0924$ GeV (MS-bar at 2 GeV) and the hadron masses of C.7 are printed. Do not attribute external values to this paper.
3. **No meson total widths $\Gamma_X$ listed.** The BR formula `eq:productionBRmesons` needs them; the paper defers to PDG 2018 (**Tanabashi:2018oca**) and only quotes $\Gamma_{B^+}/\Gamma_{B^0}\approx0.93$.
4. **$B\to K^*(1410)/K^*(1680)$:** the pole mass is $m_B$ (5.3 GeV assumed); no $m_R$, $m_{\text{fit}}$ given. The paper does not restate $m_B$ in that subsection.
5. **Tensor form factor:** the two poles use $m_X$ (the decaying meson), as printed; $F_0,a_T,b_T$ have no uncertainties.
6. **Three-body leptonic formula** (B.5) as printed lacks an explicit $\theta^2$ while the table divides by it. This is a genuine normalisation inconsistency in the source; the text supports $\text{BR}\propto\theta^2$.
7. **$X\to SS$ matrix element** (D.3 above): missing $\alpha$ in the printed prefactor.
8. **Quark-bilinear ordering** differs between `eq:hadronic-matrix-element-1` and `eq:meson-decay-scalar-matrix-element`; sign conventions for $q_\mu$ differ between `eq:trick`, the pseudoscalar subsection, and `eq:hadronic-matrix-element-1`. Only $|M|^2$ and $q^2$ are used in practice.
9. **Parity/notation looseness:** the scalar-meson weak matrix element is called $M^\mu_{X\tilde S}$ although for a parity-odd transition one needs the axial current $M^{\mu5}$; likewise the vector case uses a $\gamma_5$ bilinear obtained from the axial current. Follow the explicit expressions, not the generic labels.
10. **$K^0_S\to S\pi^0$:** assumed to vanish in this paper (CP argument, **Leutwyler:1989xj** not reproduced). The source contains a commented-out estimate $1.2\cdot10^{-5}/\theta^2$ at closing mass 0.360 GeV.
11. **Formulas referenced but not reproduced here:** in every case the paper itself *does* reproduce the formula it uses (with its own parameter fit); the only external results taken over as-is are (i) the raw $B\to K_0^*$ form-factor data of **Sun:2010nv** (refit here), (ii) the CP suppression argument of **Leutwyler:1989xj**, (iii) PDG masses/lifetimes **Tanabashi:2018oca**, (iv) scalar hadronic decay widths of **Winkler:2018qyg** which are given only as a *figure* in this paper (Fig. `fig:Winkler-Br`) and **not** as a formula, and (v) the scalar-nucleon form factor model **Batell:2020vqn, Foroughi-Abari:2021zbm**.
12. **Uncertainties:** the $\xi_{ij}$, scalar decay constants $f_X$, and tensor FF parameters are quoted **without** uncertainties; the pseudoscalar, vector, pseudovector and scalar FF parameters carry uncertainties as tabulated above.

---

## F. Implementation checklist

1. Pick the transition $X\to X'$ and read $\xi_{ij}$ from `tab:xi` (down-type: $\xi_{sb}$ for $B\to K^{(*)}$, $\xi_{ds}$ for $K\to\pi$; up-type: $\xi_{uc}$ for $D$). $v$ symbolic.
2. Evaluate the form factor at $q^2=m_S^2$ with the parameters of Sec. C (use the JSON companion).
3. Build $M_{XX'}(m_S^2)$ from the transition-specific expression in Sec. C (watch the $i$'s — they drop in $|M|^2$).
4. $|\bm p_S|$ from Sec. 0.5.
5. $\Gamma=\theta^2|\xi_{ij}|^2 m_{Q_i}^2|M_{XX'}|^2|\bm p_S|/(32\pi v^2 m_X^2)$; $\text{BR}=\Gamma/\Gamma_X$ with $\Gamma_X$ from PDG 2018 (**Tanabashi:2018oca**).
6. For $B^0$ multiply the $B^\pm$ BR by $0.93$.
7. Quartic channels: Sec. D; multiply by $\alpha^2$ (tables are at $\alpha=1$ GeV).
