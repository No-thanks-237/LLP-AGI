# Two-channel Muskhelishvili–Omnès recipe — verbatim extraction from the LaTeX source

**Source.** `SFF_v2.tex` (1344 lines) in `...\work\widths\references_thesis\arXiv-1204.6251v2\` — the LaTeX
source of M. Hoferichter, C. Ditsche, B. Kubis, U.-G. Meißner, *Dispersive analysis of the scalar form factor of
the nucleon*, JHEP **06** (2012) 063 [arXiv:1204.6251v2].

**Scope.** Every equation below is copied verbatim from the source; only the local macros `\beq` / `\eeq` are
expanded to `\begin{equation}` / `\end{equation}` and indentation is normalised. Nothing is paraphrased.
Prose quotations are marked as such; inside them the math is transliterated to ASCII, but the equations
themselves are exact.

**Equation numbers.** The source uses `\numberwithin{equation}{section}`, so equations are numbered
`<section>.<n>` (appendices: `A.n`, `B.n`). The numbers printed here were recomputed from the source by
counting numbered rows (rows carrying `\nt` = `\notag` are not counted); they are the numbers a reader sees in
the paper. Two anchors fix the mapping: the `\xi`-matrix definition is **(2.43)** and the final Omnès solution
for `\ff(t)` is **(3.8)**. Where a multi-line `align` block carries only *one* printed number, that number is
quoted once for the whole block (the `\label` sits on the block's first row, exactly as in the source).

**Macro dictionary** (source preamble, lines 33–95):

| macro | meaning |
|---|---|
| `\mpi`, `\mK`, `\meta`, `\mpii` | `M_\pi`, `M_K`, `M_\eta`, `M_{\pi^0}` |
| `\tpi`, `\tK`, `\tN`, `\tm` | `t_\pi=4M_\pi^2`, `t_K=4M_K^2`, `t_N=4m^2`, `t_m` = matching point |
| `\ff`, `\FF`, `\GG`, `\Del` | `\mathbf{f}`, `\mathbf{F}`, `\mathbf{G}`, `\boldsymbol{\Delta}` |
| `\unity` | `\mathds{1}` (unit matrix) |
| `\sigma_t^\pi`, `\sigma_t^K`, `\sigma_t^N` | two-body phase-space factors; `q_t=(\sqrt t/2)\sigma_t^\pi` etc., cf. (2.4) |
| `\ec`, `\ep` | `\;` and `\;.` |
| `\nt` | `\notag` |
| `\diff` | `\text{d}` |
| `\dashint` | principal-value integral sign |
| `\Fpi` | `F_\pi` (pion decay constant) |

---

# 1. The two-channel construction — Section 2 (`sec:MOsol_finite_tm`)

Section 2 structure: **2.1** *Formal solution*, **2.2** *Dispersive representation of the Omnès matrix*,
**2.3** *Construction of the Omnès matrix* with **2.3.1** *Infinite matching point*, **2.3.2** *Finite matching
point, single-channel case*, **2.3.3** *Finite matching point, two-channel case*.

## 1.1 Coupled-channel integral equation, unitarity, the T matrix (Section 2, opening)

(2.1) is the inhomogeneous MO equation actually solved: below `t_m` the kernel is the unitarity kernel
`T^*\Sigma`, above `t_m` the imaginary part of `\ff` is *input*. `\Del(t)` contains the left-hand cuts only
and is real for `t\ge t_\pi`.

**Eq. (2.1)** — the generic two-channel Muskhelishvili–Omnès integral equation with finite matching point

```latex
\begin{equation}
\ff(t)=\Del(t)+\frac{1}{\pi}\int\limits_{\tpi}^{\tm}\diff t'\frac{T^*(t')\Sigma(t')\ff(t')}{t'-t}
+\frac{1}{\pi}\int\limits_{\tm}^\infty\diff t'\frac{\Im \ff(t')}{t'-t}\ec
\end{equation}
```

**Eq. (2.2)** — unitarity relation, two-channel form

```latex
\begin{equation}
 \Im \ff(t)=T^*(t)\Sigma(t)\ff(t)\ec
\end{equation}
```

**Eq. (2.3)** — **T matrix and Sigma matrix** in terms of delta, psi, eta, |g| — this is *the* two-channel S-matrix parameterisation used throughout the paper  
`Tmatrix`

```latex
\begin{equation}
\label{Tmatrix}
T(t)=\begin{pmatrix}
      \frac{\eta(t) e^{2i\delta(t)}-1}{2i\sigma_t^\pi q_t^{2J}} & |g(t)|e^{i\psi(t)}\\ 
  |g(t)|e^{i\psi(t)} & \frac{\eta(t) e^{2i(\psi(t)-\delta(t))}-1}{2i\sigma_t^K k_t^{2J}}
\end{pmatrix}\ec\qquad
\Sigma(t)=\text{diag}\Big(\sigma_t^\pi q_t^{2J}\theta\big(t-\tpi\big), \sigma_t^Kk_t^{2J}\theta\big(t-\tK\big)\Big)\ec
\end{equation}
```

**Eq. (2.4)** — t-channel momenta q_t, k_t, p_t and the phase-space factors

```latex
\begin{equation}
q_t=\sqrt{\frac{t}{4}-\mpi^2}=\frac{\sqrt{t}}{2}\sigma_t^\pi\ec\qquad k_t=\sqrt{\frac{t}{4}-\mK^2}=\frac{\sqrt{t}}{2}\sigma_t^K\ec\qquad
p_t=\sqrt{\frac{t}{4}-m^2}=\frac{\sqrt{t}}{2}\sigma_t^N \ec
\end{equation}
```

**Eq. (2.5)** — inelasticity eta expressed through |g| (only delta, psi, |g| are independent input)

```latex
\begin{equation}
\eta(t)=\sqrt{1-4\sigma_t^\pi\sigma_t^K(q_tk_t)^{2J}|g(t)|^2\theta\big(t-\tK\big)}\ec
\end{equation}
```

**Eq. (2.6)** — S = 1 + 2i Sigma^{1/2} T Sigma^{1/2} (S-matrix in terms of T and the phase-space factors)

```latex
\begin{equation}
S(t)=\unity+2i\,\Sigma^{1/2}(t)T(t)\Sigma^{1/2}(t)\ep
\end{equation}
```

> *Verbatim (source lines 184–186), immediately after (2.5):* "The scattering phases `\delta(t)` and `\psi(t)`
> are required as input for `\tpi\leq t\leq\tm`, where Watson's theorem [14] demands `\psi(t)=\delta(t)` for
> `t\leq \tK`. Moreover, the modulus of the `\pi\pi\to\bar K K` `S`-wave `g(t)` is needed in the full
> range `\tpi\leq t\leq\tm`, and thus has to be analytically continued into the pseudophysical region
> `\tpi\leq t\leq \tK`."

## 1.2 Formal solution: FF = Omega G (Section 2.1)

The Omnès matrix is *defined* by (2.7) — that is the entire content of the problem: a matrix with the right
discontinuity on [t_pi, t_m] and real elsewhere, normalised Omega(0) = 1. (2.8)–(2.11) derive the jump of G, and
(2.12) is the resulting representation.

**Eq. (2.7)** — **defining property of the Omnès matrix**: Im Omega = T* Sigma Omega on [t_pi, t_m], Im Omega = 0 elsewhere; normalisation Omega(0) = 1 (stated in the text right after the equation)  
`def_finite`

```latex
\begin{equation}
\label{def_finite}
\begin{Bmatrix}\Im \Omega(t)=T^*(t)\Sigma(t)\Omega(t) \\ \Im \Omega(t)=0\end{Bmatrix}
\qquad \text{for} \qquad
\begin{Bmatrix}\tpi\leq t \leq \tm \\ \text{otherwise}\end{Bmatrix}
\end{equation}
```

**Eq. (2.8)** — FF = f - Delta = Omega G (definition of the auxiliary vector G)

```latex
\begin{equation}
\FF(t)=\ff(t)-\Del(t)=\Omega(t)\GG(t)\ec
\end{equation}
```

**Eq. (2.9)** — jump condition for G

```latex
\begin{equation}
\big(\unity-2i\, T^*(t)\Sigma(t)\big)\Omega(t_+)\big(\GG(t_+)-\GG(t_-)\big)=2i\,T^*(t)\Sigma(t)\Del(t)\ec
\end{equation}
```

**Eq. (2.10)** — unitarity in the form (1 - 2i T* Sigma)^{-1} = 1 + 2i T Sigma

```latex
\begin{equation}
\big(\unity-2i\, T^*(t)\Sigma(t)\big)^{-1}=\unity+2i\,T(t)\Sigma(t)\ec
\end{equation}
```

**Eq. (2.11)** — G(t+) - G(t-) = 2i Omega^{-1} T Sigma Delta

```latex
\begin{equation}
\GG(t_+)-\GG(t_-)=2i\,\Omega^{-1}(t)T(t)\Sigma(t)\Del(t)\ec
\end{equation}
```

**Eq. (2.12)** — **the formal MO solution**; the second line needs only Im Omega^{-1} on [t_pi, t_m]  
`omnes_fin`

```latex
\begin{align}
\label{omnes_fin}
\ff(t)&=\Del(t)+\frac{\Omega(t)}{\pi}\int\limits_{\tpi}^{\tm}\diff t'\frac{\Omega^{-1}(t')T(t')\Sigma(t')\Del(t')}{t'-t}+
\frac{\Omega(t)}{\pi}\int\limits_{\tm}^\infty\diff t'\frac{\Omega^{-1}(t')\Im \ff(t')}{t'-t}\nt\\
&=\Del(t)-\frac{\Omega(t)}{\pi}\int\limits_{\tpi}^{\tm}\diff t'\frac{\Im\Omega^{-1}(t')\Del(t')}{t'-t}+
\frac{\Omega(t)}{\pi}\int\limits_{\tm}^\infty\diff t'\frac{\Omega^{-1}(t')\Im \ff(t')}{t'-t}\ep
\end{align}
```

> *Note on notation:* in the Omnès literature one writes F = Omega . P. Here P is the vector `\GG` of (2.8)
> and the representation is (2.12); the polynomial/input part is carried explicitly by `\Del(t)`.

## 1.3 The determinant of Omega and the cusp at t_m (Sections 2.1, 2.2)

(2.13)–(2.14) give the discontinuity of the determinant, (2.15) its **closed-form solution** (the analytic result
that any numerical Omega must reproduce), and (2.16)–(2.18) the behaviour at the matching point, which defines the
exponents x and x_ij.

**Eq. (2.13)** — Omega(t+) = (1 + 2i T Sigma) Omega(t-) for t_pi <= t <= t_K

```latex
\begin{equation}
\Omega(t_+)=\big(\unity+2i\,T(t)\Sigma(t)\big)\Omega(t_-)\ep
\end{equation}
```

**Eq. (2.14)** — det Omega(t+) = exp(2i psi) det Omega(t-)

```latex
\begin{equation}
\det \Omega(t_+)=e^{2i\psi(t)}\det\Omega(t_-)\ec
\end{equation}
```

**Eq. (2.15)** — **analytic determinant** of the Omnès matrix (Moussallam's formula) — the validation target for a numerical Omega  
`det_fin`

```latex
\begin{equation}
\label{det_fin}
\det \Omega(t)=\exp\Bigg\{\frac{t}{\pi}\int\limits_{\tpi}^{\tm}\diff t'\frac{\psi(t')}{t'(t'-t)}\Bigg\}\ep
\end{equation}
```

**Eq. (2.16)** — cusp behaviour det Omega ~ |t_m - t|^x with x = psi(t_m)/pi  
`det_asym`

```latex
\begin{equation}
\label{det_asym}
\det\Omega(t)\sim |\tm-t|^x\ec\qquad x=\frac{\psi(\tm)}{\pi}\ep
\end{equation}
```

**Eq. (2.17)** — behaviour for t -> t_m **from below**: the powers x_ij and the phases delta_ij(t_m)  
`def_below`

```latex
\begin{equation}
\label{def_below}
 \det\Omega(t)=\det\bar\Omega(\tm)e^{i\pi x} |\tm-t|^x\ec\qquad \Omega_{ij}(t)=\bar\Omega_{ij}(\tm)e^{i\delta_{ij}(\tm)}|\tm-t|^{x_{ij}}\ec
\end{equation}
```

**Eq. (2.18)** — same **from above**, where Omega is real  
`def_above`

```latex
\begin{equation}
\label{def_above}
 \det\Omega(t)=\det\bar\Omega(\tm) |\tm-t|^x\ec\qquad \Omega_{ij}(t)=\bar\Omega_{ij}(\tm)|\tm-t|^{x_{ij}}\ec
\end{equation}
```

## 1.4 Dispersive representation with finite t_m and the powers x_ij (Section 2.2)

(2.19) is the **subtracted dispersive representation of Omega** with Omega(0)=1 and Omega(t_m)=0: this is the
equation a discretised solver has to solve. (2.20) are the two asymptotic integrals (principal value) that make the
cusp integrable; (2.21)–(2.22) are the resulting component-wise matching conditions.

**Eq. (2.19)** — **subtracted dispersive representation of Omega(t)** with Omega(0)=1, Omega(t_m)=0; note the kernel T* Sigma Omega

```latex
\begin{equation}
\Omega(t)=\Big(1-\frac{t}{\tm}\Big)\unity+\frac{t(t-\tm)}{\pi}\int\limits_{\tpi}^{\tm}\diff t'\frac{T^*(t')\Sigma(t')\Omega(t')}{t'(t'-\tm)(t'-t)}\ec
\end{equation}
```

**Eq. (2.20)** — the two asymptotic integrals (PV = principal value) giving -cot(pi x) and 1/sin(pi x)  
`integrals`

```latex
\begin{align}
\label{integrals}
&\frac{|\tm-t|^{x}}{\pi}\dashint[9pt]\limits_{\tm-\eps}^{\tm}\frac{\diff t'}{(t'-t)|\tm-t'|^x}
\overset{\eps\to 0,\,t\to\tm}{\longrightarrow}\frac{1}{\pi}\dashint[0.5pt]\limits_0^\infty\frac{\diff v}{v^x(1-v)}=-\cot\pi x\ec\nt\\

&\frac{|\tm-t|^{x}}{\pi}\int\limits_{\tm}^{\tm+\eps}\frac{\diff t'}{(t'-t)|\tm-t'|^x}
\overset{\eps\to 0,\,t\to\tm}{\longrightarrow}\frac{1}{\pi}\int\limits_0^\infty\frac{\diff v}{v^x(1+v)}=\frac{1}{\sin \pi x}\ec
\end{align}
```

**Eq. (2.21)** — component-wise matching conditions for t -> t_m **from below** (all four equations)  
`Omegaij_above`

```latex
\begin{align}
\label{Omegaij_above}
 \bar\Omega_{11}|\tm-t|^{x_{11}}e^{i\delta_{11}}&=\bar\Omega_{11}|\tm-t|^{x_{11}}e^{i\delta_{11}}\frac{e^{i\pi x_{11}}}{\sin\pi x_{11}}\frac{1-\eta e^{-2i\delta}}{2i}+\bar\Omega_{21}|\tm-t|^{x_{21}}e^{i\delta_{21}}\frac{e^{i\pi x_{21}}}{\sin\pi x_{21}}|g|\sigma_{\tm}^Ke^{-i\psi}\ec\nt\\
 \bar\Omega_{12}|\tm-t|^{x_{12}}e^{i\delta_{12}}&=\bar\Omega_{12}|\tm-t|^{x_{12}}e^{i\delta_{12}}\frac{e^{i\pi x_{12}}}{\sin\pi x_{12}}\frac{1-\eta e^{-2i\delta}}{2i}+\bar\Omega_{22}|\tm-t|^{x_{22}}e^{i\delta_{22}}\frac{e^{i\pi x_{22}}}{\sin\pi x_{22}}|g|\sigma_{\tm}^Ke^{-i\psi}\ec\nt\\

\bar\Omega_{21}|\tm-t|^{x_{21}}e^{i\delta_{21}}&=\bar\Omega_{11}|\tm-t|^{x_{11}}e^{i\delta_{11}}\frac{e^{i\pi x_{11}}}{\sin\pi x_{11}}|g|\sigma_{\tm}^\pi e^{-i\psi}\nt\\
&\qquad+\bar\Omega_{21}|\tm-t|^{x_{21}}e^{i\delta_{21}}\frac{e^{i\pi x_{21}}}{\sin\pi x_{21}}\frac{1-\eta e^{-2i(\psi-\delta)}}{2i}\ec\nt\\
\bar\Omega_{22}|\tm-t|^{x_{22}}e^{i\delta_{22}}&=\bar\Omega_{12}|\tm-t|^{x_{12}}e^{i\delta_{12}}\frac{e^{i\pi x_{12}}}{\sin\pi x_{12}}|g|\sigma_{\tm}^\pi e^{-i\psi}\nt\\
&\qquad+\bar\Omega_{22}|\tm-t|^{x_{22}}e^{i\delta_{22}}\frac{e^{i\pi x_{22}}}{\sin\pi x_{22}}\frac{1-\eta e^{-2i(\psi-\delta)}}{2i}\ec
\end{align}
```

**Eq. (2.22)** — component-wise matching conditions for t -> t_m **from above**  
`Omegaij_below`

```latex
\begin{align}
\label{Omegaij_below}
 \bar\Omega_{11}|\tm-t|^{x_{11}}&=\bar\Omega_{11}|\tm-t|^{x_{11}}\frac{e^{i\delta_{11}}}{\sin\pi x_{11}}\frac{1-\eta e^{-2i\delta}}{2i}+\bar\Omega_{21}|\tm-t|^{x_{21}}\frac{e^{i\delta_{21}}}{\sin\pi x_{21}}|g|\sigma_{\tm}^K e^{-i\psi}\ec\nt\\
 \bar\Omega_{12}|\tm-t|^{x_{12}}&=\bar\Omega_{12}|\tm-t|^{x_{12}}\frac{e^{i\delta_{12}}}{\sin\pi x_{12}}\frac{1-\eta e^{-2i\delta}}{2i}+\bar\Omega_{22}|\tm-t|^{x_{22}}\frac{e^{i\delta_{22}}}{\sin\pi x_{22}}|g|\sigma_{\tm}^K e^{-i\psi}\ec\nt\\

\bar\Omega_{21}|\tm-t|^{x_{21}}&=\bar\Omega_{11}|\tm-t|^{x_{11}}\frac{e^{i\delta_{11}}}{\sin\pi x_{11}}|g|\sigma_{\tm}^\pi e^{-i\psi}
+\bar\Omega_{21}|\tm-t|^{x_{21}}\frac{e^{i\delta_{21}}}{\sin\pi x_{21}}\frac{1-\eta e^{-2i(\psi-\delta)}}{2i}\ec\nt\\
\bar\Omega_{22}|\tm-t|^{x_{22}}&=\bar\Omega_{12}|\tm-t|^{x_{12}}\frac{e^{i\delta_{12}}}{\sin\pi x_{12}}|g|\sigma_{\tm}^\pi e^{-i\psi}
+\bar\Omega_{22}|\tm-t|^{x_{22}}\frac{e^{i\delta_{22}}}{\sin\pi x_{22}}\frac{1-\eta e^{-2i(\psi-\delta)}}{2i}\ec
\end{align}
```

> *Verbatim (source line 299):* "If we assume that `g(\tm)\neq 0` and `\bar \Omega_{ij}(\tm)\neq 0` ... we can
> conclude from the first line of (2.21) that `x_{21}\geq x_{11}` ... Conversely, the third line requires
> `x_{21}\leq x_{11}` by the same argument, and hence `x_{11}=x_{21}`. Similarly, we find `x_{12}=x_{22}`.
> Moreover, as the determinant behaves according to (2.16), we can conclude that `x_{11}+x_{12}=x` ..."

**Eq. (2.23)** — consistency relation from dividing the first lines of (2.21) and (2.22)

```latex
\begin{equation}
e^{i\delta_{11}}\bigg(1-\frac{e^{i\pi x_{11}}}{\sin\pi x_{11}}\frac{1-\eta e^{-2i\delta}}{2i}\bigg)
=e^{i\pi x_{21}}\bigg(1-\frac{e^{i\delta_{11}}}{\sin\pi x_{11}}\frac{1-\eta e^{-2i\delta}}{2i}\bigg)\ec
\end{equation}
```

**Eq. (2.24)** — exp(i delta_11) = exp(i pi x_11)

```latex
\begin{equation}
e^{i\delta_{11}}=e^{i\pi x_{11}}\ec
\end{equation}
```

**Eq. (2.25)** — **the exponents**: x_11 = x_21, x_12 = x_22, x_ij = delta_ij/pi, x_11 + x_12 = x  
`xij`

```latex
\begin{equation}
\label{xij}
x_{11}=x_{21}\ec\qquad x_{12}=x_{22}\ec\qquad x_{ij}=\frac{\delta_{ij}}{\pi}\ec \qquad x_{11}+x_{12}=x\ep
\end{equation}
```

**Eq. (2.26)** — the two 2x2 homogeneous systems for (Omega-bar_11, Omega-bar_21) and (Omega-bar_12, Omega-bar_22) at t_m  
`spectral_rep_const`

```latex
\begin{align}
\label{spectral_rep_const}
\renewcommand{\arraystretch}{1.3}
 \begin{pmatrix}
 \frac{e^{i\pi x_{11}}}{\sin\pi x_{11}}\frac{1-\eta e^{-2i\delta}}{2i}-1 & \frac{e^{i\pi x_{11}}}{\sin\pi x_{11}}|g|\sigma_{\tm}^K e^{-i\psi}\\
\frac{e^{i\pi x_{11}}}{\sin\pi x_{11}}|g|\sigma_{\tm}^\pi e^{-i\psi} & \frac{e^{i\pi x_{11}}}{\sin\pi x_{11}}\frac{1-\eta e^{-2i(\psi-\delta)}}{2i}-1
\end{pmatrix}
\begin{pmatrix}
 \bar\Omega_{11}\\ 
\bar\Omega_{21}
\end{pmatrix}
&=0\ec\nt\\

\renewcommand{\arraystretch}{1.3}
 \begin{pmatrix}
 \frac{e^{i\pi x_{12}}}{\sin\pi x_{12}}\frac{1-\eta e^{-2i\delta}}{2i}-1 & \frac{e^{i\pi x_{12}}}{\sin\pi x_{12}}|g|\sigma_{\tm}^K e^{-i\psi}\\
\frac{e^{i\pi x_{12}}}{\sin\pi x_{12}}|g|\sigma_{\tm}^\pi e^{-i\psi} & \frac{e^{i\pi x_{12}}}{\sin\pi x_{12}}\frac{1-\eta e^{-2i(\psi-\delta)}}{2i}-1
\end{pmatrix}
\begin{pmatrix}
 \bar\Omega_{12}\\ 
\bar\Omega_{22}
\end{pmatrix}
&=0\ep
\end{align}
```

**Eq. (2.27)** — determinant conditions: cos pi(2x_11 - x) - eta cos pi(2y - x) = 0, y = delta(t_m)/pi

```latex
\begin{equation}
\cos \pi (2x_{11}-x)-\eta\cos\pi (2y-x)=0\ec\qquad \cos \pi (2x_{12}-x)-\eta\cos\pi (2y-x)=0\ec\qquad y=\frac{\delta(\tm)}{\pi}\ep
\end{equation}
```

**Eq. (2.28)** — **closed-form solution for the exponents** x_11, x_12 in terms of x, eta, y  
`x11result`

```latex
\begin{equation}
\label{x11result}
x_{11}=\frac{1}{2}\Big\{x+\frac{1}{\pi}\arccos(\eta\cos\pi(2y-x))\Big\}\ec\qquad x_{12}=\frac{1}{2}\Big\{x-\frac{1}{\pi}\arccos(\eta\cos\pi(2y-x))\Big\}\ep
\end{equation}
```

**Eq. (2.29)** — sin and cos of pi x_11, pi x_12 in terms of z = eta cos pi(2y - x)  
`sinx11`

```latex
\begin{align}
\label{sinx11}
 \sin\pi x_{11}&=\sin\frac{\pi x}{2}\sqrt{\frac{1+z}{2}}+\cos\frac{\pi x}{2}\sqrt{\frac{1-z}{2}}\ec\qquad
 \cos\pi x_{11}=\cos\frac{\pi x}{2}\sqrt{\frac{1+z}{2}}-\sin\frac{\pi x}{2}\sqrt{\frac{1-z}{2}}\ec\nt\\
\sin\pi x_{12}&=\sin\frac{\pi x}{2}\sqrt{\frac{1+z}{2}}-\cos\frac{\pi x}{2}\sqrt{\frac{1-z}{2}}\ec\qquad
 \cos\pi x_{12}=\cos\frac{\pi x}{2}\sqrt{\frac{1+z}{2}}+\sin\frac{\pi x}{2}\sqrt{\frac{1-z}{2}}\ec\nt\\
z&=\eta\cos\pi(2y-x)\ec
\end{align}
```

**Eq. (2.30)** — **the ratios Omega-bar_21/Omega-bar_11 and Omega-bar_12/Omega-bar_22 at t_m**, with N defined here (used everywhere below)  
`ratio_Omij`

```latex
\begin{equation}
\label{ratio_Omij}
 \frac{\bar\Omega_{21}}{\bar\Omega_{11}}=\frac{N}{2|g|\sigma_{\tm}^K}\ec\qquad
 \frac{\bar\Omega_{12}}{\bar\Omega_{22}}=-\frac{N}{2|g|\sigma_{\tm}^\pi}\ec\qquad N=\sqrt{1-\eta^2\cos^2\pi(2y-x)}-\eta\sin\pi(2y-x)\ep
\end{equation}
```

## 1.5 The infinite-matching-point solution Omega^inf (Section 2.3.1)

Omega^inf satisfies the same unitarity property on [t_pi, inf) (2.31) and is computed from the **unsubtracted**
dispersion relation (2.32), discretised on Gauss–Legendre points. (2.35) rewrites unitarity as a real linear system
for the real part — that is the equation actually solved numerically; (2.37) is its simplification below t_K.

**Eq. (2.31)** — **defining property of Omega^inf(t)** (infinite matching point)  
`Omnes_def`

```latex
\begin{equation}
\label{Omnes_def}
\begin{Bmatrix}\Im \Omega^\infty(t)=T^*(t)\Sigma(t)\Omega^\infty(t) \\ \Im \Omega^\infty(t)=0\end{Bmatrix}
\qquad \text{for} \qquad
\begin{Bmatrix}t\geq \tpi \\ \text{otherwise}\end{Bmatrix}\ep
\end{equation}
```

**Eq. (2.32)** — **unsubtracted dispersion relation for Re Omega^inf** (the equation discretised on Gauss–Legendre points); dash = principal value

```latex
\begin{equation}
\Re \Omega^\infty(t)=\frac{1}{\pi}\dashint[0.5pt]\limits_{\tpi}^{\infty}\diff t'\frac{\Im \Omega^\infty(t')}{t'-t}
\end{equation}
```

**Eq. (2.33)** — closed-form single-channel Omega^inf (limiting case / sanity check)  
`sol_single_infinite`

```latex
\begin{equation}
\label{sol_single_infinite}
\Omega^\infty(t)=\exp\Bigg\{\frac{t}{\pi}\int\limits_{\tpi}^\infty\diff t'\frac{\delta(t')}{t'(t'-t)}\Bigg\}\ec
\end{equation}
```

**Eq. (2.34)** — asymptotic behaviour of Omega^inf for large t, with exponent delta(inf)/pi

```latex
\begin{equation}
\Omega^\infty(t)\sim t^{-\frac{\delta(\infty)}{\pi}}
\end{equation}
```

**Eq. (2.35)** — **unitarity rewritten as a linear equation for Re Omega^inf_i** — the 2x2 matrix multiplying the real part  
`unitarity`

```latex
\begin{equation}
\label{unitarity}
\renewcommand{\arraystretch}{1.3}
\Im\mathbf{\Omega}^\infty_i=\begin{pmatrix}
                    \frac{\eta\sin(2\delta-\psi)+\sin\psi}{\eta\cos(2\delta-\psi)+\cos\psi} & \frac{2|g|\sigma_t^K\theta(t-\tK)}{\eta\cos(2\delta-\psi)+\cos\psi}\\
\frac{2|g|\sigma_t^\pi}{\eta\cos(2\delta-\psi)+\cos\psi} & -\frac{\eta \sin(2\delta-\psi)-\sin\psi}{\eta\cos(2\delta-\psi)+\cos\psi}
                   \end{pmatrix}
\Re\mathbf{\Omega}^\infty_i\ec\qquad i\in\{1,2\}\ec
\end{equation}
```

**Eq. (2.36)** — definition of the column vectors Omega^inf_1, Omega^inf_2  
`Omnes_vector`

```latex
\begin{equation}
\label{Omnes_vector}
\mathbf{\Omega}_1^\infty=\begin{pmatrix}
                                          \Omega^\infty_{11}\\ \Omega^\infty_{21}
                                         \end{pmatrix}\ec \qquad
\mathbf{\Omega}_2^\infty=\begin{pmatrix}
                                          \Omega^\infty_{12}\\ \Omega^\infty_{22}
                                         \end{pmatrix}\ec
\end{equation}
```

**Eq. (2.37)** — the same relation below t_K (Watson regime): tan delta and |g| sigma^pi_t / cos delta  
`Omnes_below_tK`

```latex
\begin{equation}
\label{Omnes_below_tK}
\Im\mathbf{\Omega}^\infty_i=\begin{pmatrix}
                    \tan\delta & 0\\
\frac{|g|\sigma_t^\pi}{\cos\delta} & 0
                   \end{pmatrix}
\Re\mathbf{\Omega}^\infty_i\ep
\end{equation}
```

## 1.6 The xi matrix and the finite-t_m two-channel Omnès matrix (Sections 2.3.2–2.3.3)

The construction is Omega(t) = Omega^inf(t) xi(t): single channel (2.38)–(2.42), two channels (2.43)–(2.56). xi
supplies the factors ((t_m - t)/t_m)^{x_ij} and is real above t_m, so it cancels the imaginary part there while
preserving unitarity below. In the two-channel case a real matrix a(t) remains free; requiring continuity of
Omega-bar_ij at t_m and the ratios (2.30) fixes it, giving the explicit a_1(t), a_2(t) of (2.53) and the
"join two vectors" prescription below.

**Eq. (2.38)** — the single-channel xi(t) = ((t_m - t)/t_m)^{x(t)} and its phase convention

```latex
\begin{equation}
\xi(t)=\bigg(\frac{\tm-t}{\tm}\bigg)^{x(t)}=\bigg|\frac{\tm-t}{\tm}\bigg|^{x(t)}e^{-i\pi x(t)\theta(t-\tm)}
\end{equation}
```

**Eq. (2.39)** — Omega(t) = Omega^inf(t) xi(t)  
`Omnes_constr`

```latex
\begin{equation}
\label{Omnes_constr}
\Omega(t)=\Omega^\infty(t)\xi(t)\ec
\end{equation}
```

**Eq. (2.40)** — the two defining properties of this Omega (unitarity below t_m; explicit imaginary part above)

```latex
\begin{equation}
\begin{Bmatrix}\Im \Omega(t)=T^*(t)\Sigma(t)\Omega(t) \\ 
\Im \Omega(t)=\Omega^\infty(t)|\xi(t)|\Big(e^{i(\pi x(t)-\delta(t))}\sin\delta(t)-\sin\pi x(t)\Big)\end{Bmatrix}
\qquad \text{for} \qquad
\begin{Bmatrix}\tpi\leq t \leq \tm \\ t\geq\tm\end{Bmatrix}\ec
\end{equation}
```

**Eq. (2.41)** — the choice of x(t): delta(t_m)/pi for t <= t_m, delta(t)/pi above

```latex
\begin{equation}
\begin{Bmatrix}x(t)=\frac{\delta(\tm)}{\pi} \\ 
x(t)=\frac{\delta(t)}{\pi}\end{Bmatrix}
\qquad \text{for} \qquad
\begin{Bmatrix}t \leq \tm \\ t\geq\tm\end{Bmatrix}
\end{equation}
```

**Eq. (2.42)** — explicit single-channel solution with finite t_m (two exponential factors)  
`omnes_sol`

```latex
\begin{equation}
\label{omnes_sol}
\Omega(t)=\exp\Bigg\{\frac{t}{\pi}\int\limits_{\tpi}^{\tm}\diff t'\frac{\delta(t')}{t'(t'-t)}\Bigg\}
\exp\Bigg\{\frac{t}{\pi}\int\limits_{\tm}^{\infty}\diff t'\frac{\delta(t')-\pi x(t)}{t'(t'-t)}\Bigg\}\ep
\end{equation}
```

**The xi matrix itself is (2.43):**

**Eq. (2.43)** — **the xi matrix** xi_ij(t) = xi-bar_ij(t) ((t_m - t)/t_m)^{x_ij(t)}, with real functions xi-bar_ij

```latex
\begin{equation}
\xi_{ij}(t)=\bar \xi_{ij}(t)  \bigg(\frac{\tm-t}{\tm}\bigg)^{x_{ij}(t)}
\end{equation}
```

**Eq. (2.44)** — the two-branch ansatz: a Omega^inf xi below t_m and a (Omega^inf)^T xi above t_m  
`sol_fin_Omega_til`

```latex
\begin{equation}
\label{sol_fin_Omega_til}
\begin{Bmatrix}\Omega(t)=a(t)\Omega^\infty(t) \xi(t) \\ 
\Omega(t)=a(t)(\Omega^\infty)^T(t)\xi(t)\end{Bmatrix}
\qquad \text{for} \qquad
\begin{Bmatrix}t \leq \tm \\ t\geq\tm\end{Bmatrix}\ec
\end{equation}
```

**Eq. (2.45)** — the two conditions on xi and on a(t) (commutator condition for a)  
`cond_xi`

```latex
\begin{align}
\label{cond_xi}
\Im\xi^T(t)+\xi^T(t) T(t)\Sigma(t)&=0\qquad \text{for}\ t\geq\tm\ec\nt\\
\big[a(t),T^*(t)\Sigma(t)\big]&=0\qquad \text{for}\ t\leq\tm\ep
\end{align}
```

**Eq. (2.46)** — **the 2x2 systems that determine xi-bar_11, xi-bar_21 and xi-bar_12, xi-bar_22**  
`eq_xi`

```latex
\begin{align}
\label{eq_xi}
\renewcommand{\arraystretch}{1.3}
 \begin{pmatrix}
 \frac{e^{i\pi x_{11}}}{\sin\pi x_{11}}\frac{1-\eta e^{-2i\delta}}{2i}-1 & \frac{e^{i\pi x_{11}}}{\sin\pi x_{11}}|g|\sigma_{\tm}^\pi e^{-i\psi}\\
\frac{e^{i\pi x_{11}}}{\sin\pi x_{11}}|g|\sigma_{\tm}^K e^{-i\psi} & \frac{e^{i\pi x_{11}}}{\sin\pi x_{11}}\frac{1-\eta e^{-2i(\psi-\delta)}}{2i}-1
\end{pmatrix}
\begin{pmatrix}
 \bar\xi_{11}\\ 
\bar\xi_{21}
\end{pmatrix}
&=0\ec\nt\\

\renewcommand{\arraystretch}{1.3}
 \begin{pmatrix}
 \frac{e^{i\pi x_{12}}}{\sin\pi x_{12}}\frac{1-\eta e^{-2i\delta}}{2i}-1 & \frac{e^{i\pi x_{12}}}{\sin\pi x_{12}}|g|\sigma_{\tm}^\pi e^{-i\psi}\\
\frac{e^{i\pi x_{12}}}{\sin\pi x_{12}}|g|\sigma_{\tm}^K e^{-i\psi} & \frac{e^{i\pi x_{12}}}{\sin\pi x_{12}}\frac{1-\eta e^{-2i(\psi-\delta)}}{2i}-1
\end{pmatrix}
\begin{pmatrix}
 \bar\xi_{12}\\ 
\bar\xi_{22}
\end{pmatrix}
&=0\ep
\end{align}
```

**Eq. (2.47)** — ratios xi-bar_21/xi-bar_11 and xi-bar_12/xi-bar_22

```latex
\begin{equation}
 \frac{\bar\xi_{21}}{\bar\xi_{11}}=\frac{N}{2|g|\sigma_{\tm}^\pi}\ec\qquad \frac{\bar\xi_{12}}{\bar\xi_{22}}=-\frac{N}{2|g|\sigma_{\tm}^K}\ep
\end{equation}
```

**Eq. (2.48)** — the normalisation choice for xi that makes det xi equal to the single-channel factor

```latex
\begin{equation}
\bar\xi_{11}=1\ec\qquad \bar\xi_{22}=\bigg(1+\frac{N^2}{1-\eta^2}\bigg)^{-1}\ec
\end{equation}
```

**Eq. (2.49)** — det xi(t) = ((t_m - t)/t_m)^x

```latex
\begin{equation}
\det\xi(t)=\bigg(\frac{\tm-t}{\tm}\bigg)^x\ec
\end{equation}
```

**Eq. (2.50)** — the constraints on a(t): a_21 = a_12 sigma^pi_t/sigma^K_t, a_22 = a_11 - a_12 eta sin(2delta-psi)/(|g| sigma^K_t)

```latex
\begin{equation}
a_{21}=a_{12}\frac{\sigma_t^\pi}{\sigma_t^K}\ec\qquad a_{22}=a_{11}-a_{12}\frac{\eta\sin(2\delta-\psi)}{|g|\sigma_t^K}\ec
\end{equation}
```

**Eq. (2.51)** — normalisation condition Omega^inf(0) xi(0) = 1 written out

```latex
\begin{equation}
\Omega^\infty(0)\xi(0)=\begin{pmatrix}
                 \bar\xi_{11}\big(\Omega^\infty_{11}(0)+\frac{N}{2|g|\sigma_{\tm}^\pi}\Omega^\infty_{12}(0)\big) & \bar\xi_{22}\big(\Omega^\infty_{12}(0)-\frac{N}{2|g|\sigma_{\tm}^K}\Omega^\infty_{11}(0)\big)\\
		  \bar\xi_{11}\big(\Omega^\infty_{21}(0)+\frac{N}{2|g|\sigma_{\tm}^\pi}\Omega^\infty_{22}(0)\big) & \bar\xi_{22}\big(\Omega^\infty_{22}(0)-\frac{N}{2|g|\sigma_{\tm}^K}\Omega^\infty_{21}(0)\big)
                \end{pmatrix}
=\mathds{1}
\end{equation}
```

**Eq. (2.52)** — **the required normalisation of Omega^inf(0)** that enforces Omega(0) = 1

```latex
\begin{equation}
\Omega^\infty(0)=\bigg(1+\frac{N^2}{1-\eta^2}\bigg)^{-1}\begin{pmatrix}
                 {\bar\xi_{11}}^{-1} & \bar\xi_{11}^{-1}\frac{N}{2|g|\sigma_{\tm}^K}\\
		  -\bar\xi_{22}^{-1}\frac{N}{2|g|\sigma_{\tm}^\pi} & \bar\xi_{22}^{-1}
                \end{pmatrix}\ep
\end{equation}
```

**Eq. (2.53)** — **the explicit real matrices a_1(t), a_2(t)** used for t >= t_K (with f(t_K)=0, f(t_m)=1); below t_K, a_1 = a_2 = 1

```latex
\begin{align}
a_1(t)&=\begin{pmatrix}
        1 & \frac{N(t)}{2|g(t)|\sigma_t^\pi}f(t)\\
\frac{N(t)}{2|g(t)|\sigma_t^K}f(t) & 1-f(t)+\frac{N^2(t)}{1-\eta^2(t)}f(t)
       \end{pmatrix}\ec\nt\\
a_2(t)&=\begin{pmatrix}
        1 & -\frac{2|g(t)|\sigma_t^K}{N(t)}f(t)\\
-\frac{2|g(t)|\sigma_t^\pi}{N(t)}f(t) & 1-f(t)+\frac{1-\eta^2(t)}{N^2(t)}f(t)
       \end{pmatrix}\ec
\end{align}
```

**Eq. (2.54)** — N(t) = sqrt(1 - eta^2 cos^2(2delta - psi)) - eta sin(2delta - psi)

```latex
\begin{equation}
N(t)=\sqrt{1-\eta(t)^2\cos^2(2\delta(t)-\psi(t))}-\eta(t)\sin(2\delta(t)-\psi(t))\ec
\end{equation}
```

**Eq. (2.55)** — det Omega in terms of the real parts

```latex
\begin{equation}
\det\Omega=\frac{2e^{i\psi}}{\eta\cos(2\delta-\psi)+\cos\psi}\big(\Re\Omega_{11}\Re\Omega_{22}-\Re\Omega_{12}\Re\Omega_{21}\big)\ec
\end{equation}
```

**Eq. (2.56)** — **the interpolation function f(t) = ((t - t_K)/(t_m - t_K))^6** that switches the a-matrices on above t_K  
`fOmnes`

```latex
\begin{equation}
\label{fOmnes}
f(t)=\bigg(\frac{t-\tK}{\tm-\tK}\bigg)^6\ep
\end{equation}
```

> *Verbatim (source lines 550, 567, 578):* "a particular choice of a(t) can enforce continuity of either
> Omega-bar_1(t) or Omega-bar_2(t), but not of both simultaneously. ... We can thus derive an Omega-bar_1^{(1)}(t)
> from a construction with an a_1(t) tailored for this component (discarding Omega-bar_2^{(1)}(t) in this case),
> Omega-bar_2^{(2)}(t) from a different a_2(t) (discarding Omega-bar_1^{(2)}(t)), and finally join these two
> vectors into the final Omnès matrix Omega(t) = {Omega-bar_1^{(1)}(t), Omega-bar_2^{(2)}(t)}."
>
> "This construction makes sure that Omega-bar_ij(t) is continuous at t_m and that the relations (2.30) hold."
>
> "This choice of the exponent ensures that f(t) decreases rapidly below t_m, but not so fast as to cause numerical
> problems in the matching-point region. In this way, the single-channel result for the determinant is accurately
> reproduced except for the energy region close to t_m."

---

# 2. The dispersion relations actually used — Section 3 (`sec:pin_appl`) and Section 4.1

## 2.1 Unitarity relations for the pi-pi -> N-bar N and K-bar K -> N-bar N S-waves

**Eq. (3.1)** — **the coupled unitarity relations** for f^0_+ (pi pi -> N-bar N) and h^0_+ (K-bar K -> N-bar N) with the I_t = 0 partial waves t^0_0, g^0_0, r^0_0  
`unitarity_KN`

```latex
\begin{align}
\label{unitarity_KN}
 \Im f^0_+(t)&=\sigma^\pi_t\big(t^0_0(t)\big)^*f^0_+(t)\,\theta\big(t-\tpi\big)
+\frac{2}{\sqrt{3}}\,\sigma^K_t\big(g^0_0(t)\big)^*h_+^0(t)\,\theta\big(t-\tK\big)\ec\nt\\
\Im h^0_+(t)&=\sigma^K_t\big(r^0_0(t)\big)^*h^0_+(t)\,\theta\big(t-\tK\big)
+\frac{\sqrt{3}}{2}\,\sigma^\pi_t\big(g^0_0(t)\big)^*f_+^0(t)\,\theta\big(t-\tpi\big)\ec
\end{align}
```

**Eq. (3.2)** — identification of t^0_0, g^0_0, r^0_0 with eta, delta, psi, |g|

```latex
\begin{equation}
t^0_0(t)=\frac{\eta(t) e^{2i\delta(t)}-1}{2i\sigma_t^\pi}\ec\qquad g^0_0(t)=|g(t)|e^{i\psi(t)}\ec \qquad 
r^0_0(t)=\frac{\eta(t) e^{2i(\psi(t)-\delta(t))}-1}{2i\sigma_t^K}\ec
\end{equation}
```

**Eq. (3.3)** — the unitarity relation reduces to Im f = T* Sigma f with the two-component vector f = (f^0_+, (2/sqrt3) h^0_+)

```latex
\begin{equation}
\Im \ff(t)=T^*(t)\Sigma(t)\ff(t)\ec\qquad 
\ff(t)=\begin{pmatrix}
               f^0_+(t)\\
\frac{2}{\sqrt{3}}h^0_+(t)
              \end{pmatrix}\ec
\end{equation}
```

## 2.2 The Roy–Steiner dispersion relation and its inhomogeneity

Twice-subtracted RS equations for pi N give (3.4) with t_N = 4m^2; the subtraction constants a, b are related to
pi N subthreshold parameters by (3.6). The inhomogeneity Delta_1 is (3.5); Delta_2 follows from it by replacing the
nucleon pole term by the hyperon poles and M_pi by M_K (text after (3.5)).

**Eq. (3.4)** — **the twice-subtracted RS dispersion relation** for f(t) (the input equation solved with the Omnès method); t_N = 4m^2  
`f_disp_rel`

```latex
\begin{equation}
\label{f_disp_rel}
\ff(t)=\boldsymbol{\Delta}(t)+(\mathbf{a}+\mathbf{b}t)(t-\tN)+\frac{t^2(t-\tN)}{\pi}\int\limits_{\tpi}^\infty\diff t'\frac{\Im\ff(t')}{t'^2(t'-\tN)(t'-t)}\ec
\end{equation}
```

**Eq. (3.5)** — the inhomogeneity Delta_1 = N-hat^0_+ + Delta-bar^0_+ (nucleon pole term plus s-channel partial-wave integral)

```latex
\begin{align}
\Delta_1(t)&=\tilde\Delta^0_+(t)=\hat N^0_+(t)+\bar\Delta^0_+(t)\ec\nt\\

\bar\Delta^0_+(t)&=\frac{1}{\pi}\int\limits^{\infty}_{W_+}\diff W'\sum\limits^\infty_{l=0}
\Big\{\twosub{\tilde G_{0l}}(t,W')\,\Im f^+_{l+}(W')+\twosub{\tilde G_{0l}}(t,-W')\,\Im f^+_{(l+1)-}(W')\Big\}\ep
\end{align}
```

**Eq. (3.6)** — the subtraction constants a_1, b_1 in terms of the pi N subthreshold parameters g^2, d_00^+, d_01^+, b_00^+

```latex
\begin{equation}
a_1=-\frac{1}{16\pi}\Bigg(\frac{g^2}{m}+d_{00}^++b_{00}^+\frac{\mpi^2}{3}\Bigg)\ec\qquad 
b_1= -\frac{1}{16\pi}\Bigg(d_{01}^+-\frac{b_{00}^+}{12}\Bigg)
\end{equation}
```

## 2.3 Subtraction scheme and the final solution

(3.7) is the spectral representation of Omega^{-1} (once- and twice-subtracted forms), used to do the integrals
involving a and b analytically. (3.8) is **the equation the authors actually evaluate**: two subtractions, the
t^2(t - t_N) weighting, the Im Omega^{-1} Delta integral over [t_pi, t_m], and the high-energy closure integral
over [t_m, inf) containing Im f.

**Eq. (3.7)** — **spectral representation of Omega^{-1}** (two equivalent forms); dot = derivative w.r.t. t

```latex
\begin{equation}
 \Omega^{-1}(t)=\unity+\frac{t}{\pi}\int\limits_{\tpi}^{\tm}\diff t'\frac{\Im\Omega^{-1}(t')}{t'(t'-t)}
=\unity-t\,\dot\Omega(0)+\frac{t^2}{\pi}\int\limits_{\tpi}^{\tm}\diff t'\frac{\Im\Omega^{-1}(t')}{t'^2(t'-t)}
\end{equation}
```

**Eq. (3.8)** — **the final Omnès solution for f(t)** — the twice-subtracted equation actually evaluated; note the explicit a and b terms via Omega-dot(0) and the closure integral above t_m  
`final_omnes_sol`

```latex
\begin{align}
\label{final_omnes_sol}
\ff(t)&=\boldsymbol{\Delta}(t)+(t-\tN)\Omega(t)(\unity-t\,\dot\Omega(0))\mathbf{a}+t(t-\tN)\Omega(t)\mathbf{b}\nt\\
&\qquad-\frac{t^2(t-\tN)}{\pi}\Omega(t)\int\limits_{\tpi}^{\tm}\diff t'\frac{\Im\Omega^{-1}(t')\boldsymbol{\Delta}(t')}{t'^2(t'-\tN)(t'-t)}
+\frac{t^2(t-\tN)}{\pi}\Omega(t)\int\limits_{\tm}^\infty\diff t'\frac{\Omega^{-1}(t')\Im \ff(t')}{t'^2(t'-\tN)(t'-t)}\ep
\end{align}
```

## 2.4 Matching point, high-energy closure, and the RS1/RS2/RS3 variants

> *Verbatim (source lines 640, 661–666):* "In view of the remaining uncertainties in the pseudophysical region we
> refrain from including information on f^0_+(t) even above t_N from a phase-shift solution [25], but just put
> Im f = 0 above t_m. Due to continuity of the MO solution at t_m this implies that the solutions for f^0_+(t) and
> h^0_+(t) will vanish at t_m as well (unless the phases happen to pass an integer multiple of pi at t_m). We take
> advantage of the fact that for kinematic reasons f^0_+(t) and h^0_+(t) have a zero at the physical threshold t_N
> and choose the matching point accordingly as **t_m = t_N**, which should allow for a reasonably smooth matching."

> "First, we take delta and psi to be constant above t_0 ("RS1"), second, we guide delta and psi smoothly to 2pi
> above t_0 by means of [15] [eq. (3.9)], keeping the phase constant above t_m ("RS2"), and third, we modify RS1 in
> such a way that Delta_2 = 0 ("RS3"). In all three cases, |g(t)| is led smoothly to zero above t_m by a
> prescription similar to (3.9)."

**Eq. (3.9)** — **the high-energy phase extrapolation** used for RS2 (delta and psi -> 2pi above t_0)  
`phase_extrapolation`

```latex
\begin{equation}
\label{phase_extrapolation}
\delta(t)=2\pi+(\delta(t_0)-2\pi)\hat f\Big(\frac{t}{t_0}\Big)\ec\qquad \hat f(x)=\frac{2}{1+x^{3/2}}\ec 
\end{equation}
```

**Eq. (3.10)** — **the reference values of the pi N parameters** (Höhler): g^2/4pi, d_00^+, d_01^+, b_00^+  
`subthresh_coupling`

```latex
\begin{equation}
\label{subthresh_coupling}
\frac{g^2}{4\pi}=14.28\ec\qquad d_{00}^+=-1.46\,\mpi^{-1}\ec\qquad d_{01}^+=1.14\,\mpi^{-3}\ec\qquad b_{00}^+=-3.54\,\mpi^{-3}\ec
\end{equation}
```

## 2.5 The scalar pion and kaon form factors (what the Omnès matrix multiplies)

The two-channel form-factor vector obeys the same unitarity relation (4.2) and is therefore a linear combination of
the two columns of Omega^inf, (4.3). With Omega^inf(0) = 1 the coefficients are fixed by F_pi^S(0), F_K^S(0),
giving (4.4) — the formula to use for the pi-pi and K-K couplings of a light scalar.

**Eq. (4.1)** — definition of the scalar pion and kaon form factors

```latex
\begin{equation}
F^S_\pi(t)=\langle \pi(p')|\hat m(\bar u u+\bar d d)|\pi(p)\rangle\ec\qquad F^S_K(t)=\langle K(p')|\hat m(\bar u u+\bar d d)|K(p)\rangle\ec\qquad t=(p'-p)^2\ep
\end{equation}
```

**Eq. (4.2)** — unitarity relation for the form-factor vector F^S = (F_pi^S, (2/sqrt3) F_K^S)

```latex
\begin{equation}
\Im\mathbf{F}^S(t)=T^*(t)\Sigma(t)\mathbf{F}^S(t)\ec\qquad \mathbf{F}^S(t)=\begin{pmatrix}
                                                                                  F^S_\pi(t)\\ \frac{2}{\sqrt{3}}F^S_K(t)
                                                                                 \end{pmatrix}\ec
\end{equation}
```

**Eq. (4.3)** — F^S = alpha Omega^inf_1 + beta Omega^inf_2

```latex
\begin{equation}
\mathbf{F}^S(t)=\alpha \mathbf{\Omega}_1^\infty+\beta \mathbf{\Omega}_2^\infty\ec
\end{equation}
```

**Eq. (4.4)** — **the form factors in terms of the infinite-matching-point Omnès matrix** and the normalisation constants F_pi^S(0), F_K^S(0)  
`SFF_Omnes`

```latex
\begin{align}
\label{SFF_Omnes}
 F^S_\pi(t)&=F^S_\pi(0)\Omega_{11}^\infty(t)+\frac{2}{\sqrt{3}}F^S_K(0)\Omega_{12}^\infty(t)\ec\nt\\
F^S_K(t)&=\frac{\sqrt{3}}{2}F^S_\pi(0)\Omega_{21}^\infty(t)+F^S_K(0)\Omega_{22}^\infty(t)\ep
\end{align}
```

**Eq. (4.5)** — F_pi^S(0), F_K^S(0) from the Feynman–Hellmann theorem

```latex
\begin{equation}
F^S_\pi(0)=\hat m\frac{\partial}{\partial \hat m}\mpi^2\ec\qquad
F^S_K(0)=\hat m\frac{\partial}{\partial \hat m}\mK^2\ec
\end{equation}
```

**Eq. (4.6)** — F_pi^S(0) at O(p^4) — the numerical value 0.984 M_pi^2

```latex
\begin{equation}
F^S_\pi(0)=\mpi^2-\frac{\mpi^4}{32\pi^2\Fpi^2}(\bar l_3-1)=(0.984\pm0.006)\mpi^2\ec
\end{equation}
```

**Eq. (4.7)** — leading-order F_K^S(0) = M_pi^2/2

```latex
\begin{equation}
F^S_K(0)=\frac{\mpi^2}{2}
\end{equation}
```

**Eq. (4.8)** — the SU(3) corrected F_K^S(0)

```latex
\begin{align}
 F^S_K(0)&=\frac{\mpi^2}{2}\Bigg\{1+\frac{\meta^2}{32\pi^2\Fpi^2}\log\frac{\meta^2}{\mu^2}+\frac{\mK^2}{72\pi^2\Fpi^2}\bigg(\log\frac{\meta^2}{\mu^2}+1\bigg)-\frac{\mpi^2}{32\pi^2\Fpi^2}\log\frac{\mpi^2}{\mu^2}\nt\\
&\qquad+\frac{8(2\mK^2-\mpi^2)}{\Fpi^2}\big(2L_8^{\rm r}-L_5^{\rm r}\big)+\frac{32\mK^2}{\Fpi^2}\big(2L_6^{\rm r}-L_4^{\rm r}\big)\Bigg\}\ec
\end{align}
```

**Eq. (4.9)** — the ranges of the low-energy constants used, giving F_K^S(0) = 0.4...0.6 M_pi^2

```latex
\begin{equation}
2L_8^{\rm r}-L_5^{\rm r}=(-0.35\ldots +0.1)\cdot 10^{-3}\ec\qquad 2L_6^{\rm r}-L_4^{\rm r}=(0\ldots +0.2)\cdot 10^{-3}\ec
\end{equation}
```

**Eq. (4.10)** — the scalar radii in terms of Omega-dot^inf(0)

```latex
\begin{align}
\langle r^2\rangle^S_\pi&=6\bigg\{\dot\Omega_{11}^\infty(0)+\frac{2}{\sqrt{3}}\frac{F^S_K(0)}{F^S_\pi(0)}\dot\Omega_{12}^\infty(0)\bigg\}\ec\nt\\
\langle r^2\rangle^S_K&=6\bigg\{\frac{\sqrt{3}}{2}\frac{F^S_\pi(0)}{F^S_K(0)}\dot\Omega_{21}^\infty(0)+\dot\Omega_{22}^\infty(0)\bigg\}\ep
\end{align}
```

**Eq. (4.11)** — **the numerical value of Omega-dot^inf(0)** for the central solution (the only printed Omnès-matrix numbers in the paper)

```latex
\begin{equation}
\dot\Omega^\infty(0)=\begin{pmatrix}
                      2.31 & 0.32\\ 1.26 & 0.89
                     \end{pmatrix}
\,{\rm GeV}^{-2}\ec
\end{equation}
```

**Eq. (4.12)** — definition of the nucleon scalar form factor sigma(t)

```latex
\begin{equation}
\sigma(t)=\frac{1}{2m}\langle N(p')|\hat m(\bar u u+\bar d d)|N(p)\rangle\ec\qquad t=(p'-p)^2\ep
\end{equation}
```

**Eq. (4.13)** — the once-subtracted dispersion relation for sigma(t), with sigma_piN as subtraction constant

```latex
\begin{equation}
\sigma(t)=\sigma_{\pi N}+\frac{t}{\pi}\int\limits_{\tpi}^\infty\diff t'\frac{\Im\sigma(t')}{t'(t'-t)}\ec
\end{equation}
```

**Eq. (4.14)** — **the spectral function of the nucleon scalar form factor** — where the meson form factors and the meson–nucleon S-waves enter  
`unitarity_sigma_term`

```latex
\begin{equation}
\label{unitarity_sigma_term}
\Im\sigma(t)=-\frac{1}{p_t^2\sqrt{t}}\bigg\{\frac{3}{4}q_t\big(F^S_\pi(t)\big)^*f^0_+(t)\,\theta\big(t-\tpi\big)+k_t\big(F^S_K(t)\big)^*h^0_+(t)\,\theta\big(t-\tK\big)\bigg\}\ep
\end{equation}
```

**Eq. (4.15)** — **the final numerical result for Delta_sigma** and its dependence on the pi N parameters  
`DelSigRes`

```latex
\begin{align}
\label{DelSigRes}
 \Delta_\sigma&=(13.9\pm 0.3)\,{\rm MeV}\nt\\
&\qquad+ Z_1 \bigg(\frac{g^2}{4\pi}-14.28\bigg)+ Z_2\Big(d_{00}^+\,\mpi+1.46\Big)+
Z_3\Big(d_{01}^+\,\mpi^{3}-1.14\Big)+Z_4\Big(b_{00}^+\,\mpi^{3}+3.54\Big)\ec\nt\\
Z_1&=0.36\,{\rm MeV}\ec\qquad Z_2= 0.57\,{\rm MeV}\ec\qquad Z_3= 12.0\,{\rm MeV}\ec\qquad Z_4=-0.81\,{\rm MeV}\ep
\end{align}
```

**Eq. (4.16)** — the ChPT comparison value

```latex
\begin{equation}
\Delta_\sigma=14.0\,{\rm MeV}+2\mpi^4\, \bar e_2\ec
\end{equation}
```

---

# 3. Numerical method, exactly as stated in the paper

The paper gives **no algorithm in pseudocode and no discretisation parameters**. Everything it says about the
numerics is the following (verbatim, math transliterated):

1. **Two possible strategies** (source line 242): "there is in general no analytic solution for the Omnès matrix
   itself even for an infinite matching point, which therefore has to be calculated numerically, either by an
   iterative procedure [7] or a discretization method, i.e. solving a matrix equation [15,17] (for a
   mathematician's point of view see [18])."

2. **The discretisation** (source lines 380–384): "For its calculation we follow [15] and discretize the
   unsubtracted dispersion relation [eq. (2.32)] on a set of Gauß–Legendre integration points. Note that an
   unsubtracted dispersion relation converges provided that the phase-shift at infinity is positive."

3. **The numerical solution of the resulting real linear system** (source line 423): "The details of the numerical
   solution of the corresponding integral equation for Re Omega^inf_i are described in [15,17]." — the actual
   linear-algebra recipe is **not** in this paper; it is in Moussallam [15] and the Descotes-Genon thesis [17].

4. **Treatment of the cusp / principal value** (source line 242): "Similarly to the single-channel case, we expect
   a cusp at t_m, which has to be taken into account in the numerical evaluation of the integrals in [eq. (2.12)]
   (see [19,1])." The cusp is handled *analytically* through the xi factors (2.43)/(2.56), not by brute force
   (source lines 429–430): "An Omnès function with a finite matching point does not allow for an unsubtracted
   dispersion relation, since the solution [eq. (2.15)] tends to a constant for t -> inf ... Moreover, the cusp at
   the matching point renders both the discretization method and an iterative procedure involving subtracted
   dispersion relations inappropriate, as neither is able to accurately reproduce the analytic behavior around
   t_m." The principal value appears explicitly only in (2.20), (2.32) and (A.1).

5. **Number of subtractions.** Two, in the physical dispersion relation (3.4)/(3.8) (the RS equations are
   twice-subtracted); zero (unsubtracted) for Omega^inf, eq. (2.32), justified by the positive phase at infinity
   (2.34); two subtractions for Omega itself, eq. (2.19), with Omega(0)=1 and Omega(t_m)=0.

6. **High-energy closure.** Im f = 0 above t_m (source line 640) with t_m = t_N = 4m^2, so the closure integral in
   (3.8) vanishes in practice; the phases are held constant above t_0 = (1.3 GeV)^2 (RS1) or extrapolated to 2pi
   with (3.9) (RS2); |g(t)| is led smoothly to zero above t_m.

7. **Accuracy statements.** Only these: (a) source line 666: "we compare the results for det Omega(t)/(1-t/t_m)^x
   to the analytically known results in Fig. [fig:Omnesdet]. As expected, the only deviations occur in the
   proximity of t_m, where the modifications originating from f(t) according to [eq. (2.56)] set in." (b) source
   line 578: the determinant "is accurately reproduced except for the energy region close to t_m". (c) source line
   834: "the dispersion relation converges quickly and the results hardly change above 1 GeV" (this refers to the
   cutoff Lambda of the Delta_sigma integral, Fig. fig:conv_disp_int). No number of grid points, no interpolation
   order, no iteration tolerance and no digit-level accuracy is quoted anywhere.

8. **Cross-checks the authors perform (useful as validation gates for an implementation):**
   - det Omega against the analytic formula (2.15), plotted as det Omega/(1-t/t_m)^x — Fig. `fig:Omnesdet`.
   - single-channel limit: (2.42) must reduce to the standard Omnès function; (2.28) must reduce to x_11 -> x in
     the single-channel limit (source line 344).
   - continuity of f(t) at t_m, proved in Appendix A via the constraints (A.7) and the ratios (2.30).
   - the phase of F_pi^S must equal delta below t_K by unitarity (source line 749, cf. (2.37) and (4.4)).

---

# 4. Numerical inputs required, and where they come from

Reference numbers are those of the paper's bibliography (which is ordered by first citation; verified against the
source). "Given numerically here?" = does this paper print the numbers/parametrisation?

| # | Input | Symbol / equation where it enters | Source in the paper | Given numerically here? |
|---|---|---|---|---|
| 1 | pi-pi S-wave phase shift | delta(t), (2.3), (2.28), (2.53) | [21] Caprini–Colangelo–Leutwyler, EPJC 72 (2012) 1860 (Regge analysis); [22] same authors, PWA "in preparation" | **No.** Only the high-energy extrapolation (3.9) is specified; the phase itself must come from [21,22]. |
| 2 | pi-pi inelasticity | eta(t), (2.3), (2.5) | [21,22] | **No.** |
| 3 | pi-pi -> K K S-wave modulus (and phase) | abs g(t), psi(t), (2.3), (2.30), (2.53) | [19] Büttiker–Descotes-Genon–Moussallam, EPJC 33 (2004) 409 (pseudophysical region from RS equations for pi K); above t_K: [23] Cohen et al. PRD 22 (1980) 2595, [24] Etkin et al. PRD 25 (1982) 1786 | **No.** Footnote 3 (source line 638): the g(t) used is a Moussallam private communication consistent with [21,22]. |
| 4 | Matching point | t_m = t_N = 4m^2 (source line 640) | choice, with m from [13] PDG | Choice given; m not printed |
| 5 | Onset of 4-pi inelasticity | t_0, sqrt(t_0) = 1.3 GeV (source line 638) | statement in text | **Yes** (1.3 GeV) |
| 6 | Meson/nucleon masses | M_pi, M_K, m | [13] PDG | **No** (only symbolic; hyperon masses in App. B are printed) |
| 7 | Scalar form factors at t = 0 | F_pi^S(0) = 0.984 M_pi^2, F_K^S(0) = M_pi^2/2 (range 0.4–0.6 M_pi^2) | (4.6) with l3-bar = 3.2 +- 0.8 [33] FLAG, [32] GL84; (4.7)–(4.9) with LECs from [33] FLAG, [34] GL85 | **Yes** (4.6), (4.9), Table 1 |
| 8 | pi N s-channel partial waves (for Delta_1) | Im f^{I}_{l+-}(W), (3.5) | [26] GWU:2008 / [27] SAID, W <= W_a = 2.5 GeV, all partial waves up to l_max = 4; [28] Regge model above W_a | **No.** |
| 9 | K N s-channel partial waves (for Delta_2) | Im f, (3.5) with M_pi -> M_K | [27] SAID, [29] GWU:KN, l <= 4 up to 2 GeV | **No.** |
| 10 | pi N subthreshold parameters + coupling | g^2/4pi = 14.28, d_00^+ = -1.46 M_pi^{-1}, d_01^+ = 1.14 M_pi^{-3}, b_00^+ = -3.54 M_pi^{-3} | (3.10), from [8] Höhler, Landolt–Börnstein 9b2 | **Yes**, eq. (3.10) |
| 11 | K N subthreshold parameters | a_2 = b_2 = 0 (approximation) | source line 673 | stated approximation, not an input |
| 12 | Hyperon masses and couplings (Delta_2 Born terms) | m_Lambda = 1.116 GeV, m_Sigma = 1.193 GeV, g^2_{KN Lambda}/4pi = 15.55, g^2_{KN Sigma}/4pi = 0.576 | footnote 4, App. B (source line 967); [48] Holzenkamp–Holinde–Speth, NPA 500 (1989) 485 | **Yes** (footnote) |
| 13 | 4-pi, omega, rho intermediate states | — | neglected; footnote 1: "we neglect 4pi intermediate states, which are phenomenologically irrelevant at low energies" | not needed |

---

# 5. What is missing for a self-contained implementation

1. **No numerical parametrisation of the S-matrix.** The paper prints *no* table of phase-shift values, *no*
   polynomial or conformal-fit coefficients for delta(t), eta(t) or abs g(t), and *no* pole parameters for the
   f_0(980). The complete set of numbers in the paper is: F_pi^S(0) (4.6), the LEC ranges (4.9), Omega-dot^inf(0)
   (4.11), the radii (Table `tab:radii`), the pi N reference parameters (3.10), the hyperon constants (footnote 4),
   and Delta_sigma (Table `tab:DelSig`, eq. (4.15)). To evaluate Omega you must obtain delta, eta, abs g from
   [21,22] and [19] independently, or from the published parametrisations of those papers.

2. **No discretisation recipe of its own.** "Discretize ... on a set of Gauss–Legendre integration points" (2.32)
   is the whole statement. The number of points, the mapping of the infinite interval onto a finite one, the
   treatment of the integrable singularities at the two thresholds and at t_m, the interpolation between grid
   points, and the iteration/linear-solve tolerances are all deferred to [15] and [17]. Anyone re-implementing this
   must reconstruct those choices (or take them from Winkler's implementation).

3. **The private-communication input.** abs g(t) in the pseudophysical region t_pi <= t <= t_K is not reproducible
   from the paper: it is a Moussallam-provided solution consistent with the CCL phase (footnote 3). Without it the
   analytic continuation of g below t_K, which (2.3)/(2.30) require, is an open modelling choice.

4. **The a(t) construction is only specified up to conventions.** (2.53) gives a_1, a_2 with f(t) from (2.56), but
   the "join the two columns" step (which column comes from which a) is described in words only, and the paper never
   states at which t the switch to the second construction happens in practice (it is implied to be applied for
   t >= t_K). Also (2.52) fixes Omega^inf(0) after the fact: one must normalise Omega^inf by Omega^inf(0)^{-1}
   computed from (2.52), which requires N, eta, abs g at t = 0, i.e. an extrapolation of abs g to t = 0.

5. **No error budget for Omega itself.** The only accuracy statements are qualitative (see §3.7). The Delta_sigma
   uncertainty in (4.15) is a physical/parametric error, not a numerical one.

6. **What is not needed for a pi-pi / K K width calculation:** the whole pi N / K N machinery (Section 3, eqs.
   (3.1)–(3.8)) is only needed for the *nucleon* scalar form factor. The Omega matrix is fully specified by
   Section 2 plus the three functions delta, eta, abs g — that is the part relevant for the two-channel MO analysis
   of Winkler arXiv:1809.01876.

---

# 6. EPS validation targets

All five files sit in the same directory as `SFF_v2.tex`. They are Grace (xmgrace 5.1.22) EPS files; the axis
labels are embedded as bitmaps, so the ranges below were read by rendering the EPS locally (no network used).
**All five plot the *finite-matching-point* Omega (the object of (2.19), which vanishes at t_m), with
t_m = t_N = 4m^2, i.e. sqrt(t_m) = 1.877 GeV — not Omega^inf.**

| File | What it plots | Axes (read from the file) |
|---|---|---|
| `Omega11.eps` | Re and Im of Omega_11(t) for input RS1 (solid / dashed black) and RS2 (dashed / dash-dot red). Legend: `Re Omega_11^{RS1}`, `Im Omega_11^{RS1}`, `Re Omega_11^{RS2}`, `Im Omega_11^{RS2}`. | x = sqrt(t)/GeV, 0.3 … 1.9; y = -0.5 … 1.5 |
| `Omega12.eps` | Re and Im of Omega_12(t), same four curves/colours. | x = sqrt(t)/GeV, 0.3 … 1.9; y = -2.5 … 1.5 (spike at sqrt(t) ~ 0.99 clipped) |
| `Omega21.eps` | Re and Im of Omega_21(t), same four curves/colours. | x = sqrt(t)/GeV, 0.3 … 1.9; y = -2.5 … 1.0 (spike at sqrt(t) ~ 0.99 clipped) |
| `Omega22.eps` | Re and Im of Omega_22(t), same four curves/colours. | x = sqrt(t)/GeV, 0.3 … 1.9; y = -0.5 … 4.5 (spike at sqrt(t) ~ 0.99 clipped) |
| `detOmega.eps` | det Omega(t) / (1 - t/t_m)^x: numerical (solid = RS1, dashed = RS2) versus the **analytic** result of (2.15) (dotted, same colours). Legend: `RS1`, `RS1, analytical solution`, `RS2`, `RS2, analytical solution`. | x = sqrt(t)/GeV, 0.3 … 1.9; y = 0 … 7 |

These correspond to Fig. `fig:Omnesmatrix` (the four Omega_ij panels) and Fig. `fig:Omnesdet` (the determinant).

Reading guide for validation:

* All four Omega_ij approach 0 as sqrt(t) -> sqrt(t_m) = 1.877 GeV, because Omega(t_m) = 0 is imposed in (2.19). A
  solver that reproduces Omega^inf instead will *not* show this.
* The sharp structure at sqrt(t) ~ 0.99 GeV in Omega_12, Omega_21, Omega_22 is the two-kaon threshold
  (2M_K = 0.9878 GeV), i.e. the opening of the inelastic channel; it is the most sensitive test of the xi
  construction and of the a(t) matching.
* Omega_11 and Omega_22 are the "diagonal" channels; below t_K the phase of Omega_11 equals delta by unitarity
  (cf. (2.37), which gives Im Omega_11 = tan(delta) Re Omega_11 there).
* The determinant plot is the cleanest single test: (2.15) is analytic, so det Omega/(1 - t/t_m)^x should agree with
  the dotted curve everywhere except in the immediate neighbourhood of t_m, where f(t) of (2.56) switches on.
* RS1 and RS2 differ only through the high-energy behaviour of delta and psi (see (3.9) and item 6 of section 3); the
  difference between the two sets is a direct measure of the high-energy sensitivity.

---

# Appendix A. Index of every equation extracted here

| # | Eq. | label | content |
|---|---|---|---|
| 1 | **2.1** | — | the generic two-channel Muskhelishvili–Omnès integral equation with finite matching point |
| 2 | **2.2** | — | unitarity relation, two-channel form |
| 3 | **2.3** | `Tmatrix` | **T matrix and Sigma matrix** in terms of delta, psi, eta, \|g\| — this is *the* two-channel S-matrix parameterisation used throughout the paper |
| 4 | **2.4** | — | t-channel momenta q_t, k_t, p_t and the phase-space factors |
| 5 | **2.5** | — | inelasticity eta expressed through \|g\| (only delta, psi, \|g\| are independent input) |
| 6 | **2.6** | — | S = 1 + 2i Sigma^{1/2} T Sigma^{1/2} (S-matrix in terms of T and the phase-space factors) |
| 7 | **2.7** | `def_finite` | **defining property of the Omnès matrix**: Im Omega = T* Sigma Omega on [t_pi, t_m], Im Omega = 0 elsewhere; normalisation Omega(0) = 1 (stated in the text right after the equation) |
| 8 | **2.8** | — | FF = f - Delta = Omega G (definition of the auxiliary vector G) |
| 9 | **2.9** | — | jump condition for G |
| 10 | **2.10** | — | unitarity in the form (1 - 2i T* Sigma)^{-1} = 1 + 2i T Sigma |
| 11 | **2.11** | — | G(t+) - G(t-) = 2i Omega^{-1} T Sigma Delta |
| 12 | **2.12** | `omnes_fin` | **the formal MO solution**; the second line needs only Im Omega^{-1} on [t_pi, t_m] |
| 13 | **2.13** | — | Omega(t+) = (1 + 2i T Sigma) Omega(t-) for t_pi <= t <= t_K |
| 14 | **2.14** | — | det Omega(t+) = exp(2i psi) det Omega(t-) |
| 15 | **2.15** | `det_fin` | **analytic determinant** of the Omnès matrix (Moussallam's formula) — the validation target for a numerical Omega |
| 16 | **2.16** | `det_asym` | cusp behaviour det Omega ~ \|t_m - t\|^x with x = psi(t_m)/pi |
| 17 | **2.17** | `def_below` | behaviour for t -> t_m **from below**: the powers x_ij and the phases delta_ij(t_m) |
| 18 | **2.18** | `def_above` | same **from above**, where Omega is real |
| 19 | **2.19** | — | **subtracted dispersive representation of Omega(t)** with Omega(0)=1, Omega(t_m)=0; note the kernel T* Sigma Omega |
| 20 | **2.20** | `integrals` | the two asymptotic integrals (PV = principal value) giving -cot(pi x) and 1/sin(pi x) |
| 21 | **2.21** | `Omegaij_above` | component-wise matching conditions for t -> t_m **from below** (all four equations) |
| 22 | **2.22** | `Omegaij_below` | component-wise matching conditions for t -> t_m **from above** |
| 23 | **2.23** | — | consistency relation from dividing the first lines of (2.21) and (2.22) |
| 24 | **2.24** | — | exp(i delta_11) = exp(i pi x_11) |
| 25 | **2.25** | `xij` | **the exponents**: x_11 = x_21, x_12 = x_22, x_ij = delta_ij/pi, x_11 + x_12 = x |
| 26 | **2.26** | `spectral_rep_const` | the two 2x2 homogeneous systems for (Omega-bar_11, Omega-bar_21) and (Omega-bar_12, Omega-bar_22) at t_m |
| 27 | **2.27** | — | determinant conditions: cos pi(2x_11 - x) - eta cos pi(2y - x) = 0, y = delta(t_m)/pi |
| 28 | **2.28** | `x11result` | **closed-form solution for the exponents** x_11, x_12 in terms of x, eta, y |
| 29 | **2.29** | `sinx11` | sin and cos of pi x_11, pi x_12 in terms of z = eta cos pi(2y - x) |
| 30 | **2.30** | `ratio_Omij` | **the ratios Omega-bar_21/Omega-bar_11 and Omega-bar_12/Omega-bar_22 at t_m**, with N defined here (used everywhere below) |
| 31 | **2.31** | `Omnes_def` | **defining property of Omega^inf(t)** (infinite matching point) |
| 32 | **2.32** | — | **unsubtracted dispersion relation for Re Omega^inf** (the equation discretised on Gauss–Legendre points); dash = principal value |
| 33 | **2.33** | `sol_single_infinite` | closed-form single-channel Omega^inf (limiting case / sanity check) |
| 34 | **2.34** | — | asymptotic behaviour of Omega^inf for large t, with exponent delta(inf)/pi |
| 35 | **2.35** | `unitarity` | **unitarity rewritten as a linear equation for Re Omega^inf_i** — the 2x2 matrix multiplying the real part |
| 36 | **2.36** | `Omnes_vector` | definition of the column vectors Omega^inf_1, Omega^inf_2 |
| 37 | **2.37** | `Omnes_below_tK` | the same relation below t_K (Watson regime): tan delta and \|g\| sigma^pi_t / cos delta |
| 38 | **2.38** | — | the single-channel xi(t) = ((t_m - t)/t_m)^{x(t)} and its phase convention |
| 39 | **2.39** | `Omnes_constr` | Omega(t) = Omega^inf(t) xi(t) |
| 40 | **2.40** | — | the two defining properties of this Omega (unitarity below t_m; explicit imaginary part above) |
| 41 | **2.41** | — | the choice of x(t): delta(t_m)/pi for t <= t_m, delta(t)/pi above |
| 42 | **2.42** | `omnes_sol` | explicit single-channel solution with finite t_m (two exponential factors) |
| 43 | **2.43** | — | **the xi matrix** xi_ij(t) = xi-bar_ij(t) ((t_m - t)/t_m)^{x_ij(t)}, with real functions xi-bar_ij |
| 44 | **2.44** | `sol_fin_Omega_til` | the two-branch ansatz: a Omega^inf xi below t_m and a (Omega^inf)^T xi above t_m |
| 45 | **2.45** | `cond_xi` | the two conditions on xi and on a(t) (commutator condition for a) |
| 46 | **2.46** | `eq_xi` | **the 2x2 systems that determine xi-bar_11, xi-bar_21 and xi-bar_12, xi-bar_22** |
| 47 | **2.47** | — | ratios xi-bar_21/xi-bar_11 and xi-bar_12/xi-bar_22 |
| 48 | **2.48** | — | the normalisation choice for xi that makes det xi equal to the single-channel factor |
| 49 | **2.49** | — | det xi(t) = ((t_m - t)/t_m)^x |
| 50 | **2.50** | — | the constraints on a(t): a_21 = a_12 sigma^pi_t/sigma^K_t, a_22 = a_11 - a_12 eta sin(2delta-psi)/(\|g\| sigma^K_t) |
| 51 | **2.51** | — | normalisation condition Omega^inf(0) xi(0) = 1 written out |
| 52 | **2.52** | — | **the required normalisation of Omega^inf(0)** that enforces Omega(0) = 1 |
| 53 | **2.53** | — | **the explicit real matrices a_1(t), a_2(t)** used for t >= t_K (with f(t_K)=0, f(t_m)=1); below t_K, a_1 = a_2 = 1 |
| 54 | **2.54** | — | N(t) = sqrt(1 - eta^2 cos^2(2delta - psi)) - eta sin(2delta - psi) |
| 55 | **2.55** | — | det Omega in terms of the real parts |
| 56 | **2.56** | `fOmnes` | **the interpolation function f(t) = ((t - t_K)/(t_m - t_K))^6** that switches the a-matrices on above t_K |
| 57 | **3.1** | `unitarity_KN` | **the coupled unitarity relations** for f^0_+ (pi pi -> N-bar N) and h^0_+ (K-bar K -> N-bar N) with the I_t = 0 partial waves t^0_0, g^0_0, r^0_0 |
| 58 | **3.2** | — | identification of t^0_0, g^0_0, r^0_0 with eta, delta, psi, \|g\| |
| 59 | **3.3** | — | the unitarity relation reduces to Im f = T* Sigma f with the two-component vector f = (f^0_+, (2/sqrt3) h^0_+) |
| 60 | **3.4** | `f_disp_rel` | **the twice-subtracted RS dispersion relation** for f(t) (the input equation solved with the Omnès method); t_N = 4m^2 |
| 61 | **3.5** | — | the inhomogeneity Delta_1 = N-hat^0_+ + Delta-bar^0_+ (nucleon pole term plus s-channel partial-wave integral) |
| 62 | **3.6** | — | the subtraction constants a_1, b_1 in terms of the pi N subthreshold parameters g^2, d_00^+, d_01^+, b_00^+ |
| 63 | **3.7** | — | **spectral representation of Omega^{-1}** (two equivalent forms); dot = derivative w.r.t. t |
| 64 | **3.8** | `final_omnes_sol` | **the final Omnès solution for f(t)** — the twice-subtracted equation actually evaluated; note the explicit a and b terms via Omega-dot(0) and the closure integral above t_m |
| 65 | **3.9** | `phase_extrapolation` | **the high-energy phase extrapolation** used for RS2 (delta and psi -> 2pi above t_0) |
| 66 | **3.10** | `subthresh_coupling` | **the reference values of the pi N parameters** (Höhler): g^2/4pi, d_00^+, d_01^+, b_00^+ |
| 67 | **4.1** | — | definition of the scalar pion and kaon form factors |
| 68 | **4.2** | — | unitarity relation for the form-factor vector F^S = (F_pi^S, (2/sqrt3) F_K^S) |
| 69 | **4.3** | — | F^S = alpha Omega^inf_1 + beta Omega^inf_2 |
| 70 | **4.4** | `SFF_Omnes` | **the form factors in terms of the infinite-matching-point Omnès matrix** and the normalisation constants F_pi^S(0), F_K^S(0) |
| 71 | **4.5** | — | F_pi^S(0), F_K^S(0) from the Feynman–Hellmann theorem |
| 72 | **4.6** | — | F_pi^S(0) at O(p^4) — the numerical value 0.984 M_pi^2 |
| 73 | **4.7** | — | leading-order F_K^S(0) = M_pi^2/2 |
| 74 | **4.8** | — | the SU(3) corrected F_K^S(0) |
| 75 | **4.9** | — | the ranges of the low-energy constants used, giving F_K^S(0) = 0.4...0.6 M_pi^2 |
| 76 | **4.10** | — | the scalar radii in terms of Omega-dot^inf(0) |
| 77 | **4.11** | — | **the numerical value of Omega-dot^inf(0)** for the central solution (the only printed Omnès-matrix numbers in the paper) |
| 78 | **4.12** | — | definition of the nucleon scalar form factor sigma(t) |
| 79 | **4.13** | — | the once-subtracted dispersion relation for sigma(t), with sigma_piN as subtraction constant |
| 80 | **4.14** | `unitarity_sigma_term` | **the spectral function of the nucleon scalar form factor** — where the meson form factors and the meson–nucleon S-waves enter |
| 81 | **4.15** | `DelSigRes` | **the final numerical result for Delta_sigma** and its dependence on the pi N parameters |
| 82 | **4.16** | — | the ChPT comparison value |
| 83 | **A.1** | — | f(t) rewritten with a principal-value integral (starting point of the continuity proof) |
| 84 | **A.2** | `cont_cond` | the limit t -> t_m of f(t), written with the integral matrices I and I-tilde |
| 85 | **A.3** | `def_I` | the definitions of I_I ... I_IV and I-tilde_I ... I-tilde_IV |
| 86 | **A.4** | `sol_I` | the solution for I_I ... I_IV |
| 87 | **A.5** | — | A = eta exp(2i pi (x-y)) - 1, B = eta exp(2i pi y) - 1 |
| 88 | **A.6** | `sol_tilde_I` | the solution for the I-tilde matrix |
| 89 | **A.7** | `cont_conditions` | **the three conditions that guarantee continuity at t_m** (they must hold for the numerical Omega; the check involves det Omega-bar and the ratios (2.30)) |

---

# Appendix B. Source locations (for re-checking)

| Content | Source lines in `SFF_v2.tex` |
|---|---|
| Section 2 header and eqs. (2.1)–(2.6) | 156–191 |
| Section 2.1, eqs. (2.7)–(2.18) | 194–257 |
| Section 2.2, eqs. (2.19)–(2.30) | 259–365 |
| Section 2.3.1 (Omega^inf), eqs. (2.31)–(2.37) | 370–423 |
| Section 2.3.2 (single channel), eqs. (2.38)–(2.42) | 426–461 |
| Section 2.3.3 (two channels), eqs. (2.43)–(2.56) | 463–578 |
| Section 3, eqs. (3.1)–(3.10) | 581–673 |
| Section 4.1, eqs. (4.1)–(4.11) | 689–779 |
| Section 4.2, eqs. (4.12)–(4.16) | 782–846 |
| Appendix A (continuity at t_m), eqs. (A.1)–(A.7) | 878–930 |
| Appendix B (K N conventions), eqs. (B.1)–(B.16) | 933–1020 |
| Bibliography (49 entries) | 1024–1340 |

### Appendix A of the paper (continuity at the matching point) — an extra numerical check

These equations prove that the construction of §1.6 makes f(t) continuous at t_m. They are not needed to *build*
Omega, but they are a strong consistency check on a numerical implementation.

**Eq. (A.1)** — f(t) rewritten with a principal-value integral (starting point of the continuity proof)

```latex
\begin{align}
 \ff(t)&=\big(\unity+i\,T(t)\Sigma(t)\big)\Del(t)+\frac{\Omega(t)}{\pi}\dashint\limits_{\tpi}^{\tm}\diff t'\frac{\Omega^{-1}(t')T(t')\Sigma(t')\Del(t')}{t'-t}+
\frac{\Omega(t)}{\pi}\int\limits_{\tm}^\infty\diff t'\frac{\Omega^{-1}(t')\Im \ff(t')}{t'-t}\ep
\end{align}
```

**Eq. (A.2)** — the limit t -> t_m of f(t), written with the integral matrices I and I-tilde  
`cont_cond`

```latex
\begin{align}
\label{cont_cond}
 \ff(\tm)&=\big(\unity+i\,T(\tm)\Sigma(\tm)\big)\Del(\tm)+\frac{e^{-i\pi x}}{\det\bar\Omega(\tm)}\begin{pmatrix} I_{\rm I} & I_{\rm II}\\ I_{\rm III} & I_{\rm IV}  \end{pmatrix} T(\tm)\Sigma(\tm)\Del(\tm)\nt\\
&\qquad- \frac{1}{\det\bar\Omega(\tm)}\begin{pmatrix} \tilde I_{\rm I} & \tilde I_{\rm II}\\ \tilde I_{\rm III} & \tilde I_{\rm IV}  \end{pmatrix} \Im\ff(\tm)\ec                                                                                      
\end{align}
```

**Eq. (A.3)** — the definitions of I_I ... I_IV and I-tilde_I ... I-tilde_IV  
`def_I`

```latex
\begin{align}
\label{def_I}
 I_{\rm I}&=-\Big(\bar\Omega_{11}\bar\Omega_{22}\cot\pi x_{11}-\bar\Omega_{12}\bar\Omega_{21}\cot\pi x_{12}\Big)e^{i\psi}\ec
&I_{\rm II}&=\bar\Omega_{11}\bar\Omega_{12}\Big(\cot\pi x_{11}-\cot\pi x_{12}\Big)e^{i\psi}\ec\nt\\
I_{\rm III}&=-\bar\Omega_{21}\bar\Omega_{22}\Big(\cot\pi x_{11}-\cot\pi x_{12}\Big)e^{i\psi}\ec
&I_{\rm IV}&=-\Big(\bar\Omega_{11}\bar\Omega_{22}\cot\pi x_{12}-\bar\Omega_{12}\bar\Omega_{21}\cot\pi x_{11}\Big)e^{i\psi}\ec\nt\\

\tilde I_{\rm I}&=-\bar\Omega_{11}\bar\Omega_{22}\frac{e^{i\delta_{11}}}{\sin\pi x_{11}}+\bar\Omega_{12}\bar\Omega_{21}\frac{e^{i\delta_{12}}}{\sin\pi x_{12}}\ec
&\tilde I_{\rm II}&=\bar\Omega_{11}\bar\Omega_{12}\bigg\{\frac{e^{i\delta_{11}}}{\sin\pi x_{11}}-\frac{e^{i\delta_{12}}}{\sin\pi x_{12}}\bigg\}\ec\nt\\
\tilde I_{\rm III}&=-\bar\Omega_{21}\bar\Omega_{22}\bigg\{\frac{e^{i\delta_{11}}}{\sin\pi x_{11}}-\frac{e^{i\delta_{12}}}{\sin\pi x_{12}}\bigg\}\ec
&\tilde I_{\rm IV}&=-\bar\Omega_{11}\bar\Omega_{22}\frac{e^{i\delta_{12}}}{\sin\pi x_{12}}+\bar\Omega_{12}\bar\Omega_{21}\frac{e^{i\delta_{11}}}{\sin\pi x_{11}}\ep
\end{align}
```

**Eq. (A.4)** — the solution for I_I ... I_IV  
`sol_I`

```latex
\begin{align}
\label{sol_I}
 I_{\rm I}&=i\det \bar\Omega\, e^{i\pi x}\frac{A-e^{2i\pi x}A^*}{A+e^{2i\pi x}A^*}\ec\qquad 
I_{\rm IV}=i\det \bar\Omega\, e^{i\pi x}\frac{B-e^{2i\pi x}B^*}{A+e^{2i\pi x}A^*}\ec\nt\\
I_{\rm II}&=\frac{\sigma_{\tm}^K}{\sigma_{\tm}^\pi}I_{\rm III}=4\det \bar\Omega\, e^{2i\pi x}\frac{|g|\sigma_{\tm}^K}{A+e^{2i\pi x}A^*}\ec
\end{align}
```

**Eq. (A.5)** — A = eta exp(2i pi (x-y)) - 1, B = eta exp(2i pi y) - 1

```latex
\begin{equation}
A=\eta e^{2i\pi(x-y)}-1\ec\qquad B=\eta e^{2i\pi y}-1\ep
\end{equation}
```

**Eq. (A.6)** — the solution for the I-tilde matrix  
`sol_tilde_I`

```latex
\begin{equation}
\label{sol_tilde_I}
 \tilde I_{\rm I}=-2i\,\det \bar\Omega\, e^{2i\pi x}\frac{A^*}{A+e^{2i\pi x}A^*}\ec\qquad\!\! 
\tilde I_{\rm IV}=-2i\,\det \bar\Omega\, e^{2i\pi x}\frac{B^*}{A+e^{2i\pi x}A^*}\ec\qquad\!\!
\tilde I_{\rm II}=\frac{\sigma_{\tm}^K}{\sigma_{\tm}^\pi}\tilde I_{\rm III}=I_{\rm II}e^{-i\pi x}\ep
\end{equation}
```

**Eq. (A.7)** — **the three conditions that guarantee continuity at t_m** (they must hold for the numerical Omega; the check involves det Omega-bar and the ratios (2.30))  
`cont_conditions`

```latex
\begin{align}
\label{cont_conditions}
 \sqrt{1-z^2}\big(\bar\Omega_{11}\bar\Omega_{22}+\bar\Omega_{12}\bar\Omega_{21}\big)&=\eta\sin\pi(2y-x) \det\bar\Omega\ec\nt\\
\sqrt{1-z^2}\bar\Omega_{11}\bar\Omega_{12}&=-|g|\sigma_{\tm}^K\det\bar\Omega\ec\qquad\qquad \sqrt{1-z^2}\bar\Omega_{21}\bar\Omega_{22}=|g|\sigma_{\tm}^\pi\det\bar\Omega\ep
\end{align}
```
