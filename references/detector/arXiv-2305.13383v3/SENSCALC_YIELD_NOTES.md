# SensCalc yield notes — reading of arXiv:2305.13383v3 (main.tex)

**Source read:** E:\AAAwork\LLP-AGI\references\detector\arXiv-2305.13383v3\main.tex (591 lines, 95 kB), main.bbl, and the figure PDFs in the same directory.

**Secondary source (flagged, outside the SensCalc paper):** Gorkavenko et al., EPJC 84 (2024) 608 — the PDF sits in the parent directory. Used only in §4 and §7, because the SensCalc paper does not document its LHCb implementation at all.

**Conventions.** Line numbers are lines of main.tex as shipped in the arXiv e-print. Quotes are verbatim (LaTeX retained where the maths matters; [...] marks an elision). Figures were rendered with PyMuPDF at 200–420 dpi; the numbers in §5 and §7 were obtained by extracting the vector paths from the PDFs (and, for one figure from the EPJC paper, by colour-digitising a raster). The method is stated wherever a number is quoted, and every such number is approximate.

---

# 1. The number-of-events formula

## 1.1 Verbatim

Lines 76–81:

> 76: This work concentrates on FIPs produced directly at the collision point or in its immediate vicinity. In this case, the production is unaffected by the surrounding infrastructure. We calculate the number of events involving a decaying FIP using the following expression:
> 77: \begin{multline}
> 78:     N_{\text{ev}} = \sum_{i}N^{(i)}_{\text{prod}}\int dE d\theta dz \ f^{(i)}(\theta,E)\cdot \epsilon_{\text{az}}(\theta,z) \cdot \\
> 79:     \cdot \frac{dP_{\text{dec}}}{dz} \cdot\epsilon_{\text{dec}}(m,\theta,E,z) \cdot \epsilon_{\text{rec}}
> 80:     \label{eq:Nevents}
> 81: \end{multline}

that is

> N_ev = Σ_i N_prod^(i) ∫ dE dθ dz · f^(i)(θ,E) · ε_az(θ,z) · (dP_dec/dz) · ε_dec(m,θ,E,z) · ε_rec        (Eq. 1)

**Note the measure: dE dθ dz — no dφ and no sin θ Jacobian.** f is a density per unit θ, not per unit solid angle, and the azimuth is already averaged out (l. 255: "the averaging over φ is already possible at this stage since the other quantities that determine the number of events (1) do not depend on the azimuthal angle"). This matters if you compare against a generator that histogrammed dN/dΩ, dN/dη or dN/dcos θ.

Supporting equations, verbatim.

Azimuthal acceptance (Eq. 2), lines 88–92:

> 88: \item[--] $\epsilon_{\text{az}}(\theta,z)$ is the azimuthal acceptance:
> 89:     \begin{equation}
> 90:     \epsilon_{\text{az}} = \frac{\Delta \phi_{\text{decay volume}}(\theta,z)}{2\pi}
> 91:     \end{equation}
> 92: where $\Delta\phi$ is the fraction of azimuthal coverage for which FIPs decaying at $(z,\theta)$ are inside the decay volume.

Differential decay probability (Eq. 3), lines 93–97:

> 93: \item[--] $\frac{dP_{\text{dec}}}{dz}$ is the differential decay probability:
> 94: \begin{equation}
> 95: \frac{dP_{\text{dec}}}{dz} = \frac{\exp[-r(z,\theta)/l_{\text{dec}}]}{l_{\text{dec}}} \frac{dr(z,\theta)}{dz},
> 96: \end{equation}
> 97: with $r = z/\cos(\theta)$ being the modulus of the displacement of the FIP decay position from its production point, and $l_{\text{dec}}=c\tau \sqrt{\gamma^2-1}$ is the FIP decay length in the lab frame.

Decay-product acceptance (Eq. 4), lines 98–103:

> 98: \item[--] $\epsilon_{\text{dec}}(m,\theta,E,z)$ is the decay products acceptance, i.e.\ among those FIPs that are within the azimuthal acceptance, the fraction of FIPs that have at least two decay products that point to the detector and that may be reconstructed. Schematically,
> 99: \begin{equation}
> 100:     \epsilon_{\text{dec}} = \text{Br}_{\text{vis}}(m)\cdot \epsilon_{\text{dec}}^{\text{(geom)}}\cdot \epsilon_{\text{dec}}^{\text{(other cuts)}}
> 101:     \label{eq:decay-acceptance-schematic}
> 102: \end{equation}
> 103: Here, $\text{Br}_{\text{vis}}$ denotes the branching ratio of the FIP decays into final states that are detectable; [...] $\epsilon_{\text{dec}}^{\text{(geom)}}$ denotes the fraction of visible decay products that point to the end of the detector, and $\epsilon_{\text{dec}}^{\text{(other cuts)}}$ is the fraction of these decay products that additionally satisfy the remaining cuts (e.g., the energy cut, etc.).

Reconstruction efficiency, line 104:

> 104: \item[--]$\epsilon_{\text{rec}}$ is the reconstruction efficiency, i.e., the fraction of the events that pass the azimuthal and decay acceptances criteria that the detector can successfully reconstruct. It results from the non-ideal performance of the detector, which introduces a finite detection efficiency and kinematics measurement resolution.

Long-lived limit, lines 139–154:

> 139: Last but not least, this semi-analytic method allows for a simple analysis of the number of events in the limit of very long-lived FIPs with lifetimes $c\tau\langle \sqrt{\gamma^{2}-1}\rangle\gg l_{\text{experiment}}$, where $l_{\text{experiment}}$ is the length scale of the experiment. In this case, the only dependence of the number of events~\eqref{eq:Nevents} on $c\tau$ is multiplicative:
> 140: \begin{equation}
> 141:     N_{\text{ev}} \approx \sum_{i}N_{\text{prod}}^{(i)}\cdot \epsilon^{(i)},
> 142:     \label{eq:number-of-events-large-lifetimes}
> 143: \end{equation}
> 144: where $\epsilon^{(i)}$ is the total acceptance for the given production channel:
> 145: \begin{equation}
> 146:     \epsilon^{(i)} = \int d\theta dE dz\ f^{(i)}\cdot \epsilon_{\text{az}}\cdot \frac{\epsilon_{\text{dec}}}{\cos(\theta)c\tau\sqrt{\gamma^{2}-1}},
> 147:     \label{eq:total-acceptance}
> 148: \end{equation}
> 149: and the function parameters have been omitted for brevity.
> 150: This quantity may be decomposed as
> 151: \begin{equation}
> 152:     \epsilon = \langle \epsilon_{\text{FIP}}\rangle\cdot \langle (\gamma^{2}-1)^{-1/2}\rangle\cdot \frac{\Delta z}{c\tau}\cdot \langle \epsilon_{\text{decay}}\rangle,
> 153: \end{equation}
> 154: where $\langle \epsilon_{\text{FIP}}\rangle$ is the mean probability for the FIP to intersect the decay volume, $\langle (\gamma^{2}-1)^{-1/2}\rangle$ is the mean inverse $p/m$ among the FIPs meeting the azimuthal criterion, $\Delta z$ is the total longitudinal size of the decay volume, and $\langle \epsilon_{\text{decay}}\rangle$ is the mean decay products acceptance.

The only closure test the paper offers is the volume check, lines 122–129:

> 122: [...] Once $\epsilon_{\text{az}}$ is computed, a simple way to cross-check it is to verify that the integral
> 123: \begin{equation}
> 125:     \mathcal{V} &= 2\pi \int d\theta dr r^{2}(z,\theta)\sin(\theta) \epsilon_{\text{az}} \\
> 126:     &= 2\pi \int d\theta dz \frac{z^{2}}{\cos^{3}(\theta)}\sin(\theta)\epsilon_{\text{az}}
> 128: \end{equation}
> 129: matches the total volume of the decay volume.

## 1.2 Every factor, and exactly how SensCalc gets it

| Factor | Where defined | How SensCalc computes it | Closed form / table / fit / code |
|---|---|---|---|
| N_prod^(i) | l. 85: "the total number of FIPs produced by the process i, e.g., decays of mesons, direct production by proton-target collisions, etc." | Book-keeping of POT or luminosity × production cross-section × production branching ratios. "the user must also provide all the relevant quantities, such as the number of protons on target (or the integrated luminosity for LHC- and FCC-hh-based experiments), the target material, and the production cross-sections for secondary particles (mesons and W/Z/H-bosons)" (l. 248) | **No formula.** Code inputs. For h→SS **no cross-section, luminosity or BR is given anywhere in the paper.** |
| f^(i)(θ,E) | l. 87; built in FIP distribution.nb, l. 286–294 | Tabulated 2-D histogram in (θ,E) per mass and production channel, from tabulated parent spectra + boost (l. 286), or MadGraph+PYTHIA for DIS (l. 288) | **Table** {m,θ,E,f^(i)} (l. 290–294). **No analytic form anywhere.** |
| ε_az(θ,z) | Eq. (2), l. 89–91; algorithm l. 252 | MC over φ on a (θ,z) grid: "the notebook randomly generates N values of the azimuthal angle φ∈(−π,π), checks if the point {θ,z,φ} is inside the decay volume, and calculates the azimuthal acceptance ε_az(θ,z)=M/N" (l. 252) | Closed-form definition; numerical evaluation. For a body of revolution it is 0 or 1. |
| dP_dec/dz | Eq. (3), l. 94–96 | Closed form, exponential kept | **Closed form**, exact. |
| ε_dec(m,θ,E,z) | Eq. (4) (schematic), l. 98–103; algorithm l. 254–261 | Full per-event simulation of the FIP decay on the (m,θ,E,z,φ) grid, then average of a binary accept/reject | **No closed form.** Factorised only schematically into Br_vis × ε_geom × ε_other. |
| ε_rec | l. 104 | **Not computed.** l. 137: "the computation of ε_rec would require running the full simulation ... As such, it goes beyond the scope of the present semi-analytic approach." | Absent from the code; effectively 1 unless the user multiplies it in by hand. |
| The integrals | l. 323: "It then maps them to a logarithmic scale and interpolates them to obtain the functions entering Eq. (1)." | Interpolation of the tabulated ε_az, ε_dec, f on a log grid | Numerical. |

**Structural point** (l. 85, quoted above): N_prod counts *FIPs*, not events. For h→SS that means two scalars per Higgs. Eq. (1) contains no event-level "OR": the two scalars of one Higgs enter the sum twice, each with its own decay probability, azimuthal acceptance and decay-product acceptance. See §6.6.

---

# 2. f(θ,E) — the LLP kinematic distribution

## 2.1 What the paper says

l. 87:

> \item[--] $f^{(i)}(\theta,E)$ is the differential distribution of FIPs in polar angle and energy for FIPs produced through the process i.

l. 131:

> Depending on the production channel, evaluating the FIP distribution function $f^{(i)}(\theta,E)$ may require some external input. For instance, for FIPs that are produced directly in inelastic proton collisions, one needs to simulate $f^{(i)}(\theta,E)$ using, e.g., \texttt{PYTHIA~8} to account for showering and hadronization. For FIPs that are produced in the interactions of secondary particles, either in their decays or scattering with the material [...], the distribution of secondaries $f_{\text{secondary}}(\theta,E)$ is needed; nevertheless, once $f_{\text{secondary}}(\theta,E)$ has been computed, the distribution of FIPs can then be derived analytically without the need for external tools.

The recipe for any parent decay (this is the h→SS case), l. 286:

> Many production mechanisms require knowing the distributions of the parent particles at the given facility, such as mesons, heavy SM bosons, and photons --- including those produced in secondary interactions. We provide them as tabulated distributions in polar angle and energy, which we generate following the literature or just using available distributions from existing studies (see also Appendix~\ref{app:references} for a description of how we have generated the distributions of parent particles). Users may easily replace the included distributions with their own differential flux. With the distribution of parent particles at hand, we then derive the distribution of FIPs. If the FIPs are produced in decays, we compute their phase space in the rest frame of the parent particle and then boost it to the lab frame. In the case of 3-body decays, the phase space takes into account the matrix element of the process.

Output format, l. 290–294:

> The output of the second notebook is a tabulated distribution of the form
> \begin{equation}
> \{m,\theta,E,f^{(i)}\},
> \end{equation}
> where the last column is the value of the FIP distribution function for the given $(m,\theta,E)$ and the production mechanism $i$.

Where the parent spectra come from — Table IV, the row at line 545:

> 545:          $W,h,Z$   & -- & -- & \cite{Kling:2021fwx} & \cite{Kling:2021fwx} \\ \hline

with caption (l. 547): "List of the references used to generate, or directly take, the distributions of secondary particles that may produce FIPs." Kling:2021fwx is **FORESEE** (F. Kling, S. Trojanowski, PRD 104 (2021) 035012, arXiv:2105.07077 — main.bbl lines 113–119).

**So for the LHC the Higgs (θ,E) distribution is not generated by SensCalc: it is taken tabulated from FORESEE.** The paper gives no Higgs p_T distribution, no rapidity distribution, no fit, no parameters and no normalisation for it. That is a hard "**not stated in the paper**".

## 2.2 One- or two-dimensional? Factorised?

**It is integrated as a genuine two-dimensional distribution.** Nothing in the paper factorises f(θ,E):

- Eq. (1) integrates f^(i)(θ,E) over dθ dE simultaneously with ε_az(θ,z) and dP_dec/dz(θ,E).
- The stored object is a 2-D table {m,θ,E,f} in which both θ and E are grid axes (l. 290–294). The E grid is logarithmic, "starting from the FIP mass and ending with the maximal energy available at the given facility" (l. 252); the θ grid is refined where the FIP points at the detector (l. 252).
- Fig. 2 of the paper (Double-Distr-scalars-quartic.pdf, l. 316–317) is literally a 2-D histogram f_S(θ_S, E_S) with colour scale log10[f_S / (GeV⁻¹ rad⁻¹)] — per unit E and per unit θ.

The only factorisation in the paper is of the *acceptance*, and only in the long-lived limit — Eq. (7), l. 152: ε = ⟨ε_FIP⟩·⟨(γ²−1)^(−1/2)⟩·(Δz/cτ)·⟨ε_decay⟩. It is a product of *means*, equal to the integral only if the quantities are uncorrelated. The paper calls it a "representation ... particularly useful when discussing the impact of the geometry on the event rate" (l. 154); the actual yield always comes from Eq. (1).

## 2.3 Is there a functional form for the energy spectrum?

**No. Not stated in the paper.** There is no fit function, power law, log-normal, slope or parameter set for f(E) or f(θ) anywhere in the 591 lines. The only quantitative shape information is the example histogram of Fig. 2 — and that example is **not** the Higgs channel. Its own in-figure title reads "Dark scalars. Process B->X_s SS. m_S = 2 GeV", for the SPS with a molybdenum target (caption l. 317). Read off that figure (rendered at 200 dpi, values taken from the axes):

- E_S spans roughly 7 GeV to 150 GeV, with weight concentrated at 20–60 GeV;
- θ_S spans roughly 7×10⁻⁴ to 1.5×10⁻¹ rad;
- log10 f_S runs from about −2 (edges) to +2 (core);
- it is **not** separable: the E_S range shrinks as θ_S grows (the 150 GeV corner is reached only at θ_S ≲ 3×10⁻² rad).

For h→SS the paper shows no distribution at all. The left panel of the same figure is the ALP-gluon case, also at the SPS.

## 2.4 Typical boost, decay length, definition of "long-lived"

The paper never quotes a typical γ, ⟨E⟩ or decay length for any channel or any experiment. What it does say:

- l. 97: l_dec = cτ√(γ²−1), r = z/cos θ — the only kinematic ingredients of the decay probability.
- l. 133: "In the simplest case of a two-body decay into massless particles, the opening angle is Δθ_dec ≃ 2 arcsin(γ⁻¹) ... Because the detector angle is smallest at the beginning of the decay volume while the opening angle decreases as E_FIP⁻¹, ε_dec effectively imposes a cut from below on the FIP energy and the displacement of its decay position from the beginning of the decay volume."
- l. 139: the long-lived regime is defined by cτ⟨√(γ²−1)⟩ ≫ l_experiment, explicitly as a *limit* used for Eqs. (5)–(7).
- l. 152: ⟨(γ²−1)^(−1/2)⟩ is "the mean inverse p/m" — the mean of the inverse, not the inverse of the mean.

**MY INFERENCE (not in the paper).** This last point is the most likely place for a factor of order 10 to hide, and it is a property of the formula rather than of a missing input. In the long-lived limit the yield is proportional to ⟨1/(γβ)⟩ = ⟨m/p⟩, which for a broad spectrum is strictly larger than 1/⟨γβ⟩ = m/⟨p⟩ (Jensen). An implementation that computes the decay probability from a mean energy — "⟨E_S⟩ = 462 GeV, hence γ = 954, hence P_dec = Δz/(954 cτ)" — silently evaluates 1/⟨γβ⟩ instead of ⟨1/(γβ)⟩, and for a steeply falling E spectrum the two can differ by a large factor. SensCalc cannot make this mistake: Eq. (1) integrates f(θ,E) point by point with 1/(γβ) inside the integrand. Whatever else differs, the number to compare is not ⟨E_S⟩ but ⟨1/E_S⟩ taken over the same (θ,E) domain.


---

# 3. ε_dec — the decay-product acceptance

## 3.1 How it is computed

The paper states the schematic decomposition (Eq. 4, l. 100: Br_vis × ε_geom × ε_other) and then says, l. 135:

> To estimate $\epsilon_{\text{dec}}$ more accurately, by accounting for such factors as the experiment geometry, the presence of a dipole magnet, different FIP decay topologies (such as multi-body decays or decays into unstable particles), and various other selections imposed on the decay products, one can perform a \emph{separate} simulation (see details in Sec.~\ref{sec:senscalc}).

The procedure is described at lines 254–261:

> 254: Having the grid $(m,E,\{\theta,z,\phi\}_{\text{inside decay volume}})$, the notebook then simulates the FIP decays using the relevant decay channels and calculates the decay acceptance $\epsilon_{\text{dec}}(m,\theta,E,z)$ by averaging over these decays and $\phi$. The averaging over $\phi$ is already possible at this stage since the other quantities that determine the number of events~\eqref{eq:Nevents} do not depend on the azimuthal angle.

> 257: The decay channels implemented for each FIP are listed in Table~\ref{tab:implemented-models-decay}. By default, all decay channels that have at least two particles that may be reconstructed at the given experiment (it depends on the presence of the calorimeters and other limitations) are included in the computation. However, users may select their own list of decay channels. For 3-body decays, the distribution of the decay products is generated by taking into account both the phase space and the squared matrix element of the process. If the FIP decay products are short-lived, the routine decays them until only metastable particles are left. By default, those are $\gamma,e,\mu,K^{0}_{L},\pi^{\pm},K^{\pm}$.

> 261: Let us now discuss the computation of $\epsilon_{\text{dec}}$ in more detail. The main acceptance criterion is the requirement that the trajectories of at least two decay products with zero total electric charge are within the acceptance of the detector until its final plane. Decays into pure neutral final states (i.e., photons or $K^{0}_{L}$) are also included if a calorimeter is present. If the detector or decay volume includes a magnetic spectrometer, the components of the charged particles' coordinates and momenta are shifted by a kick right after the magnet to approximate the effect of the magnetic field. In addition to this geometric requirement, $\epsilon_{\text{dec}}$ may include various kinematic cuts on the visible final state particles resulting from the FIP's decays. The currently implemented cuts include cuts on the energy, transverse momentum, transverse impact parameter, and, for neutral particles in the calorimeter, their spatial separation. [...] Note that the acceptance criterion includes partially reconstructible states, i.e., the final states for which the FIP invariant mass cannot be reconstructed from the detected decay products.

So: **the actual decay products' angles and energies are used**, not a parametrisation. Each decay is generated with its matrix element (2-body isotropic, 3-body with |M|²), unstable products are decayed to metastable states, and a binary geometric+kinematic test is applied; ε_dec is the average of that binary test over the (m,θ,E,z,φ) grid.

## 3.2 "At least two charged particles" and a table of decay modes

- The criterion is "at least two decay products with **zero total electric charge**" (l. 261) — i.e. two oppositely charged tracks, not "two charged particles" and not "≥2 charged" without charge balance.
- Neutral final states (photons, K⁰_L) count only if a calorimeter is present (l. 261, and l. 103: "depending on the presence of a calorimeter ..., Br_vis may encompass only those states featuring at least two charged particles, or it may also include some neutral states such as photons and K⁰_L").
- There is **no table of "detectability" per decay mode**. What exists is Table III (l. 296–310), which lists the *implemented decay channels per FIP model*. For BC4/BC5 dark scalars (line 302):

> 302: BC4, BC5 & \cite{Boiarska:2019jym,Winkler:2018qyg} & \makecell{$ee,\mu\mu,\tau\tau$} & \makecell{ $\pi\pi,KK,4\pi, m\lesssim 2 \text{ GeV}$ \\ $c\bar{c},s\bar{s},b\bar{b},GG, m\gtrsim 2 \text{ GeV}$} \\ \hline

with the caption (l. 308) adding: "For dark scalars, we have included in their decay width into gluons the NLO correction from Ref.~\cite{Spira:1995rr}, which was previously missing in Ref.~\cite{Boiarska:2019jym}."

- The analytic short-cut for ε_dec is only a qualitative argument (l. 133, quoted in §2.4): opening angle 2 arcsin(1/γ) versus the detector's angular size, "ε_dec ≈ 0" if Δθ_dec ≳ Δθ_det, "ε_dec ≈ 1" otherwise. The paper does not use this to compute; it is offered as intuition.

## 3.3 ε_rec and ε_S/B

- ε_rec is **not computed by SensCalc at all**. l. 137:

> Finally, the computation of $\epsilon_{\text{rec}}$ would require running the full simulation, including the detector response. As such, it goes beyond the scope of the present semi-analytic approach. However, we believe that it is possible to perform an adequate pre-selection with the help of $\epsilon_{\text{dec}}^{\text{(other cuts)}}$ (for instance, by requiring the energy or $p_{\mathrm{T}}$ of the final state particles to exceed a threshold above which they are detected with high efficiency; see, e.g.,~\cite{Aberle:2839677}), such that, conditioned on this pre-selection, $\epsilon_{\text{rec}} \sim \mathcal{O}(1)$.

- **ε_S/B (the signal-preserving efficiency of the background rejection) does not appear in the SensCalc paper at all.** A full-text search of main.tex for "S/B", "background rejection" returns zero matches. It is a Gorkavenko-paper quantity (§7).
- Where ε_S/B does exist in the literature chain: Gorkavenko et al. apply ε_rec ≈ 0.4 **and** ε_S/B = 0.87 on top of SensCalc's ε_dec. If you apply only 0.4, your rate is 1/0.87 = 1.15× *higher* than theirs — this cannot explain a deficit.

---

# 4. The LHCb Downstream setup as implemented in SensCalc

## 4.1 What the SensCalc paper says: almost nothing

LHCb appears **once** in the whole paper, at line 222, in the list of implemented experiments:

> 222:         LHC   & \makecell{FASER/FASER2/FASER$\nu$/FASER$\nu 2$~\cite{FASER:2019aik,FASER:2020gpr,Feng:2022inv}\\ SND@LHC/advSND~\cite{SHiP:2020sos,Feng:2022inv}\\  FACET~\cite{Cerci:2021nlb}, MATHUSLA~\cite{MATHUSLA:2022sze}, CODEX-b~\cite{Aielli:2019ivi} \\ ANUBIS in the shaft and ceiling configurations~\cite{Bauer:2019vqk} \\ LHCb}  \\ \hline

Note that **LHCb is the only experiment in that table with no citation attached**, and there is no other sentence about it anywhere in the paper. There is no mention of the decay volume z-range, no angular range, no η range, no SciFi, no UT, no p or E cut, no magnet, no ε_rec value. The word "SciFi" does not occur in main.tex.

**Verdict: the SensCalc paper does not document its LHCb Downstream implementation. Not stated in the paper — at all.** The only public description is Appendix A of Gorkavenko et al. (see §7.1), which is consistent with your Table 2 but adds the details that matter.

## 4.2 Your specific question about the η range

From the Gorkavenko appendix (quoted in full in §7.1): the SensCalc LHCb decay volume is a **conical frustum covering 2 < η < 5**, azimuthally complete, from z = 1 m to z_max. So:

- the LLP **does** carry a 2 < η < 5 cut, implemented as ε_az (= 1 inside the frustum, 0 outside), not as an extra weight. Your 2 < η_S < 5 is therefore *not* an extra cut relative to SensCalc — good, that hypothesis is dead;
- Table 2's θ range ≈ (1.3, 260) mrad converts to η ∈ (2.03, 7.34) — because θ = 1.3 mrad corresponds to η = 7.34, not 5. So Table 2 (θ_min = 1.3 mrad) and the text (2 < η < 5, which would need θ_min = 13.5 mrad) are mutually inconsistent by a factor ~10 in θ_min. Your reading is right that "the table implies η up to ~7.3";
- the likely origin of 1.3 mrad: with z_max = 7.7 m, 1.3 mrad × 7.7 m ≈ 1 cm — i.e. a ~1 cm inner radius, whereas the SciFi beam-pipe hole quoted in the same appendix has R = 9 cm (≈ 11.7 mrad at z = 7.7 m). Whatever the provenance, the frustum that is actually implemented is the 2 < η < 5 one, so the θ_min = 1.3 mrad entry should be treated as an error or as a leftover from a different configuration.

## 4.3 Differences that would matter for your reproduction

Compared with your list (z 1–2.5 m, 2 < η_S < 5, two opposite-charge daughters reaching the SciFi, ε_rec = 0.4):

| Item | SensCalc (per Gorkavenko App. A) | Your implementation | Effect on your yield |
|---|---|---|---|
| z range | 1 m → 2.5 m for the "LHCb-downstream" option (z_max = 7.7 m exists as another option, "T-tracks-only") | 1 → 2.5 m | same |
| η range | frustum 2 < η < 5 | 2 < η_S < 5 | same |
| daughter requirement | "at least two decay products with zero total electric charge ... within the acceptance of the detector until its final plane" (l. 261), detector = SciFi parallelepiped 6.48 × 4.83 × 1.7 m with a 9 cm hole | two opposite-charge daughters reaching the SciFi envelope | approximately same |
| p/E cut | p > 5 GeV/c in Table 2, but **E > 5 GeV in the body text** ("it is required that the charged decay products have E > 5 GeV") | none quoted | you are *less* restrictive → you should get more, not fewer. For pions/kaons the two cuts differ by only a few %; for protons by ~8% |
| ε_rec | 0.4 (Gorkavenko), not from SensCalc | 0.4 | same |
| ε_S/B | 0.87 (Gorkavenko only) | — | yours 15% higher |
| magnet | kick after the dipole, ∫B dl = 4 T·m, z = 3.5–7.5 m | unclear | small |

**Conclusion for §4:** nothing in the LHCb geometry, as far as it is documented, can account for a factor of 10. Every documented difference either is neutral or works in the direction of giving *you more* events, not fewer.

---

# 5. The validation

## 5.1 What is validated and against what

Fig. 2 of the paper (files BC6-validation.pdf = left panel = HNL; BC1-validation.pdf = right panel = dark photon), caption lines 171:

> \caption{Comparison of the 90\% CL sensitivity of the SHiP experiment to heavy neutral leptons (\textbf{left panel}) and dark photons (\textbf{right panel}), obtained using Eq.~\eqref{eq:Nevents} within the framework of \texttt{SensCalc} and derived using the \texttt{FairShip} simulations~\cite{SHiP:2020vbd,SHiP:2018xqw}. The old ECN4 configuration of SHiP has been considered here.}

and the discussion, l. 164:

> In particular, Fig.~\ref{fig:validation} shows the comparison of the sensitivity of the SHiP experiment to heavy neutral leptons (HNLs) and dark photons obtained using Eq.~\eqref{eq:Nevents} with the sensitivity obtained by the SHiP collaboration using the \texttt{FairShip} simulation. In the case of dark photons, the slight differences in the sensitivity can be explained by the different elastic proton form factors used to describe the production probability. In the case of HNLs, the discrepancy at the upper bound follows from the monochromatic approximation of the HNL energy spectrum assumed when computing the sensitivity shown in the SHiP paper~\cite{SHiP:2018xqw}.

Second validation (the only one involving scalars), Fig. 3 (SHADOWS_scalar_sensitivity_CL90.pdf and MATHUSLA_scalar_sensitivity_CL90.pdf), l. 192:

> Fig.~\ref{fig:validation-off-axis} compares the sensitivity of SHADOWS and MATHUSLA to dark scalars as computed via Eq.~\eqref{eq:Nevents} and calculated independently by \texttt{SensMC}, a simple weight-based Monte-Carlo that we have implemented as described in Appendix~\ref{app:toy-mc} [...] In these calculations, we did not impose any cuts on the decay products apart from the geometric requirement $\epsilon_{\text{decay}}^{\text{geom}}$; therefore, the sensitivities shown are optimistic. The agreement between the two approaches is within 10--20\% depending on the scalar mass; the discrepancies may be explained by numeric differences in the total number of produced scalars, the sampling of the distribution of $B$ mesons (for MATHUSLA), and slightly different treatments of the decay chain of the scalar.

Third: ALPINIST, l. 383: "The predictions of \texttt{SensCalc} agree well with the results of \texttt{ALPINIST} (see Fig.~\ref{fig:comparison-alpinist}) and FORESEE (the comparison between the semi-analytic approach and FORESEE is discussed in Ref.~\cite{Ovchynnikov:2022its})."

**Claimed accuracy: 10–20%** (l. 192), repeated as "All the characteristic quantities produced by the two approaches mostly agree within 20\%" in the Fig. 3 caption (l. 199).

## 5.2 Reading the two figures

BC6-validation.pdf (left panel of Fig. 2), rendered at 200 dpi. Axes: x = m_N in GeV on a log scale from 0.1 to ~5; y = U_e² on a log scale from 10⁻¹⁰ to ~10⁻⁴. Blue solid = SensCalc; dark-red solid = the SHiP HNL paper [1811.00930]; grey filled region = excluded; grey curve = the old SHiP ECN4 reach. Reading values off the raster (tick-calibrated, ±20%):

| m_N (GeV) | SensCalc U_e² | [1811.00930] U_e² |
|---|---|---|
| 0.2 | ~1.7×10⁻⁸ | ~1.9×10⁻⁸ (curves cross near 0.15–0.2) |
| 0.5 | ~2.0×10⁻⁹ | ~2.5×10⁻⁹ |
| 1.0 | ~5×10⁻¹⁰ | ~6×10⁻¹⁰ |
| 3.0 | ~1.6×10⁻⁸ | ~1.9×10⁻⁸ |

The two curves track each other to within ~20–25% over the whole mass range, diverge visibly at the dip near m_N ≈ 2 GeV (where the HNL lifetime/decay-chain treatment differs) and above ~2.5 GeV, where the SensCalc curve turns up less steeply. This supports the paper's statement that "the discrepancy at the upper bound follows from the monochromatic approximation of the HNL energy spectrum".

BC1-validation.pdf (right panel of Fig. 2). Axes: x = m_V in GeV (log, 0.1–2); y = ε² (log, 10⁻¹⁸ to 10⁻⁵). Blue = SensCalc, dark-red = [2011.05115]. The two curves are almost on top of each other over 0.1–1 GeV (differences of a few per cent), separate slightly at the dip near m_V ≈ 0.8 GeV, and re-join above 1 GeV; the residual difference is what the paper attributes to the elastic proton form factor.

Method note: for these two plots the axes are logarithmic and the tick labels are cleanly spaced, so the calibration is reliable; the *visual* comparisons above agree with a proper vector-path extraction to within the quoted 20%.

## 5.3 Is there a validation for a Higgs-portal scalar specifically?

**No — not stated in the paper, and not present in the figures.** The paper's scalar validation (Fig. 3) is for *B-meson-produced* scalars at SHADOWS and MATHUSLA, with no cuts on the decay products beyond geometry, compared against SensMC, which itself "does not include the showering and hadronization of partonic decay products nor the dipole magnet's effect" (Fig. 3 caption, l. 199). Nothing anywhere validates the h→SS channel — not the Higgs spectrum, not the two-scalar book-keeping, not the boost. The only h→SS-related validation in the whole literature chain is Gorkavenko's Fig. 12, and that one validates ε_det(m_S, z_S) against RapidSim **for B-produced scalars** ("the event rate for the dark scalar mixed with the standard Higgs boson is analyzed", i.e. BC4, mixing only, no h→SS channel), with ε_rec set to 1.

---

# 6. Anything that could explain a factor ~10 at low mass — checklist

Each candidate you listed, with the verdict and the evidence. "Paper" = main.tex; "Gork." = Gorkavenko et al. EPJC 84 (2024) 608.

### 6.1 Energy cut versus momentum cut
**Paper:** SensCalc applies whatever the user configures; the implemented cut variables are "the energy, transverse momentum, transverse impact parameter, and, for neutral particles in the calorimeter, their spatial separation" (l. 261). No LHCb-specific value is given.
**Gork.:** internally inconsistent — Table 2 says "p > 5 GeV/c", the body says "it is required that the charged decay products have E > 5 GeV".
**Verdict: not a factor 10.** For pions/kaons at 5 GeV the two cuts differ by (E−p)/p ≈ m²/2p² ≈ 2%; for protons ≈ 8%. Only if the cut were applied as p > 5 vs E > 5 on a very soft spectrum would it matter, and the daughters of a 400-GeV scalar are nowhere near threshold.

### 6.2 Does the LLP or the daughters carry the η cut?
**Paper:** the only angular requirement on the LLP is ε_az (l. 88–92); the daughters carry ε_dec's geometric requirement (l. 261). For LHCb the LLP cut is built into the decay-volume geometry as a 2 < η < 5 frustum (§4.2), so it *is* on the LLP.
**Verdict: not a factor 10.** Your 2 < η_S < 5 is the same domain. See §4.2.

### 6.3 Are decays before z_min counted?
**Paper:** No. Eq. (3) is evaluated at r(z) inside the volume, so the survival factor exp[−r(z_min)/l_dec] is included and only decays with z ∈ (z_min, z_max) contribute. This is exact and cannot inflate a published curve.

### 6.4 Is the effective decay volume treated as flat instead of exponentially?
**Paper:** Eq. (1)+(3) keep the exponential; the flat form appears only in the long-lived limit, Eq. (6), where dP_dec/dz → 1/(cos θ cτ √(γ²−1)) (l. 146). The paper is explicit that this is a limit (l. 139).
**Gork.:** Eq. (5)–(7) of that paper are written "considering the limit when cτ⟨γ⟩ ≫ X_exp, where X_exp is the geometric size of the whole experiment", and Eq. (5) is labelled N_ev,lower.
**Verdict: a genuine candidate, and cheap to test.** The two treatments differ by the factor [1 − exp(−Δz/l_dec)]·exp(−z_min/l_dec) versus Δz/l_dec. With z_min = 1 m, Δz = 1.5 m and l_dec = γβcτ: if l_dec = 0.4 m the flat formula is ~40× the exact one; if l_dec = 4 m it is ~1.5×; if l_dec = 40 m the two agree to 2%. **Compute l_dec = (p_S/m_S)·cτ at the (m_S, θ²) points where you disagree.** If l_dec ≳ 3·z_min, the exponential is irrelevant and this is not your problem; if l_dec ~ 0.3–1 m, a code using Eq. (5) beyond its range would sit well above you, by exactly the kind of factor you are missing.

### 6.5 Is the long-lifetime approximation used for the whole curve?
**Paper:** no — it is used "in the limit of very long-lived FIPs" (l. 139) and the full Eq. (1) is what the code evaluates (l. 323).
**Gork.:** the lower branch is explicitly the long-lived limit (Eq. 5), and the upper branch is the short-lived regime (P_dec ≈ exp[−z_min m/(cτ p)], their Eq. 9 discussion).
**Verdict:** if you compare against the *lower* branch of the published BC5 curve, the two codes should agree once l_dec ≫ Δz; the branch where they would not is the upper one.

### 6.6 The two scalars from h→SS
**Paper:** N_prod is "the total number of **FIPs** produced by the process i" (l. 85). So h→SS gives **two** scalars per Higgs, each integrated independently in Eq. (1). There is no event-level OR and no "at least one scalar decays" logic anywhere.
**Verdict: a real factor 2 candidate.** If your implementation produces N_h × BR(h→SS) scalars (one per Higgs, or an event-level OR), you are short by exactly 2 relative to SensCalc. Note this is a factor that is *specific to a 2-body LLP-pair final state* and would be invisible in the B-meson channels, where one B gives one scalar.

### 6.7 Is the α² / θ² normalisation different from BR(h→SS)?
**Paper, l. 486:** "where θ is the mixing angle and α is the quartic coupling. By default, the sensitivity is evaluated assuming a constant branching ratio Br(h→SS) ∝ α²."
**Paper, l. 340:** "for dark scalars, one needs to specify the value of the branching ratio Br(h→SS), which is non-zero in the presence of the quartic coupling L ∝ hSS".
**Verdict: SensCalc never computes BR(h→SS) from α.** It takes BR(h→SS) as a free input (0.01 in your case) and scans the *mixing* θ to set the lifetime. Production through the quartic and the lifetime through the mixing are therefore decoupled by construction — there is no α² normalisation to get wrong, provided you also fix BR = 0.01 as an input and vary only the lifetime. If instead you derived BR(h→SS) from an α and then also varied α somewhere, you would introduce a factor (α/α₀)² that SensCalc does not have.

### 6.8 The one I would look at first: the low-mass BC5 curve is pure h→SS
**MY INFERENCE, with the argument written out.** In BC5 the quartic is fixed (BR(h→SS) = 0.01) while the mixing θ is what the sensitivity curve scans. Production from B mesons scales as θ² and the decay probability also scales as θ², so the B-induced rate scales as θ⁴, whereas the h→SS-induced rate scales only as θ² (production fixed, decay probability ∝ Γ ∝ θ²). At the position of the published lower branch at m_S ≈ 0.5–2 GeV (θ² ≈ 10⁻¹⁰–10⁻¹³, see §7.2) the B contribution is smaller than the Higgs one by many orders of magnitude (BR(b→sS) at θ² = 10⁻¹⁰ is ~10⁻¹⁴ against a fixed 10⁻²). **Therefore the entire low-mass end of the curve you are trying to reproduce is the h→SS channel.** Two consequences:

1. your factor-10 problem is not in the B/K channels at all — it must be in N_h × 2 × BR × ⟨ε_LLP⟩ × ⟨1/(γβ)⟩ × ⟨ε_det⟩ × ε_rec;
2. on that branch N_ev ∝ θ², so a factor 10 in rate shows up as **one full decade in θ²** — a very visible, easily diagnosed vertical shift of your curve, not a subtle one.

### 6.9 Summary ranking
1. **⟨1/(γβ)⟩ vs 1/⟨γβ⟩** (§2.4): if the decay probability was computed from a single mean energy of 462 GeV, the error is unbounded from below and a factor 10 is entirely plausible.
2. **The factor 2 for two scalars per Higgs** (§6.6): a clean, certain factor if you used an event-level OR.
3. **The exponential versus the long-lived approximation** (§6.4): check l_dec = γβcτ against z_min = 1 m.
4. **ε_det / Br_vis details** (§3): at m_S = 0.5–2 GeV Br_vis is the dominant part of ε_dec, and the modes that count change quickly with mass (ππ, KK, 4π, μμ, ee on; γγ, π⁰π⁰ off without a calorimeter).
5. **Not** the η range, **not** ε_S/B, **not** the p-vs-E cut, **not** decays before z_min.

---

# 7. Numbers to check your yield against

## 7.1 The exact LHCb implementation (from Gorkavenko App. A — outside the SensCalc paper)

Quoted because the SensCalc paper is silent (§4):

> For the decay volume, conical frustum covering pseudorapidities 2 < η < 5 and located in the longitudinal displacement z from zmin = 1 m to zmax = 7.7 m is considered, where the first SciFi layer is located. If the tracks must also intersect the UT, the size of the decay volume shrinks to zmax ≈ 2.5 m, which is the beginning of the UT. For the geometry of the SciFi layers, a parallelepiped with dimensions of 6.48 m × 4.83 m × 1.7 m with a hole of radius R = 9 cm to account for the beam pipe is used, following [11]. The magnetic field of the dipole magnet is extended from z = 3.5 m to z = 7.5 m, with the integrated field ∫B dl = 4 T·m.

> – LHCb-downstream
> – LHCb-downstream-T-tracks-only
> – LHCb-downstream-full.
> The first one corresponds to the setup considered in this paper— the one with the decay volume extending from z = 1 m to z = 2.5 m and SciFi as the detector.

and

> The vertex is reconstructed with the help of at least two tracks from decay products passing through both the UT and SciFi trackers. For the present study, only charged particles are considered detectable. Therefore, decays into solely neutral particles such as π0(→2γ), γ, K0_L are treated as invisible.

The validation of that implementation (their Fig. 12) is against RapidSim, with ε_rec = 1, for the geometric part only:
> the acceptance for the dark scalars to have 2 < η < 5 and the z-dependence of the pure geometric part of the decay product acceptance (i.e., with ε_rec = 1) is studied, which is defined as ε_det(mS, zS) ≡ ⟨f_LLP (dP_decay/dz) ε_decay⟩_{θ,E} / ⟨f_LLP (dP_decay/dz)⟩_{θ,E}, and results compared with the LHCb simulations. [...] As shown in Fig. 12, good agreement is obtained.

## 7.2 The published BC5 curve itself (Gorkavenko Fig. 7), digitised

Method: the EPJC figure is a vector plot; I rendered the panel at 420 dpi and digitised the blue pixels, calibrating x from the m_S tick marks (0.5, 1, 5, 10, 50 on a log axis) and y from the θ² ticks (10⁻⁷, 10⁻¹⁰, 10⁻¹³, 10⁻¹⁶). The solid and dashed blue curves are 25 fb⁻¹ and 300 fb⁻¹; where they were not separable I report the lower of the two. **Uncertainty ≈ ±0.3 dex**, dominated by the steepness of the curve at low mass.

| m_S (GeV) | θ², lower branch, 25 fb⁻¹ (solid) | θ², lower branch, 300 fb⁻¹ (dashed) |
|---|---|---|
| 0.5 | ≈ 1.5×10⁻¹⁰ | ≈ 8×10⁻¹¹ |
| 0.7 | ≈ 5×10⁻¹¹ | — |
| 1.0 | ≈ 4×10⁻¹¹ | ≈ 2×10⁻¹¹ |
| 1.5 | ≈ 9×10⁻¹¹ | — |
| 2.0 | ≈ 7×10⁻¹¹ | — |
| 3.0 | ≈ 4×10⁻¹¹ | — |
| 5.0 | ≈ 2×10⁻¹³ | ≈ 3×10⁻¹⁴ |

The shape to reproduce: a steep rise at m_S < 0.5, a plateau at θ² ≈ 3–10×10⁻¹¹ between 0.7 and 3 GeV, then a fall of ~2.5 decades between 3 and 5 GeV. The upper (short-lived) branch of the same band sits at θ² ≈ 10⁻⁶–10⁻⁸ over the same range. **If your curve is one decade above this one at 0.5–2 GeV, that is exactly "10× fewer events"** on a branch where N ∝ θ² (§6.8).

## 7.3 Worked decomposition of the formula, from the SensCalc paper's own figure

acceptances.pdf (Fig. 5 of the paper, l. 409–413) plots I0…I3 for SHiP, dark scalars BC4, cτ = 100 m; l. 407 defines them: I0 = total number of FIPs produced, I1 = "the number of FIPs intersecting the decay volume", I2 = "the number of FIPs decaying inside the decay volume", I3 = "the number of FIPs decays for which the decay products passed the decay acceptance". **Caption caveat (l. 412): "The overall normalization is arbitrary"** — the ratios are meaningful, the absolute values are not.

Method: vector paths extracted with PyMuPDF, colours mapped to I0 (blue), I1 (red), I2 (green), I3 (black), solid = SHiP, dashed = ANUBIS; axis calibrated from the major tick marks (x: 113.05 pt per unit m_S; y: 40.585 pt per decade, 10¹⁰ at y = 306.4).

| m_S (GeV) | I0 | I1 | I2 | I3 | I1/I0 = ⟨ε_LLP⟩ | I2/I1 | I3/I2 = ⟨ε_dec⟩ |
|---|---|---|---|---|---|---|---|
| 0.5 | 8.9×10¹⁴ | 1.6×10¹⁴ | 5.7×10¹¹ | 4.2×10¹¹ | 0.184 | 3.5×10⁻³ | 0.73 |
| 1.0 | 8.8×10¹⁴ | 1.7×10¹⁴ | 1.2×10¹² | 7.7×10¹¹ | 0.195 | 6.8×10⁻³ | 0.65 |
| 2.0 | 8.4×10¹⁴ | 2.2×10¹⁴ | 3.9×10¹² | 1.4×10¹² | 0.265 | 1.7×10⁻² | 0.35 |
| 3.0 | 7.0×10¹⁴ | 3.0×10¹⁴ | 7.6×10¹² | 4.0×10¹² | 0.426 | 2.6×10⁻² | 0.53 |
| 4.0 | 2.4×10¹⁴ | 1.4×10¹⁴ | 4.2×10¹² | 1.8×10¹² | 0.587 | 3.1×10⁻² | 0.44 |

Two things worth internalising from this table: (i) ⟨ε_dec⟩ is only 0.4–0.7 even when *nothing* but the geometric requirement is imposed — "ε_dec ≈ 1 because the opening angle is small" is not what the code actually finds; (ii) I2/I1 at cτ = 100 m is 10⁻³–10⁻², i.e. the decay probability is the dominant suppression, exactly as in your case.

## 7.4 Your own numbers, turned into a test

At m_S = 0.484 GeV with ⟨E_S⟩ = 462 GeV: γβ = p/m = √(462² − 0.484²)/0.484 = **954.5**. Over Δz = 1.5 m the decay probability in the long-lived regime is

> P_dec = Δz/(γβ·cτ) = 1.5/(954.5·cτ[m]) = **1.571×10⁻³ m / cτ[m]**,

so with cτ = 1 mm you get P_dec = 1.57 (already saturated), and with cτ = 1 m, P_dec = 1.57×10⁻³.

**A factor 10 in the yield at fixed cτ requires γβ = 95.4, i.e. ⟨E_S⟩ = 46.2 GeV.** That is the single number to keep in mind: either the effective energy entering your decay probability is ~46 GeV (not 462), or you are missing a factor 10 elsewhere in N_h·2·BR·ε.

Reference normalisation for the h→SS channel (13 TeV, LHCb Run 3 partial):

> N_h = σ(pp→h) × L = 55.45 pb × 25 fb⁻¹ = 5.545×10⁴ fb × 25 fb⁻¹ = **1.386×10⁶** Higgs bosons;
> N(h→SS) = 1.386×10⁴; N_scalars = 2× that = **2.77×10⁴** (factor 2 per l. 85);
> at 300 fb⁻¹: 3.33×10⁵ scalars.

So the yield is

> N_ev ≈ 2.77×10⁴ × ⟨ε_LLP⟩ × ⟨1/(γβ)⟩ × Δz/(cτ) × ⟨ε_det⟩ × ε_rec,

with ⟨ε_LLP⟩ the fraction of scalars inside 2 < η < 5 and 1 < z < 2.5 m. Feed in your own ⟨ε_LLP⟩, ⟨ε_det⟩ and your *acceptance-weighted* ⟨1/(γβ)⟩ (not 1/⟨γβ⟩) at m_S = 1 GeV and a cτ read off from the θ² of §7.2 — if you do not reproduce ≈ 2.3 events there, the discrepancy is in one of those four bracketed quantities, and the BC5 low-mass branch being Higgs-dominated (§6.8) tells you it is not in the B/K channels.

---

# Appendix A: method notes

- All figure values were obtained from the PDFs, not from the raster where possible. For acceptances.pdf, BC6-validation.pdf and BC1-validation.pdf the axis calibration used the major tick marks emitted as separate short vector segments, then the curves were extracted by stroke colour (and, for acceptances.pdf, by dash pattern, since solid and dashed curves share a colour). For Gorkavenko Fig. 7 only a raster was available at the needed crop, so the blue channel was thresholded in a 420-dpi render and the tick positions taken from the same page's vector content.
- The rendered PNGs and the extraction scripts are in the subfolder _figs (render.py, digitize_acc2.py, digitize_val2.py, digit_fig7.py, crop.py, gork_figs.py) if you want to redo or refine any number.
- Nothing here is taken from the web; both papers were read from the local PDFs/LaTeX.

# Appendix B: what the SensCalc paper does **not** contain (explicit list)

For the avoidance of doubt, each of these is absent from main.tex, not merely hard to find:

1. Any LHCb parameter (z range, θ/η range, SciFi, UT, magnet, p or E cut) — §4.
2. Any Higgs p_T, rapidity or energy spectrum; any fit; any parameter — §2.3.
3. Any h→SS cross-section, luminosity, or BR number — §1.2.
4. Any value of ε_rec, for any experiment — §3.3.
5. Any table of "detectability" per decay mode (only a list of implemented channels, Table III) — §3.2.
6. Any validation of the Higgs-portal h→SS channel — §5.3.
7. Any statement about event-level combinatorics when one parent produces two LLPs — §6.6.
8. Any statement of ε_S/B — §3.3.