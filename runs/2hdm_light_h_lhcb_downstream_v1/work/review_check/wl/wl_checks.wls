(* =====================================================================
   wl_checks.wls -- independent symbolic/numeric checks for the T5 audit
   of the 2HDM light-H benchmark (arXiv:2212.06186v1).
   Written by the model/audit owner.  Host Wolfram 15 (licensed).
   Checks:
     C1  App.D g0,g1,g2: x -> xH limit, x -> 0 scaling, tanb-scaling of xi_H
     C2  xi_H -> 0 as 1/tanb for tanb -> inf
     C3  g_C exact expansion in c = cos(beta-alpha) vs -mh^2 c^2/(4v)
     C4  Br(h->HH) vs 4700 (g_C/v)^2
     C5  App.A loop form factors vs the standard (Djouadi) ones, factor 1/2
     C6  SM ground-truth test of the Gamma(gaga) and Gamma(gg) normalisations
   ===================================================================== *)
$Assumptions = {};
mW = 80.385; v = 246; mHsm = 125; mHp = 600; GF = 1.1663787*^-5;
mh = mHsm;
Print["=== C1: App.D g0,g1,g2 ==="];
g0[x_, xH_, cb_] := -cb^2 (3 xH^2 - 4 xH x + x^2 - 2 x (2 xH - x) Log[xH/x])/(16 (xH - x)^3);
g1[x_, xH_, cb_] := -3/4 + cb^2 x (5 xH^2 - 8 xH x + 3 x^2 - 2 xH (2 xH - x) Log[xH/x])/(4 (xH - x)^3);
X1[x_, xH_] := -1/4 ((xH/(xH - x) - 6/(x - 1)^2 + 3) - xH (3 xH - 2 x)/(xH - x)^2 Log[xH]
     + (xH (3 xH - 2 x)/(xH - x)^2 + 3 (x + 1)/(x - 1)^3) Log[x]);
X2[x_, xH_] := x (5 xH - 3 x)/(4 (xH - x)^2) - xH x (2 xH - x)/(2 (xH - x)^3) Log[xH/x];
g2[x_, xH_, cb_] := cb X1[x, xH] + cb^3 X2[x, xH];
(* limit x -> xH : x = xH(1-eps) *)
g0lim = Limit[g0[xH (1 - e), xH, cb], e -> 0, Assumptions -> xH > 0];
g1lim = Limit[g1[xH (1 - e), xH, cb], e -> 0, Assumptions -> xH > 0];
g2lim = Limit[g2[xH (1 - e), xH, cb], e -> 0, Assumptions -> xH > 0];
Print["g0(x->xH) = ", Simplify[g0lim]];
Print["g1(x->xH) = ", Simplify[g1lim]];
Print["g2(x->xH) = ", Simplify[g2lim]];
Print["g0(x->0)  = ", Limit[g0[x, xH, cb], x -> 0, Assumptions -> xH > 0]];
Print["g1(x->0)  = ", Limit[g1[x, xH, cb], x -> 0, Assumptions -> xH > 0]];
Print["g2(x->0)  = ", Limit[g2[x, xH, cb], x -> 0, Assumptions -> xH > 0]];
Print["series g0 in 1/xH (xH->inf): ", Series[g0[x, xH, cb] /. xH -> 1/ee, {ee, 0, 1}] // Normal // Simplify];
Print["series g1 in 1/xH (xH->inf): ", Series[g1[x, xH, cb] /. xH -> 1/ee, {ee, 0, 1}] // Normal // Simplify];

Print[""];
Print["=== C2: xi_H(bs) scaling with tanb ==="];
cotbOf[tb_] := 1/tb; cOf[tb_] := 1/tb; sOf[tb_] := Sqrt[1 - cOf[tb]^2];
xCk = {2.2*^-3, 1.28, 173.1}^2/mW^2;  mQk = {2.2*^-3, 1.28, 173.1};
(* |V_us V_ub|, |V_cs V_cb|, |V_ts V_tb| : the CKM products that enter xi_H^{bs} *)
Vckm = {0.2243*0.00382, 0.97369*0.0408, 0.040122*0.99916};
lamHHHp[mH_, tbIn_] := -(2 mHp^2 - mH^2) cOf[tbIn]/v;
brk[mH_, tb_, k_] := With[{c = cOf[tb], s = sOf[tb], cb = cotbOf[tb], xH = mHp^2/mW^2},
   Module[{g0v, g1v, g2v},
    {g0v, g1v, g2v} = {g0[xCk[[k]], xH, cb], g1[xCk[[k]], xH, cb], g2[xCk[[k]], xH, cb]};
    g1v c - g2v s - g0v (2 v/mW^2) lamHHHp[mH, tb]]];
xiB[mH_, tb_] := -(4 Sqrt[2] GF)/(16 Pi^2) Sum[Vckm[[k]] mQk[[k]]^2 brk[mH, tb, k], {k, 1, 3}];
Print["xi_bs(tanb=100,mH=1)   = ", xiB[1, 100]];
Print["xi_bs*tanb at tanb=1e2,1e3,1e4 : ",
  Table[{tb, xiB[1, tb]*tb}, {tb, {100, 1000, 10000}}] // N];
Print["Ratio xi(1e3)/xi(1e2)*10 = ", N[(xiB[1, 1000]/xiB[1, 100])*10], "  (should be ~1)"];
Print["xi_bs(tanb=100,mH=1) numeric = ", N[xiB[1, 100]]];
Print["xi_bs(tanb=5,mH=1)   numeric = ", N[xiB[1, 5]]];
Print["ratio xi(5)/xi(100)*100/5 = ", N[(xiB[1,5]/xiB[1,100])*100/5], " (should be ~1)"];
Print["g0,g1,g2 at tanb->inf, x_t: ",
  N[Limit[{g0[xCk[[3]], mHp^2/mW^2, 1/tb], g1[xCk[[3]], mHp^2/mW^2, 1/tb],
      g2[xCk[[3]], mHp^2/mW^2, 1/tb]}, tb -> Infinity]]];

Print[""];
Print["=== C3: g_C exact vs approximation ==="];
gC[mH_, tb_] := With[{c = cOf[tb], s = sOf[tb], t2b = (1 - tb^2)/(2 tb)},
   s/(2 v) ((mH^2 - mh^2) (2 t2b s c - c^2 + s^2) - mH^2)];
Print["g_C(1,100)  = ", gC[1, 100] // N];
Print["-mh^2 c^2/(4v) at tanb=100 = ", (-mh^2 (1/100)^2/(4 v)) // N];
Print["ratio = ", (gC[1, 100]/(-mh^2 (1/100)^2/(4 v))) // N];
(* exact small-c expansion *)
expGc = Series[gC[mH, 1/c] /. tb -> 1/c, {c, 0, 4}] // Normal // Simplify;
Print["g_C(c) expansion = ", Collect[expGc, c]];
Print["i.e. g_C = s/(2v)[ mh^2 c^2/2 - mH^2(1+c^2/2) ] + O(c^4):  ",
  Simplify[expGc - Sqrt[1 - c^2]/(2 v) (mh^2 c^2/2 - mH^2 (1 + c^2/2))]];
Print["vanishing c vs mH (tanb=100): c_van = ", Sqrt[2] mH/Sqrt[mh^2 - mH^2] /. mH -> 1, " -> tanb = ", 1/Sqrt[2 mH^2/(mh^2 - mH^2)] /. mH -> 1];
Print["mH where paper c=1/tanb gives vanishing: mH = ", mh/(Sqrt[2] 100)];

Print[""];
Print["=== C4: Br(h->HH) vs 4700 (g/v)^2 ==="];
gamHsm = 4.07*^-3;
g3[mH_, tb_] := -2 gC[mH, tb];
gamHH[mH_, tb_] := g3[mH, tb]^2/(32 Pi mh) Sqrt[1 - 4 mH^2/mh^2];
brHH[mH_, tb_] := gamHH[mH, tb]/(gamHsm + gamHH[mH, tb]);
Print["Br(h->HH)(mH=0.021,tanb=100) = ", brHH[0.021, 100] // N];
Print["4700 (g_C/v)^2                = ", 4700 (gC[0.021, 100]/v)^2 // N];
Print["ratio (exact/approx)          = ", brHH[0.021, 100]/(4700 (gC[0.021, 100]/v)^2) // N];
Print["Br(h->HH)(mH=1,tanb=100) = ", brHH[1, 100] // N, "  (DESIGN: ~1e-8)"];

Print[""];
Print["=== C5: App.A form factors vs Djouadi (factor 1/2) ==="];
fLoop[t_] := If[t <= 1, ArcSin[Sqrt[t]]^2, -1/4 (Log[(1 + Sqrt[1 - 1/t])/(1 - Sqrt[1 - 1/t])] - I Pi)^2];
A0p[t_] := -1/2 (t - fLoop[t])/t^2;
A12p[t_] := (t + (t - 1) fLoop[t])/t^2;
A1p[t_] := -1/2 (2 t^2 + 3 t + 3 (2 t - 1) fLoop[t])/t^2;
(* Djouadi hep-ph/0503172 forms (the ones with the usual 4/3, -7 asymptotic limits) *)
A12D[t_] := 2 (t + (t - 1) fLoop[t])/t^2;
A1D[t_] := -(2 t^2 + 3 t + 3 (2 t - 1) fLoop[t])/t^2;
A0D[t_] := -(t - fLoop[t])/t^2;
Print["A12^paper/A12^Djouadi at tau=0.01,0.1,0.5 : ",
  Table[N[A12p[t]/A12D[t]], {t, {0.01, 0.1, 0.5}}]];
Print["A1^paper/A1^Djouadi  at tau=0.01,0.1,0.5 : ",
  Table[N[A1p[t]/A1D[t]], {t, {0.01, 0.1, 0.5}}]];
Print["A0^paper/A0^Djouadi  at tau=0.01,0.1,0.5 : ",
  Table[N[A0p[t]/A0D[t]], {t, {0.01, 0.1, 0.5}}]];
Print["tau->0 limits: A12p=", Limit[A12p[t], t -> 0, Direction -> "FromAbove"] // N,
      " A1p=", Limit[A1p[t], t -> 0, Direction -> "FromAbove"] // N,
      " (paper) ; A12D=", Limit[A12D[t], t -> 0, Direction -> "FromAbove"] // N,
      " A1D=", Limit[A1D[t], t -> 0, Direction -> "FromAbove"] // N];

Print[""];
Print["=== C6: SM ground-truth test of the width normalisations ==="];
alphaEm = 1/137; alphas125 = 0.1133;
tauW = mh^2/(4 mW^2);
(* paper Eq.(hgaga): Gamma = GF alpha^2 m^3/(32 Sqrt[2] Pi^3) |xi_gamma|^2,
   xi_gamma = sum Nc Q^2 xi_f A12 + xi_V A1  (paper normalisation) *)
xiGag = 4/3 A12p[mh^2/(4 1.28^2)] + 1/3 A12p[mh^2/(4 4.18^2)] + 4/3 A12p[mh^2/(4 173.1^2)]
       + 1/3 A12p[mh^2/(4 0.096^2)] + 4/3 A12p[mh^2/(4 0.0022^2)] + 1/3 A12p[mh^2/(4 0.0047^2)]
       + A12p[mh^2/(4 0.10566^2)] + A12p[mh^2/(4 1.7768^2)] + A1p[tauW];
gamGagPaper = GF alphaEm^2 mh^3/(32 Sqrt[2] Pi^3) Abs[xiGag]^2;
Print["Gamma(gaga) paper formula, xi=1 : ", gamGagPaper*1000, " MeV   (PDG SM 9.28 keV)"];
Print["   -> ratio to PDG 9.28 keV      : ", (gamGagPaper*1000)/0.00928];
xiGagStd = 4/3 A12D[mh^2/(4 1.28^2)] + 1/3 A12D[mh^2/(4 4.18^2)] + 4/3 A12D[mh^2/(4 173.1^2)]
       + 1/3 A12D[mh^2/(4 0.096^2)] + 4/3 A12D[mh^2/(4 0.0022^2)] + 1/3 A12D[mh^2/(4 0.0047^2)]
       + A12D[mh^2/(4 0.10566^2)] + A12D[mh^2/(4 1.7768^2)] + A1D[tauW];
gamGagStd = GF alphaEm^2 mh^3/(128 Sqrt[2] Pi^3) Abs[xiGagStd]^2;
Print["Gamma(gaga) standard formula, xi=1 : ", gamGagStd*1000, " MeV  -> ratio ", (gamGagStd*1000)/0.00928];
Print["xi_gamma(paper)/xi_gamma(std) = ", Abs[xiGag/xiGagStd], "  (should be 1/2)"];
Print["ratio of the two Gamma(gaga)   = ", gamGagPaper/gamGagStd, "  (should be 1)"];
Print[""];
xiG = Sum[A12p[mh^2/(4 mm^2)], {mm, {1.28, 4.18, 173.1, 0.096, 0.0022, 0.0047}}];
gamGgCode = alphas125^2 mh^3/(32 Pi^3 v^2) Abs[xiG]^2;
gamGgPaper = GF alphas125^2 mh^3/(36 Sqrt[2] Pi^3) Abs[3/2 xiG]^2;
gamGgCorrect = GF alphas125^2 mh^3/(36 Sqrt[2] Pi^3) Abs[2 xiG]^2;
Print["Gamma(gg): code form (32 Pi^3 v^2)        = ", gamGgCode*1000, " MeV"];
Print["Gamma(gg): paper form (3/2, 36 Sqrt2 Pi^3)= ", gamGgPaper*1000, " MeV"];
Print["Gamma(gg): paper form with 2 instead of 3/2 = ", gamGgCorrect*1000, " MeV"];
Print["PDG SM Gamma(h->gg) = 0.3358 MeV ; ratios: code ", (gamGgCode*1000)/0.3358,
      "  paper ", (gamGgPaper*1000)/0.3358, "  paper*16/9 ", (gamGgPaper*1000)/0.3358*16/9];
Print["Xi_g values: sum_paper(A12) = ", xiG, "  (3/2 sum) = ", 3/2 xiG, "  (2 sum) = ", 2 xiG,
      "  (4/3 limit for tau->0) = ", 4/3];
Print[""];
Print["=== DONE ==="];
