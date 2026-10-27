mW = 80.385; v = 246; mHsm = 125; mHp = 600; GF = 1.1663787*^-5; mh = mHsm;
g0[x_, xH_, cb_] := -cb^2 (3 xH^2 - 4 xH x + x^2 - 2 x (2 xH - x) Log[xH/x])/(16 (xH - x)^3);
g1[x_, xH_, cb_] := -3/4 + cb^2 x (5 xH^2 - 8 xH x + 3 x^2 - 2 xH (2 xH - x) Log[xH/x])/(4 (xH - x)^3);
X1[x_, xH_] := -1/4 ((xH/(xH - x) - 6/(x - 1)^2 + 3) - xH (3 xH - 2 x)/(xH - x)^2 Log[xH]
     + (xH (3 xH - 2 x)/(xH - x)^2 + 3 (x + 1)/(x - 1)^3) Log[x]);
X2[x_, xH_] := x (5 xH - 3 x)/(4 (xH - x)^2) - xH x (2 xH - x)/(2 (xH - x)^3) Log[xH/x];
g2[x_, xH_, cb_] := cb X1[x, xH] + cb^3 X2[x, xH];
tb = 100.; c = 1/tb; s = Sqrt[1 - c^2]; cb = 1/tb; xH = mHp^2/mW^2; xt = 173.1^2/mW^2;
Print["xt=", xt, " xH=", xH];
Print["g0,g1,g2(xt)=", {g0[xt, xH, cb], g1[xt, xH, cb], g2[xt, xH, cb]}];
lam = -(2 mHp^2 - 1)/c/v;
Print["lambda_HH+H- = ", lam];
br = g1[xt, xH, cb] c - g2[xt, xH, cb] s - g0[xt, xH, cb] (2 v/mW^2) lam;
Print["bracket_t = ", br];
Print["xi_t_term = ", -(4 Sqrt[2] GF)/(16 Pi^2) 0.0400887^2 173.1^2 br];
Print["total xi (approx, t only) = ", -(4 Sqrt[2] GF)/(16 Pi^2) 0.0400887^2 173.1^2 br];
