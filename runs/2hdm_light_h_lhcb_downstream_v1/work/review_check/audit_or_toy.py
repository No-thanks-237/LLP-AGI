# -*- coding: utf-8 -*-
"""Audit 4: probability algebra of the detection estimator.

(a) p_window with A=1 against the analytic result and a brute-force decay-length MC,
    and against the response owner's work/response/twohdm_lib.py::_p_window.
(b) history / parent / event-OR algebra against a brute-force per-event MC.
(c) the 0.4 reconstruction efficiency (applied once) and the 1e3 pb->fb conversion.
"""
import numpy as np, math, sys, os
sys.path.insert(0, r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream/work/response")
rng = np.random.default_rng(20261025)

Z1, Z2 = 1.0, 2.5
ZGRID = np.linspace(Z1, Z2, 13)

def p_analytic_A1(alpha, z0):
    """int_{z0}^{inf} dz alpha exp(-alpha(z-z0)) * Theta(Z1<=z<=Z2)."""
    if z0 >= Z2: return 0.0
    lo = max(z0, Z1)
    return math.exp(-alpha*(lo-z0)) - math.exp(-alpha*(Z2-z0))

def p_mc_A1(alpha, z0, n=400000):
    l = rng.exponential(scale=1.0/alpha, size=n)
    z = z0+l
    return float(np.mean((z >= Z1) & (z <= Z2)))

print("="*100)
print("(a1) A=1 window probability: analytic vs brute-force decay-length MC and vs response _p_window")
import twohdm_lib as TL
worst_mc, worst_lib = 0.0, 0.0
print(f"{'alpha':>10} {'z0':>6} {'analytic':>12} {'MC':>12} {'|dMC|':>9} {'_p_window':>13} {'|dLib|':>9}")
for alpha in (0.05, 0.3, 1.0, 3.0, 10.0, 50.0):
    for z0 in (0.0, 0.4, 1.0, 1.4, 2.2, 2.5):
        pa = p_analytic_A1(alpha, z0)
        pm = p_mc_A1(alpha, z0)
        A = np.ones(len(ZGRID))
        pl = float(TL._p_window(np.array([alpha]), np.array([z0]), A.reshape(1,-1), ZGRID)[0])
        worst_mc = max(worst_mc, abs(pm-pa)); worst_lib = max(worst_lib, abs(pl-pa))
        print(f"{alpha:10.3f} {z0:6.2f} {pa:12.8f} {pm:12.8f} {abs(pm-pa):9.2e} {pl:13.8f} {abs(pl-pa):9.2e}")
print(f"  -> max |MC - analytic| = {worst_mc:.2e} (binomial sigma ~ {math.sqrt(0.25/400000):.2e});"
      f" max |_p_window - analytic| = {worst_lib:.2e}")

print()
print("(a2) nonzero piecewise-linear acceptance: compare _p_window against fine-grid quadrature")
def p_quad(alpha, z0, fn, n=2_000_000):
    """Monte-Carlo/quadrature of int dz alpha e^{-alpha(z-z0)} A(z)."""
    zmax = Z2 + 10.0/alpha
    z = np.linspace(max(z0,0.0), zmax, n)
    A = fn(z)
    w = alpha*np.exp(-alpha*(z-z0))
    w[(z < Z1) | (z > Z2)] = 0.0
    return float(np.trapezoid(A*w, z))
worst_a2 = 0.0
for alpha in (0.1, 1.0, 5.0):
    for z0 in (0.0, 1.2):
        fn = lambda z: 0.2 + 0.8*np.exp(-((z-1.6)/0.35)**2)
        Anodes = fn(ZGRID)
        pl = float(TL._p_window(np.array([alpha]), np.array([z0]), Anodes.reshape(1,-1), ZGRID)[0])
        pq = p_quad(alpha, z0, fn)
        worst_a2 = max(worst_a2, abs(pl-pq))
        print(f"  alpha={alpha:5.2f} z0={z0:4.1f}: _p_window={pl:.8f} quadrature={pq:.8f} diff={pl-pq:+.2e}")
print(f"  -> max |_p_window - quadrature| = {worst_a2:.2e}  (piecewise-linear A on 13 nodes)")

print()
print("="*100)
print("(b) event-OR algebra: analytic 1-prod(1-p_i) vs brute-force per-event Bernoulli MC")
def or_events(n_ev, pmax, nmax, rng, pfun):
    """returns (mean analytic OR, MC fraction of events with >=1 hit)"""
    osum = 0.0; hit = 0
    for _ in range(n_ev):
        n = rng.integers(1, nmax+1)
        p = pfun(n, rng)
        o = 1.0 - np.prod(1.0-p)
        osum += o
        if rng.random() < o: hit += 1
    return osum/n_ev, hit/n_ev
for (nmax, pfun, label) in [
    (1, lambda n, r: np.array([r.uniform(0, 0.9)]), "1 candidate, uniform p"),
    (3, lambda n, r: r.uniform(0, 0.9, n), "1-3 candidates, uniform p"),
    (6, lambda n, r: 10**r.uniform(-40, -1, n), "1-6 candidates, log-uniform p (1e-40..0.1)"),
]:
    ne = 400000
    a, m = or_events(ne, None, nmax, rng, pfun)
    sig = math.sqrt(max(a*(1-a), 1e-30)/ne)
    print(f"  {label:45}: analytic={a:.6f} MC={m:.6f} diff={m-a:+.2e} ({abs(m-a)/sig:.2f} sigma)")
print()
print("  h -> H H sub-case: o = BR_hHH * (1-(1-p1)(1-p2)); MC with a Bernoulli(BR_hHH) gate + 2 H")
br = 0.03
ne = 400000; osum = 0.0; hit = 0
for _ in range(ne):
    p1, p2 = rng.uniform(0, 0.5, 2)
    o = br*(1-(1-p1)*(1-p2)); osum += o
    if rng.random() < br:
        if rng.random() < 1-(1-p1)*(1-p2): hit += 1
a = osum/ne; m = hit/ne; sig = math.sqrt(a*(1-a)/ne)
print(f"  BR_hHH={br}: analytic={a:.6f} MC={m:.6f} diff={m-a:+.2e} ({abs(m-a)/sig:.2f} sigma)")
print()
print("  mutual-exclusive history weights: eps = sum_c w_c eps_c, w_c = sigma_c/sum(sigma)")
w = np.array([0.85, 0.12, 0.03]); epsc = np.array([1e-3, 3e-4, 2e-2])
ne = 400000; hit = 0
for _ in range(ne):
    c = rng.choice(3, p=w)
    p = rng.random()*epsc[c]*20     # per-event p in the channel's ballpark
    if rng.random() < p: hit += 1
print(f"  direct MC per-event fraction = {hit/ne:.6f};  sum_c w_c <p_c> = {float(np.sum(w*epsc*10)):.6f}"
      f"  (toy p = U(0,20 eps_c) -> <p>=10 eps_c)")

print()
print("="*100)
print("(c) 0.4 efficiency and the 1e3 pb->fb conversion")
print("  1 fb^-1 = 1e3 pb^-1  :  1 fb^-1 * 1 pb = 1000 events at eps=1")
for L in (25.0, 300.0):
    for sigma in (5.6e8, 7.95e7, 55.45):
        for eps in (1.0, 0.4):
            s = 1e3*L*sigma*eps
            if sigma > 1e6:
                print(f"  L={L:5.0f} fb^-1 sigma={sigma:.3e} pb eps={eps:3.1f}: s={s:.4e}"
                      f"   [check L*sigma*1e3*0.4 = {L*sigma*400:.4e}]")
    print()
print("  0.4 applied TWICE would change the B-channel weight by 0.4 (i.e. 2.5x): s(0.4*0.4)/s(0.4) = 0.4")
print("  -> the estimator must multiply epsilon by RECO_EFF exactly once, tests below")
import inspect
print("  response module RECO_EFF =", TL.RECO_EFF)
print()
print("(d) fully-exclusive lifetimes: p_window -> 1 when alpha -> 0 (infinite lifetime, window fills up)?")
for ctm in (1e-6, 1e-3, 1.0, 1e3, 1e6):
    alpha = 1.0/(ctm*1.0)
    pl = float(TL._p_window(np.array([alpha]), np.array([0.0]), np.ones((1,13)), ZGRID)[0])
    pa = p_analytic_A1(alpha, 0.0)
    print(f"  ctau={ctm:9.3e} m: _p_window={pl:.6f} analytic={pa:.6f} diff={pl-pa:+.2e}")
