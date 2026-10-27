# -*- coding: utf-8 -*-
"""Audit 2: the FORESEE 2HDM_H B -> X_s H production formula vs our implementations.

FORESEE notebook (cell 7), for B+ (511), B0 (521) and their antiparticles:
  br = 9.68544e-13*mB/Gamma_B*(1-(m/mB)^2)^2
       *(-3/4/tb -0.1415/tb^3 -(1.0386/tb-0.1415/tb^3)*sqrt(1-1/tb^2) -0.608446664/tb^3)^2
with tb = tanbeta.  Compare with:
  (a) our independent implementation audit_2hdm.br_B_to_XsH  (paper Eq. after (Lfceff_phi), m_b=4.18)
  (b) the production module work/production/twohdm_production.py
"""
import math, json, sys, importlib.util, statistics
import audit_2hdm as AU

MB_FOR = 5.27934        # PDG-ish B+ mass used by FORESEE (masses('511'))
GAM_B_FOR = 6.582119569e-25/1.638e-12   # Gamma_B+ [GeV]
MB0_FOR = 5.27965
GAM_B0_FOR = 6.582119569e-25/1.519e-12

def br_foresee_B(m, tanb, meson='B+'):
    mB, GB = (MB_FOR, GAM_B_FOR) if meson == 'B+' else (MB0_FOR, GAM_B0_FOR)
    tb = tanb
    brk = (-3.0/4.0/tb - 0.1415/tb**3
           - (1.0386/tb - 0.1415/tb**3)*math.sqrt(max(0.0, 1.0-1.0/tb**2))
           - 0.608446664/tb**3)
    return 9.68544e-13*mB/GB*(1.0-(m/mB)**2)**2*brk**2

def bracket_foresee(tanb):
    tb = tanb
    return (-3.0/4.0/tb - 0.1415/tb**3
            - (1.0386/tb - 0.1415/tb**3)*math.sqrt(max(0.0, 1.0-1.0/tb**2))
            - 0.608446664/tb**3)

def bracket_ours(f, tanb):
    """our factor f(tanb) = xi_bs / (prefactor * m_t^2 * Vts Vtb)  (t term dominates);
    returns the t-quark bracket so that f ~ bracket."""
    c, s = AU.mixing(tanb); cotb = 1.0/tanb
    g0, g1, g2 = AU.g012(AU.MQ['t']**2/AU.MW**2, AU.MHC**2/AU.MW**2, cotb)
    return g1*c - g2*s - g0*(2*AU.V_VEV/AU.MW**2)*AU.lam_HHpHm(f, tanb)

# production module
spec = importlib.util.spec_from_file_location(
    "prod", r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream/work/production/twohdm_production.py")
prod = importlib.util.module_from_spec(spec); spec.loader.exec_module(prod)

print("="*100)
print(f"FORESEE constants: mB={MB_FOR} Gamma_B+={GAM_B_FOR:.6e} GeV ; coefficient 9.68544e-13*mB/Gamma_B = "
      f"{9.68544e-13*MB_FOR/GAM_B_FOR:.6f}")
print()
print(f"{'mH':>7} {'tanb':>8} | {'br_FORESEE':>12} {'br_ours':>12} {'br_prod':>12} | "
      f"{'F/ours':>9} {'prod/ours':>9} | {'brk_F':>10} {'brk_ours':>10}")
grid_m = [0.021, 0.1396, 0.27, 0.5, 1.0, 2.0, 3.0, 3.7, 4.0, 4.18, 4.5, 5.0]
grid_t = [5.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 1e4, 1e5]
rows = []
for m in grid_m:
    for t in grid_t:
        bf = br_foresee_B(m, t); bo = AU.br_B_to_XsH(m, t)[0]
        bp = prod.br_B_to_XsH(m, t, meson='B+')
        bpf = bp['total']
        bpf_f = prod.br_B_to_XsH(m, t, meson='B+', mode='foresee')['total']
        rf, rp = bf/bo, bpf/bo
        rows.append(dict(mH=m, tanb=t, brF=bf, brO=bo, brP=bpf, rF=rf, rP=rp,
                         brkF=bracket_foresee(t), brkO=bracket_ours(m, t)))
        print(f"{m:7.4f} {t:8.1f} | {bf:12.5e} {bo:12.5e} {bpf:12.5e} | {rf:9.4f} {rp:9.4f} | "
              f"{bracket_foresee(t):10.5e} {bracket_ours(m,t):10.5e}  prod_foresee_mode/ours={bpf_f/bo:.4f}")

def st(xs, key):
    v = sorted(abs(math.log(x[key])) for x in xs if x[key] > 0)
    return f"median={statistics.median(v):.4f} max={v[-1]:.4f} n={len(v)}"
print()
print("ratio statistics: prod/ours :", st(rows, 'rP'))
print("                  FORESEE/ours:", st(rows, 'rF'))
print()
print("bracket ratio (FORESEE/ours) over the grid (should be tanb-only if mH-independent):")
for t in grid_t:
    print(f"  tanb={t:8.1f}: " + " ".join(f"{bracket_foresee(t)/bracket_ours(m,t):.5f}" for m in (0.5, 1.0, 3.0, 4.0)))
print()
print("limits:")
for t in (1e3, 1e5, 1e8):
    print(f"  tanb={t:9.1e}: brk_F*tanb={bracket_foresee(t)*t:.6f}  brk_ours*tanb={bracket_ours(1.0,tanb=t)*t:.6f}"
          f"   (mH=1)     brk_ours*tanb(mH=0.021)={bracket_ours(0.021,tanb=t)*t:.6f}"
          f"   (mH=4.0)={bracket_ours(4.0,tanb=t)*t:.6f}")
print()
print("mH -> 0 and the phase-space convention (mB vs m_b):")
for m in (0.021, 1.0, 2.0, 3.0, 4.0, 4.18, 4.5):
    print(f"  mH={m:5.2f}: (1-m^2/mB^2)^2/(1-m^2/mb^2)^2 = "
          f"{((1-(m/MB_FOR)**2)/(1-(m/4.18)**2))**2:.5f}" if m < 4.18 else
          f"  mH={m:5.2f}: (1-m^2/mB^2)^2 = {(1-(m/MB_FOR)**2)**2:.5f}  [our formula closes at 4.18]")
print()
print("tanb where our B->XsH branch saturates Br=0.1 (validity):",
      {t: f"{AU.br_B_to_XsH(1.0, t)[0]:.3e}" for t in (5,10,20,30,50,100)})
json.dump(rows, open(r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream/work/review_check/foresee_B_compare.json","w"), indent=1)
print("\nwrote foresee_B_compare.json")
