import importlib.util, math, statistics
import audit_2hdm as AU
spec = importlib.util.spec_from_file_location("prod", r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream/work/production/twohdm_production.py")
prod = importlib.util.module_from_spec(spec); spec.loader.exec_module(prod)
print("DEFAULT_B_MODE =", prod.DEFAULT_B_MODE, " B_MODES =", prod.B_MODES)
print()
print("B -> Xs H: mode comparison, grid 8 mH x 6 tanb")
ms = [0.021,0.27,0.5,1.0,2.0,3.0,4.0,4.18]
ts = [5.0,10.0,30.0,100.0,1000.0,100000.0]
lp = []; lf = []; ld = []
print(f"{'mH':>7} {'tanb':>8} {'ours(paper,mb)':>15} {'prod paper':>12} {'prod foresee':>12} {'ratio_paper':>11} {'ratio_foresee':>13}")
for m in ms:
    for t in ts:
        o = AU.br_B_to_XsH(m,t)[0]
        p = prod.br_B_to_XsH(m,t,mode='paper_mb')['total']
        f = prod.br_B_to_XsH(m,t,mode='foresee')['total']
        d = prod.br_B_to_XsH(m,t)['total']
        if o>0 and p>0: lp.append(math.log(p/o))
        if o>0 and f>0: lf.append(math.log(f/o))
        if o>0 and d>0: ld.append(math.log(d/o))
        rp_ = p/o if o>0 else float('nan'); rf_ = f/o if o>0 else float('nan')
        print(f"{m:7.4f} {t:8.1f} {o:15.5e} {p:12.4e} {f:12.4e} {rp_:11.6f} {rf_:13.4f}")
def st(v):
    a=sorted(abs(x) for x in v); return f"n={len(a)} median={statistics.median(a):.3e} max={a[-1]:.3e}"
print()
print("log-ratio |prod_paper/ours|   :", st(lp))
print("log-ratio |prod_foresee/ours| :", st(lf))
print("log-ratio |prod_default/ours| :", st(ld), " (default = %s)" % prod.DEFAULT_B_MODE)
print()
print("="*95)
print("K -> pi H: production vs independent implementation")
print(f"{'mH':>7} {'tanb':>8} {'K+ ours':>13} {'K+ prod':>13} {'ratio':>9} | {'K0L ours':>12} {'K0L prod':>12} {'ratio':>9}")
for m in [0.021,0.05,0.1396,0.2,0.3,0.354]:
    for t in [5.0,100.0,1000.0]:
        ko = AU.br_K_to_piH(m,t,'K+'); kp = prod.br_K_to_piH(m,t,'K+')
        kp = kp['total'] if isinstance(kp,dict) else kp
        lo = AU.br_K_to_piH(m,t,'K0L'); lpv = prod.br_K_to_piH(m,t,'K0L')
        lpv = lpv['total'] if isinstance(lpv,dict) else lpv
        print(f"{m:7.4f} {t:8.1f} {ko:13.5e} {kp:13.5e} {kp/ko if ko>0 else float('nan'):9.6f} | "
              f"{lo:12.4e} {lpv:12.4e} {lpv/lo if lo>0 else float('nan'):9.6f}")
print()
for m in [0.0,0.05,0.1396]:
    for t in [5.0,100.0]:
        po = AU.br_semileptonic(m,t,'pi'); pp = prod.br_pi_to_Henu(m,t)
        pp = pp['total'] if isinstance(pp,dict) else pp
        ko = AU.br_semileptonic(m,t,'K'); kp = prod.br_K_to_Henu(m,t)
        kp = kp['total'] if isinstance(kp,dict) else kp
        print(f"semileptonic mH={m:6.4f} tanb={t:6.1f}: pi: ours={po:.5e} prod={pp:.5e} r={pp/po if po else float('nan'):.6f} ;"
              f" K: ours={ko:.5e} prod={kp:.5e} r={kp/ko if ko else float('nan'):.6f}")
print()
print("K->pi H four-quark/FCNC ratio (|A_4q/A_FCNC|):")
for t in [5.0,100.0,1000.0]:
    M,f4,fc = AU.M_KpiH(0.0,t)
    print(f"  tanb={t:8.1f}: |4q|={abs(f4):.4e} |FCNC|={abs(fc):.4e} ratio={abs(f4/fc):.5f}"
          f"   Br(K+->piH)={AU.br_K_to_piH(0.0,t,'K+'):.5e}")
