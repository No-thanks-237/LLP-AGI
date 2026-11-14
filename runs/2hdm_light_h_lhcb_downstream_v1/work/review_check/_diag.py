import importlib.util, math
import audit_2hdm as AU
spec = importlib.util.spec_from_file_location("prod", r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream/work/production/twohdm_production.py")
prod = importlib.util.module_from_spec(spec); spec.loader.exec_module(prod)
mH, tb = 1.0, 100.0
print("prod constants: v=", prod.V_VEV, " mW=", getattr(prod,'M_W',None), " mb=", prod.M_B_QUARK, " mc=", prod.M_C_QUARK,
      " GF=", getattr(prod,'G_F',None), " BR_BXc=", prod.BR_B_XC_ENU, " F_PHASE_B=", prod.F_PHASE_B)
print("prod S():", prod._S())
d = prod.xi_H_details(mH, tb, 's', 'b')
import json
print("prod xi_H_details:")
for k,v in d.items():
    print("   ", k, "=", v)
print()
xi, terms = AU.xi_bs(mH, tb, return_terms=True)
print("mine: xi_bs =", xi, " terms:", {k: str(v) for k,v in terms.items()})
print("prod xi_H =", prod.xi_H(mH, tb, 's','b'), " |prod|=", abs(prod.xi_H(mH,tb,'s','b')), " |mine|=", abs(xi))
print()
c,s = AU.mixing(tb); xH = AU.MHC**2/AU.MW**2
print("my g012 for u,c,t:", [AU.g012(AU.MQ[n]**2/AU.MW**2, xH, 1/tb) for n in ('u','c','t')])
try:
    print("prod g_aux_functions(tanb):", prod.g_aux_functions(tb))
except Exception as e:
    print("g_aux err", e)
print()
print("prod br dict keys:", list(prod.br_B_to_XsH(mH,tb).keys()))
print("prod ratio_to_BXcenu:", prod.br_B_to_XsH(mH,tb)['ratio_to_BXcenu'], " total:", prod.br_B_to_XsH(mH,tb)['total'])
print("mine ratio:", AU.br_B_to_XsH(mH,tb)[0]/0.104, " total:", AU.br_B_to_XsH(mH,tb)[0])
print("F_PHASE (prod) vs mine:", prod.F_PHASE_B, AU.f_phase((1.28/4.18)**2))
print("V: prod V_CKM=", prod.V_CKM)
print("my VCKM=", AU.VCKM)
