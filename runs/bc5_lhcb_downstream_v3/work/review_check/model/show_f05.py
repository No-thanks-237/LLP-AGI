import json,os
HERE=os.path.dirname(os.path.abspath(__file__))
d=json.load(open(os.path.join(HERE,'f05_results.json')))
print('%-16s %7s %8s %8s %8s %8s %8s %8s %8s'%('key','BRmulti','Abase','A_cur','A_4ch','dRR_4ch','A_rho','dRR_rho','dRR_eta'))
for k,v in d.items():
    print('%-16s %7.3f %8.3f %8.3f %8.3f %+8.1f %8.3f %+8.1f %+8.1f'%(k,v['BR_multi'],v['A_base'],v['A_multi_current'],
      v['all 4 charged pions']['A_multi'],100*v['all 4 charged pions']['delta_R_over_R'],
      v['rho rho -> 4 charged pions']['A_multi'],100*v['rho rho -> 4 charged pions']['delta_R_over_R'],
      100*v['current 3-way + eta decays']['delta_R_over_R']))
a=json.load(open(os.path.join(HERE,'f05_addendum.json'))) if os.path.exists(os.path.join(HERE,'f05_addendum.json')) else {}
for k,v in a.items():
    print('add %-12s iso4pi A=%.3f dRR=%+.1f%%  pi2pi0 dRR=%+.1f%%'%(k,v['isospin 4pi']['A'],100*v['isospin 4pi']['delta_R_over_R'],100*v['pi pi 2pi0']['delta_R_over_R']))
