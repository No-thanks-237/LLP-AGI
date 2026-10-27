import sys, os, csv, math, json
import numpy as np
RUN=r'E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3'
sys.path.insert(0, os.path.join(RUN,'work','widths'))
import scalar_widths as sw, reference_br_model as RB
sw.load_hadronic_table(os.path.join(RUN,'work','widths','hadronic_table.csv'))
p=os.path.join(RUN,'outputs','data','bc5_branching_ratios_senscalc.csv')
with open(p) as fh:
    rows=[r for r in csv.DictReader(l for l in fh if not l.startswith('#'))]
print('columns',list(rows[0].keys()))
for m in (0.3,0.5,1.0,1.5,1.99,2.5):
    i=min(range(len(rows)),key=lambda k:abs(float(rows[k]['mass_GeV'])-m))
    r=rows[i]; v={k:float(r[k]) for k in r if k!='mass_GeV'}
    s=sum(v.values())
    brmu=v['mu_plus_mu_minus']; bre=v['e_plus_e_minus']
    gtot=sw.gamma_ll(float(r['mass_GeV']),sw.M_MU)/brmu if brmu>1e-12 else sw.gamma_ll(float(r['mass_GeV']),sw.M_E)/bre
    gw=sw.widths(float(r['mass_GeV'])); totw=sum(gw.values())
    print('m=%.3f colsum=%.5f BRmu=%.4f  Gref_tot=%.4e Gours_tot=%.4e ratio=%.3f'%(
        float(r['mass_GeV']),s,brmu,gtot,totw,totw/gtot))
    print('   ref pipi=%.4f KK=%.4f 2pi2pi=%.4f ss=%.3g gg=%.3g cc=%.3g'%(
        v['pi_plus_pi_minus'],v['K_plus_K_minus'],v['2pi_plus_2pi_minus'],v['Jets_ss'],v['Jets_GG'],v['Jets_cc']))
    print('   ours pipi=%.4e KK=%.4e multi=%.4e gamma=%.4e  BRours: pipi=%.4f KK=%.4f'%(
        gw.get('pipi',0),gw.get('KK',0),gw.get('multi',0),gw.get('gamma',0),
        gw.get('pipi',0)/totw,gw.get('KK',0)/totw))
    print('   partial-width ratios: pipi ours/ref=%.3f KK ours/ref=%.3f'%(
        (gw.get('pipi',0)/(v['pi_plus_pi_minus']*gtot)) if v['pi_plus_pi_minus']>0 else float('nan'),
        (gw.get('KK',0)/(v['K_plus_K_minus']*gtot)) if v['K_plus_K_minus']>0 else float('nan')))
