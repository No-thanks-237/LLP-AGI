
import sys, math
sys.path.insert(0, r'E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3/work/analysis')
import bc5_production as bp
print('our br_K_to_S(0) K+   = %.6e'%bp.br_K_to_S(0.0,'K+'))
print('our br_K_to_S(0) K0L  = %.6e'%bp.br_K_to_S(0.0,'K0L'))
print('SensCalc tables: K+ = 1.4261876687078384e-03, K0L = 6.198330036596367e-03')
print('ratio our K0L/K+ = %.4f ; SensCalc ratio = %.4f'%(bp.br_K_to_S(0.0,'K0L')/bp.br_K_to_S(0.0,'K+'), 6.198330036596367e-03/1.4261876687078384e-03))
print()
print('our br_B_to_S(0) =', bp.br_B_to_S(0.0,'B+'), bp.br_B_to_S(0.0,'B0'))
print('SensCalc BtoXsStotal(0)=3.3095920009617386 (default/exclusive), BtoXsSinclusive(0)=6.629098459885584')
print('our br_B_to_S(2.0) =', bp.br_B_to_S(2.0,'B+'))
print('our br_B_to_S(4.0) =', bp.br_B_to_S(4.0,'B+'))
print('our br_B_to_S(4.70) =', bp.br_B_to_S(4.70,'B+'))
print('our br_B_to_S(4.75) =', bp.br_B_to_S(4.75,'B+'))
