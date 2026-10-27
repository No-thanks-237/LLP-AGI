
import sys
sys.path.insert(0, r'E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3/work/widths')
import scalar_widths as sw
print('m     Gamma_tot[GeV]   ctau[m] @theta2=1    BR(mumu)')
for m in [0.603,3.822,5.0,13.0,33.0,53.0,62.0]:
    ws=sw.widths(m)
    G=sw.gamma_total(m)
    ct=sw.ctau_m(m,1.0)
    brmumu=ws.get('mumu', ws.get('mu+mu-', None))
    print('%6.3f  %.6e  %.6e   %s'%(m,G,ct,brmumu))
print()
print('channels:', sorted(sw.widths(13.0).keys()))
