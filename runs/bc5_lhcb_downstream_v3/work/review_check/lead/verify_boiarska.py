
import sys, os, math
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\analysis")
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\widths")
import bc5_production as P
# Boiarska tab:BR at m_S = 0 (B+-), their column 2
tab = {'K1_1270': 0.91, 'K0*700': 0.76, 'K*892': 0.47, 'K': 0.43, 'K2*1430': 0.30,
       'K*1410': 0.21, 'K*1680': 0.13, 'K0*1430': 0.081, 'K1_1400': 0.016, 'pi': 0.013}
print("B+- channels at m_S -> 0: code vs Boiarska tab:BR (col 2)")
br = P.br_B_to_S(1e-6, 'B+')
tot = 0.0
for k, v in sorted(br.items()):
    ref = tab.get(k)
    print("  %-9s code = %8.4f   table = %7.3f  ratio = %.3f" % (k, v, ref if ref else float('nan'), v/ref if ref else float('nan')))
    tot += v
print("  sum = %.4f" % tot)
print("\nB0 / B+ ratio (table caption: 0.93):")
b0 = P.br_B_to_S(1e-6, 'B0')
for k in br:
    if br[k] > 1e-6:
        print("  %-9s %.4f" % (k, b0[k]/br[k]))
print("\nB_s mixing channels:", P.br_B_to_S(1e-6, 'Bs'))
print("\nK channels at m_S -> 0 (table: K+->Spi+ 1.7e-3, K0L->Spi0 7e-3):")
print("  K+   = %.4e" % P.br_K_to_S(1e-6, 'K+'))
print("  K0L  = %.4e" % P.br_K_to_S(1e-6, 'K0L'))
print("\nxi_sb used: %.4e (Boiarska Table 2: 3.6e-4)" % P.XI_SB)
print("quartic channels at 1 GeV:", [(q[0], q[1], round(q[3], 12)) for q in P.quartic_channels(1.0)][:4])