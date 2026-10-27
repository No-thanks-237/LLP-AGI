
import sys, os, math
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\widths")
import scalar_widths as sw
sw.load_hadronic_table(os.path.join(os.path.dirname(sw.__file__), 'hadronic_table.csv'))
print("SensCalc 2025 table (audit Q4) vs our strict-Winkler widths")
print("%8s %14s %14s %8s %10s" % ("m_S","Gamma_senscalc","Gamma_ours","ratio","ctau_ours(m)"))
for m, g in ((5.082255, 9.09e-6), (13.023095, 1.425e-4), (33.371213, 8.82e-4),
             (53.419632, 1.390e-3), (62.49, 1.606e-3)):
    ours = sum(sw.widths(m).values())
    print("%8.3f %14.4e %14.4e %8.3f %10.3g" % (m, g, ours, ours/g, sw.HBAR_C/ours))
