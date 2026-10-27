
import sys, os
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\widths")
import scalar_widths as sw
sw.load_hadronic_table(os.path.join(os.path.dirname(sw.__file__), 'hadronic_table.csv'))
for m in (2.0, 3.0, 4.984, 10.0, 30.0, 62.49):
    w = sw.widths(m)
    tot = sum(w.values())
    print("m=%7.3f  tot=%.4e  sum(BR)=%.4f" % (m, tot, sum(v/tot for v in w.values())))
    print("    ", {k: round(v/tot, 4) for k, v in w.items() if v/tot > 1e-3})
