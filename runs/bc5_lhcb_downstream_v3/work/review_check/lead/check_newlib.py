
import sys, os, math, json
import numpy as np
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\analysis")
import bc5_lib as L
L.load_hadron_lib(r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\magnus\hadron_lib")
masses = json.load(open(r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\outputs\data\signal_scan.json"))["masses"]
for ch, c in L.HADRON_LIB.chan.items():
    nodes = c["mass_nodes"]
    print("%-7s nodes=%3d  %.6f .. %.6f" % (ch, len(nodes), nodes[0], nodes[-1]))
bad = []
for m in masses:
    br,_ = L.channel_branching(m)
    for ch in ('gg','ss','cc','bb','tau'):
        if br.get(ch,0.0) < 1e-4: continue
        c = L.HADRON_LIB.chan[L._LIB_MAP[ch]]
        i = int(np.argmin(np.abs(c['mass_nodes']-m)))
        rel = abs(c['mass_nodes'][i]-m)/max(m, c['mass_nodes'][i])
        if rel > L.NODE_MTOL:
            bad.append((m, ch, float(c['mass_nodes'][i]), rel))
print("\npairs needing dilation: %d" % len(bad))
for b in bad[:12]: print("  m=%.6f %s node=%.6f rel=%.2f%%" % (b[0], b[1], b[2], 100*b[3]))
