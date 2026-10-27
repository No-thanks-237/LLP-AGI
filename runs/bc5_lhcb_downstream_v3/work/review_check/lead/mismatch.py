
import sys, os, math, json
import numpy as np
RUN = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3"
sys.path.insert(0, os.path.join(RUN, "work", "analysis"))
import bc5_lib as L
L.load_hadron_lib(os.path.join(RUN, "work", "magnus", "hadron_lib"))
scan = json.load(open(os.path.join(RUN, "outputs", "data", "signal_scan.json")))
masses = scan["masses"]
LIBCH = ['gg','ss','cc','bb','tau']   # 'tau' maps to tautau
print("%9s | %-28s | %-6s | %s" % ("m", "channel", "BR", "node / mismatch"))
worst=[]
for m in masses:
    br,_ = L.channel_branching(m)
    for c in LIBCH:
        b = br.get(c,0.0)
        if b < 1e-4: continue
        ch = L._LIB_MAP[c]
        cdata = L.HADRON_LIB.chan[ch]
        nodes = cdata['mass_nodes']
        i = int(np.argmin(np.abs(nodes-m)))
        node = float(nodes[i])
        rel = abs(node-m)/max(m,node)
        ok = rel <= 0.20
        flag = "" if rel < 1e-3 else ("  <-- MISMATCH %.1f%% (%s)" % (100*rel, "used" if ok else "rejected->analytic"))
        worst.append((b*rel if ok else b*0.5, m, c, b, node, rel, ok))
        if rel >= 1e-3:
            print("%9.6f | %-6s | %.4f | %.6f %s" % (m, c, b, node, flag))
print()
tot = sum(w for w,*_ in worst)
print("sum over (mass,channel) of BR*reldev = %.4f  (unit: BR-weighted fractional mass error)" % tot)
