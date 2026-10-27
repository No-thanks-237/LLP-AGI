
import sys, os, math, json, random
import numpy as np
RUN = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3"
sys.path.insert(0, os.path.join(RUN, "work", "analysis"))
sys.path.insert(0, os.path.join(RUN, "work", "widths"))
import bc5_lib as L, scalar_widths as sw
L.load_hadron_lib(os.path.join(RUN, "work", "magnus", "hadron_lib"))
scan = json.load(open(os.path.join(RUN, "outputs", "data", "signal_scan.json")))
masses = scan["masses"]
print("n scan masses", len(masses))
for ch, c in L.HADRON_LIB.chan.items():
    print(ch, "n_nodes", len(c["mass_nodes"]), "range", c["mass_nodes"][0], c["mass_nodes"][-1])
    print("   ", [round(float(x),4) for x in c["mass_nodes"]])
