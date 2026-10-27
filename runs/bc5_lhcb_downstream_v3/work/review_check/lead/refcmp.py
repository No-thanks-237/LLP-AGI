
import json, os, csv, math
import numpy as np
RUN = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3"
ref = {}
with open(os.path.join(RUN, "work", "widths", "reference_figures", "bc5_reference_curves.csv")) as f:
    for row in csv.DictReader(f):
        ref.setdefault(row.get("lumi_fb") or row.get("luminosity") or "?", {})
print("ref columns:", list(row.keys()) if 'row' in dir() else "?")
# inspect file head
print(open(os.path.join(RUN, "work", "widths", "reference_figures", "bc5_reference_curves.csv")).read()[:600])
