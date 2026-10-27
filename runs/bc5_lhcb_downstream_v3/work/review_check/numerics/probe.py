
import json, os, sys
RUN = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
p = os.path.join(RUN,"outputs/data/signal_scan.json")
d = json.load(open(p))
print("keys:", list(d.keys()))
print("n masses:", len(d["masses"]), "first/last:", d["masses"][:4], d["masses"][-3:])
print("theta2:", d["theta2"])
k = str(d["masses"][0])
print("row0 type:", type(d["signal_per_fb_per_channel"][k]), "len", len(d["signal_per_fb_per_channel"][k]))
print("row0[0]:", d["signal_per_fb_per_channel"][k][0])
print("row0[-1]:", d["signal_per_fb_per_channel"][k][-1])
for kk in d.keys():
    if kk not in ("masses","theta2","signal_per_fb_per_channel"):
        print(kk, "->", json.dumps(d[kk])[:400])
print("norm exists:", os.path.exists(os.path.join(RUN,"work/magnus/samples/normalisation.json")))
