import json
R = r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream"
d = json.load(open(R + "/work/analysis/scan_2hdm/m_1.json"))
print("top keys:", list(d.keys()))
for k,v in d.items():
    if isinstance(v, dict):
        print(k, "->", list(v.keys())[:14])
    elif isinstance(v, list):
        print(k, "-> list len", len(v), "first:", (v[0] if not isinstance(v[0],(dict,list)) else list(v[0].keys()) if isinstance(v[0],dict) else len(v[0])))
    else:
        print(k, "=", v)
