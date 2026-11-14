import json
R = r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream"
d = json.load(open(R + "/work/analysis/scan_2hdm/m_1.json"))
print("table_meta:", json.dumps(d['table_meta'], indent=1)[:1800])
print("table_counts:", d['table_counts'])
c = d['curves']['25.0']
print("curves['25.0'] keys:", list(c.keys()) if isinstance(c,dict) else type(c))
if isinstance(c, dict):
    for k,v in c.items():
        if isinstance(v, list):
            print("  ", k, "list len", len(v), "sample:", v[:3])
        else:
            print("  ", k, "=", str(v)[:120])
