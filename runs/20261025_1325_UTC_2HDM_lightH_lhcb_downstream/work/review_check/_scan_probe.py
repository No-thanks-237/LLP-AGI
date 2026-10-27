import json
p = r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream/outputs/data/signal_scan.json"
d = json.load(open(p))
print("model/version/b_mode:", d['model'], d['version'], d['b_mode'], "| s90 =", d['s_threshold'])
print("mass_grid[0:3]", d['mass_grid'][:3], "...", d['mass_grid'][-3:])
print("tanb_grid:", [round(t,4) for t in d['tanb_grid'][:6]], "...", [round(t,4) for t in d['tanb_grid'][-4:]])
print("lumis:", d['lumis'])
print("n points:", len(d['points']))
print("keys of point 0:", list(d['points'][0].keys()))
print(json.dumps(d['points'][0], indent=1)[:2500])
