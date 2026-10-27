import json
R = r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream"
d = json.load(open(R + "/outputs/data/signal_scan.json"))
mg, tg = d['mass_grid'], d['tanb_grid']
idx = lambda v, g: min(range(len(g)), key=lambda i: abs(g[i]-v))
P = {}
for p in d['points']:
    P[(idx(p['mH'], mg), idx(p['tanb'], tg))] = p
im = idx(1.0, mg); iml = idx(0.021, mg)
print("mH =", mg[im], " and ", mg[iml])
print(f"{'tanb':>9} {'ctau_m':>11} {'s_25':>11} {'B':>11} {'inel':>9} {'h':>11} {'eps_B':>11} {'eps_h':>11} {'BR_hHH':>10}")
for it in (0,1,2,3,5,8,11,14,17,20):
    p = P[(im, it)]
    print(f"{p['tanb']:9.3f} {p['ctau_m']:11.4g} {p['s']['25']:11.5g} {p['s_comp']['25']['B']:11.4g}"
          f" {p['s_comp']['25']['inel']:9.3g} {p['s_comp']['25']['h']:11.4g} {p['eps']['B']:11.4g}"
          f" {p['eps']['h']:11.4g} {p['BR_hHH']:10.3g}")
print()
print("mH = 0.021:")
for it in (0,1,2,3,4,8):
    p = P[(iml, it)]
    print(f"  tanb={p['tanb']:9.3f} ctau={p['ctau_m']:11.4g} s25={p['s']['25']:11.4g} B={p['s_comp']['25']['B']:11.4g}"
          f" inel={p['s_comp']['25']['inel']:9.3g} h={p['s_comp']['25']['h']:9.3g} eps_B={p['eps']['B']:11.4g}")
