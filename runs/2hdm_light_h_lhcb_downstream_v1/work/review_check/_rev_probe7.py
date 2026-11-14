import json, math
R = r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream"
d = json.load(open(R + "/outputs/data/signal_scan.json"))
mg, tg = d['mass_grid'], d['tanb_grid']
idx = lambda v, g: min(range(len(g)), key=lambda i: abs(g[i]-v))
P = {}
for p in d['points']:
    P[(idx(p['mH'], mg), idx(p['tanb'], tg))] = p
S90 = -math.log(0.1)
print("s(tanb) profile in the VALID region (Br(B->XsH)<0.1) for several masses, 25 fb^-1:")
print(f"{'mH':>9} {'tanb_valid':>10} {'s_max(valid)':>14} {'s_max/s90':>10} {'tanb at s_max':>13}  s at 5,10,30,100,300,1000,1e4,1e5")
for m in (0.021, 0.289535, 1.0, 3.50112, 4.1800):
    im = idx(m, mg)
    # validity boundary from the scan's own validity flag
    valid = [P[(im,it)] for it in range(len(tg)) if P[(im,it)].get('validity',{}).get('tanb_BR_gt_0.1',1e18) < tg[it]]
    smax = max((p['s']['25'], p['tanb']) for p in [P[(im,it)] for it in range(len(tg))])
    sv = []
    for t in (5,10,30,100,300,1000,10000,100000):
        it = idx(t, tg); sv.append(f"{P[(im,it)]['s']['25']:.2g}")
    print(f"{m:9.5f} {valid[0]['tanb'] if valid else float('nan'):10.2f} {smax[0]:14.4g} {smax[0]/S90:10.3g} {smax[1]:13.2f}  {', '.join(sv)}")
