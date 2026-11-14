import json, math
R = r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream"
d = json.load(open(R + "/outputs/data/signal_scan.json"))
mg, tg = d['mass_grid'], d['tanb_grid']
def idx(v, grid):
    return min(range(len(grid)), key=lambda i: abs(grid[i]-v))
pts = {}
bad = 0
for p in d['points']:
    im, it = idx(p['mH'], mg), idx(p['tanb'], tg)
    if abs(mg[im]-p['mH']) > 1e-6*max(1,mg[im]) or abs(tg[it]-p['tanb']) > 1e-6*tg[it]:
        bad += 1
    pts[(im, it)] = p
print("points:", len(pts), " unmatched:", bad, " expected 66*40 =", 66*40)
# check the tanb-key rounding issue
print("example: tg[16] =", repr(tg[16]), " point tanb =", repr([p['tanb'] for p in d['points'] if abs(p['tanb']-tg[16])<1e-9][0]))
miss = [(im,it) for im in range(len(mg)) for it in range(len(tg)) if (im,it) not in pts]
print("missing grid cells:", len(miss), miss[:10])
