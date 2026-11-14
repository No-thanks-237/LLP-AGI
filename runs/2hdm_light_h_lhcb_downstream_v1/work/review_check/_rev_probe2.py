import json, csv, math
R = r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream"
d = json.load(open(R + "/outputs/data/signal_scan.json"))
key = lambda L: str(int(L))
pts = {}
for p in d['points']:
    pts[(round(p['mH'],6), round(p['tanb'],8))] = p
tg = d['tanb_grid']
m = 1.0
print("signal_scan.json grid for mH=1.0, L=25 fb^-1 (tanb >= 100):")
print(f"{'tanb':>10} {'s':>14} {'ds':>12} {'ds/s':>9} {'ctau_m':>12} {'eps_B':>12}")
for t in tg:
    p = pts.get((round(m,6), round(t,8)))
    if p is None: continue
    if t >= 100:
        print(f"{t:10.2f} {p['s'][key(25)]:14.6g} {p['ds'][key(25)]:12.4g} {p['ds'][key(25)]/p['s'][key(25)]:9.5f} {p['ctau_m']:12.4g} {p['eps']['B']:12.5g}")
print()
rows = list(csv.DictReader(open(R + "/outputs/data/exclusion_s_curves.csv")))
print("exclusion_s_curves.csv: n rows =", len(rows), "columns:", list(rows[0].keys()))
import collections
c = collections.Counter(float(r['mH_GeV']) for r in rows)
print("masses:", len(c), " rows per mass:", sorted(set(c.values())))
print("tanb values for mH=1.0:", sorted(float(r['tanb']) for r in rows if float(r['mH_GeV'])==1.0))
