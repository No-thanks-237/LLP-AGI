import json, csv
R = r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream"
d = json.load(open(R + "/outputs/data/signal_scan.json"))
mg = d['mass_grid']
print("mass grid (66):")
for i in range(0, len(mg), 6):
    print("   ", " ".join(f"{m:9.4f}" for m in mg[i:i+6]))
print()
rows = list(csv.DictReader(open(R + "/outputs/data/exclusion_boundary.csv")))
print("boundary columns:", rows[0].keys())
print("target masses: 0.3 ->", [m for m in mg if abs(m-0.3)<0.06], " 1.0 ->", [m for m in mg if abs(m-1.0)<0.06],
      " 3.5 ->", [m for m in mg if abs(m-3.5)<0.15])
seen = set()
for r in rows:
    m = float(r['mass_GeV'])
    if m in (min(mg, key=lambda x: abs(x-0.3)), min(mg, key=lambda x: abs(x-1.0)), min(mg, key=lambda x: abs(x-3.5)),
             0.021, 10.0):
        print(r)
