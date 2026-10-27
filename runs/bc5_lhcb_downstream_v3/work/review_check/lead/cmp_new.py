
import json, os
old = {round(r['m'],6): r for r in json.load(open('impact_base.json'))['rows']}
try:
    new = {round(r['m'],6): r for r in json.load(open('impact_newbase.json'))['rows']}
except FileNotFoundError:
    raise SystemExit('newbase not ready')
print("v3.3 -> v3.4 (all review fixes except the library regeneration and the Higgs sample)")
print("%10s %12s %12s %8s" % ("m", "v3.3 lo", "v3.4 lo", "ratio"))
for m in sorted(old):
    if m not in new: continue
    a, b = old[m]['lo'], new[m]['lo']
    print("%10.5f %12.3e %12.3e %8.3f" % (m, a, b, b/a if a else float('nan')))
