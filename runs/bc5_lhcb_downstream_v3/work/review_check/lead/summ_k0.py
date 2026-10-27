
import json
d = json.load(open('impact_k0half.json'))
b = {round(r['m'],6): r for r in json.load(open('impact_base.json'))['rows']}
print("k0half (neutral kaon mixing weight x 0.5):")
for r in d['rows']:
    m = round(r['m'],6)
    lo, lo0 = r['lo'], b.get(m,{}).get('lo')
    print("  m=%8.3f  lo=%.4e  base=%.4e  ratio=%.3f" % (m, lo, lo0 if lo0 else float('nan'),
          (lo/lo0) if (lo and lo0) else float('nan')))
