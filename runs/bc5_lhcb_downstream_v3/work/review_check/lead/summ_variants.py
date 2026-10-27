
import json, os, glob
base = json.load(open('impact_base.json'))
b = {round(r['m'],6): r for r in base['rows']}
print("%9s %11s | %s" % ("m", "base lo", "  ".join("%9s" % v for v in
      ['rot','libmass','kLdead','multi4pi','gamma','bxi','bincl','bsext'])))
data = {}
for v in ['rot','libmass','kLdead','multi4pi','gamma','bxi','bincl','bsext']:
    p = 'impact_%s.json' % v
    if not os.path.exists(p):
        print('missing', v); continue
    d = json.load(open(p))
    data[v] = {round(r['m'],6): r for r in d['rows']}
for m in sorted(b):
    row = []
    for v in ['rot','libmass','kLdead','multi4pi','gamma','bxi','bincl','bsext']:
        if v not in data or m not in data[v]:
            row.append('%9s' % '-'); continue
        lo = data[v][m]['lo']; lo0 = b[m]['lo']
        row.append('%9s' % (('%.3f' % (lo/lo0)) if (lo and lo0) else '-'))
    print("%9.5f %11.3e | %s" % (m, b[m]['lo'], "  ".join(row)))
