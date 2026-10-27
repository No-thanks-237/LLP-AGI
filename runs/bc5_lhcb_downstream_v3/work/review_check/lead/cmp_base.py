
import json, os, csv
RUN = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3"
imp = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'impact_base.json')))
prod = {}
with open(os.path.join(RUN, 'outputs', 'data', 'exclusion_boundary.csv')) as f:
    for row in csv.DictReader(f):
        if float(row['lumi_fb']) != 25.0: continue
        prod[round(float(row['mass_GeV']), 6)] = (float(row['theta2_long_lived_edge']),
                                                 float(row['theta2_short_lived_edge']))
print("%10s %12s %12s %8s" % ("m", "impact lo", "prod lo", "ratio"))
for r in imp['rows']:
    m = round(r['m'], 6)
    p = prod.get(m)
    print("%10.5f %12.3e %12s %8s" % (m, r['lo'], ('%.3e' % p[0]) if p else '-',
          ('%.3f' % (r['lo']/p[0])) if p and p[0] else '-'))
