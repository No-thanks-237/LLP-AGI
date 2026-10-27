
import csv, os
CSV = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\outputs\data\bc5_branching_ratios_senscalc.csv"
rows = [r for r in csv.DictReader(l for l in open(CSV) if not l.startswith('#'))]
cols = [c for c in rows[0].keys() if c != 'mass_GeV']
print("columns:", cols)
for target in (0.3, 0.5, 0.8, 1.0, 1.5, 1.99, 2.0, 3.0, 10.0, 40.0):
    best = min(rows, key=lambda r: abs(float(r['mass_GeV'])-target))
    m = float(best['mass_GeV'])
    vals = {c: float(best[c]) for c in cols}
    s = sum(vals.values())
    pipi = vals.get('pi_plus_pi_minus', 0.0); kk = vals.get('K_plus_K_minus', 0.0)
    print("m=%8.4f sum=%.5f  1-sum=%.5f  pipi=%.5f kk=%.5f  1.5*pipi=%.5f 2*kk=%.5f" % (
        m, s, 1-s, pipi, kk, 1.5*pipi, 2*kk))
