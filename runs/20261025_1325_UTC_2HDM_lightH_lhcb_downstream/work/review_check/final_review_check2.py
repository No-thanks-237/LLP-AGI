# -*- coding: utf-8 -*-
"""Final review, check 2 (exact reproduction of the CSV error metric) + validity flags."""
import json, csv, math
R = r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream"
S90 = -math.log(0.1); L90 = math.log(S90)
d = json.load(open(R + "/outputs/data/signal_scan.json"))
mg, tg = d['mass_grid'], d['tanb_grid']
idx = lambda v, g: min(range(len(g)), key=lambda i: abs(g[i]-v))
P = {}
for p in d['points']:
    P[(idx(p['mH'], mg), idx(p['tanb'], tg))] = p
bnd = {}
for r in csv.DictReader(open(R + "/outputs/data/exclusion_boundary.csv")):
    bnd[(float(r['mass_GeV']), float(r['lumi_fb']))] = r
key = lambda L: str(int(L))
print("="*100)
print("CHECK 2: exact re-derivation of the CSV metric (make_curves.boundary_table convention)")
print("   crossing bracket index i = lower grid point; ds/s at i; slope = centred (i-1,i+1)")
print("="*100)
worst_e = 0.0; worst_t = 0.0; rows_out = []
for mreq in (0.3, 1.0, 3.5):
    im = idx(mreq, mg); m = mg[im]
    for L in (25.0, 300.0):
        sv = [P[(im, it)]['s'][key(L)] for it in range(len(tg))]
        dv = [P[(im, it)]['ds'][key(L)] for it in range(len(tg))]
        cr = []
        for it in range(len(tg)-1):
            l1, l2 = math.log(sv[it]), math.log(sv[it+1])
            if (l1-L90)*(l2-L90) <= 0:
                f = (L90-l1)/(l2-l1)
                tc = math.exp(math.log(tg[it]) + f*(math.log(tg[it+1])-math.log(tg[it])))
                i0, i1 = max(it-1,0), min(it+1, len(tg)-1)
                slope = (math.log(sv[i1])-math.log(sv[i0]))/(math.log(tg[i1])-math.log(tg[i0]))
                cr.append((tc, it, slope, dv[it]/sv[it]))
        b = bnd[(m, L)]
        for (tc, i, slope, re) in cr:
            short = slope > 0
            ct = float(b['tanb_short']) if short else float(b['tanb_long'])
            ce = float(b['tanb_short_err']) if short else float(b['tanb_long_err'])
            et = tc*re/abs(slope)
            rows_out.append((m, L, 'short' if short else 'long', tc, ct, abs(tc/ct-1), re, et, ce))
            if not short:
                worst_e = max(worst_e, abs(et/ce-1)); worst_t = max(worst_t, abs(tc/ct-1))
                print(f"  mH={m:8.5f} L={L:3.0f} LONG : tanb_grid={tc:.6f} CSV={ct:.6f} rel={abs(tc/ct-1):.2e}"
                      f" | rel_ds(ds/s at i)={re:.6f} (ds/S90 would be {re/S90:.6f}, x{S90:.3f} smaller)"
                      f" | err_formula={et:.5f} CSV={ce:.5f} rel={abs(et/ce-1):.2e}")
print()
print(f"CHECK 2 verdict: max |rel diff| of tanb_cross = {worst_t:.2e};"
      f" max |rel diff| of the tanb error = {worst_e:.2e}  (make_curves convention)")
print()
print("="*100)
print("VALIDITY of the SHORT-LIVED edges: Br(B->XsH) at the crossing vs the model-validity bound 0.1")
print("="*100)
sys.path = [r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream/work/review_check"] + sys.path
import audit_2hdm as AU
for (m, L, tag, tc, ct, rel, re, et, ce) in rows_out:
    if tag == 'short':
        br = AU.br_B_to_XsH(m, tc)[0]
        print(f"  mH={m:8.5f} L={L:3.0f} short-lived crossing tanb={tc:8.3f}: Br(B->XsH) = {br:.4g}"
              f"  -> {'INVALID (>0.1)' if br > 0.1 else 'ok'} ; MC rel err ds/s = {re:.3f}")
print()
print("  (DESIGN 3.1b: the Type-I B rate is valid only where Br(B->XsH) <~ 0.1;")
print("   above that the inclusive formula saturates and the point is outside the model's validity.)")
