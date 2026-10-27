# -*- coding: utf-8 -*-
"""Final review, checks 1-2: crossings from the scan grid vs exclusion_boundary.csv,
and the MC-error metric.  (Pass 2, with the s90 = -ln 0.1 comparison done in log space.)
"""
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
print("CHECK 1: crossings of s(tanb) = 2.302585 (ln s vs ln tanb linear interpolation)")
print("="*100)
res = []
for mreq in (0.3, 1.0, 3.5):
    im = idx(mreq, mg); m = mg[im]
    for L in (25.0, 300.0):
        s = [P[(im, it)]['s'][key(L)] for it in range(len(tg))]
        ds = [P[(im, it)]['ds'][key(L)] for it in range(len(tg))]
        cr = []
        for it in range(len(tg)-1):
            l1, l2 = math.log(s[it]), math.log(s[it+1])
            if (l1-L90)*(l2-L90) <= 0.0:
                f = (L90-l1)/(l2-l1)
                tc = math.exp(math.log(tg[it]) + f*(math.log(tg[it+1])-math.log(tg[it])))
                slope = (l2-l1)/(math.log(tg[it+1])-math.log(tg[it]))
                rel = (ds[it]/s[it]) + f*((ds[it+1]/s[it+1])-(ds[it]/s[it]))
                cr.append(dict(t=tc, slope=slope, i=it, f=f, rel=rel,
                               rel_lo=ds[it]/s[it], rel_hi=ds[it+1]/s[it+1]))
        b = bnd[(m, L)]
        print(f"\nmH = {m} GeV ({mreq}), L = {L:.0f} fb^-1 : grid crossings = {len(cr)}, CSV n_crossings = {b['n_crossings']}")
        for c in cr:
            short = c['slope'] > 0
            kind = 'short-lived (s rises with tanb)' if short else 'long-lived (s falls with tanb)'
            ct = float(b['tanb_short']) if short else float(b['tanb_long'])
            ce = float(b['tanb_short_err']) if short else float(b['tanb_long_err'])
            exp_err = c['t']*c['rel']/abs(c['slope'])
            print(f"  tanb_cross = {c['t']:.4f}   dlns/dlntanb = {c['slope']:+.4f}   [{kind}]")
            print(f"    CSV tanb = {ct:.4f} +- {ce:.4f}   -> rel diff = {abs(c['t']/ct-1):.4e}"
                  f"   (bracket tanb = {tg[c['i']]:.2f} .. {tg[c['i']+1]:.2f})")
            if not short:
                print(f"    MC error: ds/s at bracket = {c['rel_lo']:.6f}, {c['rel_hi']:.6f} (interp {c['rel']:.6f});"
                      f"  CSV rel_ds_at_long = {float(b['rel_ds_at_long']):.6f}"
                      f"  -> ratio {float(b['rel_ds_at_long'])/c['rel']:.5f}")
                print(f"              (ds/S90 instead would be {c['rel']*S90/S90:.4f}? no: ds/s90 = {c['rel']/S90:.4f}"
                      f" and s90 = {S90:.6f})")
                print(f"    tanb error: tanb*(ds/s)/|dlns/dlntanb| = {exp_err:.4f} vs CSV {ce:.4f}"
                      f"  ratio {exp_err/ce:.4f}")
                res.append((m, L, c['t'], ct, abs(c['t']/ct-1), exp_err, ce))
            else:
                print(f"    (short-lived edge: MC rel at crossing {c['rel']:.4f}; expected tanb err "
                      f"{exp_err:.4f} vs CSV {ce:.4f})")
print()
print("long-lived crossings: max |rel diff| =", f"{max(r[4] for r in res):.3e}",
      " median =", f"{sorted(r[4] for r in res)[len(res)//2]:.3e}")
print("tanb-error check: max |ratio-1| =", f"{max(abs(r[5]/r[6]-1) for r in res):.3e}")
