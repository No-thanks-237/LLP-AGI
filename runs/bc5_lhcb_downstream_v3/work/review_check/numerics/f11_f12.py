# -*- coding: utf-8 -*-
"""F11: crossing_errors level/ luminosity mismatch.  F12: gamma channel missing from CHANNELS."""
import sys, os, math, json
import numpy as np
RUN = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, os.path.join(RUN,"work","analysis"))
sys.path.insert(0, os.path.join(RUN,"work","widths"))
sys.path.insert(0, os.path.join(RUN,"work","response"))

scan = json.load(open(os.path.join(RUN,"outputs","data","signal_scan.json")))
th = np.array(scan['theta2'], dtype=float)
sig = scan['signal_per_fb_per_channel']
masses = np.array(scan['masses'], dtype=float)
nsamp = {float(k): float(v) for k, v in scan['nsamp_by_mass'].items()}

def s_rel(row):
    s = np.array([r['B']+r['K']+r['H'] for r in row], dtype=float)
    var = np.array([(r['B']*r.get('relerr_B',0))**2 + (r['K']*r.get('relerr_K',0))**2
                    + (r['H']*r.get('relerr_H',0))**2 for r in row], dtype=float)
    return s, np.where(s>0, np.sqrt(var)/np.maximum(s,1e-300), 0.0)

def crossings(x, y, level):
    lx=np.log10(x); ly=np.log10(np.maximum(y,1e-300)); lv=math.log10(level); out=[]
    for i in range(len(x)-1):
        f0=ly[i]-lv; f1=ly[i+1]-lv
        if f0*f1<0:
            t=f0/(f0-f1); out.append(10**(lx[i]+t*(lx[i+1]-lx[i])))
    return out

def rel_at(x, rel):
    return float(np.interp(math.log10(x), np.log10(th), rel))

def analyse(level, tag):
    n5=0; rels=[]; bad=[]
    for m in masses:
        s, rel = s_rel(sig[str(m)])
        cr = crossings(th, s, level)
        if not cr: continue
        lo = min(cr)
        r = rel_at(lo, rel)
        rels.append(r)
        if r > 0.05:
            n5 += 1; bad.append((m, r, nsamp.get(m, float('nan'))))
    rels=np.array(rels)
    print("  %-28s level=%.4g   n=%d  median rel=%.2f%%  max=%.2f%%  >5%%: %d/%d" % (
        tag, level, len(rels), 100*np.median(rels), 100*rels.max(), n5, len(masses)))
    return bad

print("="*90)
print("F11  crossing_errors() uses the PER-FB signal with level=2.3 (run_sensitivity.py:95,212,254)")
print("="*90)
print("  current code behaviour:")
analyse(2.3, "code check (L=1 fb^-1)")
bad25 = analyse(2.3/25.0, "L=25 fb^-1")
bad300 = analyse(2.3/300.0, "L=300 fb^-1")
print()
print("  masses whose long-lived-edge rel MC error exceeds 5%:")
for tag, bad in (("25 fb^-1", bad25), ("300 fb^-1", bad300)):
    print("   %s (%d):" % (tag, len(bad)))
    for (m, r, n) in sorted(bad):
        need = n*(r/0.05)**2
        print("     m=%8.5f  rel=%.3f%%  nsamp=%.0f  -> need ~%.0f" % (m, 100*r, n, need))
# short-lived edge: all masses with two crossings
print()
for L in (25.0, 300.0):
    hi_rels=[]; n2=0
    for m in masses:
        s, rel = s_rel(sig[str(m)])
        cr = crossings(th, s, 2.3/L)
        if len(cr)>=2:
            n2+=1; hi_rels.append(rel_at(max(cr), rel))
    hi_rels=np.array(hi_rels)
    print("  L=%3.0f fb^-1 short-lived edge: %d masses with two edges, median rel=%.1f%%, max=%.1f%%, >5%%: %d/%d" % (
        L, n2, 100*np.median(hi_rels), 100*hi_rels.max(), int((hi_rels>0.05).sum()), n2))
# also median long-lived rel for reference
for L in (1.0,25.0,300.0):
    rels=[]
    for m in masses:
        s, rel = s_rel(sig[str(m)])
        cr = crossings(th, s, 2.3/L)
        if cr: rels.append(rel_at(min(cr), rel))
    print("  long-lived edge L=%5.0f: median %.3f%% max %.3f%%" % (L, 100*np.median(rels), 100*np.max(rels)))

# ---------------------------------------------------------------- F12
print()
print("="*90)
print("F12  CHANNELS (bc5_lib.py:180) has no 'gamma'; channel_branching() renormalises over it")
print("="*90)
import bc5_lib as L
import scalar_widths as sw
w0 = sw.widths(0.1)
print("  widths(0.1) keys :", sorted(w0.keys()))
print("  CHANNELS         :", L.CHANNELS)
missing=[k for k in w0 if k not in L.CHANNELS]
print("  keys with width but not in CHANNELS:", {k: w0[k] for k in missing})
br, tot = L.channel_branching(0.1)
print("  0.1 GeV: sum(BR over CHANNELS) = %.6f  (gamma/total = %.6f)" % (sum(br.values()), w0['gamma']/tot))
print()
print("  %-10s %-12s %-12s %-12s %-12s %-14s" % ("m_S","Gamma_tot","Gamma_gamma","BR_gamma","sum_BR_ch","BR_gamma/sum"))
rows=[]
for m in masses:
    w = sw.widths(float(m)); tot = sum(w.values())
    brg = w.get('gamma',0.0)/tot
    brc = {c: w.get(c,0.0)/tot for c in L.CHANNELS}
    rows.append((float(m), tot, w.get('gamma',0.0), brg, sum(brc.values())))
for m in (0.021, 0.05, 0.1, 0.15, 0.2, 0.2118, 0.3, 0.5, 0.9884, 1.0, 1.5, 2.0, 3.5567, 10.5834, 62.49):
    r = [x for x in rows if abs(x[0]-m)<1e-9]
    if r:
        x=r[0]; print("  %-10.5f %-12.4e %-12.4e %-12.5f %-12.6f %-14.5f" % (x[0],x[1],x[2],x[3],x[4],x[3]/x[4]))
brg = np.array([x[3] for x in rows])
print("  max BR_gamma over the 67 scan masses = %.5f at m=%.5f GeV ; at m<0.2 GeV max %.5f" % (
    brg.max(), rows[int(brg.argmax())][0], brg[masses<0.2].max() if (masses<0.2).any() else 0))

# does the 'gg' channel produce charged daughters at those masses?
ok = L.load_hadron_lib(os.path.join(RUN,'work','magnus','hadron_lib'))
print("  hadron lib loaded:", ok)
if L.HADRON_LIB is not None:
    for ch in ('gg',):
        c = L.HADRON_LIB.chan.get(ch)
        print("  hadron_lib_%s mass nodes: min %.4f max %.4f (%d nodes)" % (ch, c['mass_nodes'].min(), c['mass_nodes'].max(), len(c['mass_nodes'])))
print("  charged multiplicity of the 'gg' fallback (s_decay_products) as a function of m_S:")
import random
rng = random.Random(3)
for m in (0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9884, 1.0, 1.5, 2.0, 5.0):
    nch = []; natt=0
    for _ in range(400):
        d = L.s_decay_products(m, 'gg', rng)
        nch.append(len(d)); natt += 1
    print("   m=%-8.4f  <n_charged> = %.3f   P(empty)=%.3f" % (m, float(np.mean(nch)), float(np.mean([x==0 for x in nch]))))
