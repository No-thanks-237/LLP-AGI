# -*- coding: utf-8 -*-
"""Audit 1: compare the FORESEE 2HDM_H model tables (author Light_scalar_decay output)
with our compiled author code work/decay/build/light_scalar_decay.exe.

Comparison metric: |ln ratio| per channel; the FORESEE ctau uses hbar*c = 2.0e-16 GeV m
(our DESIGN convention is 1.973269804e-16), so we divide the ratio out explicitly.
"""
import subprocess, math, json, os, statistics

RUN = r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream"
FMO = r"E:/AAAwork/LLP-AGI/references/tools/FORESEE/Models/2HDM_H/model"
EXE = RUN + r"/work/decay/build/light_scalar_decay.exe"
HBARC_OURS, HBARC_FORESEE = 1.973269804e-16, 2.0e-16
CH = ['gaga','ee','mumu','tautau','pipi','KK','4pi','cc','ss','gg']   # stdout order
FILES = {'gaga':'ga_ga.txt','ee':'e_e.txt','mumu':'mu_mu.txt','tautau':'tau_tau.txt',
         'pipi':'pi_pi.txt','KK':'k_k.txt','4pi':'4pi.txt','cc':'c_c.txt','ss':'s_s.txt','gg':'g_g.txt'}

def load_br(fn):
    d = {}
    for ln in open(os.path.join(FMO,'br',fn)):
        a = ln.split()
        if len(a) < 3: continue
        d[(float(a[0]), float(a[1]))] = float(a[2])
    return d
def load_ctau():
    d = {}
    for ln in open(os.path.join(FMO,'ctau.txt')):
        a = ln.split()
        if len(a) >= 3: d[(float(a[0]), float(a[1]))] = float(a[2])
    return d

br_tab = {ch: load_br(f) for ch, f in FILES.items()}
ctau_tab = load_ctau()
pts = sorted(set().union(*[set(t) for t in br_tab.values()], set(ctau_tab)))
print(f"FORESEE table points: {len(pts)} (union), channel table sizes:",
      {ch: len(v) for ch, v in br_tab.items()})

# ---- run our executable on all points (one process, batch stdin)
inp = "\n".join(f"{m!r} {t!r}" for m, t in pts) + "\n"
pr = subprocess.run([EXE], input=inp, capture_output=True, text=True, timeout=1800)
lines = [l for l in pr.stdout.splitlines() if l.strip()]
print("exe rc =", pr.returncode, " stdout lines =", len(lines), " (expected", len(pts), ")")
if pr.stderr.strip(): print("exe stderr head:", pr.stderr.splitlines()[:3])
ours = {}
for ln in lines:
    a = ln.split()
    if len(a) < 13: continue
    m, t = float(a[0]), float(a[1])
    v = [float(x) for x in a[2:13]]
    ours[(m, t)] = dict(mH=m, tanb=t, GammaTot=v[0],
                        **{ch: v[1+i] for i, ch in enumerate(CH)})

# ---- compare
def stats(ratios):
    r = [abs(x) for x in ratios if x == x and math.isfinite(x)]
    r.sort()
    return dict(n=len(r), max=r[-1] if r else float('nan'),
                median=statistics.median(r) if r else float('nan'),
                p90=r[int(0.9*len(r))-1] if r else float('nan'),
                mean=statistics.fmean(r) if r else float('nan'))

report = {}
for ch in CH:
    key = 'Gamma_gaga' if ch == 'gaga' else None
    ratios = []
    for (m, t), brf in br_tab[ch].items():
        if (m, t) not in ours: continue
        o = ours[(m, t)]
        if o['GammaTot'] <= 0: continue
        bro = o[ch]/o['GammaTot']
        if brf <= 0 or bro <= 0:
            continue       # below threshold: both must be ~0; handled separately
        ratios.append(math.log(bro/brf))
    report[ch] = dict(stats=stats(ratios), n_zero_foresee=sum(1 for v in br_tab[ch].values() if v <= 0))
# count "threshold" mismatches (one side zero, other not)
thresh = {}
for ch in CH:
    n = 0
    for (m, t), brf in br_tab[ch].items():
        if (m, t) not in ours: continue
        o = ours[(m, t)]
        bro = o[ch]/o['GammaTot'] if o['GammaTot'] > 0 else 0.0
        if (brf <= 1e-20) != (bro <= 1e-20): n += 1
    thresh[ch] = n
# ---- ctau
cr, cr_raw = [], []
for (m, t), cf in ctau_tab.items():
    if (m, t) not in ours: continue
    o = ours[(m, t)]
    if o['GammaTot'] <= 0 or cf <= 0: continue
    ct_ours = HBARC_OURS/o['GammaTot']
    cr.append(math.log(ct_ours/cf*HBARC_FORESEE/HBARC_OURS))   # FORESEE convention
    cr_raw.append(math.log(ct_ours/cf))
report['ctau'] = dict(stats=stats(cr), stats_after_hbarc_correction=stats(cr_raw))

print()
print(f"{'channel':8} {'N':>6} {'max|ln r|':>10} {'median':>10} {'p90':>10}   note")
for ch in CH:
    s = report[ch]['stats']
    print(f"{ch:8} {s['n']:6d} {s['max']:10.3e} {s['median']:10.3e} {s['p90']:10.3e}"
          f"   threshold-mismatch pts: {thresh[ch]}")
s = report['ctau']['stats']
print(f"{'ctau':8} {s['n']:6d} {s['max']:10.3e} {s['median']:10.3e} {s['p90']:10.3e}   "
      f"(FORESEE hbar*c=2e-16 applied); without that correction: max={report['ctau']['stats_after_hbarc_correction']['max']:.3e}")

# where are the largest deviations?
print()
print("largest |ln ratio| points per channel (mH, tanb, ours, FORESEE):")
for ch in CH:
    worst = []
    for (m, t), brf in br_tab[ch].items():
        if (m, t) not in ours: continue
        o = ours[(m, t)]
        if o['GammaTot'] <= 0 or brf <= 0: continue
        bro = o[ch]/o['GammaTot']
        if bro <= 0: continue
        worst.append((abs(math.log(bro/brf)), m, t, bro, brf))
    worst.sort(reverse=True)
    for w in worst[:2]:
        print(f"  {ch:6} mH={w[1]:7.4f} tanb={w[2]:10.5g} ours={w[3]:.6e} FORESEE={w[4]:.6e} |lnr|={w[0]:.4f}")
json.dump({'report': report, 'threshold_mismatch': thresh,
           'n_points': len(pts)}, open(RUN + '/work/review_check/foresee_vs_exe.json','w'), indent=1)
print("\nwrote work/model/audit/foresee_vs_exe.json")
