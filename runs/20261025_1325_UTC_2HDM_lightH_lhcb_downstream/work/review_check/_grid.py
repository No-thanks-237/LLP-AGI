import collections
W = r"E:/AAAwork/LLP-AGI/references/tools/FORESEE/Models/2HDM_H/model"
def load(fn):
    d = {}
    for ln in open(W+"/br/"+fn):
        a = ln.split()
        if len(a) < 3: continue
        d[(float(a[0]), float(a[1]))] = float(a[2])
    return d
for fn in ("ga_ga.txt","pi_pi.txt","k_k.txt","g_g.txt"):
    d = load(fn)
    ms = sorted(set(k[0] for k in d)); ts = sorted(set(k[1] for k in d))
    print(fn, "n=",len(d), "nmass=",len(ms), "ntanb=",len(ts), "mass", ms[0], ms[-1], "tanb", ts[0], ts[-1])
d = load("pi_pi.txt"); ms = sorted(set(k[0] for k in d)); ts=sorted(set(k[1] for k in d))
print("mass grid:", [round(m,4) for m in ms][:12], "...", [round(m,4) for m in ms][-6:])
print("tanb grid sample:", [round(t,4) for t in ts][:8], "...", [round(t,4) for t in ts][-4:])
ct = {}
for ln in open(W+"/ctau.txt"):
    a = ln.split()
    if len(a)>=3: ct[(float(a[0]),float(a[1]))] = float(a[2])
print("ctau n=", len(ct), "sample", list(ct.items())[:3])
