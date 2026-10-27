import subprocess, math, os
RUN = r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream"
FMO = r"E:/AAAwork/LLP-AGI/references/tools/FORESEE/Models/2HDM_H/model"
EXE = RUN + r"/work/decay/build/light_scalar_decay.exe"
CH = ['gaga','ee','mumu','tautau','pipi','KK','4pi','cc','ss','gg']
pts = [(0.01,1e10),(0.01,1e5),(0.01,1000.0),(0.01,1.0)]
inp = "\n".join(f"{m} {t}" for m,t in pts)+"\n"
pr = subprocess.run([EXE], input=inp, capture_output=True, text=True)
for ln in pr.stdout.splitlines():
    a = ln.split()
    if len(a)<13: continue
    m,t = float(a[0]), float(a[1]); v=[float(x) for x in a[2:13]]
    tot=v[0]
    print(f"mH={m} tanb={t:g}: GammaTot={tot:.6e}")
    for i,ch in enumerate(CH):
        print(f"    {ch:6} Gamma={v[1+i]:.6e}  BR={v[1+i]/tot:.6e}")
# FORESEE values at (0.01, 1e10)
FILES = {'gaga':'ga_ga.txt','ee':'e_e.txt','mumu':'mu_mu.txt','tautau':'tau_tau.txt',
         'pipi':'pi_pi.txt','KK':'k_k.txt','4pi':'4pi.txt','cc':'c_c.txt','ss':'s_s.txt','gg':'g_g.txt'}
print()
for ch,fn in FILES.items():
    best=None
    for ln2 in open(os.path.join(FMO,'br',fn)):
        a=ln2.split()
        if len(a)<3: continue
        if abs(float(a[0])-0.01)<1e-12 and float(a[1])>=9e9: best=float(a[2])
    print(f"  FORESEE {ch:6} at (0.01,tanb=1e10) = {best}")
