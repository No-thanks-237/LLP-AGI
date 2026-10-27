import gzip, numpy as np
R = r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream"
rows = []
with gzip.open(R + "/work/magnus/samples/bbar_heavy.txt.gz", 'rt') as f:
    for i, ln in enumerate(f):
        if i >= 2000000: break
        a = ln.split()
        if len(a) >= 10 and a[1] in ('521','-521','511','-511'):
            rows.append((float(a[2]), float(a[3]), float(a[4]), float(a[5])))
D = np.array([r[:3]+(r[3],) for r in rows]); p = np.linalg.norm(D[:,:3],axis=1); pz = D[:,3]
print("B mesons:", len(D), "(full sample scanned)")
for thr in (100, 200, 300, 500, 1000):
    print(f"  |p| >> {thr:5d} GeV: {np.sum(p>thr):7d} ({np.sum(p>thr)/len(D)*100:.4f}%)   pz>thr: {np.sum(pz>thr)}")
print("  max |p| =", p.max().round(1), " GeV; max pz =", pz.max().round(1), " GeV")
print("  median |p| =", np.median(p).round(1), " GeV")
# p_window for a 1 GeV H with ctau = 1.01e-4 m as a function of pz
import sys
sys.path.insert(0, R + "/work/response")
import twohdm_lib as L
A = np.ones(13).reshape(1,-1)
for pzv in (100, 200, 300, 400, 600, 1000, 2000):
    print(f"  pz={pzv:5d} GeV: p_window(ctau=1.01e-4 m, A=1) = {float(L.p_window(1.0, np.array([pzv]), np.array([0.0]), A, 1.01e-4)[0]):.3e}")
