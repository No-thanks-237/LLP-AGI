import sys, math, json
import numpy as np
sys.path.insert(0, r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream/work/response")
import twohdm_lib as L
Z = L.ZGRID; A = np.ones(len(Z))
print("p_window(A=1) for mH=1 GeV as a function of pz and ctau:")
for ctm in (1.147e-4, 5.09, 0.4436):
    print(f"  ctau = {ctm:.4g} m:")
    for pz in (5, 10, 50, 200, 1000):
        print(f"     pz={pz:5.0f} GeV: p_window={float(L.p_window(1.0, np.array([pz]), np.array([0.0]), A.reshape(1,-1), ctm)[0]):.4e}"
              f"   (c*tau*bg = {ctm*pz/1.0:.3g} m)")
print()
# --- B -> X_s H kinematics from the actual B parent sample
import gzip
P = r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream/work/magnus/samples/bbar_heavy.txt.gz"
data = []
with gzip.open(P, 'rt') as f:
    for i, ln in enumerate(f):
        if i >= 400000: break
        a = ln.split()
        if len(a) < 10: continue
        if a[1] in ('521','-521','511','-511'):
            data.append((float(a[2]), float(a[3]), float(a[4]), float(a[5])))
D = np.array(data)
print("B mesons read:", D.shape, " mean |p| =", np.mean(np.linalg.norm(D[:,:3],axis=1)).round(2), "GeV")
rng = np.random.default_rng(7)
px,py,pz,E = D.T
mB = 5.27934
# B -> K H two-body, isotropic in the B rest frame
def h_kin(mH):
    mK = 0.493677
    ps = math.sqrt(max((mB**2-(mK+mH)**2)*(mB**2-(mK-mH)**2),0.0))/(2*mB)
    Es = math.sqrt(ps*ps+mH*mH)
    ct = rng.uniform(-1,1,size=len(D)); ph = rng.uniform(0,2*np.pi,size=len(D))
    st = np.sqrt(1-ct*ct)
    hx, hy, hz = ps*st*np.cos(ph), ps*st*np.sin(ph), ps*ct
    b = np.linalg.norm(D[:,:3],axis=1)/E
    g = E/mB
    # boost along the B direction
    ux, uy, uz = px/np.linalg.norm(D[:,:3],axis=1), py/np.linalg.norm(D[:,:3],axis=1), pz/np.linalg.norm(D[:,:3],axis=1)
    bp = hx*ux+hy*uy+hz*uz
    gb = g*b
    k = (g-1.0)*bp + gb*Es          # p' = p + k*u   (Lorentz boost along u with beta=b)
    pxp = hx + k*ux
    pyp = hy + k*uy
    pzp = hz + k*uz
    Ep = g*(Es + b*bp)
    return pxp,pyp,pzp,Ep
for mH, edge in ((1.0, 357.03), (3.0, None), (3.50112, 273.37)):
    hx,hy,hz,En = h_kin(mH)
    pt = np.hypot(hx,hy); p = np.sqrt(pt**2+hz**2)
    eta = 0.5*np.log((p+hz)/(p-hz))
    ok = (eta>2)&(eta<5)
    bgz = hz/mH
    print(f"\nmH={mH}: H in 2<eta<5: {ok.mean()*100:.2f}% of B decays")
    for q in (10,25,50,75,90):
        print(f"   beta*gamma_z (pz/mH) p{q} = {np.percentile(bgz[ok],q):.2f}", end="")
        if q==50: print("  |  p = ", f"{np.percentile(p[ok],q):.2f} GeV")
        else: print()
