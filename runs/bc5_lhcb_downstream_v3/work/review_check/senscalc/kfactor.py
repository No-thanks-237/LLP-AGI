
import numpy as np, gzip
mh=125.21
rng=np.random.default_rng(12345)
def Kfactor(E,pz,ptx,pty,mS,nphi=6):
    # h 4-mom (lab). decay to SS isotropic in h rest frame.
    p=np.sqrt(max(E*E-mh*mh,0.0))
    uz=np.divide(pz,p,out=np.zeros_like(pz),where=p>0)
    perp=np.sqrt(np.maximum(ptx*ptx+pty*pty,0.0))
    ux=np.divide(ptx,np.where(perp>0,perp,1.0),out=np.zeros_like(ptx),where=perp>0)
    uy=np.divide(pty,np.where(perp>0,perp,1.0),out=np.zeros_like(pty),where=perp>0)
    pstar=np.sqrt(max(mh*mh/4-mS*mS,0.0))
    tot=0.0; acc=0.0; cnt=0
    for k in range(nphi):
        phi=rng.uniform(0,2*np.pi,len(E))
        cth=rng.uniform(-1,1,len(E))
        sth=np.sqrt(1-cth*cth)
        # build orthonormal basis (u, v, w)
        # choose v = (u x zhat) normalized, w = u x v
        vx=-uy; vy=ux; vz=np.zeros_like(ux)
        vn=np.sqrt(vx*vx+vy*vy); good=vn>0
        vx=np.where(good,vx/np.where(good,vn,1),1.0); vy=np.where(good,vy/np.where(good,vn,1),0.0)
        wx=uy*vz-uz*vy; wy=uz*vx-ux*vz; wz=ux*vy-uy*vx
        # momentum of S in rest frame: direction n = sth cos phi * v + sth sin phi * w + cth * u
        nx=sth*np.cos(phi)*vx+sth*np.sin(phi)*wx+cth*ux
        ny=sth*np.cos(phi)*vy+sth*np.sin(phi)*wy+cth*uy
        nz=sth*np.cos(phi)*vz+sth*np.sin(phi)*wz+cth*uz
        Es=E/2.0
        gam=E/mh; beta=p/E
        pzS=Es*(beta+gam*nz)   # lz component of momentum
        # eta from angle between S momentum and beam: pS vector = Es*beta*u_hat? careful: boost
        Es2=Es*gam*(1+beta*nz)
        pxS=Es*(gam*beta*ux+ (nx - uz*nz*0 )*0 )  # placeholder
        # full boost: p_lab = p*_par + ...
        # p* components: pstar*(n.v) etc. Use vector form:
        pSx=pstar*nx + (Es)*gam*beta*ux
        pSy=pstar*ny + (Es)*gam*beta*uy
        pSz=pstar*nz + (Es)*gam*beta*uz
        pS2=pSx*pSx+pSy*pSy+pSz*pSz
        pSm=np.sqrt(pS2)
        eta=0.5*np.log((Es2+pSz)/np.maximum(Es2-pSz,1e-12))
        ok=(eta>2)&(eta<5)&(pSm>0)
        w=np.where(ok, mS/np.maximum(np.abs(pSz),1e-12),0.0)
        tot+=w.sum(); cnt+=len(E)
    return tot/(cnt)

def ours_events(maxn=None,seed=1):
    E=[];px=[];py=[];pz=[]
    with gzip.open(r'E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3/work/magnus/samples/higgs.txt.gz','rt') as f:
        for line in f:
            if line.startswith('#'): continue
            a=line.split()
            if int(a[10])!=62: continue
            E.append(float(a[5])); px.append(float(a[2])); py.append(float(a[3])); pz.append(float(a[4]))
            if maxn and len(E)>=maxn: break
    return np.array(E),np.array(px),np.array(py),np.array(pz)

# SensCalc table sample
p=r'E:/AAAwork/LLP/code/SensCalc-main/spectra/SM particles/DoubleDistr_LHC_h.txt'
d=np.loadtxt(p); th=np.unique(d[:,0]); EE=np.unique(d[:,1]); V=d[:,2].reshape(len(th),len(EE))
dth=np.gradient(th); dE=np.gradient(EE)
W=V*dth[:,None]*dE[None,:]
ii,jj=np.nonzero(W>0)
ws=W[ii,jj]; ws=ws/ws.sum()
N=200000
idx=rng.choice(len(ws),size=N,p=ws)
TH=th[ii[idx]]; EN=EE[jj[idx]]
rng.shuffle(TH); rng.shuffle(EN)
pS=np.sqrt(np.maximum(EN*EN-mh*mh,0))
pxS=pS*np.sin(TH); pzS=pS*np.cos(TH)  # random azimuth not needed (symmetric) but for boost basis use azimuth 0
# For the table we sample theta and E independently of phi; using phi=0 is fine (rotational symmetry) but the boost basis uses pt direction; build with phi random:
phiT=rng.uniform(0,2*np.pi,N)
ptS=pS*np.sin(TH)
pxT=ptS*np.cos(phiT); pyT=ptS*np.sin(phiT)
print('table: sampling %d h'%N)
for mS in [5.08,13.02,33.37,53.42]:
    k_tab=Kfactor(EN,pzS,pxT,pyT,mS,nphi=4)
    print('  mS=%.2f  K_table = %.4e'%(mS,k_tab))
EO,PXO,PYO,PZO=ours_events(20000)
print('ours: %d h'%len(EO))
for mS in [5.08,13.02,33.37,53.42]:
    k_o=Kfactor(EO,PZO,PXO,PYO,mS,nphi=8)
    print('  mS=%.2f  K_ours  = %.4e'%(mS,k_o))
