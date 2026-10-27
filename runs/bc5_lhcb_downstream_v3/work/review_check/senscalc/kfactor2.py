
import numpy as np, gzip
mh=125.21
rng=np.random.default_rng(12345)

def Kfactor(E,pvec,mS,nphi=6):
    E=np.asarray(E,float); pvec=np.asarray(pvec,float)
    p=np.sqrt(np.maximum(E*E-mh*mh,0.0))
    u=np.divide(pvec,p[:,None],out=np.zeros_like(pvec),where=p[:,None]>0)
    beta=np.divide(p,E,out=np.zeros_like(p),where=E>0); gam=1/np.sqrt(np.maximum(1-beta*beta,1e-12))
    pstar=np.sqrt(max(mh*mh/4-mS*mS,0.0)); Estar=mh/2
    N=len(E); tot=0.0; acc=0.0
    for k in range(nphi):
        cth=rng.uniform(-1,1,N); phi=rng.uniform(0,2*np.pi,N); sth=np.sqrt(1-cth*cth)
        # random orthonormal basis perpendicular to u
        a=rng.normal(size=(N,3)); a-=np.sum(a*u,axis=1)[:,None]*u; an=np.linalg.norm(a,axis=1)
        an=np.where(an>1e-9,an,1.0); e1=a/an[:,None]
        e2=np.cross(u,e1)
        n=(sth*np.cos(phi))[:,None]*e1+(sth*np.sin(phi))[:,None]*e2+cth[:,None]*u
        mu=np.sum(n*u,axis=1)
        E_lab=gam*(Estar+beta*pstar*mu)
        p_lab=pstar*n+(((gam-1)*mu+gam*beta*Estar/pstar)*pstar)[:,None]*u
        pz=p_lab[:,2]
        eta=0.5*np.log(np.maximum(E_lab+pz,1e-12)/np.maximum(E_lab-pz,1e-12))
        ok=(eta>2.0)&(eta<5.0)
        # weight: 1_{2<etaS<5} * mS/pz  (z-integral factor dropped)
        w=np.where(ok, mS/np.maximum(np.abs(pz),1e-12),0.0)
        tot+=float(np.sum(w)); acc+=float(np.sum(ok))
    return tot/(N*nphi), acc/(N*nphi)

# --- SensCalc table ---
p=r'E:/AAAwork/LLP/code/SensCalc-main/spectra/SM particles/DoubleDistr_LHC_h.txt'
d=np.loadtxt(p); th=np.unique(d[:,0]); EE=np.unique(d[:,1]); V=d[:,2].reshape(len(th),len(EE))
dth=np.gradient(th); dE=np.gradient(EE); W=V*dth[:,None]*dE[None,:]
ii,jj=np.nonzero(W>0); ws=W[ii,jj]; ws/=ws.sum()
N=200000; idx=rng.choice(len(ws),size=N,p=ws)
TH=th[ii[idx]]; EN=EE[jj[idx]]; rng.shuffle(TH); rng.shuffle(EN)
ps=np.sqrt(np.maximum(EN*EN-mh*mh,0)); ph=rng.uniform(0,2*np.pi,N)
PV=np.stack([ps*np.sin(TH)*np.cos(ph), ps*np.sin(TH)*np.sin(ph), ps*np.cos(TH)],axis=1)
# --- ours ---
Eo=[];P=[]; st=[]
with gzip.open(r'E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3/work/magnus/samples/higgs.txt.gz','rt') as f:
    for line in f:
        if line.startswith('#'): continue
        a=line.split()
        if int(a[10])!=62: continue
        Eo.append(float(a[5])); P.append([float(a[2]),float(a[3]),float(a[4])])
        if len(Eo)>=40000: break
Eo=np.array(Eo); P=np.array(P)
print('N_table_hists=%d  N_ours=%d'%(N,len(Eo)))
for mS in [5.08,13.02,33.37,53.42]:
    kt,ft=Kfactor(EN,PV,mS,4)
    ko,fo=Kfactor(Eo,P,mS,10)
    print('mS=%6.2f  K_table=%.4e (eta cut eff %.4f)  K_ours=%.4e (%.4f)  ratio=%.3f'%(mS,kt,ft,ko,fo,kt/ko))
