
import numpy as np, json, os, math
base=r'E:/AAAwork/LLP/code/SensCalc-main/Tabulated Nevents/Scalar/LHCb-downstream-D'
files={f.split('_at_')[0].replace('Nevents_Scalar_',''):f for f in os.listdir(base)}
print(sorted(files))
def load(prod):
    d=np.loadtxt(os.path.join(base,files[prod]))
    ms=np.unique(d[:,0]); th=np.unique(d[:,1])
    V=d[:,2].reshape(len(ms),len(th))
    return ms,th,V
def interp(prod,m,th2, cache={}):
    if prod not in cache: cache[prod]=load(prod)
    ms,th,V=cache[prod]
    if not (ms[0]<=m<=ms[-1] and th[0]<=th2<=th[-1]): return 0.0
    i=np.searchsorted(ms,m); j=np.searchsorted(th,th2)
    i=min(max(i,1),len(ms)-1); j=min(max(j,1),len(th)-1)
    lm1,lm2=math.log(ms[i-1]),math.log(ms[i]); lt1,lt2=math.log(th[j-1]),math.log(th[j])
    x=(math.log(m)-lm1)/(lm2-lm1); y=(math.log(th2)-lt1)/(lt2-lt1)
    def L(a,b): return math.log10(max(V[a,b],1e-300))
    lv=(1-x)*(1-y)*L(i-1,j-1)+(1-x)*y*L(i-1,j)+x*(1-y)*L(i,j-1)+x*y*L(i,j)
    return 10**lv
jp=json.load(open(r'E:/AAAwork/LLP/code/SensCalc-main/Sensitivity domains/Scalar/LHCb-downstream-D/Sensitivity_Scalar_at_LHCb-downstream-D_BrhToSS=0.01_Nev=2.3_Npot=1.75e15.json'))
pts=np.array(jp[0]['Sensitivity domains'][0],dtype=float)
prods=['B-mixing','B-quartic','Bs-quartic','h-quartic','Upsilon-mixing','Upsilon-quartic','Bremsstrahlung-AP-Central-mixing']
print('%-14s %-14s %-10s %-12s %-10s'%('m','t1','t2','N_total','N*0.4'))
res=[]
for k in range(0,len(pts),max(1,len(pts)//40)):
    m,t=pts[k]
    N=sum(interp(p,m,t) for p in prods)
    res.append((m,t,N,N*0.4))
for m,t,N,N4 in res:
    print('%-14.5f %-14.4e %-10s %-12.4g %-10.4g'%(m,t,'',N,N4))
