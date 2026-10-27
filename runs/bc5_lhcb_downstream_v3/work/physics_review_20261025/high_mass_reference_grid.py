"""Separate the reference's tabulated response from its H->SS spectrum.
Uses the same log-coordinate, log-acceptance interpolation as BlockEpsilonDecay.
Input reference tables are diagnostic only; no production files are changed.
"""
from pathlib import Path
import json,sys,math
import numpy as np
from scipy.interpolate import RegularGridInterpolator
HERE=Path(__file__).resolve().parent;RUN=HERE.parent.parent
sys.path.insert(0,str(RUN/'work/widths'))
import reference_br_model as W,scalar_widths as sw
a=np.loadtxt(HERE/'reference_acceptance_high.csv',delimiter=',')
a=a[np.lexsort(tuple(a[:,i] for i in [3,2,1,0]))]
axes=[np.unique(a[:,i]) for i in range(4)]
eps=(a[:,4]*a[:,5]).reshape(tuple(len(x) for x in axes))
interp=RegularGridInterpolator(tuple(np.log(x) for x in axes),np.log(np.maximum(eps,1e-90)),bounds_error=False,fill_value=-207.2326583694641)
lin=RegularGridInterpolator(tuple(np.log(x) for x in axes),eps,bounds_error=False,fill_value=0)
zgrid=np.linspace(1,2.5,31)
raw=np.loadtxt(HERE/'reference_H_scalar_spectrum.csv',delimiter=',');meshes={}
for mass in np.unique(raw[:,0]):
    b=raw[raw[:,0]==mass];b=b[np.lexsort((b[:,2],b[:,1]))]
    tt=np.unique(b[:,1]);ee=np.unique(b[:,2]);ff=b[:,3].reshape(len(tt),len(ee))
    meshes[mass]=(tt,ee,ff)
massnodes=np.array(sorted(meshes))
def eps_profile(m,t,E,use_linear=False):
    ans=[];ip=lin if use_linear else interp
    for z in zgrid:
        p=np.column_stack([np.full_like(t,m),t,E,np.full_like(t,z)])
        v=ip(np.log(p));ans.append(v if use_linear else np.exp(v))
    return np.array(ans).T
def spectrum(m,t,E,log=False):
    i=min(max(np.searchsorted(massnodes,m)-1,0),len(massnodes)-2);m0,m1=massnodes[i:i+2]
    tt,ee,f0=meshes[m0];f1=meshes[m1][2]
    u=(m-m0)/(m1-m0);ff=(1-u)*f0+u*f1
    norm=np.trapezoid(np.trapezoid(ff,ee,axis=1),tt)
    if log:
        # The actual SensCalc spectrum interpolates log f in log(m,theta,E).
        u=np.log(m/m0)/np.log(m1/m0);lf=(1-u)*np.log(np.maximum(f0,1e-92))+u*np.log(np.maximum(f1,1e-92))
        ip=RegularGridInterpolator((np.log(tt),np.log(ee)),lf,bounds_error=False,fill_value=-211.8378285554522)
        return np.exp(ip(np.c_[np.log(t),np.log(E)]))
    ip=RegularGridInterpolator((tt,ee),ff,bounds_error=False,fill_value=0)
    return ip(np.c_[t,E])/norm
h=np.loadtxt(RUN/'work/magnus/samples/higgs.txt.gz');h=h[abs(h[:,-1])==62]
E=h[:,5];hp=h[:,2:5];hm=np.sqrt(np.maximum(E*E-(hp*hp).sum(1),1e-12));be=hp/E[:,None];g=E/hm;b2=(be*be).sum(1)
rng=np.random.default_rng(25220);ct=rng.uniform(-1,1,len(h));ph=rng.uniform(-np.pi,np.pi,len(h));st=np.sqrt(1-ct*ct)
out={}
for m in [5.082255,13.023095,33.371213,53.419632,62.49]:
    p=np.sqrt(np.maximum(hm*hm/4-m*m,0));v=np.c_[p*st*np.cos(ph),p*st*np.sin(ph),p*ct];bp=(be*v).sum(1)
    cnat=0;cnatlin=0
    for sign in [1,-1]:
        kk=(g-1)*sign*bp/np.maximum(b2,1e-30)+g*hm/2;mom=sign*v+kk[:,None]*be;en=g*(hm/2+sign*bp)
        theta=np.arctan2(np.hypot(mom[:,0],mom[:,1]),mom[:,2]);sel=(theta>=axes[1][0])&(theta<=axes[1][-1]);t=theta[sel];ee=en[sel]
        prof=eps_profile(m,t,ee);pl=eps_profile(m,t,ee,True)
        cnat+=float(np.sum(m/mom[sel,2]*np.trapezoid(prof,zgrid,axis=1)))/(2*200000)
        cnatlin+=float(np.sum(m/mom[sel,2]*np.trapezoid(pl,zgrid,axis=1)))/(2*200000)
    theta=np.geomspace(axes[1][0],axes[1][-1],220);ee=np.geomspace(m+1e-5,6500,450)
    tt,en=np.meshgrid(theta,ee,indexing='ij');t=tt.ravel();enflat=en.ravel();prof=eps_profile(m,t,enflat)
    geom=m/(np.sqrt(enflat*enflat-m*m)*np.cos(t))
    response=geom*np.trapezoid(prof,zgrid,axis=1)
    row={'coefficient_native_spectrum_reference_acceptance':cnat,
         'coefficient_native_spectrum_linear_epsilon':cnatlin}
    for label,log in [('linear_spectrum',False),('log_spectrum',True)]:
        f=spectrum(m,t,enflat,log);vv=(f*response).reshape(tt.shape)
        row['coefficient_reference_spectrum_and_acceptance_'+label]=float(np.trapezoid(np.trapezoid(vv,ee,axis=1),theta))
    out[str(m)]=row
    (HERE/'high_mass_reference_grid.json').write_text(json.dumps(out,indent=2),encoding='utf8')
    print(m,json.dumps(row),flush=True)
print('DONE',flush=True)
