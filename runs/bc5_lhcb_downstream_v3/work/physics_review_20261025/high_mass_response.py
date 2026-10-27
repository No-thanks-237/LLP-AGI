"""Common-response comparison of native and reference H->SS spectra.

Reference spectra are diagnostic inputs only. Production sources/results are untouched.
The reported coefficient is <1_eta int dz A(z) m/pz> per produced scalar.
All comparisons use the SAME current conditional decay library and reference BRs.
"""
from pathlib import Path
import json, math, random, sys, time, argparse
import numpy as np
from scipy.interpolate import RegularGridInterpolator

HERE=Path(__file__).resolve().parent; RUN=HERE.parent.parent
sys.path.insert(0,str(RUN/'work/analysis'));sys.path.insert(0,str(RUN/'work/widths'))
import run_scan as R, bc5_lib as L, reference_br_model as W
L.sw.widths=W.widths; L.sw.gamma_total=W.gamma_total
# W uses gamma_ll from the same module; its widths implementation remains independent.
L._CH_CACHE.clear()
ap=argparse.ArgumentParser();ap.add_argument('--samples',type=int,default=3500)
ap.add_argument('--decays',type=int,default=6);args=ap.parse_args()
L.load_hadron_lib(str(RUN/'work/magnus/hadron_lib'))
sm=R.load_sample(R.SAMPLE_H);ix=np.flatnonzero(abs(sm['status'])==62)
hE=sm['e'][ix]; hp=np.column_stack([sm[c][ix] for c in ['px','py','pz']])
hm=np.sqrt(np.maximum(hE*hE-(hp*hp).sum(axis=1),1e-12))
raw=np.loadtxt(HERE/'reference_H_scalar_spectrum.csv',delimiter=',')
massnodes=np.unique(raw[:,0]); meshes={}
for mass in massnodes:
    a=raw[raw[:,0]==mass];a=a[np.lexsort((a[:,2],a[:,1]))]
    t=np.unique(a[:,1]);en=np.unique(a[:,2]);f=a[:,3].reshape(len(t),len(en))
    meshes[mass]=(t,en,f)
out={'samples':args.samples,'ndecay':args.decays,'comparison':'long-lifetime common decay response, per scalar',
     'reference_interpolation':'linear density in mass/theta/E; normalized to one scalar over full grid',
     'masses':{}}
def save(): (HERE/'high_mass_response.json').write_text(json.dumps(out,indent=2),encoding='utf8')

def scalar_mesh(m):
    i=min(max(np.searchsorted(massnodes,m)-1,0),len(massnodes)-2);m0,m1=massnodes[i:i+2]
    t,en,f0=meshes[m0];f1=meshes[m1][2];u=(m-m0)/(m1-m0)
    f=np.maximum((1-u)*f0+u*f1,0)
    norm=np.trapezoid(np.trapezoid(f,en,axis=1),t)
    tlo,thi=2*np.arctan(np.exp(-5)),2*np.arctan(np.exp(-2))
    tc=np.r_[tlo,t[(t>tlo)&(t<thi)],thi];ec=np.r_[m,en[en>m]]
    interp=RegularGridInterpolator((t,en),f,bounds_error=False,fill_value=0)
    tt,ee=np.meshgrid(tc,ec,indexing='ij');fc=interp(np.column_stack([tt.ravel(),ee.ravel()])).reshape(tt.shape)
    area=np.diff(tc)[:,None]*np.diff(ec)[None,:]
    avg=(fc[:-1,:-1]+fc[1:,:-1]+fc[:-1,1:]+fc[1:,1:])/4
    weights=area*avg;fraction=weights.sum()/norm
    # Exact rejection sampler for the bilinear density on each rectangle.
    flat=weights.ravel();flat/=flat.sum(); ni,nj=weights.shape
    pts=[];rng=np.random.default_rng(25200+round(m*100))
    while len(pts)<args.samples:
        k=rng.choice(len(flat),p=flat);i,j=divmod(k,nj)
        ceiling=max(fc[i,j],fc[i+1,j],fc[i,j+1],fc[i+1,j+1])
        while True:
            a,b=rng.random(2)
            dens=(1-a)*(1-b)*fc[i,j]+a*(1-b)*fc[i+1,j]+(1-a)*b*fc[i,j+1]+a*b*fc[i+1,j+1]
            if rng.random()*ceiling<=dens:break
        theta=tc[i]+a*(tc[i+1]-tc[i]);energy=ec[j]+b*(ec[j+1]-ec[j]);ph=rng.uniform(-math.pi,math.pi)
        p=math.sqrt(max(energy*energy-m*m,0));pts.append((energy,p*math.sin(theta)*math.cos(ph),p*math.sin(theta)*math.sin(ph),p*math.cos(theta)))
    return np.array(pts),fraction,norm

def native(m):
    rng=np.random.default_rng(25201);ct=rng.uniform(-1,1,len(ix));ph=rng.uniform(-math.pi,math.pi,len(ix))
    st=np.sqrt(1-ct*ct);p=np.sqrt(np.maximum(hm*hm/4-m*m,0));v=np.column_stack([p*st*np.cos(ph),p*st*np.sin(ph),p*ct])
    beta=hp/hE[:,None];g=hE/hm;b2=(beta*beta).sum(axis=1);bp=(beta*v).sum(axis=1)
    allpts=[]
    for sign in [1,-1]:
        energy=g*(hm/2+sign*bp);k=(g-1)*sign*bp/np.maximum(b2,1e-30)+g*hm/2
        mom=sign*v+k[:,None]*beta;eta=np.arcsinh(mom[:,2]/np.maximum(np.hypot(mom[:,0],mom[:,1]),1e-30))
        sel=(eta>2)&(eta<5)&(hm>2*m)
        allpts.append(np.column_stack([energy[sel],mom[sel]]))
    pts=np.vstack(allpts);fraction=len(pts)/(2*R.N_GEN_HIGGS)
    chosen=rng.choice(len(pts),args.samples,replace=False)
    return pts[chosen],fraction

original=L._passes_pair
for m in [5.082255,13.023095,33.371213,53.419632,62.49]:
    started=time.time(); row={}
    pref,fref,norm=scalar_mesh(m);pnative,fnative=native(m)
    for label,pts,frac in [('native',pnative,fnative),('reference_spectrum',pref,fref)]:
        vals=[];geom=[];variants=[]
        # Capture a paired profile with the dipole position from the paper, leaving all
        # other physics and every random draw identical. No UT cut is added.
        def paired(daughters,x,y,z):
            value=original(daughters,x,y,z)
            old=L.dk.DEFAULT_CONFIG.copy()
            L.dk.DEFAULT_CONFIG.update(fieldZMin=3.5,fieldZMax=7.5)
            try: alt=original(daughters,x,y,z)
            finally:L.dk.DEFAULT_CONFIG.clear();L.dk.DEFAULT_CONFIG.update(old)
            variants.append(alt)
            return value
        L._passes_pair=paired;rng=random.Random(25202)
        for q,p4 in enumerate(pts):
            variants.clear(); a=L.acceptance_profile(m,p4,(0,0,0),rng,ndecay=args.decays)
            aa=np.array(variants).reshape(-1,len(L.ZGRID)).sum(axis=0)/args.decays if variants else np.zeros(len(L.ZGRID))
            vals.append([m/p4[3]*np.trapezoid(a,L.ZGRID),m/p4[3]*np.trapezoid(aa,L.ZGRID)])
            geom.append(m/p4[3]*(L.ZGRID[-1]-L.ZGRID[0]))
        v=np.array(vals);g=np.array(geom)
        row[label]={'fraction_in_eta':frac,'geometric_lifetime_coefficient':float(frac*g.mean()),
                    'coefficient_current_kernel':float(frac*v[:,0].mean()),
                    'coefficient_paper_magnet_position':float(frac*v[:,1].mean()),
                    'relative_sampling_error':float(v[:,0].std(ddof=1)/math.sqrt(len(v))/v[:,0].mean()),
                    'mean_E_conditional':float(pts[:,0].mean()),'mean_inv_pz_conditional':float(np.mean(1/pts[:,3]))}
    L._passes_pair=original
    row['reference_mesh_norm']=float(norm)
    row['native_over_reference_common_response']=row['native']['coefficient_current_kernel']/row['reference_spectrum']['coefficient_current_kernel']
    row['paper_magnet_over_current_native']=row['native']['coefficient_paper_magnet_position']/row['native']['coefficient_current_kernel']
    out['masses'][str(m)]=row;save()
    print(m,json.dumps(row),'seconds',time.time()-started,flush=True)
print('DONE',flush=True)
