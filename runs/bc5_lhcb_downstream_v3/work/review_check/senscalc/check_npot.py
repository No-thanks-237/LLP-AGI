
import numpy as np, json, os, math
base=r'E:/AAAwork/LLP/code/SensCalc-main/Tabulated Nevents/Scalar/LHCb-downstream-D'
# build lookup for a given production file
def load(prod):
    fn=[f for f in os.listdir(base) if f.startswith('Nevents_Scalar_'+prod+'_')][0]
    d=np.loadtxt(os.path.join(base,fn))
    return d, fn
for prod, masses in [('h-quartic',[13.004,30.053,33.389,53.410]),('B-mixing',[2.0,3.0,4.5])]:
    d,fn=load(prod)
    print('===',prod,fn, d.shape)
    ms=np.unique(d[:,0]); th=np.unique(d[:,1])
    print('  m grid n=%d range %.4f-%.4f; coupling grid n=%d range %.2e-%.2e'%(len(ms),ms[0],ms[-1],len(th),th[0],th[-1]))
    # json boundary points
    jp=json.load(open(r'E:/AAAwork/LLP/code/SensCalc-main/Sensitivity domains/Scalar/LHCb-downstream-D/Sensitivity_Scalar_at_LHCb-downstream-D_BrhToSS=0.01_Nev=2.3_Npot=1.75e15.json'))
    pts=np.array(jp[0]['Sensitivity domains'][0],dtype=float)
    for m0 in masses:
        i=np.argmin(np.abs(pts[:,0]-m0)); m,th2v=pts[i]
        # interpolate log-log bilinear manually: find bracket
        def val(mm,tt):
            im=np.searchsorted(ms,mm); it=np.searchsorted(th,tt)
            im=min(max(im,1),len(ms)-1); it=min(max(it,1),len(th)-1)
            m1,m2=ms[im-1],ms[im]; t1,t2=th[it-1],th[it]
            # build dict
            idx={(round(row[0],10),round(row[1],20)):row[2] for row in d}
            try:
                v11=idx[(round(m1,10),round(t1,20))]; v21=idx[(round(m2,10),round(t1,20))]
                v12=idx[(round(m1,10),round(t2,20))]; v22=idx[(round(m2,10),round(t2,20))]
            except KeyError:
                return None
            # bilinear in log space
            x=(math.log(mm)-math.log(m1))/(math.log(m2)-math.log(m1))
            y=(math.log(tt)-math.log(t1))/(math.log(t2)-math.log(t1))
            lv=(1-x)*(1-y)*math.log10(max(v11,1e-300))+(1-x)*y*math.log10(max(v12,1e-300))+x*(1-y)*math.log10(max(v21,1e-300))+x*y*math.log10(max(v22,1e-300))
            return 10**lv
        v=val(m,th2v)
        print('  m=%.4f th2=%.4e  N_table(m,th2)=%s   -> N*eps0.4=%s'%(m,th2v,v,None if v is None else v*0.4))
