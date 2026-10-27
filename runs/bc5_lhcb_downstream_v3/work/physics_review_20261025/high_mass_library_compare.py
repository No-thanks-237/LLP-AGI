"""Exact common-mass common-kinematics comparison: current bb library vs SensCalc bb.
Only audit files are written. This isolates the conditional final-state model.
"""
from pathlib import Path
import json,random,math,sys,time,argparse
import numpy as np
HERE=Path(__file__).resolve().parent;RUN=HERE.parent.parent
sys.path.insert(0,str(RUN/'work/analysis'))
import bc5_lib as L,run_scan as R
ms=62.49;N=2200;ndecay=5
ap=argparse.ArgumentParser();ap.add_argument('--old',action='store_true');args=ap.parse_args()
outfile=HERE/('high_mass_library_compare_2023.json' if args.old else 'high_mass_library_compare.json')
L.load_hadron_lib(str(RUN/'work/magnus/hadron_lib'))
aa=np.loadtxt(HERE/('reference_2023_bb_62p49.csv' if args.old else 'reference_bb_62p49.csv'),delimiter=',')
nref=int(aa[:,0].max())+1
groups=[]
for i in range(nref):
    d=aa[aa[:,0]==i]
    groups.append([(float(x[5]),float(x[4]),float(x[1]),float(x[2]),float(x[3]),int(x[7])) for x in d])
sm=R.load_sample(R.SAMPLE_H);ix=np.flatnonzero(abs(sm['status'])==62)
pool=[];rng=random.Random(25210)
for i in ix:
    E=sm['e'][i];hp=tuple(sm[c][i] for c in ['px','py','pz'])
    mh=math.sqrt(max(E*E-sum(x*x for x in hp),1e-12))
    if mh<2*ms:continue
    d=L.boost(L.two_body(mh,ms,ms,rng),E,hp)
    for x in d:
        if x[4]<=0:continue
        eta=math.asinh(x[4]/max(math.hypot(x[2],x[3]),1e-30))
        if 2<eta<5:pool.append(tuple(x[1:5]))
    if len(pool)>N*5:break
rng=random.Random(25211);p4s=rng.sample(pool,N)
out={'mass':ms,'samples':N,'decays':ndecay,'reference_attempts':nref,'reference_version':'2023-10-25' if args.old else 'current local tables','results':{}}
for name in ['current','SensCalc_conditional_bb']:
    vals=[];rng=random.Random(25212);t=time.time()
    for p4 in p4s:
        a=np.zeros(len(L.ZGRID))
        for _ in range(ndecay):
            rest=L.HADRON_LIB.sample(ms,'bb',rng) if name=='current' else groups[rng.randrange(nref)]
            rest=L._rotate_all(rest,rng)
            b=L.boost([tuple(x[:5]) for x in rest],p4[0],p4[1:4])
            lab=[tuple(x)+(r[5],None) for x,r in zip(b,rest)]
            for k,z in enumerate(L.ZGRID):
                a[k]+=L._passes_pair(lab,z*p4[1]/p4[3],z*p4[2]/p4[3],z)
        vals.append(ms/p4[3]*np.trapezoid(a/ndecay,L.ZGRID))
    v=np.array(vals)
    out['results'][name]={'coefficient_conditional_eta':float(v.mean()),
        'relerr':float(v.std(ddof=1)/math.sqrt(len(v))/v.mean()),'seconds':time.time()-t}
    outfile.write_text(json.dumps(out,indent=2),encoding='utf8')
    print(name,out['results'][name],flush=True)
out['current_over_reference_bb_response']=out['results']['current']['coefficient_conditional_eta']/out['results']['SensCalc_conditional_bb']['coefficient_conditional_eta']
out['reference_bb_mean_Nch']=len(aa)/nref
outfile.write_text(json.dumps(out,indent=2),encoding='utf8')
print('ratio',out['current_over_reference_bb_response'],flush=True)
