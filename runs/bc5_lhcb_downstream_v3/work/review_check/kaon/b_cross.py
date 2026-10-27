
import sys, os, math, json
import numpy as np
R = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, R+'/work/analysis'); sys.path.insert(0, R+'/work/response'); sys.path.insert(0, R+'/work/widths')
import bc5_lib as L, bc5_production as P, run_scan as RS
scan = json.load(open(R+'/outputs/data/signal_scan.json'))
THE=np.array(scan['theta2'])
def scan_val(ms,t2,c):
    k=min(scan['signal_per_fb_per_channel'].keys(), key=lambda x: abs(float(x)-ms))
    y=np.maximum(np.array([r.get(c,0.0) for r in scan['signal_per_fb_per_channel'][k]]),1e-300)
    return float(np.exp(np.interp(math.log(t2), np.log(THE), np.log(y))))
# composition of the B sample
smp=RS.load_sample(RS.SAMPLE_B)
sel=np.isin(smp['pdg'],[521,-521,511,-511,531,-531]); idxs=np.where(sel)[0]
pdg=smp['pdg'][idxs]
frac={}
for code,bm in ((521,'B+'),(-521,'B+'),(511,'B0'),(-511,'B0'),(531,'Bs'),(-531,'Bs')):
    frac[bm]=frac.get(bm,0.0)+float((pdg==code).sum())/len(pdg)
print('B sample composition', {k:round(v,4) for k,v in frac.items()})
for ms,t2 in ((1.0,4.453e-11),(1.5,1.169e-10),(2.0,1.262e-10)):
    Sb=sum(frac[bm]*sum(P.br_B_to_S(ms,bm).values()) for bm in ('B+','B0','Bs'))
    epsB=scan_val(ms,t2,'B')/(1e3*5.6e8)
    byc=1.918/ (float(sel.sum())/800000)
    p_uncond = epsB/(0.4*1.8855*byc*Sb*t2)
    print('m=%.2f t2=%.3e  SumBR=%.4f  B_YIELD_CORR=%.4f  scan p_uncond=%.4e'%(ms,t2,Sb,byc,p_uncond))
