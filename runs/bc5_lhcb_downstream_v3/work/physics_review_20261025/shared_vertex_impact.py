"""Isolate per-daughter vs per-kaon flight; no production files are modified."""
from pathlib import Path
import json,math,random,sys
import numpy as np
HERE=Path(__file__).resolve().parent;RUN=HERE.parent.parent
sys.path.insert(0,str(RUN/'work/analysis'))
import run_scan as R,bc5_lib as L,run_sensitivity as RS
scan=json.loads((RUN/'outputs/data/signal_scan.json').read_text())
cr25=RS.crossing_errors(scan,2.3,25); cr300=RS.crossing_errors(scan,2.3,300)
smp=R.load_sample(R.SAMPLE_B)
R.B_YIELD_CORR=R.REF_B_MESON_YIELD/(np.isin(smp['pdg'],[521,-521,511,-511,531,-531]).sum()/R.N_GEN_BBAR)
original=L.s_decay_products
def shared(ms,ch,rng,s_p4=None):
    d=original(ms,ch,rng,s_p4=s_p4)
    if ch!='KK' or not d or d[0][6] is None:return d
    # Each neutral-kaon decay returns either zero or exactly two charged products here.
    # Keep the first product's already drawn Exp(1), using it for that kaon's whole pair;
    # discard (but still consume) the extra draw, to keep the two runs paired in all RNG calls.
    out=[]
    for i in range(0,len(d),2):
        off=d[i][6]
        out.extend([tuple(a[:6])+(off,) for a in d[i:i+2]])
    return out
results={}
for ms in [1.0,1.44937,2.0]:
    row={}
    for mode in ['current','shared_vertex']:
        L.s_decay_products=original if mode=='current' else shared
        R._np_rng=np.random.default_rng(25107)
        tab=R.run_mass(ms,5000,random.Random(25106),ndecay=8,verbose=False)
        row[mode]={}
        for lum,cr in [(25,cr25),(300,cr300)]:
            theta=cr[ms][0]
            d=R.signal_per_lumi(ms,tab,theta,lum,R.SIGMA_BB_PB,R.SIGMA_INEL_PB,
                R.N_GEN_BBAR,R.N_GEN_MINBIAS,R.N_GEN_HIGGS)
            row[mode][str(lum)]={'theta2':theta,'signal':sum(d[c] for c in ['B','K','H']),
                                 'B':d['B'],'K':d['K'],'H':d['H']}
    results[str(ms)]=row
    print(ms,row,flush=True)
L.s_decay_products=original
(HERE/'shared_vertex_impact.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
print('DONE',flush=True)
