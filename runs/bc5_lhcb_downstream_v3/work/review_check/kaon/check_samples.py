
import sys, gzip, numpy as np, json
RUN = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
SAM = RUN + "/work/magnus/samples"
MAG = RUN + "/work/magnus/outputs"

def load(path):
    ev, pdg, px, py, pz, e = [],[],[],[],[],[]
    op = gzip.open if path.endswith('.gz') else open
    with op(path,'rt') as f:
        for line in f:
            if not line or line[0]=='#': continue
            v=line.split()
            if len(v)<10: continue
            ev.append(int(v[0])); pdg.append(int(v[1]))
            px.append(float(v[2])); py.append(float(v[3])); pz.append(float(v[4])); e.append(float(v[5]))
    return dict(ev=np.array(ev), pdg=np.array(pdg), px=np.array(px), py=np.array(py), pz=np.array(pz), e=np.array(e))

for name,path in [
    ("PROD minbias_light (used by production)", SAM+"/minbias_light.txt"),
    ("nocut2 minbias_light", MAG+"/nocut2/out/minbias_light.txt.gz"),
    ("v3fix3 minbias_light", MAG+"/v3fix3/out/minbias_light.txt.gz"),
    ("v3fix minbias_light", MAG+"/v3fix/out/minbias_light.txt.gz"),
]:
    d = load(path)
    p3 = np.sqrt(d['px']**2+d['py']**2+d['pz']**2)
    eta = 0.5*np.log((p3+d['pz'])/np.maximum(p3-d['pz'],1e-300))
    m = np.sqrt(np.maximum(d['e']**2-p3**2,0))
    print("###", name, "rows", len(d['pdg']), "events", len(np.unique(d['ev'])))
    for p in np.unique(d['pdg']):
        s = d['pdg']==p
        print("  pdg %5d : n=%7d  eta[%.4f,%.4f]  pT[%.3f,%.3f]  p[%.2f,%.2f]"%(p,s.sum(),eta[s].min(),eta[s].max(),
              np.min(np.hypot(d['px'][s],d['py'][s])), np.max(np.hypot(d['px'][s],d['py'][s])), p3[s].min(), p3[s].max()))
    print("  all eta range [%.5f, %.5f]"%(eta.min(),eta.max()))
