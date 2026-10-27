
import sys, os, math, random, json
import numpy as np
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\analysis")
os.environ['BC5_SAMPLE_H'] = os.environ.get('H_SAMPLE', 'higgs.txt.gz')
import run_scan as R, bc5_lib as L
L.load_hadron_lib(r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\magnus\hadron_lib")
smpB = R.load_sample(R.SAMPLE_B)
sel = np.isin(smpB['pdg'], [521,-521,511,-511,531,-531])
R.B_YIELD_CORR = R.REF_B_MESON_YIELD/(float(sel.sum())/float(R.N_GEN_BBAR))
print("sample:", R.SAMPLE_H)
out = {}
for ms in (5.082255,):
    rng = random.Random(20261016)
    tab = R.run_mass(ms, 12000, rng, ndecay=8, verbose=False)
    for t2 in (1e-12,):
        d = R.signal_per_lumi(ms, tab, t2, 1.0, R.SIGMA_BB_PB, R.SIGMA_INEL_PB,
                              R.N_GEN_BBAR, R.N_GEN_MINBIAS, R.N_GEN_HIGGS)
        out['%.3f_%.0e' % (ms, t2)] = (d['H'], d['relerr_H'])
print(json.dumps(out))