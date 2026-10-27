
import json, os, sys, math, random
import numpy as np
RUN = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3"
sys.path.insert(0, os.path.join(RUN, "work", "analysis"))
import run_scan as R
smpB = R.load_sample(R.SAMPLE_B)
sel = np.isin(smpB['pdg'], [521, -521, 511, -511, 531, -531])
R.B_YIELD_CORR = R.REF_B_MESON_YIELD/(float(sel.sum())/float(R.N_GEN_BBAR))
print('corr', R.B_YIELD_CORR)
scan = json.load(open(os.path.join(RUN, "outputs", "data", "signal_scan.json")))
m = 5.082255
for nsamp in (600, 1900):
    rng = random.Random(20261016)
    tab = R.run_mass(m, nsamp, rng, ndecay=8, verbose=True)
    print("nsamp", nsamp)
    for i in [0, 30, 60, 90, 120, 150, 200, 239]:
        t2 = scan["theta2"][i]
        d = R.signal_per_lumi(m, tab, t2, 1.0, R.SIGMA_BB_PB, R.SIGMA_INEL_PB,
                              R.N_GEN_BBAR, R.N_GEN_MINBIAS, R.N_GEN_HIGGS)
        row = scan["signal_per_fb_per_channel"]["5.082255"][i]
        print("  t2=%.4e mine B=%.3g K=%.3g H=%.4g | stored B=%.3g K=%.3g H=%.4g  ratioH=%.3f" % (
            t2, d['B'], d['K'], d['H'], row['B'], row['K'], row['H'],
            d['H']/row['H'] if row['H'] else float('nan')))
