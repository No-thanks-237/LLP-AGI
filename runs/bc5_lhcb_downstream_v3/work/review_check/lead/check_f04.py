
import sys, os, math, random
import numpy as np
sys.path.insert(0, r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\analysis")
import run_scan as R, bc5_lib as L
L.load_hadron_lib(r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\magnus\hadron_lib")
smpB = R.load_sample(R.SAMPLE_B)
sel = np.isin(smpB['pdg'], [521,-521,511,-511,531,-531])
R.B_YIELD_CORR = R.REF_B_MESON_YIELD/(float(sel.sum())/float(R.N_GEN_BBAR))
th = np.exp(np.linspace(math.log(R.THETA2_MIN), math.log(R.THETA2_MAX), 240))
def crossing(ms, nsamp=2000, seed=4242):
    rng = random.Random(seed)
    tab = R.run_mass(ms, nsamp, rng, ndecay=8, verbose=False)
    s = np.array([R.signal_per_lumi(ms, tab, float(t), 1.0, R.SIGMA_BB_PB, R.SIGMA_INEL_PB,
                 R.N_GEN_BBAR, R.N_GEN_MINBIAS, R.N_GEN_HIGGS)['B'] +
                 R.signal_per_lumi(ms, tab, float(t), 1.0, R.SIGMA_BB_PB, R.SIGMA_INEL_PB,
                 R.N_GEN_BBAR, R.N_GEN_MINBIAS, R.N_GEN_HIGGS)['K'] +
                 R.signal_per_lumi(ms, tab, float(t), 1.0, R.SIGMA_BB_PB, R.SIGMA_INEL_PB,
                 R.N_GEN_BBAR, R.N_GEN_MINBIAS, R.N_GEN_HIGGS)['H'] for t in th])
    f = np.log10(np.maximum(s,1e-300))-math.log10(2.3/25.0)
    cr = []
    for i in range(len(th)-1):
        if f[i]*f[i+1] < 0:
            cr.append(10**(np.log10(th[i])+f[i]/(f[i]-f[i+1])*(np.log10(th[i+1])-np.log10(th[i]))))
    return min(cr) if cr else float('nan')
ms = 1.5
on = crossing(ms)
# switch the kaon flight off by stripping the offset (reproduces the old behaviour)
orig = L._passes_pair
def no_offspring(parts, x, y, z):
    return orig([(a[0],a[1],a[2],a[3],a[4],a[5],None) for a in parts], x, y, z)
L._passes_pair = no_offspring
off = crossing(ms)
L._passes_pair = orig
print("m_S=1.5: with K flight %.4e, without %.4e, ratio %.4f" % (on, off, on/off))
print("audit-kaon independent implementation: crossing +8.6% (ratio 1.086)")
