
# -*- coding: utf-8 -*-
"""Cross-check: recompute the H-channel per-Higgs visible probability with the production
functions and compare with the value implied by signal_scan.json."""
import sys, os, math, json
import numpy as np
R = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
sys.path.insert(0, R+'/work/analysis'); sys.path.insert(0, R+'/work/response'); sys.path.insert(0, R+'/work/widths')
import bc5_lib as L, bc5_production as P, run_scan as RS
import random

scan = json.load(open(R+'/outputs/data/signal_scan.json'))
THE = np.array(scan['theta2'])
def scan_val(ms, t2, c):
    k = min(scan['signal_per_fb_per_channel'].keys(), key=lambda x: abs(float(x)-ms))
    y = np.maximum(np.array([r.get(c,0.0) for r in scan['signal_per_fb_per_channel'][k]]), 1e-300)
    return float(np.exp(np.interp(math.log(t2), np.log(THE), np.log(y))))

def per_scalar(ms, t2, p4, z0m, A):
    e,px,py,pz = p4
    ctau = L.sw.HBAR_C/(L.sw.gamma_total(ms)*t2)
    z = L.ZGRID
    alpha = ms/(ctau*pz) if pz>0 else 0.0
    w = np.zeros(len(z)); m = z>z0m
    w[m] = alpha*np.exp(-(z[m]-z0m)*alpha)
    return float(np.trapezoid(w*A, z))

smp = RS.load_sample(RS.SAMPLE_H)
ev, first = np.unique(smp['ev'], return_index=True); idxs = np.sort(first)
rng = random.Random(2024)
for ms, t2 in ((1.0, 4.453e-11), (1.5, 1.169e-10), (2.0, 1.262e-10)):
    o_sum = 0.0; n = 0; p_sum = 0.0; ns = 0
    for _ in range(3000):
        i = int(idxs[rng.randrange(len(idxs))])
        d = L.two_body(P.M_H, ms, ms, rng)
        if d is None: continue
        pair = L.boost(d, smp['e'][i], (smp['px'][i], smp['py'][i], smp['pz'][i]))
        probs = []
        for s4 in pair:
            s = (s4[1], s4[2], s4[3], s4[4], smp['x'][i], smp['y'][i], smp['z'][i])
            r = RS.decay_vertex_and_selection(ms, s, rng, ndecay=8)
            if r is None: continue
            p4, v0, A = r
            probs.append(per_scalar(ms, t2, p4, v0[2]*1e-3, A))
            p_sum += probs[-1]; ns += 1
        o = 1.0 - np.prod([1.0-min(max(pp,0.0),1.0) for pp in probs])
        o_sum += o; n += 1
    eps_scan = scan_val(ms, t2, 'H')/(1e3*P.SIGMA_H_TOTAL*P.B_H)
    print('m=%.2f t2=%.3e  <o_pair>=%.4e  eps_H=0.4<o>=%.4e   scan eps_H=%.4e  ratio=%.3f  (ns S/Higgs=%.2f)' % (
        ms, t2, o_sum/n, 0.4*o_sum/n, eps_scan, (0.4*o_sum/n)/eps_scan, ns/n))
