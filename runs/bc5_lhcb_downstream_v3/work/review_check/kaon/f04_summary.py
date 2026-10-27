
import json, numpy as np, math
R = r"E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3"
res = json.load(open(R+'/work/review_check/kaon/f04_results.json'))
scan = json.load(open(R+'/outputs/data/signal_scan.json'))
THE = np.array(scan['theta2']); SC = scan['signal_per_fb_per_channel']
print('mass  brKK   r_KK_comb   Ds/s(B)   Ds/s(H)   Ds/s(tot)   slope   dth2/th2   (short-edge Ds/s)')
for key in sorted(res, key=float):
    a = res[key]; ms = a['mass']; brkk = a['br_KK']
    e = a['long']; fB = e['B']['f']; fH = e['H']['f']; rB = e['B']['r']; rH = e['H']['r']
    wB = e['B']['dss']*0; wBf = e['wBfrac']; wHf = 1-wBf
    dssB = fB*(rB-1); dssH = fH*(rH-1); tot = wBf*dssB + wHf*dssH
    rcomb = (wBf*fB*rB + wHf*fH*rH)/(wBf*fB + wHf*fH)
    _key = min(SC.keys(), key=lambda k: abs(float(k)-ms))
    s_all = np.array([r.get('B',0)+r.get('K',0)+r.get('H',0) for r in SC[_key]])
    i = int(np.argmin(np.abs(np.log(THE/a['edges']['long'])))); lo=max(0,i-2); hi=min(len(THE)-1,i+2)
    slope = (math.log(max(s_all[hi],1e-300))-math.log(max(s_all[lo],1e-300)))/(math.log(THE[hi])-math.log(THE[lo]))
    es = a['short']; dssB_s = es['B']['f']*(es['B']['r']-1); dssH_s = es['H']['f']*(es['H']['r']-1)
    tot_s = es['wBfrac']*dssB_s + (1-es['wBfrac'])*dssH_s
    print('%.2f  %.3f  %.4f      %+.4f   %+.4f   %+.4f     %.2f    %+.3f     %+.3f' % (
        ms, brkk, rcomb, dssB, dssH, tot, slope, -tot/slope, tot_s))
    print('       (per-source dss raw: B %.4f H %.4f ; Bq r=%.3f f=%.3f)' % (e['B']['dss'], e['H']['dss'], e['Bq']['r'], e['Bq']['f']))