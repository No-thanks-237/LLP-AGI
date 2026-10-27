import json, numpy as np, os
RUN=r'E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3'
sc=json.load(open(os.path.join(RUN,'outputs','data','signal_scan.json')))
print('keys',list(sc.keys()))
m=np.array(sc['masses']); th=np.array(sc['theta2'])
print('nmass',len(m),'masses',np.round(m,4))
print('theta2 range',th[0],th[-1],len(th))
sig=sc['signal_per_fb_per_channel']
print('nsamp_by_mass sample',list(sc.get('nsamp_by_mass',{}).items())[:5])
print('model',sc.get('model_tag'),'sigma_bb',sc.get('sigma_bb_pb'),'byc',sc.get('b_yield_corr'))
# at some masses print B/K/H at theta 1e-8 and the total curve near crossing
for mm in [1.0,2.0,4.0,4.3,4.4,5.0,30.0,62.0]:
    i=np.argmin(abs(m-mm)); row=sig[str(m[i])]
    j=np.argmin(abs(np.log(th)-np.log(1e-8)))
    r=row[j]
    print('m=%.4f th=%.2e B=%.4g K=%.4g H=%.4g epsB=%.3g epsK=%.3g epsH=%.3g relB=%.3g relH=%.3g'%(
        m[i],th[j],r['B'],r['K'],r['H'],r['eps_B'],r['eps_K'],r['eps_H'],r.get('relerr_B',0),r.get('relerr_H',0)))
