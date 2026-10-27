import json, numpy as np, math, os
RUN=r'E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3'
sc=json.load(open(os.path.join(RUN,'outputs','data','signal_scan.json')))
m=np.array(sc['masses']); th=np.array(sc['theta2']); sig=sc['signal_per_fb_per_channel']
def crossings(x,y,level):
    out=[]; lx=np.log10(x); ly=np.log10(np.maximum(y,1e-300))
    for i in range(len(x)-1):
        f0=ly[i]-math.log10(level); f1=ly[i+1]-math.log10(level)
        if f0*f1<0:
            t=f0/(f0-f1); out.append(10**(lx[i]+t*(lx[i+1]-lx[i])))
    return out
print('%9s %6s %10s %10s %8s %8s %8s %8s'%('m','ncr','th_lo','th_hi','Bfrac','Kfrac','Hfrac','s@th_lo/L'))
for mm in [1.0,1.3,1.6,1.99,4.0,4.3446,5.0823,6.9547,8.1355,2.0]:
    i=np.argmin(abs(m-mm)); row=sig[str(m[i])]
    s=np.array([r['B']+r['K']+r['H'] for r in row])
    for L in (25.0,300.0):
        cr=crossings(th,s,2.3/L)
        if not cr: 
            print('%9.4f %6d L=%3.0f NO CROSSING s(1e-8)=%.3g'%(m[i],len(cr),L,s[np.argmin(abs(np.log(th)-np.log(1e-8)))]))
            continue
        lo=min(cr)
        j=np.argmin(abs(np.log(th)-np.log(lo)))
        r=row[j]; tot=r['B']+r['K']+r['H']
        print('%9.4f %6d L=%3.0f lo=%.3e hi=%.3e B=%.3f K=%.3f H=%.3f s=%.3f'%(
            m[i],len(cr),L,lo,max(cr),r['B']/tot,r['K']/tot,r['H']/tot,tot*L))
