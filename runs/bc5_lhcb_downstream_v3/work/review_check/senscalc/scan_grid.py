
import struct, numpy as np
p=r'E:/AAAwork/LLP/code/SensCalc-main/phenomenology/Scalar/decay widths/ctauS.mx'
d=open(p,'rb').read()
# align at every offset 0..7 relative to file start, decode all doubles
best=[]
for shift in range(8):
    n=(len(d)-shift)//8
    a=np.frombuffer(d, dtype='<f8', count=n, offset=shift)
    # look for increasing runs within [0.02,63.5]
    good=np.isfinite(a) & (a>0.019) & (a<63.6)
    # find runs of >=25 consecutive good & strictly increasing & ratio<3
    i=0
    N=len(a)
    while i<N-1:
        if good[i]:
            j=i+1
            while j<N and good[j] and a[j]>a[j-1] and a[j]/a[j-1]<3.0:
                j+=1
            if j-i>=25:
                best.append((shift,i,j,float(a[i]),float(a[j-1]),j-i))
            i=j
        else:
            i+=1
best.sort(key=lambda t:-t[5])
for b in best[:8]:
    print('shift',b[0],'idx',b[1],'len',b[5],'from %.5g to %.5g'%(b[3],b[4]))
    shift,i0,j,n0,n1,L=best[0] if False else b
    n=(len(d)-shift)//8
    a=np.frombuffer(d, dtype='<f8', count=n, offset=shift)
    print('   first 8:', ['%.6g'%x for x in a[i0:i0+8]])
    print('   last 8:', ['%.6g'%x for x in a[j-8:j]])
