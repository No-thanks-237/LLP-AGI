
import struct
p=r'E:/AAAwork/LLP/code/SensCalc-main/phenomenology/Scalar/decay widths/ctauS.mx'
d=open(p,'rb').read()
i=d.find(b'Stage 1')+8
best=None
for shift in range(0,12):
    off=i+shift
    vals=[]
    while off+9<=len(d) and d[off]==1:
        v=struct.unpack('<d', d[off+1:off+9])[0]
        if not (abs(v)<1e14): break
        vals.append(v); off+=9
    if best is None or len(vals)>len(best[1]):
        best=(shift,vals)
    if vals: print('shift',shift,'n',len(vals),'first', ['%.8g'%x for x in vals[:6]])
print('BEST shift',best[0],'n',len(best[1]))
vals=best[1]
for k in range(0,min(30,len(vals))):
    print(k,'%.12g'%vals[k])
