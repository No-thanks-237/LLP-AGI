
import struct
p=r'E:/AAAwork/LLP/code/SensCalc-main/phenomenology/Scalar/decay widths/ctauS.mx'
d=open(p,'rb').read()
i=d.find(b'Stage 1')+8
off=i+3
vals=[]
while off+9<=len(d) and d[off]==1:
    vals.append(struct.unpack('<d', d[off+1:off+9])[0]); off+=9
print('n',len(vals))
for k,v in enumerate(vals[:80]):
    print(k,'%.10g'%v)
print('...last 20:')
for k,v in enumerate(vals[-20:]):
    print(len(vals)-20+k,'%.10g'%v)
