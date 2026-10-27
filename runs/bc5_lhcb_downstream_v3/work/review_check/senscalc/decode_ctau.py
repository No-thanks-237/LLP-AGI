
import struct
p=r'E:/AAAwork/LLP/code/SensCalc-main/phenomenology/Scalar/decay widths/ctauS.mx'
d=open(p,'rb').read()
i=d.find(b'Stage 1')
print('Stage1 at',i)
start=i+8
print('bytes after:', d[start:start+6].hex())
# try: sequence of 0x01 tag + 8 bytes double
off=start
vals=[]
while off+9<=len(d) and d[off]==1:
    v=struct.unpack('<d', d[off+1:off+9])[0]
    vals.append(v); off+=9
print('nvals',len(vals))
for k in range(0,min(40,len(vals))):
    print(k, '%.10g'%vals[k])
