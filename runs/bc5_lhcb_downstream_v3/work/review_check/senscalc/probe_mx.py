
import struct, re, os

def runs(d, start, n=12):
    out=[]
    for k in range(n):
        off=start+8*k
        if off+8>len(d): break
        out.append(struct.unpack('<d', d[off:off+8])[0])
    return out

def probe(path, limit=6):
    d=open(path,'rb').read()
    print('===', os.path.basename(path), len(d), 'bytes')
    cnt=0
    for m in re.finditer(b'\x05\x06\x02', d):
        off=m.start()+3
        vals=runs(d, off, 10)
        ok=sum(1 for v in vals if v==v and abs(v)<1e12 and (v==0 or 1e-12<abs(v)))
        if ok>=8:
            print(' at',m.start(), ['%.6g'%v for v in vals])
            cnt+=1
            if cnt>=limit: break
    print(' found',cnt)

base=r'E:/AAAwork/LLP/code/SensCalc-main/'
probe(base+'Auxiliary data/Scalar/LHCb-downstream-D/AcceptanceAveraged-Scalar.m', 3)
probe(base+'phenomenology/Scalar/decay widths/ctauS.mx', 6)
