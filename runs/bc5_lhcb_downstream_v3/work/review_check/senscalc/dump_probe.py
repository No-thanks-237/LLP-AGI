
import struct
p=r'E:/AAAwork/LLP/code/SensCalc-main/Acceptances/LHCb-downstream-D/Acceptance_LHCb-downstream-D_for_Scalar_decayDescription_Default.m'
d=open(p,'rb').read()
i=d.find(b'Npot')
print('Npot at',i)
print(d[i-8:i+60])
print(d[i-8:i+60].hex())
# search for 'Ph' followed by type tag and real
import re
for m in re.finditer(rb'Ph', d):
    off=m.start()
    seg=d[off:off+24]
    print(off, seg.hex())
