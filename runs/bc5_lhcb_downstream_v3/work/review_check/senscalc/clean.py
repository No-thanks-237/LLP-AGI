
import re
d=open(r'E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3/work/review_check/senscalc/manual.txt','rb').read()
s=d.decode('latin1')
s=''.join(c if (32<=ord(c)<127 or c=='\n') else ' ' for c in s)
s=re.sub(r'[ \t]+',' ',s)
open(r'E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3/work/review_check/senscalc/manual_ascii.txt','w',encoding='ascii').write(s)
print(len(s.splitlines()))
