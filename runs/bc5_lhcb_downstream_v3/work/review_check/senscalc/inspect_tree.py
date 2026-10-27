import json
import numpy as np
d=r'E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3/work/physics_review_20261025/'
t=json.load(open(d+'reference_2023_tree.json'))
print('tree type',type(t))
print('keys' , list(t.keys()) if isinstance(t,dict) else len(t))
