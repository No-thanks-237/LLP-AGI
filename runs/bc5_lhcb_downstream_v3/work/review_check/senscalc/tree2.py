import json
d=r'E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3/work/physics_review_20261025/'
t=json.load(open(d+'reference_2023_tree.json'))
print('sha',t['sha'])
print('url',t['url'])
print('truncated',t['truncated'])
print('n entries',len(t['tree']))
for e in t['tree'][:12]: print(e)
