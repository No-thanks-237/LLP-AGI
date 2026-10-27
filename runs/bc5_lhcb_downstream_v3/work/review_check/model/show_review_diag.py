import json,os
p=r'E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\physics_review_20261024\additional_results.json'
d=json.load(open(p))
print(json.dumps(d.get('higgs_acceptance_diagnostic'),indent=1)[:6000])
