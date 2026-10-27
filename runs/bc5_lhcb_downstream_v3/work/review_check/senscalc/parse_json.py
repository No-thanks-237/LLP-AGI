
import json
p=r'E:/AAAwork/LLP/code/SensCalc-main/Sensitivity domains/Scalar/LHCb-downstream-D/Sensitivity_Scalar_at_LHCb-downstream-D_BrhToSS=0.01_Nev=2.3_Npot=1.75e15.json'
d=json.load(open(p))
print(type(d), len(d))
for i,item in enumerate(d):
    print('--- entry',i,'keys:',list(item.keys()))
    for k,v in item.items():
        if isinstance(v,list):
            print('   ',k,'list len',len(v), 'first:', (v[0] if len(v)<3 else (v[0], '...')))
        else:
            print('   ',k,'=',repr(v)[:100])
