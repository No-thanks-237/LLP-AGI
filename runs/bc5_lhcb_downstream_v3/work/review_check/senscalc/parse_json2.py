
import json
p=r'E:/AAAwork/LLP/code/SensCalc-main/Sensitivity domains/Scalar/LHCb-downstream-D/Sensitivity_Scalar_at_LHCb-downstream-D_BrhToSS=0.01_Nev=2.3_Npot=1.75e15.json'
d=json.load(open(p))
print(type(d), len(d))
for i,item in enumerate(d):
    print('--- entry',i)
    for k,v in item.items():
        if isinstance(v,list):
            print('   %-28s list len %d'%(k,len(v)))
        else:
            print('   %-28s = %r'%(k,v))
    # for sensitivity domains: number of sublists and mass ranges
    for k,v in item.items():
        if isinstance(v,list) and k!='Production modes':
            try:
                print('     ',k,'nsub',len(v),'npts',len(v[0]),'m range', v[0][0][0], v[0][-1][0], 'th2 range', v[0][0][1], v[0][-1][1])
            except Exception as e:
                print('     ',k,'?',e)
