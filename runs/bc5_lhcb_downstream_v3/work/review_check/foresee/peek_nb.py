
import json, glob, re
for f in ['temp/foresee/Models/DarkHiggs/DarkHiggs.ipynb','temp/foresee/Tutorial/DarkPhoton/DarkPhoton.ipynb','temp/foresee/Models/UpPhilic/UpPhilic.ipynb','temp/foresee/Models/2HDM_H/2HDM_H.ipynb']:
    try:
        nb=json.load(open(f))
    except Exception as e:
        print('skip',f,e); continue
    print('\n########',f)
    for cell in nb.get('cells',[]):
        src=''.join(cell.get('source',[]))
        if 'add_production_2bodydecay' in src or 'add_production_3bodydecay' in src:
            print('--- cell ---')
            print(src[:1500])
