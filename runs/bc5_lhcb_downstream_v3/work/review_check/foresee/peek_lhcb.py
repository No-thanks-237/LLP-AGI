
import json, glob, re
for f in glob.glob('temp/foresee/**/*.ipynb', recursive=True):
    try: nb=json.load(open(f))
    except Exception: continue
    hits=[]
    for ci,cell in enumerate(nb.get('cells',[])):
        src=''.join(cell.get('source',[]))
        if re.search(r'LHCb', src, re.I):
            for line in src.split('\n'):
                if re.search(r'LHCb', line, re.I):
                    hits.append((ci,line.strip()[:300]))
    if hits:
        print('#'*8, f)
        for ci,l in hits[:8]:
            print(' cell',ci,':',l)
