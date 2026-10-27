
import json, glob, re
seen=set()
for f in glob.glob('temp/foresee/**/*.ipynb', recursive=True):
    try: nb=json.load(open(f))
    except Exception: continue
    for cell in nb.get('cells',[]):
        src=''.join(cell.get('source',[]))
        if 'set_detector' in src:
            # print the blocks
            for m in re.finditer(r'set_detector\(', src):
                s=src[m.start():m.start()+700]
                key=s[:200]
                if key in seen: continue
                seen.add(key)
                print('#'*10, f)
                print(s)
                print()
