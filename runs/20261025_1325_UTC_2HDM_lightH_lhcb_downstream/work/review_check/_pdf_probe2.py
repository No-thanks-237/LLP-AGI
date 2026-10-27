import fitz
p = r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream/inputs/references/LightH.pdf"
doc = fitz.open(p); pg = doc[0]
d = pg.get_text("dict")
print("spans:")
for b in d['blocks']:
    for l in b.get('lines', []):
        for s in l['spans']:
            t = s['text'].strip()
            if t:
                print(f"  {[round(x,1) for x in s['bbox']]} size={s['size']:.1f} text={t!r}")
print()
dr = pg.get_drawings()
print("n drawings:", len(dr))
from collections import Counter
cols = Counter()
for it in dr:
    for sub in it['items']:
        pass
    cols[(it.get('color'), it.get('fill'), it.get('dashes'), it.get('width'))] += 1
for k, v in list(cols.items())[:20]:
    print("  ", k, v)
