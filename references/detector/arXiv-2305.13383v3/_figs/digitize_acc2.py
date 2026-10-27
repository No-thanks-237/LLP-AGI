import fitz, os, math
d = r"E:\AAAwork\LLP-AGI\references\detector\arXiv-2305.13383v3"
doc = fitz.open(os.path.join(d,"acceptances.pdf")); p = doc[0]
dr = p.get_drawings()
X0, SX, M0 = 154.7, 113.05, 1.0
Y0, SY, LOG0 = 306.4, 40.585, 10.0
names = {(0.0,0.0,1.0):"I0", (0.67,0.0,0.0):"I1", (0.0,0.67,0.0):"I2", (0.0,0.0,0.0):"I3"}
out = {}
for g in dr:
    if len(g["items"]) < 50: continue
    col = g.get("color")
    if col is None: continue
    col = tuple(round(c,2) for c in col)
    dashes = g.get("dashes")
    style = "dashed" if (dashes and dashes not in ("[] 0","() 0")) else "solid"
    pts = [(it[1].x, it[1].y) for it in g["items"] if it[0]=="l"]
    pts += [(g["items"][-1][2].x, g["items"][-1][2].y)]
    out.setdefault((names.get(col,"?"+str(col)), style), []).append(pts)
for k,v in out.items():
    print(k, [len(x) for x in v], [ (round(min(q[0] for q in x),1), round(max(q[0] for q in x),1)) for x in v])
print()
for k,v in sorted(out.items()):
    pts = v[0]
    pts = [(x,y) for x,y in pts if 47.3<=x<=573.0 and 22.1<=y<=347.0]
    pts.sort()
    row=[]
    for mm in [0.5,1.0,1.5,2.0,3.0,4.0]:
        t = X0+(mm-M0)*SX
        best = min(pts, key=lambda q: abs(q[0]-t))
        row.append((mm, 10**(LOG0+(Y0-best[1])/SY)))
    print(k, " ".join(f"m={m}:{v:.2e}" for m,v in row))
