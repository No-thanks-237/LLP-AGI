import fitz, sys, math, statistics
sys.stdout.reconfigure(encoding="utf-8")
p_ = r"E:\AAAwork\LLP-AGI\references\detector\Gorkavenko 等 - 2024 - LHCb potential to discover long-lived new physics particles with lifetimes above 100 ps.pdf"
doc = fitz.open(p_); p = doc[9]
cands=[]
for k,g in enumerate(p.get_drawings()):
    r = g["rect"]; col = g.get("color")
    if col is None: continue
    col = tuple(round(c,2) for c in col)
    if r.y0>338 and r.y1<480 and r.x0>80 and r.x1<300 and len(g["items"])>15:
        cands.append((k,col,g.get("dashes"),len(g["items"]),r))
for k,col,d,n,r in cands:
    print(k, col, d, n, [round(v,1) for v in (r.x0,r.y0,r.x1,r.y1)])
print("total", len(cands))
X1, DX, Y1, DY, LY = 144.6, 77.1, 366.8, 9.523, -7.0
def val(y): return 10**(LY + (y-Y1)/DY)
for k,col,d,n,r in cands:
    if col[2] < 0.5: continue   # keep blue-ish
    g = p.get_drawings()[k]
    pts=[(it[1].x,it[1].y) for it in g["items"] if it[0]=="l"]+[(g["items"][-1][2].x,g["items"][-1][2].y)]
    pts=[(x,y) for x,y in pts if 90<x<287 and 347<y<469]
    if not pts: continue
    row=[]
    for mm in [0.5,0.7,1.0,1.5,2.0,3.0,5.0]:
        t = X1+DX*math.log10(mm)
        ys=[y for x,y in pts if abs(x-t)<1.0]
        row.append((mm, val(statistics.median(ys)) if ys else float("nan")))
    print(f"curve {k} n={n} dash={d}: " + " ".join(f"{m}:{v:.2e}" for m,v in row))
