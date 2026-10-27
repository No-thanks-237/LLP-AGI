import fitz, os, math, statistics
d = r"E:\AAAwork\LLP-AGI\references\detector\arXiv-2305.13383v3"

def grab(name, colors, X1, DXPX, Y1, DYPX, LOGY1, masses):
    doc = fitz.open(os.path.join(d,name)); p = doc[0]
    out = {}
    for g in p.get_drawings():
        col = g.get("color")
        if col is None or len(g["items"]) < 40: continue
        col = tuple(round(c,2) for c in col)
        if col not in colors: continue
        pts = [(it[1].x, it[1].y) for it in g["items"] if it[0]=="l"] + [(g["items"][-1][2].x, g["items"][-1][2].y)]
        out.setdefault(colors[col], []).append(pts)
    res = {}
    for lab, lst in out.items():
        pts = max(lst, key=len)
        row = []
        for mm in masses:
            t = X1 + math.log10(mm)*DXPX
            ys = [y for x,y in pts if abs(x-t) < 1.2]
            if not ys: row.append((mm, float("nan"))); continue
            row.append((mm, 10**(LOGY1 + (Y1-statistics.median(ys))/DYPX)))
        res[lab] = row
    return res

m6=[0.2,0.3,0.5,1.0,2.0,3.0,5.0]
r6 = grab("BC6-validation.pdf", {(0.0,0.0,1.0):"SensCalc", (0.67,0.0,0.0):"ref 1811.00930"}, 256.1, 192.4, 228.1, 34.9, -10.0, m6)
print("BC6: U_e^2")
for k,v in r6.items(): print("  ", k, " | ".join(f"{m:.2f}:{val:.2e}" for m,val in v))
m1=[0.1,0.2,0.3,0.5,0.8,1.0,1.5,2.0,3.0]
r1 = grab("BC1-validation.pdf", {(0.0,0.0,1.0):"SensCalc", (0.67,0.0,0.0):"ref 2011.05115"}, 306.5, 186.7, 230.5, 14.64, -18.0, m1)
print("BC1: eps^2")
for k,v in r1.items(): print("  ", k, " | ".join(f"{m:.2f}:{val:.2e}" for m,val in v))
print("ratio (ref/SensCalc) BC6:", [f"{a/b:.2f}" for (m,a),(_,b) in zip(r6['ref 1811.00930'], r6['SensCalc'])])
print("ratio (ref/SensCalc) BC1:", [f"{a/b:.2f}" for (m,a),(_,b) in zip(r1['ref 2011.05115'], r1['SensCalc'])])
