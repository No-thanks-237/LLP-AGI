import fitz, os, math, sys
d = r"E:\AAAwork\LLP-AGI\references\detector\arXiv-2305.13383v3"

def grab(name, colors, X1, DXPX, Y1, DYPX, LOGY1):
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
        pts = max(lst, key=len); pts.sort()
        row = []
        for mm in [0.15,0.3,0.5,1.0,2.0,3.0,5.0]:
            t = X1 + math.log10(mm)*DXPX
            b = min(pts, key=lambda q: abs(q[0]-t))
            row.append((mm, 10**(LOGY1 + (Y1-b[1])/DYPX), abs(b[0]-t)))
        res[lab] = row
    return res

r6 = grab("BC6-validation.pdf", {(0.0,0.0,1.0):"SensCalc", (0.67,0.0,0.0):"[1811.00930], 0.0,0.0)"}, 256.1, 192.4, 228.1, 34.9, -10.0)
print("BC6 (U_e^2 vs m_N):")
for k,v in r6.items(): print("  ", k, " | ".join(f"{m:.2f}:{val:.2e}(dx{a:.0f})" for m,val,a in v))
r1 = grab("BC1-validation.pdf", {(0.0,0.0,1.0):"SensCalc", (0.67,0.0,0.0):"[2011.05115]"}, 306.5, 186.7, 230.5, 14.64, -18.0)
print("BC1 (eps^2 vs m_V):")
for k,v in r1.items(): print("  ", k, " | ".join(f"{m:.2f}:{val:.2e}(dx{a:.0f})" for m,val,a in v))
