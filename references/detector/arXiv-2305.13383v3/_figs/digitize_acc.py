import fitz, os, math
d = r"E:\AAAwork\LLP-AGI\references\detector\arXiv-2305.13383v3"
doc = fitz.open(os.path.join(d,"acceptances.pdf")); p = doc[0]
dr = p.get_drawings()

X0, SCALE_X, M0 = 154.7, 113.05, 1.0     # x pixel of m=1, px per unit mass
Y0, SCALE_Y, LOG0 = 306.4, 40.585, 10.0  # y pixel of 1e10, px per decade

def m_of(x): return M0 + (x - X0)/SCALE_X
def logv_of(y): return LOG0 + (Y0 - y)/SCALE_Y

curves = {}
for g in dr:
    pts = [it[1] for it in g["items"]] + [g["items"][-1][2]]
    if len(pts) < 50: continue
    col = g.get("color")
    if col is None: continue
    col = tuple(round(c,2) for c in col)
    # continuity: fraction of consecutive segments sharing an endpoint
    items = [(it[1], it[2]) for it in g["items"] if it[0]=="l"]
    cont = sum(1 for i in range(len(items)-1) if abs(items[i][1].x-items[i+1][0].x)<0.3 and abs(items[i][1].y-items[i+1][0].y)<0.3)
    frac = cont/max(1,len(items)-1)
    key = (col, "solid" if frac>0.5 else "dashed")
    curves.setdefault(key, []).append((pts, frac))

names = {(0.0,0.0,1.0):"I0", (0.67,0.0,0.0):"I1", (0.0,0.67,0.0):"I2", (0.0,0.0,0.0):"I3"}
print("m     " + "  ".join(f"{names.get(c,'?')}-{st}" for (c,st) in sorted(curves)))
table = {}
for (col,st), lst in sorted(curves.items()):
    pts = max(lst, key=lambda t: len(t[0]))[0]
    # keep only points inside frame
    pts = [(q.x,q.y) for q in pts if 47.3<=q.x<=573.0 and 22.1<=q.y<=347.0]
    pts.sort()
    row = []
    for mm in [1.0,1.5,2.0,2.5,3.0,3.5,4.0]:
        target = X0 + (mm-M0)*SCALE_X
        best = min(pts, key=lambda q: abs(q[0]-target))
        row.append(10**logv_of(best[1]) if abs(best[0]-target)<6 else float('nan'))
    table[(names.get(col,'?'),st)] = row
    print(f"{names.get(col,'?')}-{st}", "  ".join(f"{v:.2e}" for v in row))
print("m =", [1.0,1.5,2.0,2.5,3.0,3.5,4.0])
print()
for st in ["solid","dashed"]:
    if ("I0",st) in table and ("I2",st) in table:
        print(st, "I2/I0:", [f"{a/b:.2e}" for a,b in zip(table[("I2",st)], table[("I0",st)])])
        print(st, "I3/I2:", [f"{a/b:.2f}" for a,b in zip(table[("I3",st)], table[("I2",st)])])
