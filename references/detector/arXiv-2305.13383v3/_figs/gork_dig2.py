import fitz, sys
sys.stdout.reconfigure(encoding="utf-8")
p_ = r"E:\AAAwork\LLP-AGI\references\detector\Gorkavenko 等 - 2024 - LHCb potential to discover long-lived new physics particles with lifetimes above 100 ps.pdf"
doc = fitz.open(p_); p = doc[9]
t = p.get_text("dict")
print("--- spans y in 240..360, x<300")
for b in t["blocks"]:
    if b["type"]!=0: continue
    for l in b["lines"]:
        for s in l["spans"]:
            x0,y0,x1,y1 = s["bbox"]
            if 40<x0<300 and 240<y0<360:
                print(f"'{s['text']}' ({x0:.1f},{y0:.1f},{x1:.1f},{y1:.1f})")
print("--- drawings in that region")
for k,g in enumerate(p.get_drawings()):
    r = g["rect"]
    if r.x0>40 and r.x1<310 and r.y0>235 and r.y1<365:
        print(k, g.get("type"), "col", g.get("color"), "fill", g.get("fill"), "n", len(g["items"]), [round(v,1) for v in (r.x0,r.y0,r.x1,r.y1)])
