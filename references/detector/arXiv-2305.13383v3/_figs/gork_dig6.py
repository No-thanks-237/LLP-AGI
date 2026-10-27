import fitz, sys
sys.stdout.reconfigure(encoding="utf-8")
p_ = r"E:\AAAwork\LLP-AGI\references\detector\Gorkavenko 等 - 2024 - LHCb potential to discover long-lived new physics particles with lifetimes above 100 ps.pdf"
doc = fitz.open(p_); p = doc[9]
for k,g in enumerate(p.get_drawings()):
    r = g["rect"]
    if r.y0>330 and len(g["items"])>15 and r.width>10:
        print(k, g.get("type"), "col", g.get("color"), "fill", g.get("fill"), "n", len(g["items"]), [round(v,1) for v in (r.x0,r.y0,r.x1,r.y1)])
