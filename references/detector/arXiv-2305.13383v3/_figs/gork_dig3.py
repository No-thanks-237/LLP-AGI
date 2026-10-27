import fitz, sys, math, statistics
sys.stdout.reconfigure(encoding="utf-8")
p_ = r"E:\AAAwork\LLP-AGI\references\detector\Gorkavenko 等 - 2024 - LHCb potential to discover long-lived new physics particles with lifetimes above 100 ps.pdf"
doc = fitz.open(p_); p = doc[9]
print("--- y ticks (vertical segments at x~90.7)")
for k,g in enumerate(p.get_drawings()):
    r = g["rect"]
    if 85 < r.x0 < 92 and r.x1 < 95 and 180 < r.y0 < 350:
        print("  ", k, g.get("type"), [round(v,1) for v in (r.x0,r.y0,r.x1,r.y1)])
print("--- curves")
for k,g in enumerate(p.get_drawings()):
    r = g["rect"]
    col = g.get("color")
    if col is None: continue
    col = tuple(round(c,2) for c in col)
    if not (88 < r.x0 < 300 and 180 < r.y0 < 350 and r.x1 < 300): continue
    if len(g["items"]) < 20: continue
    print("  ", k, col, g.get("dashes"), "n", len(g["items"]), [round(v,1) for v in (r.x0,r.y0,r.x1,r.y1)])
