import fitz, sys, math, statistics
sys.stdout.reconfigure(encoding="utf-8")
p_ = r"E:\AAAwork\LLP-AGI\references\detector\Gorkavenko 等 - 2024 - LHCb potential to discover long-lived new physics particles with lifetimes above 100 ps.pdf"
doc = fitz.open(p_); p = doc[9]
d = p.get_drawings()
print("--- frame lines in right panel")
for k,g in enumerate(d):
    r=g["rect"]
    if 315<r.x0<480 and 430<r.y0<580 and (r.width>100 or r.height>100) and len(g["items"])<=2:
        print("  ",k,[round(v,1) for v in (r.x0,r.y0,r.x1,r.y1)], len(g["items"]))
print("--- y ticks near x=326")
for k,g in enumerate(d):
    r=g["rect"]
    if 320<r.x0<330 and 430<r.y0<580 and r.width<6 and r.height<6:
        print("  ",k,[round(v,1) for v in (r.x0,r.y0,r.x1,r.y1)])
print("--- x ticks along bottom")
for k,g in enumerate(d):
    r=g["rect"]
    if 550<r.y0<575 and 320<r.x0<480 and r.width<6 and r.height<6:
        print("  ",k,[round(v,1) for v in (r.x0,r.y0,r.x1,r.y1)])
print("--- x labels")
for b in p.get_text("dict")["blocks"]:
    if b["type"]!=0: continue
    for l in b["lines"]:
        for s in l["spans"]:
            x0,y0,x1,y1=s["bbox"]
            if 320<x0<480 and 565<y0<585:
                print(f"   '{s['text']}' ({x0:.1f},{y0:.1f},{x1:.1f},{y1:.1f})")
