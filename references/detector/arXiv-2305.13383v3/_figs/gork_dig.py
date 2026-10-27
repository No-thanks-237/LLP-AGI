import fitz, sys, math, statistics
sys.stdout.reconfigure(encoding="utf-8")
p_ = r"E:\AAAwork\LLP-AGI\references\detector\Gorkavenko 等 - 2024 - LHCb potential to discover long-lived new physics particles with lifetimes above 100 ps.pdf"
doc = fitz.open(p_); p = doc[9]
print("page rect", p.rect)
t = p.get_text("dict")
for b in t["blocks"]:
    if b["type"]!=0: continue
    for l in b["lines"]:
        for s in l["spans"]:
            x0,y0,x1,y1 = s["bbox"]
            if 40<x0<640 and 600<y0<900:
                print(f"'{s['text']}' bbox=({x0:.1f},{y0:.1f},{x1:.1f},{y1:.1f})")
