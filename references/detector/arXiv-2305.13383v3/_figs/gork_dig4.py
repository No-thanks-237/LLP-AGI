import fitz, sys
sys.stdout.reconfigure(encoding="utf-8")
p_ = r"E:\AAAwork\LLP-AGI\references\detector\Gorkavenko 等 - 2024 - LHCb potential to discover long-lived new physics particles with lifetimes above 100 ps.pdf"
doc = fitz.open(p_); p = doc[9]
print("--- texts x1<92, y in 330..500")
for b in p.get_text("dict")["blocks"]:
    if b["type"]!=0: continue
    for l in b["lines"]:
        for s in l["spans"]:
            x0,y0,x1,y1 = s["bbox"]
            if x1<95 and 330<y0<500:
                print(f"  '{s['text']}' ({x0:.1f},{y0:.1f},{x1:.1f},{y1:.1f}) yc={(y0+y1)/2:.1f}")
print("--- short segments x0 in 84..96, y 330..500")
for k,g in enumerate(p.get_drawings()):
    r = g["rect"]
    if 84<=r.x0<=96 and 330<r.y0<500 and (r.x1-r.x0)<6 and (r.y1-r.y0)<6:
        print("  ", k, [round(v,1) for v in (r.x0,r.y0,r.x1,r.y1)])
