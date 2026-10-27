import fitz, os
d = r"E:\AAAwork\LLP-AGI\references\detector\arXiv-2305.13383v3"
for n in ["BC6-validation.pdf","BC1-validation.pdf"]:
    doc = fitz.open(os.path.join(d,n)); p = doc[0]
    print("=====", n)
    t = p.get_text("dict")
    for b in t["blocks"]:
        if b["type"]!=0: continue
        for l in b["lines"]:
            for s in l["spans"]:
                x0,y0,x1,y1 = s["bbox"]
                if x1 < 70:
                    print(f"   '{s['text']}' x=[{x0:.1f},{x1:.1f}] yc={(y0+y1)/2:.1f}")
