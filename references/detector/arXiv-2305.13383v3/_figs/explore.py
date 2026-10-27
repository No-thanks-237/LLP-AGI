import fitz, os
d = r"E:\AAAwork\LLP-AGI\references\detector\arXiv-2305.13383v3"
for n in ["acceptances.pdf","PlotDiffNevents.pdf","BC6-validation.pdf","BC1-validation.pdf"]:
    doc = fitz.open(os.path.join(d,n)); p = doc[0]
    dr = p.get_drawings()
    print("=====", n, "drawings:", len(dr))
    for k,g in enumerate(dr[:8]):
        types = {}
        for it in g["items"]:
            types[it[0]] = types.get(it[0],0)+1
        print("  ", k, "type", g.get("type"), "color", g.get("color"), "fill", g.get("fill"), "w", g.get("width"), "items", types, "rect", g["rect"])
    t = p.get_text("dict")
    print("  text spans:", sum(len(l["spans"]) for b in t["blocks"] if b["type"]==0 for l in b["lines"]))
    print("  first spans:", [(s["text"], [round(v,1) for v in s["bbox"]]) for b in t["blocks"] if b["type"]==0 for l in b["lines"] for s in l["spans"]][:14])
