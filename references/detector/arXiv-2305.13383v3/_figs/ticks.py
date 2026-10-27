import fitz, os
d = r"E:\AAAwork\LLP-AGI\references\detector\arXiv-2305.13383v3"
doc = fitz.open(os.path.join(d,"acceptances.pdf")); p = doc[0]
dr = p.get_drawings()
for k,g in enumerate(dr):
    if len(g["items"])<=6:
        print(k, "color", g.get("color"), "n", len(g["items"]), [ (round(it[1].x,1), round(it[1].y,1), round(it[2].x,1), round(it[2].y,1)) if it[0]=='l' else it[0] for it in g["items"] ][:6])
