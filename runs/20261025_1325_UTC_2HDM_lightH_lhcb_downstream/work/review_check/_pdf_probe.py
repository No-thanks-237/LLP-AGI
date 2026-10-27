import fitz
p = r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream/inputs/references/LightH.pdf"
doc = fitz.open(p)
print("pages:", len(doc))
for i, pg in enumerate(doc):
    print("page", i, "size", pg.rect)
    words = pg.get_text("words")
    print("  n words:", len(words))
    for w in words[:80]:
        print("   ", [round(x,1) for x in w[:4]], repr(w[4]))
