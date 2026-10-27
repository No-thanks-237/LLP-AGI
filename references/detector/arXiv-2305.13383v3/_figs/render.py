import fitz, os, sys
d = r"E:\AAAwork\LLP-AGI\references\detector\arXiv-2305.13383v3"
out = os.path.join(d, "_figs")
os.makedirs(out, exist_ok=True)
names = ["BC6-validation.pdf","BC1-validation.pdf","Double-Distr-scalars-quartic.pdf",
         "PlotDiffNevents.pdf","acceptances.pdf","NeventsDensityPlotExample.pdf",
         "SensitivityScalar.pdf","FIP-production.pdf","acceptance-illustration.pdf",
         "Comparison-ALPINIST.pdf"]
for n in names:
    doc = fitz.open(os.path.join(d, n))
    p = doc[0]
    pix = p.get_pixmap(dpi=200)
    base = os.path.join(out, n.replace(".pdf", ".png"))
    pix.save(base)
    print(n, "pages", len(doc), "rect", p.rect, "->", pix.width, pix.height)
