import fitz, os
p_ = r"E:\AAAwork\LLP-AGI\references\detector\Gorkavenko 等 - 2024 - LHCb potential to discover long-lived new physics particles with lifetimes above 100 ps.pdf"
out = r"E:\AAAwork\LLP-AGI\references\detector\arXiv-2305.13383v3\_figs"
doc = fitz.open(p_)
for i in [8,9,10,12]:
    pix = doc[i].get_pixmap(dpi=170)
    f = os.path.join(out, f"gork_p{i+1}.png")
    pix.save(f); print(f, pix.width, pix.height)
