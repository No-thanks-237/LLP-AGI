import fitz, os
p_ = r"E:\AAAwork\LLP-AGI\references\detector\Gorkavenko 等 - 2024 - LHCb potential to discover long-lived new physics particles with lifetimes above 100 ps.pdf"
out = r"E:\AAAwork\LLP-AGI\references\detector\arXiv-2305.13383v3\_figs"
doc = fitz.open(p_); p = doc[9]
clips = {"fig7_left": fitz.Rect(40, 330, 300, 490), "fig7_right": fitz.Rect(315, 430, 560, 580)}
for name, c in clips.items():
    pix = p.get_pixmap(dpi=420, clip=c)
    f = os.path.join(out, name+".png"); pix.save(f); print(f, pix.width, pix.height)
