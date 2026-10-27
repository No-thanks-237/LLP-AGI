
import fitz, re
p = r"E:\AAAwork\LLP-AGI\references\detector\Gorkavenko 等 - 2024 - LHCb potential to discover long-lived new physics particles with lifetimes above 100 ps.pdf"
doc = fitz.open(p)
print("pages:", len(doc))
full = "\n".join("=== page %d ===\n" % (i+1) + pg.get_text() for i, pg in enumerate(doc))
open(r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\review_check\lead\Gorkavenko_text.txt", "w", encoding="utf-8").write(full)
print("chars:", len(full))
for pat in ("0.87", "efficien", "N_ev", "2.3"):
    idxs = [m.start() for m in re.finditer(re.escape(pat), full)]
    print(pat, "->", len(idxs), "hits")
i = full.find("0.87")
if i > 0:
    print("\n=== context of 0.87 ===")
    print(full[max(0,i-2500):i+1200])
