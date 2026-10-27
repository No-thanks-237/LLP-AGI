
import fitz, re
p = r"E:\AAAwork\LLP-AGI\references\model\BC5.pdf"
doc = fitz.open(p)
print("pages:", len(doc))
txt = []
for i, pg in enumerate(doc):
    txt.append("=== page %d ===\n" % (i+1) + pg.get_text())
full = "\n".join(txt)
open(r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\review_check\lead\BC5_text.txt", "w", encoding="utf-8").write(full)
print("chars:", len(full))
# find the efficiency discussion
for m_ in re.finditer(r"signal|efficien|0\.87|0\.4|epsilon", full, re.I):
    pass
i = full.lower().find("0.87")
while i != -1 and i < len(full):
    print("--- context around 0.87 (char %d) ---" % i)
    print(full[max(0,i-1200):i+600].replace("\n"," "))
    break
