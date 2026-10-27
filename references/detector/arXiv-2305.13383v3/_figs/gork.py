import fitz, os, re, io, sys
sys.stdout.reconfigure(encoding="utf-8")
p_ = r"E:\AAAwork\LLP-AGI\references\detector\Gorkavenko 等 - 2024 - LHCb potential to discover long-lived new physics particles with lifetimes above 100 ps.pdf"
doc = fitz.open(p_)
print("pages", len(doc))
txt = []
for i,pg in enumerate(doc):
    txt.append(f"\n===== PAGE {i+1} =====\n" + pg.get_text())
full = "".join(txt)
open(r"E:\AAAwork\LLP-AGI\references\detector\arXiv-2305.13383v3\_figs\gork.txt","w",encoding="utf-8").write(full)
print(len(full))
for kw in ["Table 2","SciFi","SensCalc","Higgs","h -> SS","h \\to SS","eta","acceptance","0.4","downstream","Downstream"]:
    idxs = [m.start() for m in re.finditer(kw, full)]
    print(kw, len(idxs))
