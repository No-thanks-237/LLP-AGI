
import re
full = open(r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\review_check\lead\Gorkavenko_text.txt", encoding="utf-8").read()
for pat in ("reconstruct", "acceptance", "background", "S/B", "0.97", "0.9", "N_{", "signal"):
    print(pat, len(re.findall(re.escape(pat), full, re.I)))
i = full.find("reconstruct")
print(full[max(0,i-800):i+1500])
