
import re, sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
full = open(r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\review_check\lead\Gorkavenko_text.txt", encoding="utf-8").read()
i = full.find("S/B")
print(full[max(0,i-2000):i+1500])
