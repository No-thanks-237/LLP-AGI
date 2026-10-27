import re, sys
sys.stdout.reconfigure(encoding="utf-8")
t = open(r"E:\AAAwork\LLP-AGI\references\detector\arXiv-2305.13383v3\_figs\gork.txt",encoding="utf-8").read()
for m in re.finditer("Higgs", t):
    s = max(0, m.start()-600); e = min(len(t), m.start()+900)
    print("\n\n############ at", m.start())
    print(t[s:e])
