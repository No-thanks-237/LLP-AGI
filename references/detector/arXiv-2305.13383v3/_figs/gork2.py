import re, sys
sys.stdout.reconfigure(encoding="utf-8")
t = open(r"E:\AAAwork\LLP-AGI\references\detector\arXiv-2305.13383v3\_figs\gork.txt",encoding="utf-8").read()
for kw in ["Table 2", "Table 3", "SensCalc"]:
    for m in re.finditer(kw, t):
        s = max(0, m.start()-700); e = min(len(t), m.start()+1200)
        print("\n\n############", kw, "at", m.start())
        print(t[s:e])
