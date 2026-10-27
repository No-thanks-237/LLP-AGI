import sys
sys.stdout.reconfigure(encoding="utf-8")
t = open(r"E:\AAAwork\LLP-AGI\references\detector\arXiv-2305.13383v3\_figs\gork.txt",encoding="utf-8").read()
i = t.find("===== PAGE 13")
print(t[i:i+9000])
