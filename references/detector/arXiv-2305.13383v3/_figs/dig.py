import fitz, os, math, json
d = r"E:\AAAwork\LLP-AGI\references\detector\arXiv-2305.13383v3"

def frame_and_ticks(p):
    dr = p.get_drawings()
    segs = []
    for g in dr:
        if len(g["items"]) == 1 and g["items"][0][0] == "l":
            a, b = g["items"][0][1], g["items"][0][2]
            segs.append((a.x, a.y, b.x, b.y))
    return segs

def report(name):
    doc = fitz.open(os.path.join(d, name)); p = doc[0]
    segs = frame_and_ticks(p)
    print("=====", name, "nseg", len(segs))
    # frame: long segments
    longs = [s for s in segs if abs(s[2]-s[0]) > 100 or abs(s[3]-s[1]) > 100]
    print(" frame:", [(round(a,1),round(b,1),round(c,1),round(e,1)) for a,b,c,e in longs])
    xs = [s for s in segs if abs(s[3]-s[1]) < 1e-6 and abs(s[2]-s[0]) < 10 and abs(s[2]-s[0]) > 3]
    ys = [s for s in segs if abs(s[2]-s[0]) < 1e-6 and abs(s[3]-s[1]) < 10 and abs(s[3]-s[1]) > 3]
    ybottom = max(s[1] for s in longs) if longs else None
    print(" x ticks (vertical-tick rows grouped):")
    from collections import defaultdict
    rows = defaultdict(list)
    for s in xs: rows[round(s[1],1)].append(round(s[0],1))
    for k in sorted(rows): print("   y=",k, sorted(rows[k]))
    colrows = defaultdict(list)
    for s in ys: colrows[round(s[0],1)].append(round(s[1],1))
    print(" y ticks:")
    for k in sorted(colrows)[:4]: print("   x=",k, sorted(colrows[k]))

for n in ["BC6-validation.pdf","BC1-validation.pdf","PlotDiffNevents.pdf"]:
    report(n)
