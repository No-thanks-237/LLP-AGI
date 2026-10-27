import numpy as np, math
from PIL import Image
f = r"E:\AAAwork\LLP-AGI\references\detector\arXiv-2305.13383v3\_figs\fig7_left.png"
im = np.asarray(Image.open(f).convert("RGB")).astype(int)
H,W,_ = im.shape
R,G,B = im[:,:,0],im[:,:,1],im[:,:,2]
blue = (B>170)&(R<120)&(G<120)
print("shape",im.shape,"blue px",blue.sum())
S = 5.835   # px per pt
def x_of_m(m): return (144.6 + 77.1*math.log10(m) - 40.0)*S
def val_of_y(py): return 10**(-7 + (214.7-py)/55.57)
for m in [0.4,0.5,0.7,1.0,1.5,2.0,3.0,5.0,10.0]:
    c = int(round(x_of_m(m)))
    if c<0 or c>=W: continue
    col = np.where(blue[:,max(0,c-1):c+2].any(axis=1))[0]
    if len(col)==0:
        print(f"m={m}: none"); continue
    # cluster
    groups=[]; cur=[col[0]]
    for y in col[1:]:
        if y-cur[-1]<=6: cur.append(y)
        else: groups.append(cur); cur=[y]
    groups.append(cur)
    vals=[val_of_y(float(np.mean(g))) for g in groups if len(g)>=1]
    print(f"m={m:5.2f} x={c:4d}: " + "  ".join(f"{v:.2e}" for v in vals))
