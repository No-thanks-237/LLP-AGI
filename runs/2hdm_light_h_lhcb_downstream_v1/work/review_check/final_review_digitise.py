# -*- coding: utf-8 -*-
"""Optional item 5 (pixel version): digitise the paper FASER/FASER2 reach from LightH.pdf and
overlay with our contour.  Calibration from the PDF text spans, expressed in PDF points."""
import fitz, json, csv
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = r"E:/AAAwork/LLP-AGI/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream"
doc = fitz.open(R + "/inputs/references/LightH.pdf"); pg = doc[0]
ZC = 4
pix = pg.get_pixmap(matrix=fitz.Matrix(ZC, ZC), alpha=False)
img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)[:, :, :3].astype(int)
print("render:", img.shape)
# calibration in PDF points -> pixels
X0, XD = 313.35, 244.85
PY = np.array([87.6, 171.1, 254.4, 337.4]); LG = np.array([4.0, 3.0, 2.0, 1.0])
def x_of(px_pt): return 10.0**((np.asarray(px_pt, float) - X0)/XD)
def y_of(py_pt): return 10.0**np.interp(np.asarray(py_pt, float), PY, LG)
R_, G_, B_ = img[:,:,0], img[:,:,1], img[:,:,2]
red  = (R_ > 150) & (G_ < 110) & (B_ < 110)
blue = (B_ > 150) & (R_ < 120) & (G_ < 120)
# mask out the legend swatches (PDF x 470-576, y 45-100)
Yg, Xg = np.mgrid[0:img.shape[0], 0:img.shape[1]]
leg = ((Xg/ZC) > 468) & ((Yg/ZC) > 44) & ((Yg/ZC) < 102)
red &= ~leg; blue &= ~leg
# also mask the plot frame border lines (thin black, not coloured) - not needed
print("red px:", red.sum(), " blue px:", blue.sum())
def trace(mask, ncols=None):
    xs, ys = [], []
    for cx in range(mask.shape[1]):
        rows = np.flatnonzero(mask[:, cx])
        if len(rows) == 0: continue
        # cluster consecutive rows, take the largest cluster (avoid stray marks)
        groups = np.split(rows, np.flatnonzero(np.diff(rows) > 3)+1)
        g = max(groups, key=len)
        xs.append(cx/ZC); ys.append(np.median(g)/ZC)
    return np.array(xs), np.array(ys)
fig, ax = plt.subplots(figsize=(6.8, 5.2))
for mask, name, colr in ((red, 'FASER 2', 'crimson'), (blue, 'FASER', 'royalblue')):
    s = np.zeros_like(mask); s[::2, ::2] = mask[::2, ::2]     # subsample to keep the plot light
    yy, xx = np.nonzero(s)
    mH = x_of(xx/ZC); tb = y_of(yy/ZC)
    keep = (mH > 0.02) & (mH < 10.5) & (tb > 1) & (tb < 2e5)
    print(f"{name}: {keep.sum()} digitised pixels; mH {mH[keep].min():.3f}-{mH[keep].max():.3f},"
          f" tanb {tb[keep].min():.1f}-{tb[keep].max():.1f}")
    ax.scatter(mH[keep], tb[keep], s=0.4, color=colr, marker='.', lw=0,
               label=f'paper {name} reach (digitised)')
    for xq in (0.1, 0.3, 1.0, 3.0):
        sel = keep & (np.abs(np.log(mH/mH_target) if False else np.abs(mH/xq - 1)) < 0.02) if False else None
    for xq in (0.1, 0.3, 1.0, 3.0):
        sel = keep & (np.abs(mH - xq) < 0.02*xq)
        if sel.sum():
            print(f"   at mH~{xq}: tanb = {np.median(tb[sel]):.0f} ({sel.sum()} px)")
rows = list(csv.DictReader(open(R + "/outputs/data/exclusion_boundary.csv")))
for L, col in ((25.0, 'k'), (300.0, '0.45')):
    ms, ts = [], []
    for r in rows:
        if float(r['lumi_fb']) != L or r['tanb_long'] == '': continue
        ms.append(float(r['mass_GeV'])); ts.append(float(r['tanb_long']))
    o = np.argsort(ms)
    ax.plot(np.array(ms)[o], np.array(ts)[o], color=col, lw=2.2,
            label=f'this work: LHCb downstream {L:.0f} fb$^{{-1}}$')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(0.02, 10); ax.set_ylim(5, 1e5)
ax.set_xlabel(r'$m_H$ [GeV]'); ax.set_ylabel(r'$\tan\beta$')
ax.legend(fontsize=8, loc='center left')
ax.set_title('Type-I 2HDM light $H$, $\\cos(\\beta-\\alpha)=1/\\tan\\beta$, $m_{A}=m_{H^\\pm}=600$ GeV\n'
             'paper FASER/FASER2 reach (digitised from LightH.pdf) vs this work', fontsize=9)
fig.tight_layout()
out = R + "/work/review_check/paper_reach_overlay.png"
fig.savefig(out, dpi=150); print("wrote", out)
