"""digitise_paper.py -- colour-based digitisation of the two light-H decay figures
of Li et al. arXiv:2212.06186v1 (inputs/references/*.png), used by the decay-module
validation (task-1, DESIGN.md section 2).

The two reference figures are
  A) mphi_t1Br_tanb102_cba=1tanb.png : Type-I, cos(beta-alpha)=1/tanb, tanb=10
     left panel  = BR(H->XX) vs mH  (log-log, 0.1..10 GeV)
     right panel = Gamma(H->XX) [GeV] vs mH (log-log)
  B) mphi_Type1_cba=1tanb-2.png : total width Gamma (left axis) and c*tau (right
     axis) vs mH for tanb = 5, 10, 70, 200, 1000.

Method (no physics input):
  * panel boxes and log-axis calibrations are measured from the PNG itself
    (major-tick pixel positions are detected by tick length, tick values are read
    off the axis labels and hard-coded in FIG);
  * each channel is identified by its matplotlib colour + solid/dashed style, read
    off the figure legend (top-to-bottom order)
    pi pi: blue solid, K K: red solid, 4pi: green solid, mu mu: orange solid,
    e e: cyan solid, gamma gamma: purple solid, c c: blue dashed, s s: red dashed,
    g g: green dashed, tau tau: orange dashed;
  * a colour mask is split into its solid and dashed curve by connected-component
    analysis after a horizontal dilation (a continuous curve covers every image
    column of its support, a dashed curve covers only ~2/3);
  * the legend box (found from its unique light-grey frame, edgecolour 214,214,214)
    is blanked, so legend swatches are never mistaken for data.
"""
import numpy as np, os, json
from PIL import Image
from scipy import ndimage

# ------------------------------------------------------------------ geometry
FIG = {
 "A": dict(file="mphi_t1Br_tanb102_cba=1tanb.png",
    panels=dict(
      br =dict(x0=350, x1=2292, y0=163, y1=1833,
               xmaj=[(351,-1.0),(1321.5,0.0),(2292,1.0)],
               ymaj=[(357,0.0),(998,-1.0),(1640,-2.0)]),
      gam=dict(x0=2729, x1=4672, y0=163, y1=1833,
               xmaj=[(2730,-1.0),(3701,0.0),(4672,1.0)],
               ymaj=[(372,-12.0),(615,-14.0),(859,-16.0),(1102,-18.0),
                     (1346,-20.0),(1589,-22.0),(1833,-24.0)]))),
 "B": dict(file="mphi_Type1_cba=1tanb-2.png",
    panels=dict(
      gam=dict(x0=395, x1=2045, y0=170, y1=1506,
               xmaj=[(396,-1.0),(1221,0.0),(2045,1.0)],
               ymaj=[(228,-9.0),(398,-11.0),(569,-13.0),(739,-15.0),
                     (909,-17.0),(1079,-19.0),(1249,-21.0),(1420,-23.0)]))),
}

PALETTE = {"blue":(0,0,255),"red":(255,0,0),"green":(0,128,0),
           "orange":(255,165,0),"cyan":(0,255,255),"purple":(128,0,128)}

CHANNELS_A = {"pipi":("blue","solid"), "KK":("red","solid"), "4pi":("green","solid"),
              "mumu":("orange","solid"), "ee":("cyan","solid"), "gaga":("purple","solid"),
              "cc":("blue","dashed"), "ss":("red","dashed"), "gg":("green","dashed"),
              "tautau":("orange","dashed")}
CHANNELS_B = {"tanb5":"blue","tanb10":"green","tanb70":"orange",
              "tanb200":"purple","tanb1000":"red"}


def _fit(tickpts):
    px=np.array([p for p,v in tickpts],float); val=np.array([v for p,v in tickpts],float)
    A=np.vstack([px,np.ones_like(px)]).T
    (m,b),*_=np.linalg.lstsq(A,val,rcond=None)
    return float(m),float(b),float(np.max(np.abs(m*px+b-val)))


def _cmask(a,rgb,tol):
    d=a.astype(np.int32)-np.array(rgb,np.int32)
    return (d*d).sum(axis=2) < tol*tol


def _legend_bbox(a,x0,x1,y0,y1):
    """bbox of the legend box from its frame colour 214,214,214 (unique in the png)."""
    d=np.abs(a.astype(np.int32)-np.array([214,214,214],np.int32)).max(axis=2)
    m=(d<=3)
    m[:y0+2,:]=False; m[y1-2:,:]=False; m[:,:x0+2]=False; m[:,x1-2:]=False
    lab,n=ndimage.label(ndimage.binary_dilation(m,structure=np.ones((3,3),bool)))
    # the legend frame can be split into several connected components where curves
    # touch its edge, so collect every interior frame-like component and take the
    # union of their bounding boxes (a data set has at most one legend per axes).
    boxes=[]
    for i in range(1,n+1):
        ys,xs=np.where(lab==i)
        if len(ys)<200: continue
        xa,xb,ya,yb=int(xs.min()),int(xs.max()),int(ys.min()),int(ys.max())
        if xb-xa<80 or yb-ya<80: continue
        if xa<=x0+8 or xb>=x1-8 or ya<=y0+8 or yb>=y1-8: continue   # spines touch the edge
        boxes.append((xa,xb,ya,yb))
    if not boxes: return None
    xa=min(b[0] for b in boxes); xb=max(b[1] for b in boxes)
    ya=min(b[2] for b in boxes); yb=max(b[3] for b in boxes)
    m=8
    return xa-m,xb+m,ya-m,yb+m


def _curves(sub):
    """split a colour sub-mask into (solid_mask, dashed_mask, components).
    A curve is found as a connected component after a horizontal (1x15) dilation
    that bridges the dash gaps; a continuous curve then covers every column of its
    support (~1.0) while a dashed curve covers only the 'on' part (~0.6-0.85)."""
    dil=ndimage.binary_dilation(sub,structure=np.ones((1,15),bool))
    lab,n=ndimage.label(dil)
    solid=np.zeros_like(sub,bool); dashed=np.zeros_like(sub,bool)
    comps=[]
    for i in range(1,n+1):
        comp=(lab==i)
        if comp.sum()<200: continue
        ys,xs=np.where(comp)
        cx0,cx1,cy0,cy1=int(xs.min()),int(xs.max()),int(ys.min()),int(ys.max())
        if cx1-cx0<20: continue
        orig=comp&sub
        cols=np.unique(np.where(orig)[1])
        cover=len(cols)/(cx1-cx0+1)
        comps.append(dict(cx0=cx0,cx1=cx1,npx=int(orig.sum()),cover=float(cover)))
        if 0.35<=cover<0.90:        dashed |= orig
        elif cover>=0.93:           solid   |= orig
    return solid,dashed,comps


def _runs(col):
    idx=np.where(col)[0]
    if len(idx)==0: return []
    o=[];s=idx[0];p=idx[0]
    for y in idx[1:]:
        if y-p>2: o.append((s,p)); s=y
        p=y
    o.append((s,p))
    return o


def _sample(mask, xs, xm, bx, x0, y1_, y0_, ym, by, maxrun=30):
    """return [(X,value)] sampled on a solid/dashed mask."""
    out=[]
    for X in xs:
        xc=int(round((np.log10(X)-bx)/xm))-x0
        if not (3<=xc<=(mask.shape[1]-3)): continue
        rs=_runs(mask[:,xc])
        if len(rs)!=1: continue
        s,e=rs[0]; L=e-s+1
        if L>maxrun: continue
        yc=(s+e)/2.0+y0_
        out.append((float(X),float(10**(ym*yc+by))))
    return out


def digitise(base, figkey, panelkey, xs):
    cfg=FIG[figkey]; p=cfg["panels"][panelkey]
    img=np.asarray(Image.open(os.path.join(base,"inputs","references",cfg["file"])).convert("RGB"))
    xm,bx,rx=_fit(p["xmaj"]); ym,by,ry=_fit(p["ymaj"])
    x0,x1,y0,y1=p["x0"],p["x1"],p["y0"],p["y1"]
    leg=_legend_bbox(img,x0,x1,y0,y1)
    img=img.copy()
    if leg:
        xa,xb,ya,yb=leg
        img[ya:yb+1,xa:xb+1]=255
    info=dict(figure=cfg["file"],panel=panelkey,
              xfit_resid=rx,yfit_resid=ry,legend_bbox=leg,
              xmap=[xm,bx],ymap=[ym,by],panel_box=[x0,x1,y0,y1])
    colour={c:_cmask(img,rgb,60) for c,rgb in PALETTE.items()}
    sub={c:colour[c][y0:y1,x0:x1] for c in colour}
    split={}; info["components"]={}
    for c in colour:
        solid,dashed,comps=_curves(sub[c])
        split[c]=(solid,dashed)
        info["components"][c]=dict(solid_bbox=None,dashed_bbox=None,ncomps=len(comps))
        for kind,m in (("solid",solid),("dashed",dashed)):
            if m.any():
                ys,xs_=np.where(m)
                info["components"][c][kind+"_bbox"]=[int(xs_.min()),int(xs_.max()),int(ys.min()),int(ys.max())]
                info["components"][c][kind+"_npx"]=int(m.sum())
        info["components"][c]["components"]=comps
    # ---- sampling -----------------------------------------------------
    # The same colour is used for a solid+dashed pair, but in these figures the two
    # channels of a pair have disjoint (or decade-separated) mass ranges by
    # construction of the model matching (hadronic channels for mH<2, partonic
    # above; mu mu vs tau tau separated by (m_tau/m_mu)^2).  We therefore assign a
    # colour run to its channel by mass range, which needs no style guess:
    def _sample_ch(col, ch, xlo, xhi, pick=None):
        """pick=None : accept only a column with exactly one run (unambiguous);
           pick='low'/'high' : of 2+ runs take the lower/higher curve (used only
           where the two same-colour channels coexist: mu mu vs tau tau)."""
        m=sub[col]; out=[]
        for X in xs:
            if not (xlo<X<xhi): continue
            xc=int(round((np.log10(X)-bx)/xm))-x0
            if not (3<=xc<=(m.shape[1]-3)): continue
            rs=[r for r in _runs(m[:,xc]) if (r[1]-r[0]+1)<=30]   # skip vertical jumps
            if len(rs)==0: continue
            if pick is None:
                if len(rs)!=1: continue
                s,e=rs[0]
            else:
                if len(rs)<2: continue
                # run[0] is a pixel row: larger value <=> smaller row.
                s,e=max(rs,key=lambda r:r[0]) if pick=="low" else min(rs,key=lambda r:r[0])
            out.append((float(X),float(10**(ym*((s+e)/2.0+y0)+by))))
        return out
    res={}
    if figkey=="A":
        RULES={"pipi":("blue",0.0,1.98,None),  "cc":("blue",3.55,11.0,None),
               "KK":("red",0.0,2.02,None),     "ss":("red",2.6,11.0,None),
               "4pi":("green",0.0,2.0,None),   "gg":("green",2.0,11.0,None),
               "mumu":("orange",0.0,3.55,None),"tautau":("orange",3.60,11.0,"high"),
               "ee":("cyan",0.0,11.0,None),    "gaga":("purple",0.0,11.0,None)}
        for ch,(col,xlo,xhi,pick) in RULES.items():
            res[ch]=_sample_ch(col,ch,xlo,xhi,pick)
    else:
        for ch,col in CHANNELS_B.items():
            res[ch]=_sample_ch(col,ch,0.0,11.0,None)
    return res,info
