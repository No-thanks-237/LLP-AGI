#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
digitize.py -- digitize three published Mathematica vector-PDF physics figures into CSV
validation targets.

Run:  D:/py/python.exe digitize.py

Inputs (read-only):
  references/model/arXiv-1809.01876v1/decayrates.pdf
  references/model/arXiv-1904.10447v4/Scalar-decay-br.pdf
  references/model/arXiv-1904.10447v4/tau-scalar.pdf

Outputs (written next to this script):
  decayrates_curves.csv
  scalardecaybr_curves.csv
  tauscalar_curves.csv
  overlay_decayrates.png / overlay_scalardecaybr.png / overlay_tauscalar.png
  digitization_report.md

Method
------
* text  : page.get_text("dict") -> line-level spans with bboxes (tick + legend labels)
* paths : page.get_drawings()  -> axis/frame/ticks (short or very long straight lines)
                                  vs data polylines (>= 5 line segments)
* legend: short horizontal single-segment strokes that do not touch the frame are
          legend swatches; each is matched to the text line to its right.
* axes  : every axis here is log10.  Calibration coord = a*log10(value) + b is fitted
          by least squares to (label value, detected major-tick coordinate) pairs.
          Independent check: every *minor* tick coordinate is mapped back to a value
          and compared with log-nice numbers n*10^k (n = 1..9, 1.5, 2.5, 3.5, 4.5).
* curves: data paths are merged per (colour, dash, width) by x-adjacency into curves,
          sorted by x, decimated to <= 2000 points, converted with the calibration.
"""
import os, re, csv, math
import numpy as np
import fitz
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
# The source PDFs (Winkler arXiv:1809.01876 decayrates.pdf, Boiarska arXiv:1904.10447) are NOT
# shipped inside the run directory.  The DIGITISED tables they produce ARE (decayrates_curves.csv,
# scalardecaybr_curves.csv, hadronic_table.csv), which is what the analysis uses.  Point BC5_REFS
# at a folder holding the papers to re-digitise.
REF = os.environ.get('BC5_REFS', os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             '..', '..', '..', 'inputs', 'references'))

DPI = 200
MAXPTS = 2000
COLORTOL = 0.02

# ----------------------------------------------------------------------------- helpers
def fmt(v, n=6):
    return ("%." + str(n) + "g") % v

def colkey(c):
    if c is None:
        return None
    return tuple(round(float(x), 3) for x in c)

DASH_NUM = re.compile(r"[-+]?\d*\.?\d+")
def dashkey(d):
    """PyMuPDF may return the dash pattern as a list or as the raw PDF string
    '[] 0' / '[ .247884 4.957673 ] 0'; normalise both to a tuple of dash lengths."""
    if not d:
        return ()
    if isinstance(d, (list, tuple)):
        vals = [float(x) for x in d]
    else:
        s = str(d)
        m = re.search(r"\[(.*?)\]", s)
        inner = m.group(1) if m else ""
        vals = [float(x) for x in DASH_NUM.findall(inner)]
    if not vals:
        return ()
    return tuple(round(v, 3) for v in vals)

def path_points(d):
    """Ordered point list of a path; Bezier segments are sampled."""
    pts = []
    for it in d["items"]:
        if it[0] == "l":
            p0, p1 = it[1], it[2]
            if not pts:
                pts.append((p0.x, p0.y))
            pts.append((p1.x, p1.y))
        elif it[0] == "c":
            p0, p1, p2, p3 = it[1], it[2], it[3], it[4]
            if not pts:
                pts.append((p0.x, p0.y))
            for t in np.linspace(0.1, 1.0, 9):
                mt = 1.0 - t
                x = (mt**3)*p0.x + 3*(mt**2)*t*p1.x + 3*mt*(t**2)*p2.x + (t**3)*p3.x
                y = (mt**3)*p0.y + 3*(mt**2)*t*p1.y + 3*mt*(t**2)*p2.y + (t**3)*p3.y
                pts.append((x, y))
        elif it[0] == "re":
            r = it[1]
            pts += [(r.x0, r.y0), (r.x1, r.y0), (r.x1, r.y1), (r.x0, r.y1), (r.x0, r.y0)]
    return pts

def nsegs(d):
    return sum(1 for it in d["items"] if it[0] == "l")

def is_axis_aligned(d, tol=0.35):
    p = path_points(d)
    if len(p) < 2:
        return False
    xs = [q[0] for q in p]; ys = [q[1] for q in p]
    return (max(xs) - min(xs) < tol) or (max(ys) - min(ys) < tol)

def text_lines(page):
    out = []
    for b in page.get_text("dict")["blocks"]:
        if b.get("type") != 0:
            continue
        for l in b["lines"]:
            txt = "".join(s["text"] for s in l["spans"])
            x0, y0, x1, y1 = l["bbox"]
            out.append(dict(text=txt, bbox=(x0, y0, x1, y1), x0=x0, x1=x1, y0=y0, y1=y1,
                            cx=0.5*(x0+x1), cy=0.5*(y0+y1)))
    return out

NUMRE = re.compile(r"^10-(\d+)$")
def parse_label(txt):
    """'10-9' -> 1e-9 ; '0.5' -> 0.5 ; '10' -> 10 ; otherwise None."""
    t = txt.strip().replace("\u2212", "-")
    m = NUMRE.match(t)
    if m:
        return 10.0 ** (-int(m.group(1)))
    try:
        return float(t)
    except ValueError:
        return None

# ----------------------------------------------------------------------------- figure specs
FIGURE_SPECS = [
    dict(
        key="decayrates",
        pdf=os.path.join(REF, "arXiv-1809.01876v1", "decayrates.pdf"),
        title="Hadronic and leptonic decay rates of a light scalar mixing with the Higgs "
              "(s_theta^2 = 1); scalar decay RATES vs scalar mass.",
        units="GeV",
        slug_map={
            "\u03bc\u03bc": "muons",
            "\u03c4\u03c4": "taus",
            "\u03c0\u03c0": "pions",
            "KK": "kaons",
            "4\u03c0,\u03b7\u03b7,\u03c1\u03c1,...": "multimeson_4pi_etaeta_rhorho",
            "gg": "gluons_gg",
            "ss": "strange_quarks_ss",
            "cc": "charm_quarks_cc",
            "hadrons": "hadrons",
        },
        manual={},
        slug_notes={
            "muons": "legend swatch dark teal; mu+mu- channel",
            "taus": "legend swatch yellow-green; tau+tau- channel",
            "pions": "legend swatch blue; pi+pi- channel (dispersive / data driven, m < 2 GeV)",
            "kaons": "legend swatch orange; K Kbar channel (dispersive / data driven, m < 2 GeV)",
            "multimeson_4pi_etaeta_rhorho": "legend swatch grey-blue; 4pi, eta eta, rho rho, ... (m < 2 GeV)",
            "gluons_gg": "legend swatch indigo; perturbative g g channel (drawn only above 2 GeV)",
            "strange_quarks_ss": "legend swatch dark orange; perturbative s sbar channel (drawn only above 2 GeV)",
            "charm_quarks_cc": "legend swatch magenta; perturbative c cbar channel (drawn only above 2 GeV)",
            "hadrons": "legend swatch black dotted; total hadronic width, dispersive below 2 GeV "
                       "and perturbative above (continuous at 2 GeV - verified)",
        },
    ),
    dict(
        key="scalardecaybr",
        pdf=os.path.join(REF, "arXiv-1904.10447v4", "Scalar-decay-br.pdf"),
        title="Branching ratios of decays of a scalar S as a function of its mass "
              "(left panel of the figure in Boiarska et al.).",
        units="dimensionless",
        slug_map={
            "e+e-": "BR_ee",
            "\u03bc+\u03bc-": "BR_mumu",
            "\u03c4+\u03c4-": "BR_tautau",
            "\u03c0\u03c0": "BR_pions",
            "KK": "BR_kaons",
            "cc": "BR_cc",
            "GG": "BR_gg",
            "ss": "BR_ss",
            "4\u03c0": "BR_4pi",
        },
        manual={},
        slug_notes={
            "BR_ee": "legend swatch red",
            "BR_mumu": "legend swatch blue",
            "BR_tautau": "legend swatch cyan",
            "BR_pions": "legend swatch green",
            "BR_kaons": "legend swatch solid black",
            "BR_cc": "legend swatch magenta",
            "BR_gg": "legend swatch orange",
            "BR_ss": "legend swatch long-dashed black",
            "BR_4pi": "legend swatch dotted black",
        },
    ),
    dict(
        key="tauscalar",
        pdf=os.path.join(REF, "arXiv-1904.10447v4", "tau-scalar.pdf"),
        title="Lifetime of a scalar S as a function of its mass with mixing angle theta^2 = 1. "
              "Solid blue: decays into leptons, kaons and pions (Winkler) plus fictitious "
              "multi-meson channel.  Solid red: perturbative QCD decays into quarks and gluons.",
        units="s",
        slug_map={},
        manual={
            (0.0, 0.0, 1.0): "lifetime_blue_winkler_plus_multimeson",
            (1.0, 0.0, 0.0): "lifetime_red_perturbative_qcd",
        },
        slug_notes={
            "lifetime_blue_winkler_plus_multimeson": "solid blue, m = 0.21 .. 2.0 GeV "
                "(drawn as two PDF paths joined at m ~ 0.693 GeV; visually continuous)",
            "lifetime_red_perturbative_qcd": "solid red, m = 2.0 .. 5.3 GeV",
        },
    ),
]

# ----------------------------------------------------------------------------- geometry
def find_frame(drawings, page):
    """Frame = extreme long black axis-aligned lines."""
    W, H = page.rect.width, page.rect.height
    vx, hy = [], []
    for d in drawings:
        if d["type"] != "s":
            continue
        c = colkey(d.get("color"))
        if c is None or max(c) > 0.25:          # black-ish only (drop the grey grid)
            continue
        if nsegs(d) != 1 or not is_axis_aligned(d):
            continue
        p = path_points(d)
        (x0, y0), (x1, y1) = p[0], p[-1]
        if abs(x1 - x0) < 0.35:                  # vertical
            if abs(y1 - y0) > 0.45 * H:
                vx.append((0.5*(x0+x1), abs(y1-y0)))
        elif abs(y1 - y0) < 0.35:                # horizontal
            if abs(x1 - x0) > 0.45 * W:
                hy.append((0.5*(y0+y1), abs(x1-x0)))
    if not vx or not hy:
        return None
    xs = sorted(v[0] for v in vx); ys = sorted(h[0] for h in hy)
    return dict(x0=xs[0], x1=xs[-1], y0=ys[0], y1=ys[-1], vlines=vx, hlines=hy)

def near_frame(v, frame, tol=1.2):
    return (abs(v - frame["x0"]) < tol or abs(v - frame["x1"]) < tol or
            abs(v - frame["y0"]) < tol or abs(v - frame["y1"]) < tol)

def find_ticks(drawings, frame):
    """Short axis-aligned black strokes attached to the frame."""
    xt, yt = [], []
    for i, d in enumerate(drawings):
        if d["type"] != "s" or nsegs(d) != 1 or not is_axis_aligned(d):
            continue
        c = colkey(d.get("color"))
        if c is None or max(c) > 0.25:
            continue
        p = path_points(d)
        (x0, y0), (x1, y1) = p[0], p[-1]
        if abs(x1 - x0) < 0.35:                              # vertical stroke -> x tick
            L = abs(y1 - y0)
            if 0.8 < L < 12 and (near_frame(min(y0, y1), frame, 2.5) or near_frame(max(y0, y1), frame, 2.5)):
                xt.append(dict(idx=i, coord=0.5*(x0+x1), len=L))
        elif abs(y1 - y0) < 0.35:                            # horizontal stroke -> y tick
            L = abs(x1 - x0)
            if 0.5 < L < 12 and (near_frame(min(x0, x1), frame, 2.5) or near_frame(max(x0, x1), frame, 2.5)):
                yt.append(dict(idx=i, coord=0.5*(y0+y1), len=L))
    def uniq(lst):
        out = []
        for t in sorted(lst, key=lambda z: z["coord"]):
            if out and abs(t["coord"] - out[-1]["coord"]) < 0.3:
                if t["len"] > out[-1]["len"]:
                    out[-1] = t
            else:
                out.append(t)
        return out
    return uniq(xt), uniq(yt)

# ----------------------------------------------------------------------------- calibration
def fit_axis(pairs, scale):
    """pairs = [(value, coord)] ; coord = a*f(value)+b."""
    vals = np.array([p[0] for p in pairs], float)
    coords = np.array([p[1] for p in pairs], float)
    f = np.log10(vals) if scale == "log" else vals
    A = np.vstack([f, np.ones_like(f)]).T
    (a, b), res, rank, sv = np.linalg.lstsq(A, coords, rcond=None)
    pred = a * f + b
    r = coords - pred
    return dict(scale=scale, a=float(a), b=float(b), n=len(pairs),
                res_max=float(np.max(np.abs(r))), res_rms=float(np.sqrt(np.mean(r**2))),
                res=r)

def value_from_coord(cal, coord):
    if cal["scale"] == "log":
        return 10.0 ** ((coord - cal["b"]) / cal["a"])
    return (coord - cal["b"]) / cal["a"]

def coord_from_value(cal, value):
    f = math.log10(value) if cal["scale"] == "log" else value
    return cal["a"] * f + cal["b"]

def calibrate_axis(ticks, labels, axis, frame, tol_px=4.0):
    warn = []
    if not ticks:
        return None, [], ["no ticks detected on the %s axis" % axis]
    maxlen = max(t["len"] for t in ticks)
    majors = [t for t in ticks if t["len"] >= 0.7 * maxlen]
    pairs = []
    for l in labels:
        v = parse_label(l["text"])
        if v is None or v <= 0:
            continue
        c = l["cx"] if axis == "x" else l["cy"]
        cand = min(majors, key=lambda t: abs(t["coord"] - c))
        if abs(cand["coord"] - c) > tol_px:
            warn.append("label %r at coord %.2f has no major tick within %.1f px" % (l["text"], c, tol_px))
            continue
        pairs.append((v, cand["coord"], l["text"]))
    if len(pairs) < 2:
        return None, pairs, warn + ["fewer than 2 usable tick labels on the %s axis" % axis]
    pv = [p[0] for p in pairs]; pc = [p[1] for p in pairs]
    cal_log = fit_axis(list(zip(pv, pc)), "log")
    cal_lin = fit_axis(list(zip(pv, pc)), "linear")
    cal = cal_log if cal_log["res_max"] <= cal_lin["res_max"] else cal_lin
    if cal_log["res_max"] < 0.5 and cal_lin["res_max"] > 0.5:
        cal = cal_log
    return cal, pairs, warn

NICE_MULT = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9]
def minor_tick_check(ticks, cal, used_coords, vmin=1e-30, vmax=1e30):
    """Map every tick coord back to a value; residual to nearest n*10^k in px."""
    errs = []
    k0 = int(math.floor(math.log10(vmin))) - 1
    k1 = int(math.ceil(math.log10(vmax))) + 1
    for t in ticks:
        if any(abs(t["coord"] - u) < 0.3 for u in used_coords):
            continue
        best = None
        for k in range(k0, k1 + 1):
            for n in NICE_MULT:
                cand = n * (10.0 ** k)
                d = abs(coord_from_value(cal, cand) - t["coord"])
                if best is None or d < best[0]:
                    best = (d, cand)
        errs.append((best[0], t["coord"], best[1], value_from_coord(cal, t["coord"])))
    errs.sort(reverse=True)
    return errs

# ----------------------------------------------------------------------------- paths
def classify_paths(drawings, frame, page):
    data, framep, fills, swatches, refs = [], [], [], [], []
    for i, d in enumerate(drawings):
        c = colkey(d.get("color"))
        p = path_points(d)
        if d["type"] == "f":
            fills.append(dict(idx=i, bbox=d["rect"], fill=colkey(d.get("fill"))))
            continue
        if len(p) < 2:
            continue
        x0 = min(q[0] for q in p); x1 = max(q[0] for q in p)
        y0 = min(q[1] for q in p); y1 = max(q[1] for q in p)
        length = max(x1 - x0, y1 - y0)
        seg = nsegs(d)
        if seg == 1 and abs(y1 - y0) < 0.35 and 6 < (x1 - x0) < 60:
            if not (near_frame(y0, frame, 3) or near_frame(x0, frame, 3) or near_frame(x1, frame, 3)):
                swatches.append(dict(idx=i, y=0.5*(y0+y1), x0=x0, x1=x1,
                                     color=c, dash=dashkey(d.get("dashes")), width=d.get("width"),
                                     bbox=d["rect"]))
                continue
        if seg <= 2:
            if length > 0.3 * max(page.rect.width, page.rect.height) or is_axis_aligned(d):
                framep.append(dict(idx=i, bbox=d["rect"], len=length, color=c, seg=seg,
                                   dash=dashkey(d.get("dashes"))))
                if seg == 1 and length > 0.3 * max(page.rect.width, page.rect.height):
                    refs.append(dict(idx=i, bbox=d["rect"], color=c, seg=seg))
                continue
        if seg >= 4 and length > 1.5:
            data.append(dict(idx=i, bbox=d["rect"], color=c, width=d.get("width"),
                             dash=dashkey(d.get("dashes")), nseg=seg, pts=p))
        else:
            framep.append(dict(idx=i, bbox=d["rect"], len=length, color=c, seg=seg,
                               dash=dashkey(d.get("dashes"))))
    return data, framep, fills, swatches, refs

def merge_paths(data):
    """Group by (colour, dash, width) and chain fragments whose x-ranges touch."""
    groups = {}
    for d in data:
        groups.setdefault((d["color"], d["dash"], round(d["width"] or 0, 2)), []).append(d)
    curves = []
    for key, frags in groups.items():
        frags = sorted(frags, key=lambda f: f["bbox"].x0)
        chains = []
        for f in frags:
            placed = False
            for ch in chains:
                if f["bbox"].x0 <= ch["x1"] + 2.0 and f["bbox"].x0 >= ch["x0"] - 2.0:
                    ch["frags"].append(f); ch["x1"] = max(ch["x1"], f["bbox"].x1)
                    placed = True
                    break
            if not placed:
                chains.append(dict(frags=[f], x0=f["bbox"].x0, x1=f["bbox"].x1))
        for ch in chains:
            pts = []