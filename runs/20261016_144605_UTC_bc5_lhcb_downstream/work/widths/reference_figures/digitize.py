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
REF = r"E:\AAAwork\LLP-AGI\references\model"

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
            for f in ch["frags"]:
                pp = list(f["pts"])
                if pp and pts and abs(pp[0][0] - pts[-1][0]) > abs(pp[-1][0] - pts[-1][0]):
                    pp = pp[::-1]
                pts.extend(pp)
            seen = set(); ordered = []
            for q in pts:
                k = (round(q[0], 4), round(q[1], 4))
                if k in seen:
                    continue
                seen.add(k); ordered.append(q)
            curves.append(dict(key=key, color=key[0], dash=key[1], width=key[2],
                               frags=[f["idx"] for f in ch["frags"]],
                               frag_bbox=[tuple(round(v, 2) for v in f["bbox"]) for f in ch["frags"]],
                               pts=ordered,
                               bbox=(min(q[0] for q in ordered), min(q[1] for q in ordered),
                                     max(q[0] for q in ordered), max(q[1] for q in ordered))))
    return curves

def dash_str(dash):
    if not dash:
        return "solid"
    if len(dash) >= 2 and dash[0] < 1.5:
        return "dotted"
    return "dashed"

def pick_overlay_color(curve_colors):
    cands = [(1,0,1), (0,1,1), (1,0,0), (0,1,0), (1,1,0), (0.5,0,1), (0,0,0)]
    def d2(a, b):
        return sum((x-y)**2 for x, y in zip(a, b))
    best = max(cands, key=lambda c: min(d2(c, k) for k in curve_colors) if curve_colors else 1)
    return "#%02x%02x%02x" % tuple(int(round(255*x)) for x in best)

# ----------------------------------------------------------------------------- main
def match_label(lines, sw):
    """The legend label is the text line immediately to the right of the swatch
    (all legends in these figures are horizontal, swatch then text)."""
    cands = [l for l in lines
             if (sw["x1"] - 3.0) <= l["x0"] <= (sw["x1"] + 15.0)
             and abs(l["cy"] - sw["y"]) <= 8.0
             and re.sub(r"[\u2013\u2014\u2212\-]", "", l["text"]).strip() != ""]
    if not cands:
        return None
    return min(cands, key=lambda l: (abs(l["cy"] - sw["y"]), abs(l["x0"] - sw["x1"])))["text"]

def process(spec):
    rep = []
    doc = fitz.open(spec["pdf"])
    page = doc[0]
    W, H = page.rect.width, page.rect.height
    lines = text_lines(page)
    drawings = page.get_drawings()
    frame = find_frame(drawings, page)
    xt, yt = find_ticks(drawings, frame)

    rep.append("### page")
    rep.append("")
    rep.append("* file: " + spec["pdf"])
    rep.append("* producer: " + str(doc.metadata.get("creator", "?")))
    rep.append("* pages: %d ; page size: %.1f x %.1f pt (landscape=%s)" % (doc.page_count, W, H, W > H))
    rep.append("* plot frame (PDF user space, y down): x0=%.2f x1=%.2f y_top=%.2f y_bottom=%.2f"
               % (frame["x0"], frame["x1"], frame["y0"], frame["y1"]))
    rep.append("* total vector drawings on page: %d" % len(drawings))
    rep.append("")

    rep.append("### raw extracted text (line level, with bboxes)")
    rep.append("")
    rep.append("| text | x0 | y0 | x1 | y1 | centre |")
    rep.append("|---|---|---|---|---|---|")
    for l in sorted(lines, key=lambda z: (round(z["cy"], 1), z["x0"])):
        t = l["text"].replace("|", "/")
        rep.append("| %s | %.2f | %.2f | %.2f | %.2f | (%.2f, %.2f) |"
                   % (t, l["x0"], l["y0"], l["x1"], l["y1"], l["cx"], l["cy"]))
    rep.append("")

    xlab = [l for l in lines if l["cy"] > frame["y1"] + 0.5 and l["cx"] > frame["x0"] - 60
            and l["cx"] < frame["x1"] + 60]
    ylab = [l for l in lines if l["cx"] < frame["x0"] - 0.5 and l["cy"] > frame["y0"] - 30
            and l["cy"] < frame["y1"] + 30]
    calx, pairsx, warnx = calibrate_axis(xt, xlab, "x", frame)
    caly, pairsy, warny = calibrate_axis(yt, ylab, "y", frame)

    rep.append("### axis scale determination and calibration")
    rep.append("")
    for nm, cal, pairs, ticks, axis in (("x", calx, pairsx, xt, "x"), ("y", caly, pairsy, yt, "y")):
        rep.append("**%s axis**" % nm)
        rep.append("")
        maxlen = max(t["len"] for t in ticks)
        rep.append("* tick strokes detected: %d (major: %d)"
                   % (len(ticks), len([t for t in ticks if t["len"] >= 0.7*maxlen])))
        rep.append("* label/tick pairs used for the fit:")
        for v, c, txt in sorted(pairs, key=lambda z: z[1]):
            rep.append("    * label %s = %s  <->  tick at %.2f" % (txt, fmt(v), c))
        if cal is None:
            rep.append("* **CALIBRATION FAILED**: " + "; ".join(pairs))
            continue
        rep.append("* chosen model: **%s**  coord = a * %s + b"
                   % (cal["scale"], "log10(value)" if cal["scale"] == "log" else "value"))
        rep.append("* a = %.6f , b = %.6f  (px per %s)"
                   % (cal["a"], cal["b"], "decade" if cal["scale"] == "log" else "unit"))
        rep.append("* residual of the fit over the %d labels: max %.4f px, rms %.4f px"
                   % (cal["n"], cal["res_max"], cal["res_rms"]))
        rep.append("* fitted -> measured tick coordinates: " +
                   ", ".join("%.2f->%.2f" % (coord_from_value(cal, v), c)
                             for v, c, _ in sorted(pairs, key=lambda z: z[1])))
        used = [c for _, c, _ in pairs]
        if axis == "x":
            lo = min(value_from_coord(cal, frame["x0"]), value_from_coord(cal, frame["x1"]))
            hi = max(value_from_coord(cal, frame["x0"]), value_from_coord(cal, frame["x1"]))
        else:
            lo = min(value_from_coord(cal, frame["y0"]), value_from_coord(cal, frame["y1"]))
            hi = max(value_from_coord(cal, frame["y0"]), value_from_coord(cal, frame["y1"]))
        chk = minor_tick_check(ticks, cal, used, lo, hi)
        if chk:
            worst = chk[0]
            rep.append("* independent check on the %d tick strokes NOT used in the fit: "
                       "max deviation from the nearest log-nice value = **%.3f px** "
                       "(worst at coord %.2f, nearest nice value %s)"
                       % (len(chk), worst[0], worst[1], fmt(worst[2])))
            rep.append("* value implied by every tick stroke NOT used in the fit "
                       "(coord -> implied value, nearest log-nice value): " + ", ".join(
                "%.2f->%s (%s)" % (c, fmt(impl), fmt(v))
                for _, c, v, impl in sorted(chk, key=lambda z: z[1])))
        if cal["scale"] == "log":
            rep.append("* frame edges in data units: %s = %s at the low edge, %s at the high edge"
                       % (("mass" if axis == "x" else "value"),
                          fmt(value_from_coord(cal, frame["x0"] if axis == "x" else frame["y1"])),
                          fmt(value_from_coord(cal, frame["x1"] if axis == "x" else frame["y0"]))))
        rep.append("")

    data, framep, fills, swatches, refs = classify_paths(drawings, frame, page)
    rep.append("### vector paths")
    rep.append("")
    rep.append("* candidate **data paths** (>= 4 straight segments, real extent): %d" % len(data))
    rep.append("")
    rep.append("| path idx | colour | line width | dash | #segments | bbox (x0,y0,x1,y1) |")
    rep.append("|---|---|---|---|---|---|")
    for d in data:
        rep.append("| %d | %s | %.3f | %s | %d | %s |"
                   % (d["idx"], d["color"], d["width"] or 0, dash_str(d["dash"]), d["nseg"],
                      ", ".join("%.2f" % v for v in d["bbox"])))
    rep.append("")
    rep.append("* frame / axis / tick / reference-line paths: %d (indices %s)"
               % (len(framep), ", ".join(str(f["idx"]) for f in framep)))
    rep.append("* filled areas (not data curves): %s"
               % ("; ".join("#%d fill=%s bbox=%s" % (f["idx"], f["fill"],
                          ", ".join("%.1f" % v for v in f["bbox"])) for f in fills) or "none"))
    rep.append("* long straight reference lines (single stroke, not data): %s"
               % ("; ".join("#%d bbox=%s" % (r["idx"], ", ".join("%.1f" % v for v in r["bbox"]))
                            for r in refs) or "none"))
    rep.append("")

    rep.append("### legend swatches and curve labels")
    rep.append("")
    if swatches:
        rep.append("| swatch idx | colour | dash | x-range | y | nearest label to the right |")
        rep.append("|---|---|---|---|---|---|")
        for s in swatches:
            lab = match_label(lines, s)
            rep.append("| %d | %s | %s | %.2f..%.2f | %.2f | %s |"
                       % (s["idx"], s["color"], dash_str(s["dash"]), s["x0"], s["x1"], s["y"],
                          lab if lab else "(none)"))
    else:
        rep.append("_no legend in this figure_")
    rep.append("")

    curves = merge_paths(data)
    for c in curves:
        best = None
        for s in swatches:
            dc = max(abs(a - b) for a, b in zip(s["color"] or (0, 0, 0), c["color"] or (0, 0, 0)))
            dd = 0.0 if dash_str(s["dash"]) == dash_str(c["dash"]) else 1.0
            score = dc + 0.5 * dd
            if best is None or score < best[0]:
                best = (score, s)
        label = None
        if best and best[0] < 0.1:
            label = match_label(lines, best[1])
        c["label"] = label
        c["slug"] = spec["slug_map"].get(label) if label is not None else None

    if spec["manual"]:
        for c in curves:
            if c["slug"] is None and c["color"] is not None:
                for col, slug in spec["manual"].items():
                    if max(abs(a-b) for a, b in zip(col, c["color"])) < COLORTOL:
                        c["slug"] = slug

    matched_sw = set()
    for c in curves:
        for s in swatches:
            if c["label"] is not None and match_label(lines, s) == c["label"] and \
               max(abs(a - b) for a, b in zip(s["color"] or (0, 0, 0), c["color"] or (0, 0, 0))) < 0.1:
                matched_sw.add(s["idx"])
    unmatched_sw = [s for s in swatches if s["idx"] not in matched_sw]

    rep.append("### curve_id -> path assignment")
    rep.append("")
    if unmatched_sw:
        rep.append("Legend swatches that match NO curve (reported, not digitized): " +
                   "; ".join("#%d colour=%s dash=%s y=%.2f label=%s"
                             % (s["idx"], s["color"], dash_str(s["dash"]), s["y"],
                                match_label(lines, s) or "(no text label)") for s in unmatched_sw))
        rep.append("")
    rep.append("| curve_id | path indices | colour | dash | n. points | x range (px) | y range (px) | source of the name |")
    rep.append("|---|---|---|---|---|---|---|---|")
    for c in curves:
        slug = c["slug"] or ("UNIDENTIFIED_" + "_".join("%.2f" % v for v in c["color"]) +
                             "_" + dash_str(c["dash"]))
        if c["label"]:
            src = "legend label " + c["label"]
        elif c["slug"]:
            src = spec["slug_notes"].get(c["slug"], "caption / line colour")
        else:
            src = "no legend match - NOT identified"
        rep.append("| %s | %s | %s | %s | %d | %.2f..%.2f | %.2f..%.2f | %s |"
                   % (slug, ", ".join(str(i) for i in c["frags"]), c["color"],
                      dash_str(c["dash"]), len(c["pts"]), c["bbox"][0], c["bbox"][2],
                      c["bbox"][1], c["bbox"][3], src))
    rep.append("")

    for c in curves:
        xs = np.array([p[0] for p in c["pts"]])
        ys = np.array([p[1] for p in c["pts"]])
        order = np.argsort(xs, kind="stable")
        xs, ys = xs[order], ys[order]
        if len(xs) > MAXPTS:
            keep = np.unique(np.round(np.linspace(0, len(xs) - 1, MAXPTS)).astype(int))
            xs, ys = xs[keep], ys[keep]
        c["mass"] = np.array([value_from_coord(calx, v) for v in xs])
        c["value"] = np.array([value_from_coord(caly, v) for v in ys])
        c["n_out"] = len(xs)

    name = {"decayrates": "decayrates_curves.csv",
            "scalardecaybr": "scalardecaybr_curves.csv",
            "tauscalar": "tauscalar_curves.csv"}[spec["key"]]
    out_csv = os.path.join(HERE, name)
    rows = []
    for c in sorted(curves, key=lambda z: (z["slug"] or "zz", z["bbox"][0])):
        slug = c["slug"] or ("UNIDENTIFIED_" + "_".join("%.2f" % v for v in c["color"]))
        for m, v in zip(c["mass"], c["value"]):
            rows.append((slug, m, v, spec["units"]))
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["curve_id", "mass_GeV", "value", "value_units"])
        for slug, m, v, u in rows:
            w.writerow([slug, fmt(m), fmt(v), u])

    pix = page.get_pixmap(dpi=DPI)
    img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
    fig, ax = plt.subplots(figsize=(pix.width / 100.0, pix.height / 100.0), dpi=100)
    ax.imshow(img, extent=[0, W, H, 0], interpolation="nearest")
    ov = pick_overlay_color([c["color"] for c in curves if c["color"]])
    for c in curves:
        xs = np.array([p[0] for p in c["pts"]])
        ys = np.array([p[1] for p in c["pts"]])
        ax.plot(xs, ys, marker="o", ms=1.6, mew=0, ls="none", color=ov, alpha=0.9)
    ax.set_xlim(0, W); ax.set_ylim(H, 0)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color("black")
    fig.subplots_adjust(0, 0, 1, 1)
    out_png = os.path.join(HERE, "overlay_%s.png" % spec["key"])
    fig.savefig(out_png, dpi=100)
    plt.close(fig)

    doc.close()
    return rep, curves, dict(calx=calx, caly=caly, frame=frame, W=W, H=H,
                             n_draw=len(drawings), n_data=len(data), n_frame=len(framep),
                             fills=fills, swatches=swatches, refs=refs,
                             lines=lines, out_csv=out_csv, out_png=out_png,
                             warnx=warnx, warny=warny, overlay=ov,
                             swatches_unmatched=unmatched_sw)

def main():
    report = []
    report.append("# Digitization report - reference figures for scalar decay validation")
    report.append("")
    report.append("Generated by digitize.py (PyMuPDF vector extraction, no image tracing).")
    report.append("")
    for spec in FIGURE_SPECS:
        rep, curves, info = process(spec)
        report.append("---")
        report.append("")
        report.append("## Figure %s - %s" % (spec["key"], os.path.basename(spec["pdf"])))
        report.append("")
        report.append("* %s" % spec["title"])
        report.append("* value units assigned: **%s**" % spec["units"])
        report.append("* CSV written: %s" % os.path.basename(info["out_csv"]))
        report.append("* overlay written: %s (digitized points drawn as %s dots on the 200 dpi page render)"
                      % (os.path.basename(info["out_png"]), info["overlay"]))
        report.append("")
        report.extend(rep)

        report.append("### axis ranges recovered")
        report.append("")
        cx, cy, f = info["calx"], info["caly"], info["frame"]
        if cx:
            report.append("* mass axis: %.4g GeV (left frame edge) .. %.4g GeV (right frame edge)"
                          % (value_from_coord(cx, f["x0"]), value_from_coord(cx, f["x1"])))
        if cy:
            report.append("* value axis: %.4g (top frame edge) .. %.4g (bottom frame edge)"
                          % (value_from_coord(cy, f["y0"]), value_from_coord(cy, f["y1"])))
        report.append("")
        for c in curves:
            slug = c["slug"] or "UNIDENTIFIED"
            d = np.diff(c["mass"])
            report.append("* %s: %d output points, mass %.4g .. %.4g GeV, value %.4g .. %.4g, "
                          "x monotonic = %s" % (slug, c["n_out"], c["mass"].min(), c["mass"].max(),
                                                c["value"].min(), c["value"].max(),
                                                bool(np.all(d >= -1e-12))))
        report.append("")
        span = [c for c in curves if c["mass"].min() < 2.0 < c["mass"].max()]
        if span:
            report.append("### continuity across the 2 GeV evaluation boundary")
            report.append("")
            report.append("| curve_id | last value below 2 GeV | first value above 2 GeV | ratio |")
            report.append("|---|---|---|---|")
            for c in span:
                lo = c["value"][c["mass"] <= 2.0]
                hi = c["value"][c["mass"] >= 2.0]
                a = float(lo[-1]); b = float(hi[0])
                report.append("| %s | %.6g | %.6g | %.4f |"
                              % (c["slug"] or "UNIDENTIFIED", a, b, (b / a) if a else float("nan")))
            report.append("")
        report.append("### warnings / uncertainty for this figure")
        report.append("")
        ws = list(info["warnx"]) + list(info["warny"])
        if info["fills"]:
            ws.append("filled grey area(s) present and NOT digitized as curves: %s"
                      % "; ".join("#%d bbox=%s" % (fi["idx"], ", ".join("%.2f" % v for v in fi["bbox"]))
                                  for fi in info["fills"]))
        for c in curves:
            if c["slug"] and c["value"].max() > 1.0 and spec["units"] == "dimensionless":
                ws.append("curve %s reaches %.4g > 1 (above the 1.0 line); this is what the "
                          "published vector figure draws - it is not a digitization error."
                          % (c["slug"], c["value"].max()))
        for r in info["refs"]:
            if abs(r["bbox"].x0 - 193.0) < 0.6 and abs(r["bbox"].x1 - 193.0) < 0.6:
                ws.append("path #%d is the vertical line at m = 2.0 GeV (PDF x = 193.0); it is "
                          "the boundary between the dispersive (m < 2 GeV) and the perturbative "
                          "quark/gluon (m > 2 GeV) evaluation. It is excluded from the data."
                          % r["idx"])
        if info["refs"]:
            ws.append("long straight reference line(s) excluded from the data: %s"
                      % "; ".join("#%d bbox=%s" % (r["idx"], ", ".join("%.2f" % v for v in r["bbox"]))
                                  for r in info["refs"]))
        for s in info["swatches_unmatched"]:
            ws.append("legend swatch #%d (colour %s, %s, y=%.2f) matches NO data curve and has "
                      "the text label %r - it is an unexplained legend entry, not digitized."
                      % (s["idx"], s["color"], dash_str(s["dash"]), s["y"],
                         match_label(info["lines"], s)))
        for c in curves:
            if c["slug"] is None:
                ws.append("path group %s (colour %s, %s, bbox %s) could NOT be matched to a "
                          "legend label / caption: emitted as UNIDENTIFIED."
                          % (c["frags"], c["color"], dash_str(c["dash"]),
                             ", ".join("%.2f" % v for v in c["bbox"])))
        if not ws:
            ws.append("none beyond the general caveats below")
        for w in ws:
            report.append("* %s" % w)
        report.append("")

        print("=" * 78)
        print(spec["key"], "->", os.path.basename(info["out_csv"]))
        for c in curves:
            print("   %-42s paths %-22s %5d pts  mass %.4g..%.4g  val %.4g..%.4g"
                  % (c["slug"] or "UNIDENTIFIED", str(c["frags"]), c["n_out"],
                     c["mass"].min(), c["mass"].max(), c["value"].min(), c["value"].max()))

    report.append("---")
    report.append("")
    report.append("## General caveats")
    report.append("")
    report.append("* All coordinates are PDF user-space points of the single page of each file "
                  "(origin top-left, y increasing downwards, 72 pt = 1 inch); extraction uses the "
                  "vector path centre lines, so stroke width does not bias the values.")
    report.append("* Curves clipped by the axis frame keep the clipped (frame) value at their end "
                  "points; those points are true plotted values but are only bounds on the "
                  "underlying model.")
    report.append("* Near-vertical parts of a curve give several values at (almost) the same mass; "
                  "the CSVs are sorted by mass but are NOT single-valued functions there.")
    report.append("* Ticks are read from vector strokes and tick labels from the text layer; the "
                  "quoted residual is the max |fitted - measured| tick position in PDF points.")
    report.append("")
    with open(os.path.join(HERE, "digitization_report.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(report))
    print("\nwrote", os.path.join(HERE, "digitization_report.md"))

if __name__ == "__main__":
    main()
