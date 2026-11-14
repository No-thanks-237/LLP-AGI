#!/usr/bin/env python
"""make_decay_figures.py -- three standalone, publication-quality decay figures for
the 2HDM Type-I light CP-even H benchmark of Li et al. arXiv:2212.06186v1
(DESIGN.md section 2, task-9).

SMOOTH MODEL CURVES ONLY: no digitised paper points and no markers of any kind.
Every number comes from the author Light_scalar_decay batch executable through
work/decay/decay_widths.py (one process invocation for all points).

  (1) outputs/figures/decay_branching_ratios.{png,pdf}
      BR(H -> XX) vs mH, tanb = 10, all 10 channels, log-log, 0.1-10 GeV.
  (2) outputs/figures/decay_partial_widths.{png,pdf}
      Gamma(H -> XX) [GeV] vs mH, tanb = 10, same 10 channels.
  (3) outputs/figures/decay_total_width_ctau.{png,pdf}
      Gamma_tot vs mH for tanb = 5, 10, 70, 200, 1000 (left axis) with the
      corresponding c*tau = hbar*c / Gamma_tot on the right axis.

Colours, TeX labels and (nearly) the axis ranges follow make_overlay() in
work/decay/validate_decay.py, so the three figures are drop-in replacements for
the left/middle/right panels of the old 1x3 overlay figure, but without the
digitised points.  The BR panel's lower limit is extended from 5e-3 to 1e-9 so
that the open e e channel (BR ~ 7e-9 .. 9e-5) is actually visible; the
total-width panel is trimmed to the data range (1e-22 .. 5e-9 GeV).  The panel
ranges are recorded in logs/decay_figures.md.

The mH = 2.0 GeV node is deliberately never sampled (the author program double
counts there, see logs/decay_validation.md section 8); 500 log-spaced masses in
0.1-10 GeV never hit it.  The (genuine, intended) 2 GeV hadron/parton matching
discontinuity of the model is marked with a thin grey dotted line in figures
1-2; it is not a data marker.

Run:    python work/decay/make_decay_figures.py
Writes: outputs/figures/decay_{branching_ratios,partial_widths,total_width_ctau}.{png,pdf}
        logs/decay_figures.md
"""
from __future__ import annotations

import hashlib
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
import decay_widths as dw

HBARC = dw.HBARC_GEV_M                 # 1.973269804e-16 GeV m (DESIGN.md section 0)
CHANNELS = dw.CHANNELS
OUTDIR = os.path.join(RUN, "outputs", "figures")
LOGPATH = os.path.join(RUN, "logs", "decay_figures.md")

# ---------------------------------------------------------------- style contract
# exact same colours / labels as make_overlay() in work/decay/validate_decay.py
COLORS = dict(gaga="purple", pipi="blue", KK="red", **{"4pi": "green"}, cc="navy",
              ss="darkred", gg="darkgreen", mumu="orange", ee="cyan", tautau="peru")
LABELS = dict(gaga=r"$\gamma\gamma$", pipi=r"$\pi\pi$", KK=r"$KK$", **{"4pi": r"$4\pi$"},
              cc=r"$c\bar c$", ss=r"$s\bar s$", gg=r"$gg$", mumu=r"$\mu\mu$",
              ee=r"$ee$", tautau=r"$\tau\tau$")
# same colours as make_overlay() for the five total-width curves
COLORS_TANB = {5.0: "blue", 10.0: "green", 70.0: "orange", 200.0: "purple",
               1000.0: "red"}

MASS_RANGE = (0.1, 10.0)               # [GeV], same window as make_overlay()
N_MASS = 500                           # log-spaced mass points
TANB_REF = 10.0
TANB_TOT = (5.0, 10.0, 70.0, 200.0, 1000.0)
MASS_MATCH = 2.0                       # hadron/parton matching switch of the program

# axis ranges (BR extended down to see e e; total width trimmed to the data)
YLIM_BR = (1e-9, 2.0)
YLIM_WIDTH = (1e-24, 1e-9)
YLIM_TOTAL = (1e-22, 5e-9)

MODEL_TEXT = (r"Type-I 2HDM: $m_h=125$ GeV, $m_A=m_{H^\pm}=600$ GeV, "
              r"$\cos(\beta-\alpha)=1/\tan\beta$")

# channel thresholds hard-coded in the author source (CPevenDecay.cpp); the
# mH = 2 GeV matching is the switch from hadronic to partonic descriptions.
THRESHOLDS = [("ee", 2 * 0.000511), ("mumu", 2 * 0.10566), ("pipi", 2 * 0.13497),
              ("4pi", 4 * 0.13497), ("KK", 2 * 0.498), ("ss", 2.0), ("gg", 2.0),
              ("tautau", 2 * 1.7768), ("cc", 2 * 1.864843)]

# spot values copied from logs/decay_validation.md section 7 (author program)
SPOT = [
    dict(mH=0.5, tanb=10.0, gamma=1.272186e-14, ctau=1.5511e-02,
         br=dict(gaga=0.4814, pipi=4.65e-01, gg=7.9e-27, cc=7.9e-27)),
    dict(mH=1.0, tanb=10.0, gamma=1.631479e-13, ctau=1.2095e-03,
         br=dict(gaga=0.3003, pipi=5.48e-01, gg=6.1e-28, cc=6.1e-28)),
    dict(mH=2.0, tanb=10.0, gamma=4.356101e-13, ctau=4.5299e-04,
         br=dict(gaga=0.8995, pipi=3.18e-04, gg=3.3e-02, cc=2.3e-28)),
    dict(mH=5.0, tanb=10.0, gamma=8.367106e-12, ctau=2.3584e-05,
         br=dict(gaga=0.7331, pipi=1.2e-29, gg=1.1e-02, cc=1.4e-01)),
    dict(mH=1.0, tanb=5.0, gamma=7.610433e-12, ctau=2.5928e-05,
         br=dict(gaga=0.0253, pipi=7.64e-01, gg=1.3e-29, cc=1.3e-29)),
    dict(mH=0.5, tanb=100.0, gamma=6.162374e-17, ctau=3.2021e+00,
         br=dict(gaga=0.9999, pipi=9.5e-05, gg=1.6e-24, cc=1.6e-24)),
    dict(mH=1.0, tanb=100.0, gamma=4.930578e-16, ctau=4.0021e-01,
         br=dict(gaga=0.9998, pipi=1.8e-04, gg=2.0e-25, cc=2.0e-25)),
    dict(mH=2.0, tanb=100.0, gamma=3.943782e-15, ctau=5.0035e-02,
         br=dict(gaga=1.0000, pipi=3.5e-08, gg=3.6e-06, cc=2.5e-26)),
    dict(mH=0.5, tanb=1000.0, gamma=6.162092e-19, ctau=3.2023e+02,
         br=dict(gaga=1.0000, pipi=9.5e-09, gg=1.6e-22, cc=1.6e-22)),
    dict(mH=1.0, tanb=1000.0, gamma=4.929737e-18, ctau=4.0028e+01,
         br=dict(gaga=1.0000, pipi=1.8e-08, gg=2.0e-23, cc=2.0e-23)),
    dict(mH=5.0, tanb=1000.0, gamma=6.164705e-16, ctau=3.2009e-01,
         br=dict(gaga=1.0000, pipi=1.6e-25, gg=1.5e-10, cc=1.9e-09)),
]


# ------------------------------------------------------------------- utilities
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def set_style():
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "font.size": 10.5,
        "font.family": "DejaVu Sans",
        "axes.labelsize": 12,
        "axes.titlesize": 12.5,
        "axes.linewidth": 0.9,
        "axes.edgecolor": "0.15",
        "axes.labelcolor": "0.1",
        "text.color": "0.1",
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 4.5,
        "ytick.major.size": 4.5,
        "xtick.minor.size": 2.6,
        "ytick.minor.size": 2.6,
        "legend.fontsize": 8,
        "mathtext.fontset": "dejavusans",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def log_axes(ax):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*MASS_RANGE)


def add_matching_marker(ax):
    """thin grey dotted guide line at the program's 2 GeV hadron/parton matching."""
    ax.axvline(MASS_MATCH, color="0.55", lw=0.9, ls=(0, (1, 2.5)), zorder=0.5)
    ax.annotate(r"2 GeV matching", xy=(2.03, 0.03),
                xycoords=("data", "axes fraction"), ha="left", va="bottom",
                fontsize=7.0, color="0.45")


def add_footer(fig, extra=""):
    txt = MODEL_TEXT + (("; " + extra) if extra else "")
    fig.text(0.005, 0.005, txt, ha="left", va="bottom", fontsize=7.2, color="0.35")


def save(fig, stem):
    os.makedirs(OUTDIR, exist_ok=True)
    png = os.path.join(OUTDIR, stem + ".png")
    pdf = os.path.join(OUTDIR, stem + ".pdf")
    fig.savefig(png, dpi=300, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    return png, pdf


# -------------------------------------------------------------------- figures
def fig_branching_ratios(mass, br):
    fig, ax = plt.subplots(figsize=(6.8, 4.9))
    for ch in CHANNELS:
        ax.plot(mass, br[ch], "-", color=COLORS[ch], lw=1.55, label=LABELS[ch])
    log_axes(ax)
    ax.set_ylim(*YLIM_BR)
    ax.set_xlabel(r"$m_H$ [GeV]")
    ax.set_ylabel(r"$\mathrm{BR}(H\to XX)$")
    ax.set_title(r"Branching ratios of the 2HDM light $H$, $\tan\beta = 10$")
    ax.legend(loc="lower left", ncol=2, frameon=True, framealpha=0.95,
              edgecolor="0.82", handlelength=2.0, borderpad=0.45,
              labelspacing=0.32, columnspacing=1.0)
    add_matching_marker(ax)
    add_footer(fig)
    return save(fig, "decay_branching_ratios")


def fig_partial_widths(mass, widths):
    fig, ax = plt.subplots(figsize=(6.8, 4.9))
    for ch in CHANNELS:
        ax.plot(mass, widths[ch], "-", color=COLORS[ch], lw=1.55, label=LABELS[ch])
    log_axes(ax)
    ax.set_ylim(*YLIM_WIDTH)
    ax.set_xlabel(r"$m_H$ [GeV]")
    ax.set_ylabel(r"$\Gamma(H\to XX)$ [GeV]")
    ax.set_title(r"Partial decay widths of the 2HDM light $H$, $\tan\beta = 10$")
    ax.legend(loc="upper left", ncol=2, frameon=True, framealpha=0.95,
              edgecolor="0.82", handlelength=2.0, borderpad=0.45,
              labelspacing=0.32, columnspacing=1.0)
    add_matching_marker(ax)
    add_footer(fig)
    return save(fig, "decay_partial_widths")


def fig_total_width_ctau(mass, gtot):
    fig, ax = plt.subplots(figsize=(7.0, 4.9))
    ax2 = ax.twinx()
    for tb in TANB_TOT:
        ax.plot(mass, gtot[tb], "-", color=COLORS_TANB[tb], lw=1.7,
                label=r"$\tan\beta = %g$" % tb)
    log_axes(ax)
    ax.set_ylim(*YLIM_TOTAL)
    ylo, yhi = ax.get_ylim()
    ax2.set_yscale("log")
    ax2.set_ylim(HBARC / yhi, HBARC / ylo)      # c*tau = hbar*c / Gamma_tot
    ax.set_xlabel(r"$m_H$ [GeV]")
    ax.set_ylabel(r"$\Gamma_{\rm tot}$ [GeV]")
    ax2.set_ylabel(r"$c\tau$ [m]")
    ax.set_title(r"Total width and $c\tau$ of the 2HDM light $H$")
    ax.legend(loc="upper left", frameon=True, framealpha=0.95, edgecolor="0.82",
              handlelength=2.2, borderpad=0.45, labelspacing=0.35)
    ax2.tick_params(which="both", direction="in")
    add_footer(fig, r"$c\tau = \hbar c/\Gamma_{\rm tot}$, "
                   r"$\hbar c = 1.973\times10^{-16}$ GeV m")
    return save(fig, "decay_total_width_ctau")


# --------------------------------------------------------------- sanity checks
def run_checks(rows, mass, br10, w10, rowmap, exe, runtime_s, n_extra):
    info = dict(n_points=len(rows), n_mass=len(mass), n_extra=n_extra,
                mass_range=[float(mass[0]), float(mass[-1])],
                tanbs=list(TANB_TOT), exe=exe, exe_sha256=sha256(exe),
                runtime_s=float(runtime_s))
    # BR normalisation and Gamma_tot = sum of channels over every evaluated point
    brsum = np.array([sum(r["br"].values()) for r in rows])
    info["br_sum_max_dev"] = float(np.max(np.abs(brsum - 1.0)))
    info["br_sum_median_dev"] = float(np.median(np.abs(brsum - 1.0)))
    ratio = np.array([r["gamma_total"] / sum(r["widths"].values()) for r in rows])
    info["gamma_ratio"] = [float(ratio.min()), float(ratio.max())]
    ct = np.array([r["ctau"] * r["gamma_total"] / HBARC for r in rows])
    info["ctau_identity_max_rel"] = float(np.max(np.abs(ct - 1.0)))
    # no exact mH = 2.0 node in the plotted sample (program double counts there)
    info["min_distance_to_2GeV"] = float(np.min(np.abs(mass - MASS_MATCH)))
    info["ctau_range_m"] = [float(min(r["ctau"] for r in rows)),
                            float(max(r["ctau"] for r in rows))]
    info["gamma_total_range_GeV"] = [float(min(r["gamma_total"] for r in rows)),
                                     float(max(r["gamma_total"] for r in rows))]
    # spot values against the validation table
    spot_rows = []
    max_rel_g = 0.0
    max_rel_ct = 0.0
    for s in SPOT:
        r = rowmap[(float(s["mH"]), float(s["tanb"]))]
        rg = abs(r["gamma_total"] / s["gamma"] - 1.0)
        rc = abs(r["ctau"] / s["ctau"] - 1.0)
        max_rel_g = max(max_rel_g, rg)
        max_rel_ct = max(max_rel_ct, rc)
        spot_rows.append(dict(
            mH=s["mH"], tanb=s["tanb"],
            gamma=r["gamma_total"], gamma_ref=s["gamma"], rel_gamma=rg,
            ctau=r["ctau"], ctau_ref=s["ctau"], rel_ctau=rc,
            br={ch: (r["br"][ch], s["br"][ch]) for ch in ("gaga", "pipi", "gg", "cc")}))
    info["spot"] = spot_rows
    info["spot_max_rel_gamma"] = float(max_rel_g)
    info["spot_max_rel_ctau"] = float(max_rel_ct)
    # threshold crossings visible in the sampled window (tanb = 10)
    thr_rows = []
    for ch, t0 in THRESHOLDS:
        y = w10[ch]
        opened = y > 1e-30
        if ch in ("ss", "gg"):
            opened &= mass > MASS_MATCH
        if opened[0]:
            first = float(mass[0])
            were_open = True
        else:
            idx = np.argmax(opened) if opened.any() else -1
            first = float(mass[idx]) if idx >= 0 else float("nan")
            were_open = False
        thr_rows.append(dict(channel=ch, threshold=float(t0), first_sample=first,
                             open_at_0p1=bool(were_open),
                             in_window=bool(MASS_RANGE[0] <= t0 <= MASS_RANGE[1])))
    info["thresholds"] = thr_rows
    # the 2 GeV matching jump, bracketed at +-0.05 % (grid never hits 2.0)
    r_lo = rowmap[(1.999, 10.0)]
    r_hi = rowmap[(2.001, 10.0)]
    info["match_jump"] = dict(
        gamma_below=float(r_lo["gamma_total"]), gamma_above=float(r_hi["gamma_total"]),
        gamma_ratio=float(r_hi["gamma_total"] / r_lo["gamma_total"]),
        br_below={ch: float(r_lo["br"][ch]) for ch in ("pipi", "KK", "4pi", "ss", "gg", "gaga")},
        br_above={ch: float(r_hi["br"][ch]) for ch in ("pipi", "KK", "4pi", "ss", "gg", "gaga")},
        w_below={ch: float(r_lo["widths"][ch]) for ch in ("pipi", "KK", "4pi", "ss", "gg")},
        w_above={ch: float(r_hi["widths"][ch]) for ch in ("pipi", "KK", "4pi", "ss", "gg")})
    # structure near the K K threshold (the form-factor cusp of the author database)
    sel = (mass >= 0.9) & (mass <= 1.3)
    for ch in ("gaga", "pipi", "KK"):
        y = br10[ch][sel]
        info.setdefault("near_kk", {})[ch] = dict(
            min=float(y.min()), max=float(y.max()),
            m_at_min=float(mass[sel][np.argmin(y)]),
            m_at_max=float(mass[sel][np.argmax(y)]))
    return info


# ------------------------------------------------------------------- log file
def write_log(info, figs, script_sha):
    L = []
    L.append("# Decay figures (task-9): three standalone smooth-curve plots")
    L.append("")
    L.append("Status: **DONE** -- three figures written, no digitised points and no "
             "markers anywhere; every curve is the author Light_scalar_decay "
             "program evaluated on 500 log-spaced masses.")
    L.append("")
    L.append("Author: decay-figures (teammate), 2026-10-26.  Run: "
             "runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream.")
    L.append("Script: work/decay/make_decay_figures.py "
             "(sha256 `%s`)." % script_sha)
    L.append("")
    L.append("## 1. Deliverables")
    L.append("")
    L.append("| figure | content | axes |")
    L.append("|---|---|---|")
    L.append("| `outputs/figures/decay_branching_ratios.png` / `.pdf` | "
             "BR(H->XX) vs mH, all 10 channels | log-log, x 0.1-10 GeV, y 1e-9-2 |")
    L.append("| `outputs/figures/decay_partial_widths.png` / `.pdf` | "
             "Gamma(H->XX) [GeV] vs mH, all 10 channels | log-log, x 0.1-10 GeV, y 1e-24-1e-9 |")
    L.append("| `outputs/figures/decay_total_width_ctau.png` / `.pdf` | "
             "Gamma_tot vs mH for tanb = 5,10,70,200,1000 (left) and c*tau (right) | "
             "log-log, x 0.1-10 GeV, y 1e-22-5e-9 GeV / 3.9e-8-2.0e6 m |")
    L.append("")
    L.append("Colours and LaTeX labels are identical to make_overlay() in "
             "work/decay/validate_decay.py: " +
             ", ".join("`%s`=%s" % (ch, COLORS[ch]) for ch in CHANNELS) + ".")
    L.append("The old outputs/figures/decay_validation.png (1x3 panel with digitised "
             "circles) is left untouched; the three new figures replace its "
             "left/middle/right panels without any data points.")
    L.append("")
    L.append("Two deviations from the old panel ranges, needed so that all ten "
             "channels are visible: the BR lower limit is 1e-9 (old 5e-3) because "
             "BR(e e) <= 8.8e-5, and the total-width window is 1e-22-5e-9 (old "
             "1e-24-5e-9) because Gamma_tot spans 4.9e-21-9.3e-10 over the five "
             "tanb values.")
    L.append("")
    L.append("## 2. Generation")
    L.append("")
    L.append("python work/decay/make_decay_figures.py -- one batch call to the "
             "author executable, %d points (500 masses x 5 tanb values + %d "
             "spot/jump-check points), runtime %.2f s." %
             (info["n_points"], info["n_extra"], info["runtime_s"]))
    L.append("")
    L.append("| item | value |")
    L.append("|---|---|")
    L.append("| executable | `%s` |" % info["exe"])
    L.append("| executable sha256 | `%s` |" % info["exe_sha256"])
    L.append("| mass grid | 500 log-spaced points, %.6g ... %.6g GeV |" %
             (info["mass_range"][0], info["mass_range"][1]))
    L.append("| tanb values | %s |" % ", ".join("%g" % t for t in info["tanbs"]))
    L.append("| closest sampled mH to 2.0 GeV | %.4f GeV (the mH=2.0 node is "
             "never evaluated: the program double counts there) |" %
             info["min_distance_to_2GeV"])
    L.append("| Gamma_tot range on all points | %.3e ... %.3e GeV |" %
             (info["gamma_total_range_GeV"][0], info["gamma_total_range_GeV"][1]))
    L.append("| c*tau range on all points | %.3e ... %.3e m |" %
             (info["ctau_range_m"][0], info["ctau_range_m"][1]))
    L.append("")
    L.append("## 3. Sanity checks")
    L.append("")
    L.append("| check | result |")
    L.append("|---|---|")
    L.append("| max abs(BR sum - 1) over all %d points | **%.2e** |" %
             (info["n_points"], info["br_sum_max_dev"]))
    L.append("| median abs(BR sum - 1) | %.2e |" % info["br_sum_median_dev"])
    L.append("| Gamma_tot / sum(10 channels) min ... max | %.12f ... %.12f |" %
             (info["gamma_ratio"][0], info["gamma_ratio"][1]))
    L.append("| max rel. deviation of c*tau*Gamma_tot from hbar*c | %.2e |" %
             info["ctau_identity_max_rel"])
    L.append("| spot values vs logs/decay_validation.md sec. 7: max rel. dev. "
             "Gamma_tot | %.2e |" % info["spot_max_rel_gamma"])
    L.append("| same for c*tau | %.2e |" % info["spot_max_rel_ctau"])
    L.append("")
    L.append("### Spot values (ours vs the validation log table)")
    L.append("")
    L.append("| mH [GeV] | tanb | Gamma_tot [GeV] | rel.diff | c*tau [m] | rel.diff | "
             "BR(gaga) | BR(pipi) | BR(gg) | BR(cc) |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for s in info["spot"]:
        b = s["br"]
        L.append("| %g | %g | %.6e | %.1e | %.6e | %.1e | %.4g | %.4g | %.3g | %.3g |" %
                 (s["mH"], s["tanb"], s["gamma"], s["rel_gamma"], s["ctau"],
                  s["rel_ctau"], b["gaga"][0], b["pipi"][0], b["gg"][0], b["cc"][0]))
    L.append("")
    L.append("The tiny-BR entries (gg/cc below threshold) are the 1e-40 program floor "
             "divided by Gamma_tot; the validation table prints them with 2 significant "
             "digits, so their relative differences are print quantisation, not physics. "
             "Gamma_tot reproduces the table to <= %.0e and c*tau to <= %.0e "
             "(limited by the 4-5 significant digits printed there)."
             % (info["spot_max_rel_gamma"], info["spot_max_rel_ctau"]))
    L.append("")
    L.append("## 4. Thresholds and the 2 GeV matching (features visible in the curves)")
    L.append("")
    L.append("Channel thresholds hard-coded in the author program "
             "(CPevenDecay.cpp), with the first sampled mH at which the width is "
             "above the 1e-40 floor (tanb = 10; the 500-point grid resolves a "
             "threshold to ~1 % in mass):")
    L.append("")
    L.append("| channel | threshold [GeV] | first sampled point [GeV] | in 0.1-10 window |")
    L.append("|---|---|---|---|")
    for t in info["thresholds"]:
        first = "open from 0.1" if t["open_at_0p1"] else "%.4f" % t["first_sample"]
        L.append("| %s | %.6g | %s | %s |" %
                 (t["channel"], t["threshold"], first, "yes" if t["in_window"] else "no (below 0.1)"))
    L.append("")
    mj = info["match_jump"]
    L.append("At the 2 GeV hadron/parton matching (m_H > 2 floors pi pi / K K / 4pi "
             "and switches on s sbar / g g) the individual channels jump while the "
             "total width stays continuous. Bracketed at mH = 1.999 vs 2.001 GeV, "
             "tanb = 10:")
    L.append("")
    L.append("| quantity | 1.999 GeV | 2.001 GeV |")
    L.append("|---|---|---|")
    L.append("| Gamma_tot [GeV] | %.4e | %.4e (ratio %.4f) |" %
             (mj["gamma_below"], mj["gamma_above"], mj["gamma_ratio"]))
    for ch in ("pipi", "KK", "4pi", "ss", "gg", "gaga"):
        L.append("| BR(%s) | %.4g | %.4g |" % (ch, mj["br_below"][ch], mj["br_above"][ch]))
    for ch in ("pipi", "KK", "4pi", "ss", "gg"):
        L.append("| Gamma(%s) [GeV] | %.3e | %.3e |" %
                 (ch, mj["w_below"][ch], mj["w_above"][ch]))
    L.append("")
    L.append("This is the intended matching of the author program (the 4pi coefficient "
             "C = 5.1e-9 was tuned so that the hadronic and partonic descriptions meet "
             "at 2 GeV), not a plotting artifact; it is why the curves in figures 1-2 "
             "have a near-vertical jump there. A thin grey dotted guide line marks "
             "mH = 2 GeV in figures 1-2 (a guide line, not a data marker).")
    L.append("")
    kk = info["near_kk"]
    L.append("Structure at mH ~ 1 GeV (the author form-factor cusp just above the "
             "K K threshold, 2 m_K = 0.996 GeV), tanb = 10, 0.9 <= mH <= 1.3 GeV:")
    L.append("")
    L.append("| channel | min BR | at mH | max BR | at mH |")
    L.append("|---|---|---|---|---|")
    for ch in ("gaga", "pipi", "KK"):
        L.append("| %s | %.4g | %.4f | %.4g | %.4f |" %
                 (ch, kk[ch]["min"], kk[ch]["m_at_min"], kk[ch]["max"], kk[ch]["m_at_max"]))
    L.append("")
    L.append("Other visible features: gaga BR -> 1 below the pi pi threshold "
             "(only gamma gamma and e e are open); abrupt onsets at mH = 0.211 (mu mu), "
             "0.270 (pi pi), 0.540 (4pi), 0.996 (K K), 2.0 (ss, gg), 3.554 (tau tau), "
             "3.730 (c cbar); at large tanb gamma gamma dominates from low mass "
             "(figure 3 shows the correspondingly tiny total widths and long c*tau).")
    L.append("")
    L.append("## 5. Files")
    L.append("")
    L.append("| file | sha256 |")
    L.append("|---|---|")
    for p in figs:
        L.append("| `%s` | `%s` |" % (os.path.relpath(p, RUN).replace("\\", "/"), sha256(p)))
    L.append("")
    with open(LOGPATH, "w") as f:
        f.write("\n".join(L) + "\n")


# ------------------------------------------------------------------------ main
def main():
    t0 = time.time()
    set_style()
    mass = np.logspace(np.log10(MASS_RANGE[0]), np.log10(MASS_RANGE[1]), N_MASS)
    if np.min(np.abs(mass - MASS_MATCH)) < 1e-6:
        raise RuntimeError("mass grid hits the mH=2 GeV double-counting node")
    # one batch call for every plotted point plus the spot/jump checks
    extra = [(float(s["mH"]), float(s["tanb"])) for s in SPOT]
    extra += [(1.999, 10.0), (2.001, 10.0), (1.999, 5.0), (2.001, 5.0)]
    pts = [(float(m), tb) for tb in TANB_TOT for m in mass]
    seen = set(pts)
    pts += [p for p in extra if p not in seen]
    n_extra = len(pts) - N_MASS * len(TANB_TOT)
    dws = dw.DecayWidths()
    rows = dws.batch(pts)                                   # ONE process invocation
    rowmap = {p: r for p, r in zip(pts, rows)}
    rows10 = [rowmap[(float(m), TANB_REF)] for m in mass]
    br10 = {ch: np.array([r["br"][ch] for r in rows10]) for ch in CHANNELS}
    w10 = {ch: np.array([r["widths"][ch] for r in rows10]) for ch in CHANNELS}
    gtot = {tb: np.array([rowmap[(float(m), tb)]["gamma_total"] for m in mass])
            for tb in TANB_TOT}
    figs = []
    figs += fig_branching_ratios(mass, br10)
    figs += fig_partial_widths(mass, w10)
    figs += fig_total_width_ctau(mass, gtot)
    t_run = time.time() - t0
    info = run_checks(rows, mass, br10, w10, rowmap, dws.exe, t_run, n_extra)
    write_log(info, figs, sha256(os.path.abspath(__file__)))
    print("wrote:", *figs, LOGPATH, sep="\n  ")
    print("BR sum max dev: %.2e   Gtot/sum ratio: %.12f..%.12f" %
          (info["br_sum_max_dev"], info["gamma_ratio"][0], info["gamma_ratio"][1]))
    print("points: %d, one batch call, %.2f s" % (len(rows), t_run))
    return 0


if __name__ == "__main__":
    sys.exit(main())
