#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""Production parameter-space figures for the BC5 LHCb-downstream analysis.

Each panel is a 2-D histogram of parent particles in the (theta, p) plane:

    x = theta = angle of the parent with respect to the beam axis [rad], log 1e-5 .. 1
    y = p     = parent momentum [GeV],                                 log 1 .. 1e4
    colour    = dsigma per 2-D bin [pb/bin], logarithmic colour scale

The style follows the reference figure of the FORESEE paper (Kling & Trojanowski,
arXiv:2105.07077, fig. 1, right panel: "B^0 at 14 TeV LHC [pb/bin], Pythia 8"):
20 bins per decade on both axes, cyclic 'hsv' colour map, a dashed p_T = m_X line and
vertical dotted lines marking the detector angular window.  The reference is 14 TeV with
the FASER / FASER 2 windows; OURS IS 13 TeV WITH THE LHCb DOWNSTREAM WINDOW 2 < eta < 5,
i.e. theta(eta) = 2 atan(exp(-eta)), theta(eta=2) = 0.269036 rad and
theta(eta=5) = 0.0134757 rad.

Normalisation:  dsigma_bin = sigma_ref * N_bin / n_gen  [pb per bin], with sigma_ref and
n_gen read from work/magnus/samples/normalisation.json.  Charge-conjugate pairs are
combined into one species.

Outputs (all under outputs/figures/production/):
    production_Bplus.{png,pdf}   production_B0.{png,pdf}   production_Bs.{png,pdf}
    production_Kplus.{png,pdf}   production_K0L.{png,pdf}  production_Higgs.{png,pdf}
    production_all_particles.png   (2 x 3 grid of the same six panels)

Usage (from the run directory, E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3):
    D:\py\python.exe work\analysis\plot_production.py

Deterministic: fixed bin edges, colour ranges derived from the data (one shared scale per
parent sample: B mesons, K mesons, Higgs), no RNG.  Does not read or write any other
analysis product.
"""
from __future__ import annotations

import json
import os
import time

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator, NullFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

# --------------------------------------------------------------------------------------
# paths
# --------------------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, "..", ".."))
SAMPLES = os.path.join(RUN, "work", "magnus", "samples")
OUTDIR = os.path.join(RUN, "outputs", "figures", "production")
NORM_JSON = os.path.join(SAMPLES, "normalisation.json")

# --------------------------------------------------------------------------------------
# binning: 20 bins per decade on both axes, as in the reference figure
# --------------------------------------------------------------------------------------
N_PER_DECADE = 20
THETA_LO, THETA_HI = 1e-5, 1.0
P_LO, P_HI = 1.0, 1.0e4
THETA_EDGES = np.logspace(np.log10(THETA_LO), np.log10(THETA_HI),
                          5 * N_PER_DECADE + 1)          # 101 edges -> 100 bins
P_EDGES = np.logspace(np.log10(P_LO), np.log10(P_HI),
                      4 * N_PER_DECADE + 1)              # 81 edges  -> 80 bins

# LHCb downstream angular window, 2 < eta_S < 5
ETA_WINDOW = (2.0, 5.0)
THETA_WINDOW = 2.0 * np.arctan(np.exp(-np.asarray(ETA_WINDOW)))   # [0.269036, 0.0134757]

# reference colour map: cyclic 'hsv'; empty bins are white, as in the reference
CMAP = plt.get_cmap("hsv").copy()
CMAP.set_bad("white")

# --------------------------------------------------------------------------------------
# species: one entry per panel
# --------------------------------------------------------------------------------------
SPECIES = [
    dict(key="Bplus", label=r"$B^+$",      sample="bbar_heavy", group="B",
         pdg=(521, -521),   mass=5.27934),
    dict(key="B0",    label=r"$B^0$",      sample="bbar_heavy", group="B",
         pdg=(511, -511),   mass=5.27965),
    dict(key="Bs",    label=r"$B_s^0$",    sample="bbar_heavy", group="B",
         pdg=(531, -531),   mass=5.36688),
    dict(key="Kplus", label=r"$K^+$",      sample="minbias_light", group="K",
         pdg=(321, -321),   mass=0.493677),
    dict(key="K0L",   label=r"$K^0_L$",    sample="minbias_light", group="K",
         pdg=(311, -311),   mass=0.497611),
    dict(key="Higgs", label=r"$h$",        sample="higgs", group="H",
         pdg=None,          mass=125.25),
]

# samples that must be read (once) and the two normalisation groups they belong to
SAMPLE_FILES = {
    "bbar_heavy": ("B", "pp -> b bbar parent sample"),
    "minbias_light": ("K", "inelastic pp (SoftQCD) sample"),
    "higgs": ("H", "SM Higgs sample (several copies per event)"),
}


# --------------------------------------------------------------------------------------
# loading and kinematics
# --------------------------------------------------------------------------------------
def load_sample(name):
    """Read one parent sample (plain .txt or .txt.gz); keep one row per physical parent.

    Two things changed in v3.6 and are handled here:
      * the samples are gzipped (the B file alone is 1.5 M rows), so both extensions are tried;
      * the Higgs file carries an 11th column with the Pythia STATUS and holds EVERY copy of the
        resonance in the event record (~7 per event).  Only the pre-decay copy, |status| = 62,
        is a physical Higgs; without that cut the Higgs panel would be ~6x too high.  The B and
        K samples have 10 columns and no such cut.
    """
    path = None
    for ext in (".txt", ".txt.gz"):
        cand = os.path.join(SAMPLES, name + ext)
        if os.path.isfile(cand):
            path = cand
            break
    if path is None:
        raise FileNotFoundError(os.path.join(SAMPLES, name + ".txt[.gz]"))
    t0 = time.time()
    raw = pd.read_csv(path, sep=r"\s+", comment="#", header=None)
    nrows_in = len(raw)
    if raw.shape[1] >= 11:
        raw = raw[np.abs(raw[10].to_numpy()) == 62]
    df = pd.DataFrame({
        "event": raw[0].to_numpy(np.int64), "pdg": raw[1].to_numpy(np.int32),
        "px": raw[2].to_numpy(np.float64), "py": raw[3].to_numpy(np.float64),
        "pz": raw[4].to_numpy(np.float64), "E": raw[5].to_numpy(np.float64)})
    print("  read %-22s %8d of %8d rows  (%5.1f MB)  %.2f s"
          % (os.path.basename(path), len(df), nrows_in,
             os.path.getsize(path) / 1e6, time.time() - t0))
    return df


def kinematics(df):
    """p = |p|, theta = angle w.r.t. the beam axis (0, pi/2), m from E and p."""
    px = df["px"].to_numpy()
    py = df["py"].to_numpy()
    pz = df["pz"].to_numpy()
    e = df["E"].to_numpy()
    pt = np.hypot(px, py)
    p = np.sqrt(pt * pt + pz * pz)
    theta = np.arctan2(pt, np.abs(pz))       # |pz|: fold the two beam directions
    m = np.sqrt(np.maximum(e * e - p * p, 0.0))
    return p, theta, m


def load_all():
    """Read each sample once and compute the kinematics on the filtered columns."""
    data = {}
    for name in SAMPLE_FILES:
        df = load_sample(name)
        p, theta, m = kinematics(df)
        data[name] = dict(df=df, p=p, theta=theta, m=m)
        print("    distinct events %d" % df["event"].nunique())
    return data


# --------------------------------------------------------------------------------------
# sanity checks on the samples
# --------------------------------------------------------------------------------------
def sanity_checks(data, norm):
    print("\n" + "=" * 78)
    print("SANITY CHECKS")
    print("=" * 78)

    for name, d in data.items():
        df = d["df"]
        vc = pd.Series(df["pdg"]).value_counts()
        print("\n%s: %d rows, %d distinct event ids"
              % (name + ".txt", len(df), df["event"].nunique()))
        print("  PDG codes present (code: rows):")
        for code, n in vc.items():
            print("    %5d : %7d" % (code, n))
        print("  pz < 0 fraction: %.4f" % np.mean(df["pz"].to_numpy() < 0.0))

    # ---- hand checks: p from the momentum components vs E -------------------------
    print("\nhand check A: B+ candidate (first row of bbar_heavy.txt)")
    bb = data["bbar_heavy"]
    row = bb["df"].iloc[0]
    px, py, pz, e = (float(row["px"]), float(row["py"]),
                     float(row["pz"]), float(row["E"]))
    p = np.sqrt(px * px + py * py + pz * pz)
    print("  px=%.7f py=%.7f pz=%.7f E=%.7f GeV" % (px, py, pz, e))
    print("  p = sqrt(px^2+py^2+pz^2)              = %.7f GeV" % p)
    print("  sqrt(p^2 + m_B+^2) (m_B+ = 5.27934)   = %.7f GeV" % np.sqrt(p * p + 5.27934 ** 2))
    print("  E_file - sqrt(p^2+m^2)                = %.3e GeV" % (e - np.sqrt(p * p + 5.27934 ** 2)))
    print("  E^2 - p^2 = %.5f GeV^2   (m_B+^2 = %.5f GeV^2)"
          % (e * e - p * p, 5.27934 ** 2))

    print("\nhand check B: a B+ candidate with negative pz (theta uses |pz|)")
    bbdf = bb["df"]
    sel = np.where(bbdf["pdg"].to_numpy() == 521)[0]
    neg = sel[bbdf["pz"].to_numpy()[sel] < 0.0][0]
    px, py, pz, e = (float(bbdf["px"].iloc[neg]), float(bbdf["py"].iloc[neg]),
                     float(bbdf["pz"].iloc[neg]), float(bbdf["E"].iloc[neg]))
    pt = np.hypot(px, py)
    p = np.sqrt(pt * pt + pz * pz)
    print("  event %d: px=%.6f py=%.6f pz=%.6f E=%.6f GeV"
          % (bbdf["event"].iloc[neg], px, py, pz, e))
    print("  pT=%.6f GeV  p=%.6f GeV  theta=atan(pT/|pz|)=%.6f rad (in (0, pi/2): %s)"
          % (pt, p, np.arctan2(pt, abs(pz)), 0.0 < np.arctan2(pt, abs(pz)) < np.pi / 2))

    # ---- mass consistency ---------------------------------------------------------
    print("\nmass consistency  sqrt(E^2-p^2) [GeV]:")
    for name, pdgs, mnom in (("bbar_heavy", (521, -521), 5.27934),
                             ("bbar_heavy", (511, -511), 5.27965),
                             ("bbar_heavy", (531, -531), 5.36688),
                             ("minbias_light", (321, -321), 0.493677),
                             ("minbias_light", (311, -311), 0.497611)):
        d = data[name]
        mask = np.isin(d["df"]["pdg"].to_numpy(), np.asarray(pdgs))
        mm = d["m"][mask]
        print("  %-12s pdg %s  nominal %.4f  median %.6f  0.1-99.9%% [%.6f, %.6f]"
              % (name, str(pdgs), mnom, np.median(mm),
                 np.percentile(mm, 0.1), np.percentile(mm, 99.9)))
    h = data["higgs"]
    ev = h["df"]["event"].to_numpy()
    print("  higgs        %d rows, %d distinct events; several copies per event,"
          " only the FIRST row per event is kept" % (len(ev), len(np.unique(ev))))
    u, first = np.unique(ev, return_index=True)
    print("    unique events %d  (event ids 0 .. %d)  rows kept %d"
          % (len(u), ev.max(), len(first)))
    mh = h["m"][first]
    print("    mass of the first-row Higgs: median %.4f GeV, 1-99%% [%.4f, %.4f] GeV"
          % (np.median(mh), np.percentile(mh, 1), np.percentile(mh, 99)))
    print("    -> the Higgs mass is generated over a finite range (not a delta);"
          " this does not affect the (theta, p) spectra.")

    # ---- kinematic coverage -------------------------------------------------------
    print("\nkinematic coverage (all rows of each sample, except Higgs: first row per event):")
    for name, d in data.items():
        if name == "higgs":
            ev = d["df"]["event"].to_numpy()
            _, first = np.unique(ev, return_index=True)
            th = d["theta"][first]
            p = d["p"][first]
        else:
            th = d["theta"]
            p = d["p"]
        eta = -np.log(np.tan(th / 2.0))
        print("  %-14s p [%.4g, %.4g] GeV   theta [%.4g, %.4g] rad   eta(|pz|) [%.4g, %.4g]"
              % (name, p.min(), p.max(), th.min(), th.max(),
                 np.nanmin(eta), np.nanmax(eta)))
        assert np.all(th > 0.0) and np.all(th <= np.pi / 2.0 + 1e-12)

    # ---- window edges ------------------------------------------------------------
    print("\nLHCb downstream window: theta(eta) = 2 atan(exp(-eta))")
    for eta, th in zip(ETA_WINDOW, THETA_WINDOW):
        print("  eta = %.1f  ->  theta = %.7f rad" % (eta, th))


# --------------------------------------------------------------------------------------
# histogram + normalisation + closure
# --------------------------------------------------------------------------------------
def build_species(data, sp, norm):
    """Select the charge-combined species, histogram it and normalise."""
    d = data[sp["sample"]]
    if sp["group"] == "H":
        ev = d["df"]["event"].to_numpy()
        _, first = np.unique(ev, return_index=True)     # first row of each event id
        mask = np.zeros(len(ev), dtype=bool)
        mask[first] = True
        print("  %-5s: kept first row of %d distinct events (of %d raw rows)"
              % (sp["key"], mask.sum(), len(ev)))
    else:
        mask = np.isin(d["df"]["pdg"].to_numpy(), np.asarray(sp["pdg"]))
        print("  %-5s: pdg %s -> %d candidates in the sample"
              % (sp["key"], sp["pdg"], mask.sum()))

    p = d["p"][mask]
    theta = d["theta"][mask]

    sigma = float(norm[sp["group"]]["sigma_pb"])
    n_gen = int(norm[sp["group"]]["n_gen"])
    n_species = int(mask.sum())

    H, _, _ = np.histogram2d(theta, p, bins=[THETA_EDGES, P_EDGES])
    n_cand = int(round(H.sum()))
    dsigma = sigma * H / n_gen                      # pb/bin, shape (100, 80)

    sigma_total = sigma * n_species / n_gen
    sigma_exp = sigma * n_cand / n_gen
    sum_bins = float(dsigma.sum())

    sp.update(dict(H=H, dsigma=dsigma, n_species=n_species, n_cand=n_cand,
                   n_gen=n_gen, sigma_ref=sigma,
                   sigma_total=sigma_total, sigma_exp=sigma_exp,
                   sum_bins=sum_bins,
                   rel_dev=abs(sum_bins - sigma_exp) / sigma_exp if sigma_exp else 0.0,
                   title=sp["label"] + r" at 13 TeV LHC [pb/bin]"))
    return sp


def print_closure(species):
    print("\n" + "=" * 100)
    print("CLOSURE CHECK   dsigma_bin = sigma_ref * N_bin / n_gen ;  N_cand = candidates")
    print("inside the plotted window (1e-5 <= theta <= 1 rad, 1 <= p <= 1e4 GeV)")
    print("=" * 100)
    hdr = ("%-6s %-12s %10s %10s %9s %13s %15s %15s %10s"
           % ("spec", "PDG", "N_species", "N_cand", "n_gen", "sigma_ref[pb]",
              "sigma_exp[pb]", "sum_bins[pb]", "rel.dev"))
    print(hdr)
    print("-" * len(hdr))
    for sp in species:
        pdg_lab = str(sp["pdg"]) if sp["pdg"] is not None else "901-909"
        print("%-6s %-12s %10d %10d %9d %13.4e %15.6e %15.6e %10.2e"
              % (sp["key"], pdg_lab, sp["n_species"], sp["n_cand"],
                 sp["n_gen"], sp["sigma_ref"], sp["sigma_exp"], sp["sum_bins"],
                 sp["rel_dev"]))
    mx = max(sp["rel_dev"] for sp in species)
    print("-" * len(hdr))
    print("max |sum_bins - sigma_expected| / sigma_expected = %.3e" % mx)
    print("tolerance 1e-2 -> %s (the binned and analytic normalisations are the same"
          " sum, so they agree to floating point)" % ("PASS" if mx < 1e-2 else "FAIL"))
    print("\nfull-sample normalisation (all candidates of the species, not only the"
          " plotted window):")
    for sp in species:
        frac = sp["n_cand"] / sp["n_species"] if sp["n_species"] else 0.0
        print("  %-6s N_species %8d  sigma_ref*N_species/n_gen = %12.6e pb"
              "   (%.2f%% inside the plotted window)"
              % (sp["key"], sp["n_species"], sp["sigma_total"], 100.0 * frac))
    print("=" * 100)


# --------------------------------------------------------------------------------------
# drawing
# --------------------------------------------------------------------------------------
def plot_panel(ax, cax, sp, fs=1.0):
    """Draw one (theta, p) production panel into ax, with its colour bar in cax."""
    H = sp["H"]
    dsigma = sp["dsigma"]
    vmin, vmax = sp["vmin"], sp["vmax"]

    masked = np.ma.masked_where(H <= 0, dsigma)
    mesh = ax.pcolormesh(THETA_EDGES, P_EDGES, masked.T, cmap=CMAP,
                         norm=LogNorm(vmin=vmin, vmax=vmax),
                         shading="flat", rasterized=True)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(THETA_LO, THETA_HI)
    ax.set_ylim(P_LO, P_HI)
    ax.set_xlabel(r"angle wrt. beam axis $\theta$ [rad]", fontsize=13.0 * fs)
    ax.set_ylabel(r"momentum $p$ [GeV]", fontsize=13.0 * fs)
    ax.tick_params(which="major", labelsize=11.2 * fs, length=4.5 * fs, width=0.8)
    ax.tick_params(which="minor", length=2.3 * fs, width=0.7)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_minor_formatter(NullFormatter())

    # ---- diagonal dashed line: p_T = p sin(theta) = m_X ---------------------------
    m_x = sp["mass"]
    th_line = np.logspace(np.log10(THETA_LO), np.log10(THETA_HI), 600)
    p_line = m_x / np.sin(th_line)
    vis = (p_line >= P_LO) & (p_line <= P_HI)
    ax.plot(th_line[vis], p_line[vis], color="black", lw=1.35,
            ls=(0, (5.5, 2.4)), zorder=6, solid_capstyle="butt")

    # label, drawn on the line and rotated along it (contour-label style)
    t_lo, t_hi = th_line[vis].min(), th_line[vis].max()
    t_lab = t_lo * (t_hi / t_lo) ** 0.42
    t1, t2 = t_lab * 0.85, t_lab * 1.18
    xy1 = ax.transData.transform((t1, m_x / np.sin(t1)))
    xy2 = ax.transData.transform((t2, m_x / np.sin(t2)))
    angle = np.degrees(np.arctan2(xy2[1] - xy1[1], xy2[0] - xy1[0]))
    p_lab = m_x / np.sin(t_lab)
    ax.text(t_lab, p_lab, r"$p_T = m_X$", rotation=angle, rotation_mode="anchor",
            ha="center", va="center", fontsize=9.6 * fs, color="black", zorder=7,
            bbox=dict(boxstyle="square,pad=0.18", facecolor="white",
                      edgecolor="none", alpha=0.85))

    # ---- vertical dotted lines: LHCb downstream 2 < eta < 5 -----------------------
    for th_w, lab in zip(THETA_WINDOW, [r"$\eta = 2$", r"$\eta = 5$"]):
        ax.axvline(th_w, color="black", lw=1.35, ls=(0, (1.5, 2.475)), zorder=6)
        ax.text(th_w, 1.45, lab, rotation=90, ha="center", va="bottom",
                fontsize=9.0 * fs, color="black", zorder=7,
                bbox=dict(boxstyle="square,pad=0.12", facecolor="white",
                          edgecolor="none", alpha=0.80))
    th_mid = float(np.sqrt(THETA_WINDOW[0] * THETA_WINDOW[1]))
    ax.text(th_mid, 1.45, r"LHCb downstream $2<\eta<5$", rotation=90,
            ha="center", va="bottom", fontsize=8.2 * fs, color="black", zorder=7,
            bbox=dict(boxstyle="square,pad=0.12", facecolor="white",
                      edgecolor="none", alpha=0.80))

    # ---- titles, in the style of the reference (right-aligned, inside the axes) ---
    ax.text(0.985, 0.978, sp["title"], transform=ax.transAxes, ha="right",
            va="top", fontsize=11.7 * fs, color="black", zorder=8,
            bbox=dict(boxstyle="square,pad=0.15", facecolor="white",
                      edgecolor="none", alpha=0.95))
    ax.text(0.985, 0.930, "Pythia 8.316", transform=ax.transAxes, ha="right",
            va="top", fontsize=9.9 * fs, color="black", zorder=8)
    ax.text(0.985, 0.888, r"$N_{\rm cand}$ = " + format(sp["n_cand"], ","),
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8.8 * fs, color="black", zorder=8,
            bbox=dict(boxstyle="square,pad=0.12", facecolor="white",
                      edgecolor="none", alpha=0.95))

    # ---- colour bar ---------------------------------------------------------------
    cb = plt.colorbar(mesh, cax=cax)
    cb.set_label("pb/bin", fontsize=10.6 * fs)
    cb.ax.tick_params(labelsize=9.8 * fs, length=3.5 * fs, width=0.8)
    cb.ax.yaxis.set_minor_locator(LogLocator(base=10.0,
                                             subs=np.arange(2.0, 10.0) * 0.1,
                                             numticks=100))
    cb.ax.yaxis.set_minor_formatter(NullFormatter())
    cb.ax.tick_params(which="minor", length=1.8 * fs, width=0.6)
    return mesh


def colour_range(nz):
    """Round the colour scale outwards to whole decades (deterministic)."""
    vmin = 10.0 ** np.floor(np.log10(nz.min()))
    vmax = 10.0 ** np.ceil(np.log10(nz.max()))
    if vmax / vmin < 100.0:            # keep at least two decades of dynamic range
        vmax = vmin * 100.0
    return float(vmin), float(vmax)


def plot_single(sp):
    fig = plt.figure(figsize=(8.0, 5.5))
    ax = fig.add_axes([0.105, 0.110, 0.755, 0.855])
    cax = fig.add_axes([0.905, 0.110, 0.030, 0.855])
    plot_panel(ax, cax, sp, fs=1.0)
    base = os.path.join(OUTDIR, "production_" + sp["key"])
    fig.savefig(base + ".png", dpi=150)
    fig.savefig(base + ".pdf")
    plt.close(fig)
    print("  wrote %s.png / .pdf" % base)


def plot_grid(species):
    fig = plt.figure(figsize=(19.5, 9.0))
    gs = fig.add_gridspec(2, 3, left=0.052, right=0.987, bottom=0.055,
                          top=0.978, wspace=0.30, hspace=0.26)
    for i, sp in enumerate(species):
        ax = fig.add_subplot(gs[i // 3, i % 3])
        cax = make_axes_locatable(ax).append_axes("right", size="4.2%", pad=0.10)
        plot_panel(ax, cax, sp, fs=0.80)
    out = os.path.join(OUTDIR, "production_all_particles.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("  wrote %s" % out)


# --------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------
def main():
    t_start = time.time()
    os.makedirs(OUTDIR, exist_ok=True)
    print("=" * 78)
    print("production parameter-space figures  (13 TeV, LHCb downstream 2 < eta < 5)")
    print("output directory: %s" % OUTDIR)
    print("=" * 78)

    with open(NORM_JSON) as f:
        norm = json.load(f)
    print("\nnormalisation from %s:" % NORM_JSON)
    for g in ("B", "K", "H"):
        print("  %s: sigma_ref = %.6g pb   n_gen = %d   file = %s"
              % (g, norm[g]["sigma_pb"], norm[g]["n_gen"], norm[g]["file"]))

    print("\nreading samples ...")
    data = load_all()

    sanity_checks(data, norm)

    print("\n" + "=" * 78)
    print("HISTOGRAMS  (%d x %d bins, %.0f bins per decade)"
          % (len(THETA_EDGES) - 1, len(P_EDGES) - 1, N_PER_DECADE))
    print("=" * 78)
    species = []
    for sp in SPECIES:
        sp = dict(sp)
        sp = build_species(data, sp, norm)
        print("       sigma_ref * N_cand / n_gen = %.6e pb" % sp["sigma_exp"])
        species.append(sp)

    # one colour scale per parent sample, so the panels of a family are directly comparable
    for grp, grp_name in (("B", "B mesons"), ("K", "K mesons"), ("H", "Higgs")):
        members = [sp for sp in species if sp["group"] == grp]
        nz = np.concatenate([sp["dsigma"][sp["dsigma"] > 0.0].ravel() for sp in members])
        vmin, vmax = colour_range(nz)
        for sp in members:
            sp["vmin"], sp["vmax"] = vmin, vmax
        print("  colour scale for the %-8s panels: %.0e .. %.0e pb/bin (shared)"
              % (grp_name, vmin, vmax))

    print_closure(species)

    print("\n" + "=" * 78)
    print("FIGURES")
    print("=" * 78)
    for sp in species:
        plot_single(sp)
    plot_grid(species)

    # keep only what is useful for the README
    print("\nDONE in %.1f s" % (time.time() - t_start))
    print("regeneration command (from the run directory):"
          "  D:\\py\\python.exe work\\analysis\\plot_production.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
