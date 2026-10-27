#!/usr/bin/env python
"""
spectrum_foresee.py -- our pp-collision particle spectra vs FORESEE's published
forward spectra (files/hadrons/13.6TeV/Pythia8/*.txt).

Inputs (all in the workspace)
-----------------------------
Our samples (normalisation: work/magnus/samples/normalisation.json):
  * B   : work/magnus/samples/bbar_heavy.txt.gz
          Pythia 8.316, pp 13 TeV, HardQCD gg/qqbar -> b bbar, 800000 events,
          sigma(bbbar) = 560 ub = 5.6e8 pb (BC5 v3.6 value); rows are
          B+/B0/Bs/Bc fragmentation hadrons.
  * inel: work/magnus/samples/minbias_pik.txt.gz
          Pythia 8.316, pp 13 TeV, SoftQCD:inelastic, 500000 events,
          sigma_inel = 79.5 mb = 7.95e10 pb; rows are pi+-,K+-,K0/K0bar
          fragmentation hadrons at status 91/92.
Columns: event pdg px py pz E x y z t  [px..E in GeV, x..t in mm].

FORESEE (arXiv:2105.07077, https://github.com/KlingFelix/FORESEE):
  files/hadrons/<E>TeV/<gen>/<gen>_<E>TeV_<pid>.txt, three columns
  (log10 theta/rad, log10 p/GeV, cross section [pb/bin]), 20 bins/decade,
  forward hemisphere only.  The '# Normalization:' header equals the sum of the
  bin weights, i.e. the forward-hemisphere cross section of that species.
  Main reference: 13.6 TeV Pythia8 (Monash tune).  Cross-checks: the NLO-P8
  beauty family (used by the authors' own 2HDM_H model for B -> X_s H) with its
  Max/Min variants, and the EPOSLHC / SIBYLL 13.6 TeV light-hadron files.

Outputs
-------
  outputs/figures/production/spectrum_vs_foresee_overview.png / .pdf
  (per-particle figures in the same directory are made by spectrum_particle_figures.py)
  a quantitative report on stdout (quoted in logs/spectrum_vs_foresee.md)

Missing FORESEE files are downloaded from the GitHub raw URL
https://raw.githubusercontent.com/KlingFelix/FORESEE/main/files/hadrons

Usage
-----
    python work/analysis/spectrum_foresee.py
    python work/analysis/spectrum_foresee.py --quick       # 100k-event pi/K cache
    python work/analysis/spectrum_foresee.py --no-download
"""

from __future__ import annotations

import argparse
import gzip
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm

# ----------------------------------------------------------------------------
# paths / constants
# ----------------------------------------------------------------------------
HERE = Path(__file__).resolve()
RUN = HERE.parents[2]                 # runs/<run>/
REPO = HERE.parents[4]                # workspace root
FORESEE = REPO / "references" / "tools" / "FORESEE"
FOR_HADRONS = FORESEE / "files" / "hadrons"
FIG_DIR = RUN / "outputs" / "figures" / "production"
B_SAMPLE = RUN / "work" / "magnus" / "samples" / "bbar_heavy.txt.gz"
LIGHT_SAMPLE = RUN / "work" / "magnus" / "samples" / "minbias_pik.txt.gz"
LIGHT_CACHE = (RUN / "work" / "analysis" / "cache" /
               "sample_minbias_pik_525049850_evsub100000_seed20261025.npz")
FORESEE_RAW = "https://raw.githubusercontent.com/KlingFelix/FORESEE/main/files/hadrons"

SIGMA_BB_PB = 5.6e8         # pp -> b bbar at 13 TeV (BC5 v3.6 value)
N_GEN_B = 800_000
SIGMA_INEL_PB = 7.95e10     # 79.5 mb measured inelastic pp cross section
N_GEN_INEL = 500_000

B_PIDS = [511, -511, 521, -521, 531, -531, 541, -541]
B4_PIDS = [511, -511, 521, -521]          # B0/B0bar/B+/B- = used by the analysis
LIGHT_PIDS = [211, -211, 321, -321]

# FORESEE binning (log10 units, 20 bins/decade)
TE = np.arange(-7.3, 0.2001, 0.05)        # 151 edges -> 150 theta bins
PE_B = np.arange(-2.0, 5.0001, 0.05)      # 141 edges -> 140 p bins (B)
PE_L = np.arange(-2.0, 4.0001, 0.05)      # 121 edges -> 120 p bins (light)

# the NLO-P8 beauty family covers theta in [1e-5, 1] rad, p in [1,1e4] GeV,
# i.e. on the Pythia8 grid theta bins 46..145 and p bins 40..119
NLO_SLICE = (slice(46, 146), slice(40, 120))

DATASETS = {
    ("13.6", "Pythia8"): B_PIDS,
    ("14", "Pythia8"): B_PIDS,
    ("13.6", "NLO-P8"): B_PIDS,
    ("14", "NLO-P8"): B_PIDS,
    ("13.6", "NLO-P8-Max"): B4_PIDS,
    ("13.6", "NLO-P8-Min"): B4_PIDS,
}
LIGHT_GENS = ("Pythia8", "EPOSLHC", "SIBYLL")

COARSE = 5            # merge 5x5 FORESEE bins (0.25 dex) for maps/ratios


# ----------------------------------------------------------------------------
# FORESEE files
# ----------------------------------------------------------------------------
def fors_path(energy, generator, pid):
    return (FOR_HADRONS / f"{energy}TeV" / generator /
            f"{generator}_{energy}TeV_{pid}.txt")


def fors_grid(generator, pid):
    """(n_theta, n_p) of a FORESEE spectrum file."""
    if generator.startswith("Pythia8"):
        return (150, 140) if abs(pid) in (511, 521, 531, 541, 5122) else (150, 120)
    if generator in ("EPOSLHC", "SIBYLL", "QGSJET", "DPMJET"):
        return (150, 120)
    return (100, 80)          # NLO-P8 and its Max/Min variants


def ensure_foresee_files(allow_download=True):
    needed = [(e, g, p) for (e, g), pids in DATASETS.items() for p in pids]
    needed += [("13.6", "Pythia8", p) for p in LIGHT_PIDS]
    needed += [("13.6", g, p) for g in ("EPOSLHC", "SIBYLL") for p in LIGHT_PIDS]
    missing = [item for item in needed if not fors_path(*item).exists()]
    if not missing:
        return
    if not allow_download:
        raise SystemExit("missing FORESEE files:\n  " +
                         "\n  ".join(str(fors_path(*m)) for m in missing))
    for energy, generator, pid in missing:
        path = fors_path(energy, generator, pid)
        path.parent.mkdir(parents=True, exist_ok=True)
        url = f"{FORESEE_RAW}/{energy}TeV/{generator}/{generator}_{energy}TeV_{pid}.txt"
        print(f"[download] {url}", flush=True)
        for attempt in range(4):
            try:
                urllib.request.urlretrieve(url, path)
                if path.stat().st_size > 1000:
                    break
            except Exception as exc:              # noqa: BLE001 - network
                if attempt == 3:
                    raise SystemExit(f"cannot download {url}: {exc}")
                time.sleep(2.0)


def read_forsee(energy, generator, pid):
    """Return (weights[n_theta, n_p], per-particle weight quantum eps).

    The files store rounded counts x (sigma/n_gen); the smallest positive weight
    is the per-particle quantum eps, so a bin of weight w holds ~w/eps entries
    and has Poisson variance eps*w.
    """
    data = np.loadtxt(fors_path(energy, generator, pid), comments="#")
    nth, npb = fors_grid(generator, pid)
    weights = data[:, 2].reshape(nth, npb)
    positive = weights[weights > 0]
    eps = float(positive.min()) if positive.size else 0.0
    return weights, eps


def forsee_species(energy, generator, pids):
    """{pid: (weights, Poisson variance)}."""
    out = {}
    for pid in pids:
        w, eps = read_forsee(energy, generator, pid)
        out[pid] = (w, eps * w)
    return out


def stack(species, pids):
    """Sum weights and variances over pids."""
    w = sum(species[p][0] for p in pids)
    v = sum(species[p][1] for p in pids)
    return np.asarray(w, float), np.asarray(v, float)


# ----------------------------------------------------------------------------
# our samples
# ----------------------------------------------------------------------------
_B_CACHE = None


def _log_theta_p(px, py, pz):
    pt = np.hypot(px, py)
    pmom = np.sqrt(px * px + py * py + pz * pz)
    return np.log10(np.arctan2(pt, pz)), np.log10(pmom)


def _b_rows():
    """Forward B rows loaded once: (pdg, log10 theta, log10 p)."""
    global _B_CACHE
    if _B_CACHE is None:
        with gzip.open(B_SAMPLE, "rt") as handle:
            data = np.loadtxt(handle, usecols=(1, 2, 3, 4))
        pdg = data[:, 0].astype(np.int32)
        px, py, pz = data[:, 1], data[:, 2], data[:, 3]
        forward = pz > 0
        log_t, log_p = _log_theta_p(px[forward], py[forward], pz[forward])
        _B_CACHE = (pdg[forward], log_t, log_p)
    return _B_CACHE


def our_b_counts(pids):
    """B counts on the FORESEE B grid (forward hemisphere)."""
    pdg, log_t, log_p = _b_rows()
    mask = np.isin(pdg, np.asarray(pids))
    counts, _, _ = np.histogram2d(log_t[mask], log_p[mask], bins=(TE, PE_B))
    return counts


def our_light_counts(use_cache=False):
    """{pid: counts} for pi/K on the FORESEE light grid (forward hemisphere)."""
    counts = {pid: np.zeros((len(TE) - 1, len(PE_L) - 1)) for pid in LIGHT_PIDS}
    if use_cache and LIGHT_CACHE.exists():
        data = np.load(LIGHT_CACHE)
        n_sub = len(np.unique(data["ev"]))
        fwd = data["pz"] > 0
        log_t, log_p = _log_theta_p(data["px"].astype(np.float64)[fwd],
                                    data["py"].astype(np.float64)[fwd],
                                    data["pz"].astype(np.float64)[fwd])
        q = data["pdg"][fwd]
        scale = N_GEN_INEL / n_sub        # counts -> full-sample equivalents
        for pid in LIGHT_PIDS:
            mask = q == pid
            hist, _, _ = np.histogram2d(log_t[mask], log_p[mask], bins=(TE, PE_L))
            counts[pid] = hist * scale
        return counts

    import pandas as pd

    reader = pd.read_csv(
        LIGHT_SAMPLE, sep=r"\s+", header=None, comment="#",
        usecols=[0, 1, 2, 3, 4], names=["ev", "pdg", "px", "py", "pz"],
        dtype={"ev": "int32", "pdg": "int16",
               "px": "float32", "py": "float32", "pz": "float32"},
        chunksize=4_000_000,
    )
    for chunk in reader:
        fwd = chunk["pz"].to_numpy(np.float64) > 0
        log_t, log_p = _log_theta_p(chunk["px"].to_numpy(np.float64)[fwd],
                                    chunk["py"].to_numpy(np.float64)[fwd],
                                    chunk["pz"].to_numpy(np.float64)[fwd])
        q = chunk["pdg"].to_numpy()[fwd]
        for pid in LIGHT_PIDS:
            mask = q == pid
            hist, _, _ = np.histogram2d(log_t[mask], log_p[mask], bins=(TE, PE_L))
            counts[pid] += hist
    return counts


# ----------------------------------------------------------------------------
# comparison helpers
# ----------------------------------------------------------------------------
def tot(hist, slc=None):
    return float(hist[slc].sum() if slc is not None else hist.sum())


def eta_of_rows(n_theta):
    """eta at the centre of each theta row (full or coarsened FORESEE grid)."""
    if n_theta == len(TE) - 1:
        centres = 0.5 * (TE[:-1] + TE[1:])
    elif n_theta == 100:                       # NLO-P8 grid
        centres = -4.975 + 0.05 * np.arange(100)
    else:                                      # coarsened uniform grid
        edges = TE[0] + (TE[-1] - TE[0]) / n_theta * np.arange(n_theta + 1)
        centres = 0.5 * (edges[:-1] + edges[1:])
    return -np.log(np.tan(0.5 * 10.0 ** centres))


def p_centres(n_p):
    """log10 p bin centres for full or coarsened FORESEE p grids."""
    if n_p == len(PE_B) - 1:
        return 0.5 * (PE_B[:-1] + PE_B[1:])
    if n_p == len(PE_L) - 1:
        return 0.5 * (PE_L[:-1] + PE_L[1:])
    if n_p == 80:                              # NLO-P8 grid
        return 0.025 + 0.05 * np.arange(n_p)
    plo, phi = (-2.0, 5.0) if n_p == len(PE_B) // COARSE - 1 or n_p == 28 else (-2.0, 4.0)
    edges = plo + (phi - plo) / n_p * np.arange(n_p + 1)
    return 0.5 * (edges[:-1] + edges[1:])


def eta_fractions(hist):
    """Fractions of the forward cross section in eta bands."""
    eta = eta_of_rows(hist.shape[0])
    total = hist.sum()
    bands = {"eta<1.5": eta < 1.5,
             "1.5<eta<5.5": (eta >= 1.5) & (eta < 5.5),
             "eta>5.5": eta >= 5.5}
    return {name: float(hist[mask, :].sum() / total) for name, mask in bands.items()}


def shape_mask(n_ours, f_hist, nmin=100, eta_range=None, logp_range=None):
    """Boolean mask of bins usable for a shape comparison."""
    mask = (n_ours >= nmin) & (f_hist > 0)
    if eta_range is not None:
        eta = eta_of_rows(n_ours.shape[0])
        mask &= ((eta >= eta_range[0]) & (eta < eta_range[1]))[:, None]
    if logp_range is not None:
        pc = p_centres(n_ours.shape[1])
        mask &= ((pc >= logp_range[0]) & (pc < logp_range[1]))[None, :]
    return mask


def ratio_metrics(n_ours, w_row, f_hist, f_var, slc=None, nmin=100,
                  eta_range=None, logp_range=None):
    """Bin-by-bin ratio of our spectrum to a FORESEE spectrum."""
    if slc is not None and n_ours.shape != f_hist.shape:
        n_o = n_ours[slc]
    else:
        n_o = n_ours
    ours = n_o.astype(float) * w_row
    normalisation = float(ours.sum() / f_hist.sum())
    mask = shape_mask(n_o, f_hist, nmin, eta_range, logp_range)
    if mask.sum() < 10:
        return None
    r = ours[mask] / f_hist[mask]
    relvar = 1.0 / n_o[mask] + f_var[mask] / f_hist[mask] ** 2
    shape = r / normalisation - 1.0
    return dict(n_bins=int(mask.sum()), total_ratio=normalisation,
                median=float(np.median(r)),
                p16=float(np.percentile(r, 16)), p84=float(np.percentile(r, 84)),
                shape_median=float(np.median(r / normalisation)),
                shape_rms=float(np.sqrt(np.mean(shape ** 2))),
                shape_iqr=(float(np.percentile(r / normalisation, 16)),
                           float(np.percentile(r / normalisation, 84))),
                frac_10=float(np.mean(np.abs(shape) < 0.10)),
                frac_20=float(np.mean(np.abs(shape) < 0.20)),
                rms_over_stat=float(np.sqrt(np.mean(shape ** 2 / relvar))))


def coarsen(hist, factor):
    """Sum factor x factor blocks (0.25 dex for factor=5)."""
    nth, npb = hist.shape
    return hist.reshape(nth // factor, factor, npb // factor, factor).sum(axis=(1, 3))


def shape_values(n_ours, w_row, f_hist, nmin=100, eta_range=None,
                 logp_range=None, factor=COARSE):
    """(bin ratio / global ratio) in coarsened bins with enough statistics."""
    n_coarse = coarsen(n_ours.astype(float), factor)
    f_coarse = coarsen(f_hist, factor)
    ours = n_coarse * w_row
    mask = shape_mask(n_coarse, f_coarse, nmin, eta_range, logp_range)
    normalisation = ours.sum() / f_coarse.sum()
    return (ours[mask] / f_coarse[mask]) / normalisation


def shape_values_err(n_ours, w_row, f_hist, f_var, nmin=100,
                     eta_range=None, logp_range=None, factor=COARSE):
    """Shape ratios and their relative statistical errors in coarsened bins."""
    n_c = coarsen(n_ours.astype(float), factor)
    f_c = coarsen(f_hist, factor)
    v_c = coarsen(f_var, factor)
    ours = n_c * w_row
    mask = shape_mask(n_c, f_c, nmin, eta_range, logp_range)
    normalisation = ours.sum() / f_c.sum()
    vals = (ours[mask] / f_c[mask]) / normalisation
    rel = np.sqrt(1.0 / n_c[mask] + v_c[mask] / f_c[mask] ** 2)
    return vals, rel


def band_profile(n_ours, w_row, f_hist, axis, edges):
    """Shape ratio (ours/FOR / global ratio) integrated in 1-decade bands."""
    ours = n_ours.astype(float) * w_row
    normalisation = ours.sum() / f_hist.sum()
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        if axis == 0:
            sel = slice(*np.searchsorted(TE, [lo, hi]))
            o = ours[sel, :].sum()
            f = f_hist[sel, :].sum()
        else:
            pc = p_centres(ours.shape[1])
            sel = np.where((pc >= lo) & (pc < hi))[0]
            o = ours[:, sel].sum()
            f = f_hist[:, sel].sum()
        out.append((lo, hi, float((o / f) / normalisation) if f > 0 else np.nan))
    return out


def peak(hist, te, pe):
    i, j = np.unravel_index(np.argmax(hist), hist.shape)
    return float(te[i]), float(pe[j]), float(hist[i, j])


# ----------------------------------------------------------------------------
# figure
# ----------------------------------------------------------------------------
def make_figure(out_png, out_pdf, data):
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.3,
                         "axes.labelsize": 9, "figure.titlesize": 12,
                         "legend.fontsize": 8, "xtick.labelsize": 8,
                         "ytick.labelsize": 8})
    fig, axes = plt.subplots(3, 3, figsize=(15.5, 12.8))
    fig.subplots_adjust(left=0.055, right=0.985, bottom=0.055, top=0.925,
                        wspace=0.28, hspace=0.35)

    def map_panel(ax, hist, pe, title, norm, label="d$\\sigma$/bin [pb]"):
        im = ax.pcolormesh(TE, pe, np.ma.masked_less_equal(hist, 0).T,
                           norm=norm, cmap="viridis", shading="flat",
                           rasterized=True)
        ax.set_xlim(TE[0], TE[-1])
        ax.set_ylim(pe[0], pe[-1])
        cb = plt.colorbar(im, ax=ax, pad=0.012, fraction=0.046)
        cb.set_label(label, fontsize=8)
        ax.set_title(title)
        ax.set_xlabel(r"$\log_{10}\theta$ [rad]")
        ax.set_ylabel(r"$\log_{10}(p/$GeV$)$")

    def ratio_panel(ax, coarse_hist, pe, title):
        te_c = 0.5 * (TE[:-COARSE] + TE[COARSE:])[::COARSE][:coarse_hist.shape[0]]
        pe_c = 0.5 * (pe[:-COARSE] + pe[COARSE:])[::COARSE][:coarse_hist.shape[1]]
        im = ax.pcolormesh(te_c, pe_c, np.ma.masked_invalid(coarse_hist).T,
                           norm=TwoSlopeNorm(vcenter=1.0, vmin=0.7, vmax=1.3),
                           cmap="RdBu_r", shading="nearest", rasterized=True)
        cb = plt.colorbar(im, ax=ax, pad=0.012, fraction=0.046)
        cb.set_label("ours / FORESEE (normalised)", fontsize=8)
        ax.set_title(title)
        ax.set_xlabel(r"$\log_{10}\theta$ [rad]")
        ax.set_ylabel(r"$\log_{10}(p/$GeV$)$")

    b, l = data["b"], data["light"]

    # row 0 -- B mesons
    map_panel(axes[0, 0], b["ours_hist"], PE_B,
              "ours: B0, B0bar, B+, B- (Pythia 8.316, 13 TeV)", LogNorm(1e2, 2e6))
    map_panel(axes[0, 1], b["p8_hist"], PE_B,
              "FORESEE Pythia8 (13.6 TeV, forward hemisphere)", LogNorm(1e2, 2e6))
    ratio_panel(axes[0, 2], b["ratio_map"], PE_B,
                "B mesons: ours / FORESEE, normalised\n" + b["ratio_title"])

    # row 1 -- pi+
    map_panel(axes[1, 0], l["ours_pip_hist"], PE_L, r"ours: $\pi^+$ (13 TeV)",
              LogNorm(1e4, 5e9))
    map_panel(axes[1, 1], l["p8_pip_hist"], PE_L,
              r"FORESEE Pythia8 $\pi^+$ (13.6 TeV)", LogNorm(1e4, 5e9))
    ratio_panel(axes[1, 2], l["ratio_map"], PE_L,
                r"$\pi^+$: ours / FORESEE, normalised" + "\n" + l["ratio_title"])

    # row 2 -- distributions and totals
    def shape_kde(ax, vals_rel, color, title, text, xlim):
        vals, rel = vals_rel
        # kernel-density estimate with a bandwidth of a few per-cent, i.e. of
        # the order of the per-bin statistical resolution: the integer bin
        # counts would otherwise produce a comb-like histogram
        grid = np.linspace(xlim[0], xlim[1], 401)
        bw = max(0.04, 1.5 * float(np.median(rel)))
        dens = np.exp(-0.5 * ((grid[:, None] - vals[None, :]) / bw) ** 2).sum(1)
        dens /= vals.size * bw * np.sqrt(2.0 * np.pi)
        ax.fill_between(grid, 0.0, dens, color=color, alpha=0.45, lw=0)
        ax.plot(grid, dens, color=color, lw=1.4)
        ax.axvline(1.0, color="k", lw=1.2)
        ax.axvline(float(np.median(vals)), color=color, lw=1.0, ls="--")
        ax.set_xlabel("bin ratio / global ratio  (0.05 dex bins, KDE)")
        ax.set_ylabel("density")
        ax.set_xlim(*xlim)
        ax.set_ylim(bottom=0.0)
        ax.set_title(title)
        ax.text(0.97, 0.95, text, transform=ax.transAxes, ha="right",
                va="top", fontsize=7.0)

    shape_kde(axes[2, 0], b["shape_vals"], "C0",
              "B mesons: shape-only bin-ratio distribution", b["shape_text"],
              (0.5, 1.5))
    shape_kde(axes[2, 1], l["shape_vals"], "C2",
              r"$\pi^\pm K^\pm$: shape-only bin-ratio distribution",
              l["shape_text"], (0.5, 1.5))

    ax = axes[2, 2]
    keys = ["ours", "p8", "nlo", "max", "min"]
    names = ["ours\n13 TeV", "FORESEE\nPythia8", "FORESEE\nNLO-P8",
             "NLO-P8\nMax", "NLO-P8\nMin"]
    vals = np.array([data["totals"][k] for k in keys]) / 1e6
    colors = ["C0", "C1", "C3", "C3", "C3"]
    bars = ax.bar(np.arange(len(vals)), vals, color=colors, alpha=0.85, width=0.65)
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, val * 1.02, f"{val:.0f}",
                ha="center", va="bottom", fontsize=8)
    expected = vals[0] * data["energy_factor"]
    ax.axhline(expected, color="k", ls="--", lw=1.0)
    ax.text(4.45, expected * 1.02,
            f"ours x {data['energy_factor']:.3f}\n(same $\\sigma(b\\bar b)$ "
            "at 13.6 TeV)", ha="right", va="bottom", fontsize=7)
    ax.set_xticks(np.arange(len(vals)))
    ax.set_xticklabels(names, fontsize=8)
    ax.set_ylabel(r"forward B cross section [$\mu$b]")
    ax.set_ylim(0, max(vals) * 1.28)
    ax.set_title("B0/B0bar/B+/B- totals, common FORESEE region\n"
                 r"($\theta<1$ rad, $1<p<10^4$ GeV)")
    ax.text(0.03, 0.97,
            "ours/other: " + "  ".join(
                f"{n} {data['totals']['ours'] / data['totals'][k]:.2f}"
                for n, k in zip(["P8", "NLO", "Max", "Min"],
                                ["p8", "nlo", "max", "min"])),
            transform=ax.transAxes, ha="left", va="top", fontsize=7.5)

    for ax, letter in zip(axes.ravel(), "abcdefghi"):
        ax.text(-0.16, 1.03, f"({letter})", transform=ax.transAxes,
                fontsize=10.5, fontweight="bold", va="bottom", ha="left")

    fig.suptitle("Our pp-collision spectra (13 TeV) vs FORESEE forward spectra "
                 "(13.6 TeV Pythia 8): shape agreement and normalisation",
                 y=0.975)
    fig.savefig(out_png, dpi=150)
    fig.savefig(out_pdf)
    plt.close(fig)


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-download", action="store_true")
    parser.add_argument("--quick", action="store_true",
                        help="use the 100k-event pi/K event-subsample cache")
    args = parser.parse_args(argv)

    start = time.time()
    ensure_foresee_files(allow_download=not args.no_download)

    print("=" * 78)
    print("our pp-collision spectra vs FORESEE forward spectra")
    print("=" * 78)
    print(f"run dir       : {RUN}")
    print(f"FORESEE files : {FOR_HADRONS}")

    # ---------------- our samples
    print("\n-- our samples ----------------------------------------------------")
    n_b4 = our_b_counts(B4_PIDS)
    n_b8 = our_b_counts(B_PIDS)
    w_b = SIGMA_BB_PB / N_GEN_B
    ours_b4_hist = n_b4 * w_b
    ours_b8_hist = n_b8 * w_b
    print(f"B   : {int(n_b8.sum()):d} forward B rows (pz>0), weight {w_b:.0f} pb/row"
          "  (sigma(bbbar)=560 ub at 13 TeV, n_gen=800000)")
    print(f"      total forward: B0/B0bar/B+/B- = {tot(ours_b4_hist) / 1e6:.1f} ub,"
          f" all 8 B-meson states = {tot(ours_b8_hist) / 1e6:.1f} ub")
    light = our_light_counts(use_cache=args.quick)
    w_l = SIGMA_INEL_PB / N_GEN_INEL
    n_light = sum(int(round(c.sum())) for c in light.values())
    print(f"inel: {n_light:d} forward pi/K rows, weight {w_l:.4g} pb/row"
          "  (sigma_inel=79.5 mb, n_gen=500000)")
    light_hist = sum(light.values()) * w_l

    # ---------------- FORESEE spectra
    p8 = forsee_species("13.6", "Pythia8", B_PIDS)
    p8_14 = forsee_species("14", "Pythia8", B_PIDS)
    nlo = forsee_species("13.6", "NLO-P8", B_PIDS)
    nlo_14 = forsee_species("14", "NLO-P8", B_PIDS)
    nmax = forsee_species("13.6", "NLO-P8-Max", B4_PIDS)
    nmin = forsee_species("13.6", "NLO-P8-Min", B4_PIDS)

    p8_b4, v_p8_b4 = stack(p8, B4_PIDS)
    p8_b8, _ = stack(p8, B_PIDS)
    p8_14_b4, _ = stack(p8_14, B4_PIDS)
    nlo_b4, v_nlo_b4 = stack(nlo, B4_PIDS)
    nlo_14_b4, _ = stack(nlo_14, B4_PIDS)
    nmax_b4, _ = stack(nmax, B4_PIDS)
    nmin_b4, _ = stack(nmin, B4_PIDS)

    print("\n-- FORESEE B-meson spectra ---------------------------------------")
    print(f"Pythia8 13.6: B0/B0bar/B+/B- {p8_b4.sum() / 1e6:7.2f} ub, "
          f"all 8 states {p8_b8.sum() / 1e6:7.2f} ub")
    for label, hist in [("Pythia8 14  ", p8_14_b4), ("NLO-P8 13.6 ", nlo_b4),
                        ("NLO-Max 13.6", nmax_b4), ("NLO-Min 13.6", nmin_b4)]:
        print(f"{label}: B0/B0bar/B+/B- {hist.sum() / 1e6:7.2f} ub")

    # ---------------- energy scaling 13 -> 13.6 TeV
    def energy_scaling(f136, f14):
        local = f14.sum() / f136.sum()
        power = np.log(local) / np.log(14.0 / 13.6)
        return float((13.6 / 13.0) ** power), float(local)

    e_p8, loc_p8 = energy_scaling(p8_b4, p8_14_b4)
    e_nlo, loc_nlo = energy_scaling(nlo_b4, nlo_14_b4)
    e_avg = 0.5 * (e_p8 + e_nlo)
    print("\n-- known energy scaling of the bbbar cross section, 13 -> 13.6 TeV --")
    print(f"Pythia8 : sigma(14)/sigma(13.6) = {loc_p8:.4f} -> "
          f"sigma(13.6)/sigma(13) = {e_p8:.4f}")
    print(f"NLO-P8  : sigma(14)/sigma(13.6) = {loc_nlo:.4f} -> "
          f"sigma(13.6)/sigma(13) = {e_nlo:.4f}")
    print(f"adopted : {e_avg:.4f}  (+{100 * (e_avg - 1):.1f} % from 13 to 13.6 TeV)")

    # ---------------- totals, normalisation
    ours_common = tot(ours_b4_hist, NLO_SLICE)
    p8_common = p8_b4[NLO_SLICE].sum()
    nlo_common = nlo_b4.sum()
    print("\n-- total forward B-meson cross sections (B0/B0bar/B+/B-) ----------")
    print(f"ours              13 TeV   : {tot(ours_b4_hist) / 1e6:7.1f} ub"
          f"   (all 8 states {tot(ours_b8_hist) / 1e6:7.1f} ub)")
    print(f"FORESEE Pythia8  13.6 TeV : {p8_b4.sum() / 1e6:7.1f} ub"
          f"   (all 8 states {p8_b8.sum() / 1e6:7.1f} ub)"
          f"   ratio ours/FOR = {tot(ours_b4_hist) / p8_b4.sum():.4f}")
    print("\ncommon FORESEE region (theta<1 rad, 1 GeV<p<1e4 GeV):")
    print(f"  ours             {ours_common / 1e6:7.1f} ub")
    print(f"  FORESEE Pythia8  {p8_common / 1e6:7.1f} ub  ratio {ours_common / p8_common:.4f}")
    print(f"  FORESEE NLO-P8   {nlo_common / 1e6:7.1f} ub  ratio {ours_common / nlo_common:.4f}")
    print(f"  FORESEE NLO-Max  {nmax_b4.sum() / 1e6:7.1f} ub  ratio "
          f"{ours_common / nmax_b4.sum():.4f}")
    print(f"  FORESEE NLO-Min  {nmin_b4.sum() / 1e6:7.1f} ub  ratio "
          f"{ours_common / nmin_b4.sum():.4f}")
    print(f"  Pythia8 -> NLO-P8 normalisation factor: {nlo_common / p8_common:.4f}")

    s_p8 = ours_common / p8_common
    s_nlo = ours_common / nlo_common
    print("\nnormalisation decomposition (common region):")
    print(f"  observed ratios              : ours/Pythia8 {s_p8:.4f}, "
          f"ours/NLO-P8 {s_nlo:.4f}")
    print(f"  energy 13->13.6 TeV would give: {e_avg:.4f} for identical sigma(bbbar)")
    print(f"  same-energy normalisation     : Pythia8 {s_p8 * e_avg:.4f}, "
          f"NLO-P8 {s_nlo * e_avg:.4f}")
    print("  implied sigma(bbbar) of the FORESEE variants at 13.6 TeV:")
    for label, forsum in [("Pythia8", p8_common), ("NLO-P8 ", nlo_common),
                          ("NLO-Max", nmax_b4.sum()), ("NLO-Min", nmin_b4.sum())]:
        sigma = SIGMA_BB_PB / (ours_common / forsum)
        print(f"    {label}: {sigma / 1e6:6.1f} ub   (ours: 560 ub at 13 TeV = "
              f"{560 * e_avg:.0f} ub at 13.6 TeV)")

    # ---------------- shape comparison, B
    print("\n-- shape comparison, B mesons ------------------------------------")
    b_p8 = ratio_metrics(n_b4, w_b, p8_b4, v_p8_b4)
    b_nlo = ratio_metrics(n_b4, w_b, nlo_b4, v_nlo_b4, slc=NLO_SLICE)
    b_core = ratio_metrics(n_b4, w_b, p8_b4, v_p8_b4,
                           eta_range=(1.5, 5.5), logp_range=(0.7, 2.7))
    for label, met in [("vs Pythia8 13.6, full forward hemisphere", b_p8),
                       ("vs NLO-P8 13.6, full forward hemisphere", b_nlo),
                       ("vs Pythia8 13.6, analysis core (eta 1.5-5.5, "
                        "5<p<500 GeV)", b_core)]:
        print(f"{label}: {met['n_bins']} bins with our counts>=100")
        print(f"  total ratio {met['total_ratio']:.4f}; bin ratio median "
              f"{met['median']:.4f}, 68% [{met['p16']:.3f}, {met['p84']:.3f}]")
        print(f"  shape-only (bin ratio / global): median {met['shape_median']:.4f}, "
              f"68% [{met['shape_iqr'][0]:.3f}, {met['shape_iqr'][1]:.3f}], "
              f"RMS {100 * met['shape_rms']:.1f}%")
        print(f"  |shape|<10%/20%: {100 * met['frac_10']:.0f}% / "
              f"{100 * met['frac_20']:.0f}%")
    ef_ours = eta_fractions(ours_b4_hist)
    ef_p8 = eta_fractions(p8_b4)
    print("  eta bands [fraction of the forward cross section] and shape ratio:")
    for band in ef_ours:
        print(f"    {band:12s}: ours {100 * ef_ours[band]:5.1f}%  "
              f"Pythia8 {100 * ef_p8[band]:5.1f}%  "
              f"ratio {ef_ours[band] / ef_p8[band]:.3f}")
    print("  same eta-band fractions at 14 vs 13.6 TeV (FORESEE energy check):")
    for label, h136, h14 in [("Pythia8", p8_b4, p8_14_b4),
                             ("NLO-P8 ", nlo_b4, nlo_14_b4)]:
        f136, f14 = eta_fractions(h136), eta_fractions(h14)
        print(f"    {label}: " + ", ".join(
            f"{band}: {f14[band] / f136[band]:.3f}" for band in f136))
    print("  shape ratio (ours/FOR after global normalisation) vs p:")
    for lo, hi, val in band_profile(n_b4, w_b, p8_b4, 1, [-1, 0, 1, 2, 3, 4]):
        print(f"    log10p {lo:+.0f}..{hi:+.0f}: {val:.3f}")
    print("  shape ratio vs log10 theta:")
    for lo, hi, val in band_profile(n_b4, w_b, p8_b4, 0, [-5, -4, -3, -2, -1, 0]):
        print(f"    log10theta {lo:+.0f}..{hi:+.0f}: {val:.3f}")
    tp = peak(ours_b4_hist, TE, PE_B)
    tpf = peak(p8_b4, TE, PE_B)
    print(f"  peak: ours (log10theta,log10p)=({tp[0]:+.2f},{tp[1]:+.2f}) "
          f"FORESEE=({tpf[0]:+.2f},{tpf[1]:+.2f})")

    # ---------------- pi/K comparison (secondary)
    print("\n-- pi/K comparison (secondary; ours = prompt fragmentation hadrons) --")
    f_light = {}
    for gen in LIGHT_GENS:
        try:
            f_light[gen] = forsee_species("13.6", gen, LIGHT_PIDS)
        except SystemExit:
            print(f"[warn] {gen} light-hadron files missing; skipped")
    print("total forward cross sections [mb]:")
    print(f"  {'species':6s} {'ours':>12s}" +
          "".join(f" {g:>12s}" for g in f_light))
    for pid, name in [(211, "pi+"), (-211, "pi-"), (321, "K+"), (-321, "K-")]:
        line = f"  {name:6s} {light[pid].sum() * w_l / 1e9:12.1f}"
        for gen in f_light:
            line += f" {f_light[gen][pid][0].sum() / 1e9:12.1f}"
        print(line)
    light_p8_hist, v_light = stack(f_light["Pythia8"], LIGHT_PIDS)
    for gen in f_light:
        forsum, _ = stack(f_light[gen], LIGHT_PIDS)
        print(f"  pi/K total ours/FOR({gen}) = {light_hist.sum() / forsum.sum():.4f}")
    light_met = ratio_metrics(sum(light.values()), w_l, light_p8_hist, v_light)
    light_core = ratio_metrics(sum(light.values()), w_l, light_p8_hist, v_light,
                               eta_range=(2.0, 5.0), logp_range=(1.0, 2.7))
    print(f"  pi/K bins with our counts>=100: {light_met['n_bins']}")
    print(f"  full forward: total ratio {light_met['total_ratio']:.4f}; bin ratio "
          f"median {light_met['median']:.3f}, 68% [{light_met['p16']:.3f}, "
          f"{light_met['p84']:.3f}]; shape RMS {100 * light_met['shape_rms']:.1f}%")
    if light_core is not None:
        print(f"  core (eta 2-5, 10<p<500 GeV): bin ratio median "
              f"{light_core['median']:.3f}, 68% [{light_core['p16']:.3f}, "
              f"{light_core['p84']:.3f}], shape RMS "
              f"{100 * light_core['shape_rms']:.1f}%")
    print("  FORESEE light-hadron generator spread (pi+ total, mb): " +
          ", ".join(f"{g} {f_light[g][211][0].sum() / 1e9:.0f}" for g in f_light))
    kp_ours = light[321].sum() + light[-321].sum()
    pp_ours = light[211].sum() + light[-211].sum()
    spread = ", ".join(
        f"{g} {(f_light[g][321][0].sum() + f_light[g][-321][0].sum()) / (f_light[g][211][0].sum() + f_light[g][-211][0].sum()):.3f}"
        for g in f_light)
    print(f"  K/pi: ours {kp_ours / pp_ours:.3f}; FORESEE {spread}")

    # ---------------- figure
    def coarse_ratio_map(ours_hist, f_hist, min_counts):
        o = coarsen(ours_hist, COARSE)
        f = coarsen(f_hist, COARSE)
        normalisation = o.sum() / f.sum()
        mask = (o >= min_counts) & (f > 0)
        ratio = np.full(o.shape, np.nan)
        ratio[mask] = (o[mask] / f[mask]) / normalisation
        return ratio

    data = {
        "b": {
            "ours_hist": ours_b4_hist,
            "p8_hist": p8_b4,
            "ratio_map": coarse_ratio_map(ours_b4_hist, p8_b4, min_counts=10),
            "shape_vals": shape_values_err(n_b4, w_b, p8_b4, v_p8_b4,
                                           eta_range=(1.5, 5.5),
                                           logp_range=(0.7, 2.7), factor=1),
            "ratio_title": (f"full hemisphere: global ratio "
                            f"{b_p8['total_ratio']:.2f}, shape RMS "
                            f"{100 * b_p8['shape_rms']:.0f}%"),
            "shape_text": (f"eta 1.5-5.5, 5<p<500 GeV:\n"
                           f"raw bin ratio median {b_core['median']:.2f}, 68% "
                           f"[{b_core['p16']:.2f}, {b_core['p84']:.2f}]\n"
                           f"shape ratio (plotted): median "
                           f"{b_core['shape_median']:.2f}, RMS "
                           f"{100 * b_core['shape_rms']:.0f}%\n"
                           f"{100 * b_core['frac_20']:.0f}% of bins within 20%"),
        },
        "light": {
            "ours_pip_hist": light[211] * w_l,
            "p8_pip_hist": f_light["Pythia8"][211][0],
            "ratio_map": coarse_ratio_map(light[211] * w_l,
                                          f_light["Pythia8"][211][0],
                                          min_counts=50),
            "shape_vals": shape_values_err(sum(light.values()), w_l,
                                           light_p8_hist, v_light,
                                           eta_range=(2.0, 5.0),
                                           logp_range=(1.0, 2.7), factor=1),
            "ratio_title": (f"full hemisphere: global ratio "
                            f"{light_met['total_ratio']:.2f}, shape RMS "
                            f"{100 * light_met['shape_rms']:.0f}%"),
            "shape_text": (f"eta 2-5, 10<p<500 GeV:\n"
                           f"raw bin ratio median {light_core['median']:.2f}, 68% "
                           f"[{light_core['p16']:.2f}, {light_core['p84']:.2f}]\n"
                           f"shape ratio (plotted): RMS "
                           f"{100 * light_core['shape_rms']:.0f}%\n"
                           "prompt vs inclusive, generator-dependent"),
        },
        "totals": {"ours": ours_common, "p8": p8_common, "nlo": nlo_common,
                   "max": nmax_b4.sum(), "min": nmin_b4.sum()},
        "energy_factor": e_avg,
    }
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    out_png = FIG_DIR / "spectrum_vs_foresee_overview.png"
    out_pdf = FIG_DIR / "spectrum_vs_foresee_overview.pdf"
    make_figure(out_png, out_pdf, data)
    print(f"\n[figure] {out_png}")
    print(f"[figure] {out_pdf}")
    print(f"\ndone in {time.time() - start:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
