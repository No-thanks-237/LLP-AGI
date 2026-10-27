#!/usr/bin/env python
"""Per-particle pp production spectra: ours vs FORESEE, one figure per particle.

Produces, for each particle species used as a parent in the analysis,
  outputs/figures/production/spectrum_<species>.{png,pdf}
with three panels: our spectrum (Pythia 8.316, pp 13 TeV), the FORESEE 13.6 TeV
Pythia8 spectrum on the identical (log10 theta, log10 p) grid, and the ratio map
(bin ratio divided by the global ratio, i.e. shape only).  Each figure also states
the numbers: forward cross sections [ub], ratio, same-energy ratio (B mesons),
bin-level shape median/RMS in the analysis acceptance, and the deviation relative
to the combined statistical error.

It also writes the machine-readable table
  outputs/data/spectrum_vs_foresee_numbers.csv
(one row per species) and a coarsened bin-by-bin table
  outputs/data/spectrum_vs_foresee_bins.csv.

Usage:  python work/analysis/spectrum_particle_figures.py [--quick]
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm

HERE = Path(__file__).resolve()
RUN = HERE.parents[2]
sys.path.insert(0, str(HERE.parent))
import spectrum_foresee as S          # noqa: E402

PROD = RUN / "outputs" / "figures" / "production"
DATA = RUN / "outputs" / "data"

#: particle species used as parents in the analysis (pid, label, kind)
SPECIES = [
    (511, "B0", "B"), (-511, "B0bar", "B"),
    (521, "Bp", "B"), (-521, "Bm", "B"),
    (211, "pip", "L"), (-211, "pim", "L"),
    (321, "Kp", "L"), (-321, "Km", "L"),
]

#: analysis acceptance used for the shape numbers (log10 p range and eta range)
CORE = {"B": ((1.5, 5.5), (0.7, 2.7)), "L": ((2.0, 5.0), (1.0, 2.7))}


def eta_to_logtheta(eta):
    """log10 of the polar angle (rad) for a pseudorapidity eta."""
    return np.log10(2.0 * np.arctan(np.exp(-eta)))


def ratio_map(ours_hist, for_hist, factor=S.COARSE, min_counts=10):
    o = S.coarsen(ours_hist, factor)
    f = S.coarsen(for_hist, factor)
    norm = o.sum() / f.sum()
    mask = (o > min_counts) & (f > 0)
    out = np.full(o.shape, np.nan)
    out[mask] = (o[mask] / f[mask]) / norm
    return out, norm


def plot_species(pid, label, kind, n_ours, w_row, for_hist, for_var, outdir):
    pe = S.PE_B if kind == "B" else S.PE_L
    ours = n_ours.astype(float) * w_row
    eta_range, logp_range = CORE[kind]
    met = S.ratio_metrics(n_ours, w_row, for_hist, for_var,
                          eta_range=eta_range, logp_range=logp_range)
    # acceptance-region integrated numbers
    m = S.shape_mask(n_ours, for_hist, nmin=0, eta_range=eta_range,
                     logp_range=logp_range)
    ours_core = float(ours[m].sum())
    for_core = float(for_hist[m].sum())

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.6),
                             gridspec_kw={"width_ratios": [1, 1, 1]})

    def panel_map(ax, hist, title, cmap="viridis"):
        pos = hist[hist > 0]
        vmin = max(float(pos.min()) if pos.size else 1e-9, 1e-6 * float(hist.sum() or 1))
        vmax = float(np.percentile(pos, 99.9)) if pos.size else 1.0
        mesh = ax.pcolormesh(S.TE, pe, hist.T, norm=LogNorm(vmin=vmin, vmax=vmax),
                             cmap=cmap, shading="auto")
        fig.colorbar(mesh, ax=ax, pad=0.01, label="d$\\sigma$/bin [pb]")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("log10 theta [rad]")
        ax.set_ylabel("log10(p/GeV)")
        ax.set_xlim(S.TE[0], S.TE[-1])
        ax.set_ylim(pe[0], pe[-1])

    panel_map(axes[0], ours, f"{label}: ours, Pythia 8.316, pp 13 TeV")
    panel_map(axes[1], for_hist, f"{label}: FORESEE Pythia8, pp 13.6 TeV")

    rmap, norm = ratio_map(ours, for_hist)
    ax = axes[2]
    mesh = ax.pcolormesh(S.TE[::5], pe[::5], rmap.T,
                         cmap="RdBu_r", shading="auto",
                         norm=TwoSlopeNorm(vcenter=1.0, vmin=0.5, vmax=1.5))
    fig.colorbar(mesh, ax=ax, pad=0.01,
                 label="bin ratio / global ratio (shape)")
    ax.set_title(f"{label}: shape ratio ours/FORESEE (global ratio {norm:.2f})",
                 fontsize=10)
    ax.set_xlabel("log10 theta [rad]")
    ax.set_ylabel("log10(p/GeV)")
    ax.set_xlim(S.TE[0], S.TE[-1])
    ax.set_ylim(pe[0], pe[-1])

    # acceptance rectangle
    t_hi = eta_to_logtheta(eta_range[0])
    t_lo = eta_to_logtheta(eta_range[1])
    for ax in axes:
        ax.add_patch(plt.Rectangle((t_lo, logp_range[0]), t_hi - t_lo,
                                   logp_range[1] - logp_range[0], fill=False,
                                   edgecolor="w", lw=1.0, ls="--", alpha=0.9))

    e_avg = 1.0463          # sigma(13.6)/sigma(13) for identical sigma(bbbar)
    lines = [f"forward total: ours {ours.sum() / 1e6:.2f} ub, FORESEE "
             f"{for_hist.sum() / 1e6:.2f} ub",
             f"  ratio ours/FOR = {ours.sum() / for_hist.sum():.3f}"
             + (f", same-energy {ours.sum() / for_hist.sum() / e_avg:.3f}"
                if kind == "B" else "")]
    if kind == "B":
        lines.append(f"  (13 -> 13.6 TeV energy factor {e_avg:.3f} applied)")
    lines += [f"acceptance region: ours {ours_core / 1e6:.2f} ub, FORESEE "
              f"{for_core / 1e6:.2f} ub, ratio {ours_core / for_core:.3f}"]
    if met:
        lines += [f"shape: median {met['shape_median']:.3f}, RMS "
                  f"{100 * met['shape_rms']:.1f}%, "
                  f"{100 * met['frac_20']:.0f}% of {met['n_bins']} bins within 20%",
                  f"shape RMS / stat. error = {met['rms_over_stat']:.1f}"]
    text = "\n".join(lines)
    fig.text(0.005, 0.010, text, fontsize=8.2, va="bottom", ha="left",
             family="monospace", linespacing=1.35)
    fig.tight_layout(rect=(0, 0.185, 1, 1))
    png = outdir / f"spectrum_{label}.png"
    pdf = outdir / f"spectrum_{label}.pdf"
    fig.savefig(png, dpi=150)
    fig.savefig(pdf)
    plt.close(fig)
    return dict(pid=pid, name=label, kind=kind,
                ours_ub=ours.sum() / 1e6, for_ub=for_hist.sum() / 1e6,
                ratio=ours.sum() / for_hist.sum(),
                ratio_same_energy=(ours.sum() / for_hist.sum() / e_avg
                                   if kind == "B" else float("nan")),
                ours_core_ub=ours_core / 1e6, for_core_ub=for_core / 1e6,
                core_ratio=ours_core / for_core,
                shape_median=(met["shape_median"] if met else float("nan")),
                shape_rms=(met["shape_rms"] if met else float("nan")),
                frac_20=(met["frac_20"] if met else float("nan")),
                rms_over_stat=(met["rms_over_stat"] if met else float("nan")),
                n_bins=(met["n_bins"] if met else 0))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quick", action="store_true",
                    help="use the 100k-event pi/K event-subsample cache")
    args = ap.parse_args(argv)
    t0 = time.time()
    PROD.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    rows = []
    bins_rows = []
    light = None
    for pid, label, kind in SPECIES:
        if kind == "B":
            n_ours = S.our_b_counts([pid])
            w_row = S.SIGMA_BB_PB / S.N_GEN_B
        else:
            if light is None:
                light = S.our_light_counts(use_cache=args.quick or True)
            n_ours = light[pid]
            w_row = S.SIGMA_INEL_PB / S.N_GEN_INEL
        w, eps = S.read_forsee("13.6", "Pythia8", pid)
        for_var = eps * w
        row = plot_species(pid, label, kind, n_ours, w_row, w, for_var, PROD)
        rows.append(row)
        # coarsened bin-by-bin table
        pe = S.PE_B if kind == "B" else S.PE_L
        o = S.coarsen(n_ours.astype(float) * w_row, S.COARSE)
        f = S.coarsen(w, S.COARSE)
        tc = S.TE[::5]
        pc = pe[::5]
        norm = o.sum() / f.sum()
        for i in range(o.shape[0]):
            for j in range(o.shape[1]):
                if o[i, j] <= 0 and f[i, j] <= 0:
                    continue
                bins_rows.append(dict(pid=pid, name=label,
                                      log10theta=0.5 * (tc[i] + tc[i + 1]),
                                      log10p=0.5 * (pc[j] + pc[j + 1]),
                                      ours_pb=o[i, j], foresee_pb=f[i, j],
                                      ratio=(o[i, j] / f[i, j] if f[i, j] > 0 else ""),
                                      shape_ratio=(o[i, j] / f[i, j] / norm
                                                   if f[i, j] > 0 else "")))
        print(f"{label:6s} ours {row['ours_ub']:8.2f} ub  FORESEE "
              f"{row['for_ub']:8.2f} ub  ratio {row['ratio']:.3f}  "
              f"shape median {row['shape_median']:.3f} RMS "
              f"{100 * row['shape_rms']:.1f}%", flush=True)

    keys = list(rows[0].keys())
    with open(DATA / "spectrum_vs_foresee_numbers.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=keys)
        wr.writeheader()
        for r in rows:
            wr.writerow(r)
    with open(DATA / "spectrum_vs_foresee_bins.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(bins_rows[0].keys()))
        wr.writeheader()
        wr.writerows(bins_rows)
    print(f"\n[table] {DATA / 'spectrum_vs_foresee_numbers.csv'}")
    print(f"[table] {DATA / 'spectrum_vs_foresee_bins.csv'} ({len(bins_rows)} rows)")
    print(f"[figures] {PROD}")
    print(f"done in {time.time() - t0:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
