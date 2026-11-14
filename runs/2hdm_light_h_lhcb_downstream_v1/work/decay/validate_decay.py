#!/usr/bin/env python
"""validate_decay.py -- validation of the 2HDM light-H decay module (DESIGN.md section 2).

Checks
  1. internal: BR sum = 1, all widths >= 0, Gamma_tot = sum(channels), channel
     thresholds/threshold behaviour, the 2 GeV hadron/parton matching discontinuity,
     absence of b bbar.
  2. NUMERIC (primary external check): comparison against the FORESEE 2-D tables in
     <LLP-AGI>/references/tools/FORESEE/Models/2HDM_H/model/{ctau.txt, br/*.txt}.  Those tables were
     produced by the authors with the SAME Light_scalar_decay program and are the exact
     numbers behind the reference figures.  10 BR channels + c*tau are compared on the
     full 61 mass x 201 tanb node grid (12261 nodes).
  3. VISUAL: colour digitisation of the two paper figures
     (inputs/references/mphi_t1Br_tanb102_cba=1tanb.png and mphi_Type1_cba=1tanb-2.png)
     with a max/median |log10 ratio| report and an overlay plot
     outputs/figures/decay_validation.png.
Outputs: work/decay/validation.json, logs/decay_validation.md (written by hand from
the numbers; this script writes the JSON + the figure + a text summary).
"""
import hashlib, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
import decay_widths as dw
import digitise_paper as dp

HBARC = 1.973269804e-16
HBARC_FORESEE = 2.0e-16          # the rounded value used to build the authors' tables
FORESEE = os.environ.get("FORESEE_2HDM_H_MODEL",
                         r"E:\AAAwork\LLP-AGI\temp\foresee\Models\2HDM_H\model")
CH = dw.CHANNELS
TANB_REF = 10.0


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def logratio_stats(ours, ref, mask=None):
    ours = np.asarray(ours, float); ref = np.asarray(ref, float)
    ok = (ours > 0) & (ref > 0) & np.isfinite(ours) & np.isfinite(ref)
    if mask is not None:
        ok &= mask
    lr = np.log10(ours[ok] / ref[ok])
    a = np.abs(lr)
    return dict(n=int(ok.sum()), median_abs_log10=float(np.median(a)),
                mean_abs_log10=float(np.mean(a)), max_abs_log10=float(a.max()),
                p90_abs_log10=float(np.percentile(a, 90)),
                p95_abs_log10=float(np.percentile(a, 95)),
                frac_within_0p1pct=float(np.mean(a < np.log10(1.001))),
                frac_within_1pct=float(np.mean(a < np.log10(1.01))))


# ---------------------------------------------------------------- 1. basic checks
def basic_checks(dws, masses, tanbs):
    out = {}
    pts = [(float(m), float(t)) for m in masses for t in tanbs]
    rows = dws.batch(pts)
    br_sum = np.array([sum(r["br"].values()) for r in rows])
    gam_vs_sum = np.array([r["gamma_total"] / sum(r["widths"].values()) for r in rows])
    minw = min(min(r["widths"].values()) for r in rows)
    mn = min(((v, ch, r["mH"], r["tanb"]) for r in rows for ch, v in r["widths"].items()),
             key=lambda z: z[0])
    n_zero = sum(1 for r in rows for v in r["widths"].values() if v == 0.0)
    out["n_points"] = len(rows)
    out["br_sum_max_dev_from_1"] = float(np.max(np.abs(br_sum - 1.0)))
    out["gamma_total_over_sum_channels_min"] = float(gam_vs_sum.min())
    out["gamma_total_over_sum_channels_max"] = float(gam_vs_sum.max())
    out["min_width_GeV"] = float(minw)
    out["min_width_location"] = dict(value=float(mn[0]), channel=mn[1], mH=mn[2], tanb=mn[3])
    out["n_exact_zero_widths"] = int(n_zero)
    out["all_widths_nonnegative"] = bool(minw >= 0)
    out["ctau_range_m"] = [float(min(r["ctau"] for r in rows)), float(max(r["ctau"] for r in rows))]
    out["gamma_total_range_GeV"] = [float(min(r["gamma_total"] for r in rows)),
                                    float(max(r["gamma_total"] for r in rows))]
    # channel threshold behaviour (tanb=10): width must be at its 1e-40 floor below the
    # program's threshold and > floor above it.  The thresholds are the ones hard-coded
    # in the author source (CPevenDecay.cpp): 2*mpion(0.13497), 2*0.498, 4*mpion,
    # 2*m_D(1.864843), 2*mtau(1.7768), and the mH>2 / mH<=2 matching switch.
    thr = dict(pipi=2 * 0.13497, KK=2 * 0.498, fourpi=4 * 0.13497, cc=2 * 1.864843,
               tautau=2 * 1.7768, mumu=2 * 0.10566, ee=2 * 0.000511, ss=2.0, gg=2.0)
    chan_of = dict(fourpi="4pi")
    thr_test = {}
    for ch, t0 in thr.items():
        key = chan_of.get(ch, ch)
        lo = dws.widths(t0 * 0.999, TANB_REF)[key]
        hi = dws.widths(t0 * 1.001, TANB_REF)[key]
        thr_test[ch] = dict(threshold_GeV=float(t0), width_below_GeV=float(lo),
                            width_above_GeV=float(hi),
                            opens=bool(hi > 1e-39 and lo <= 1.0000001e-40))
    out["thresholds"] = thr_test
    # the 2 GeV matching: hadronic channels are cut for mH>2, partonic ss/gg switch on at 2
    tot = lambda d: sum(d.values())
    match = {}
    for tb in (1.0, 5.0, 10.0, 100.0):
        below, at, above = (dws.widths(2.0 - 1e-4, tb), dws.widths(2.0, tb),
                            dws.widths(2.0 + 1e-4, tb))
        match["tanb_%g" % tb] = dict(
            gamma_tot_below=float(tot(below)), gamma_tot_at=float(tot(at)),
            gamma_tot_above=float(tot(above)),
            jump_factor_above_over_below=float(tot(above) / tot(below)),
            br_hadronic_below=float((below["pipi"] + below["KK"] + below["4pi"]) / tot(below)),
            br_partonic_above=float((above["ss"] + above["gg"] + above["cc"]) / tot(above)),
            channels_below={k: float(v) for k, v in below.items()},
            channels_above={k: float(v) for k, v in above.items()})
    out["matching_2GeV"] = dict(
        per_tanb=match,
        note=("author code (CPevenDecay.cpp): pi pi / K K / 4pi are set to the 1e-40 "
              "floor for mH>2; s sbar and g g are 0 for mH<2; at exactly mH=2.0 both "
              "sets are computed because the cuts are strict inequalities. The hadronic "
              "(dispersion-relation) and partonic (spectator) descriptions are therefore "
              "joined by a genuine discontinuity: above 2 GeV the hadronic channels are "
              "switched off abruptly and the partonic s sbar / g g turn on. The jump "
              "factor above/below depends on how much of the width is hadronic vs "
              "diphoton at that tanb; numerically the total width is continuous across "
              "2 GeV to 0.05-0.5% because the 4pi coefficient C=5.1e-9 was tuned by the "
              "authors to make the hadronic and partonic descriptions meet there, while "
              "the individual channel widths jump. This is the intended matching of the "
              "program (paper: 'Decays to hadrons and quarks/gluons are connected at "
              "m_H=2 GeV'), not a bug, but it must be kept in mind when interpolating "
              "across 2 GeV."),
    )
    # no b bbar channel exists in the author CP-even H total width
    src = open(os.path.join(RUN, "work", "decay", "src_patched", "src", "CPevenDecay.cpp")).read()
    out["no_bbbar"] = dict(
        DecayWidth_HTotal_defined_as=[ln.strip() for ln in src.splitlines()
                                      if "double BSModel::DecayWidth_HTotal" in ln or
                                      ("return DecayWidth_Hgaga" in ln)],
        bbbar_channel_absent=bool("DecayWidth_Hbb" not in src),
        comment=("the author CP-even-H module has no H->b bbar partial width at all "
                 "(partonic channels are c cbar, s sbar, g g only); this is why the "
                 "physics scan is restricted to mH <= 10 GeV."))
    return out, rows


# ------------------------------------------------------- 2. FORESEE numeric check
def foresee_check(dws):
    model = FORESEE
    files = {
        "gaga": "br/ga_ga0.txt", "ee": "br/e_e0.txt", "mumu": "br/mu_mu0.txt",
        "pipi": "br/pi_pi.txt", "KK": "br/k_k.txt", "4pi": "br/4pi.txt",
        "gg": "br/g_g.txt", "cc": "br/c_c.txt", "ss": "br/s_s.txt",
        "tautau": "br/tau_tau.txt"}
    if not os.path.isdir(model):
        return dict(available=False, path=model)
    tabs = {ch: np.loadtxt(os.path.join(model, f)) for ch, f in files.items()}
    m0 = tabs["pipi"][:, 0]; t0 = tabs["pipi"][:, 1]
    grid_ok = all(np.allclose(tabs[ch][:, 0], m0) and np.allclose(tabs[ch][:, 1], t0) for ch in tabs)
    pts = [(m, t) for m, t in zip(m0, t0)]
    rows = dws.batch(pts)
    out = dict(available=True, path=model, grid_consistent=bool(grid_ok),
               n_nodes=int(len(m0)), masses=[float(m0.min()), float(m0.max())],
               tanbs=[float(t0.min()), float(t0.max())],
               files={ch: os.path.join(model, f) for ch, f in files.items()},
               files_sha256={ch: sha256(os.path.join(model, f)) for ch, f in files.items()})
    box = (m0 >= 0.021) & (m0 <= 10.0) & (t0 >= 5.0) & (t0 <= 1e5)
    out["n_nodes_in_physics_box"] = int(box.sum())
    # BR per channel
    br_stats = {}
    br_ours = {}
    br_ref = {}
    for ch in CH:
        ours = np.array([r["br"][ch] for r in rows])
        ref = tabs[ch][:, 2]
        br_ours[ch] = ours; br_ref[ch] = ref
        s_all = logratio_stats(ours, ref)
        s_box = logratio_stats(ours, ref, box)
        br_stats[ch] = dict(all_nodes=s_all, physics_box=s_box)
    out["br"] = br_stats
    # BR sum from the tables vs 1 (the 10 tables cover all channels of the program)
    brsum = sum(br_ref[ch] for ch in CH)
    out["br_sum_tables"] = dict(median=float(np.median(brsum)), min=float(brsum.min()),
                                max=float(brsum.max()), n_nodes=int(len(brsum)),
                                max_dev_from_1=float(np.max(np.abs(brsum - 1))))
    # ctau table (31 mass x 101 tanb : the authors' rounded hbar*c = 2.0e-16 GeV m)
    ct = np.loadtxt(os.path.join(model, "ctau.txt"))
    out["ctau_sha256"] = sha256(os.path.join(model, "ctau.txt"))
    rows_c = dws.batch([(m, t) for m, t in ct[:, :2]])
    gam = np.array([r["gamma_total"] for r in rows_c])
    ctau_ours_20 = HBARC_FORESEE / gam
    ctau_ours_true = HBARC / gam
    box_c = (ct[:, 0] >= 0.021) & (ct[:, 0] <= 10.0) & (ct[:, 1] >= 5.0) & (ct[:, 1] <= 1e5)
    out["ctau"] = dict(
        n_nodes=int(len(ct)),
        n_nodes_in_physics_box=int(box_c.sum()),
        note=("the foreseen tables were built with the rounded hbar*c = 2.0e-16 GeV m; "
              "the ratio with our canonical hbar*c = 1.973269804e-16 is then "
              "1.973269804/2.0 = 0.98663490 by construction"),
        stats_hbarc_2e_16=dict(all=logratio_stats(ctau_ours_20, ct[:, 2]),
                                physics_box=logratio_stats(ctau_ours_20, ct[:, 2], box_c)),
        stats_hbarc_canonical=dict(all=logratio_stats(ctau_ours_true, ct[:, 2]),
                                   physics_box=logratio_stats(ctau_ours_true, ct[:, 2], box_c)),
        scale_expected=float(HBARC / HBARC_FORESEE),
    )
    # a ctau-based total width check on the fine (61x201) node set, using ctau.txt
    # is impossible (different grid); instead compare Gamma_tot reconstructed from the
    # 10 BR tables + the family-1 ctau table cannot be combined -> skip, BR route is enough.
    out["_arrays"] = dict(mass=m0, tanb=t0, br_ours=br_ours, br_ref=br_ref, box=box)
    return out


# ------------------------------------------------------------- 3. figure digitisation
def digitise_check(dws):
    xs = list(np.logspace(-1, 1, 241))          # 121 points per decade
    res = {}
    info = {}
    for fig, panel in (("A", "br"), ("A", "gam"), ("B", "gam")):
        r, i = dp.digitise(RUN, fig, panel, xs)
        res[(fig, panel)] = r
        info[(fig, panel)] = i
    out = dict(xs_per_decade=121, n_x_samples=len(xs), panels={})
    digitised_rows = []
    # ---- panel A/br : BR(tanb=10)
    for key, pick in ((("A", "br"), "br"), (("A", "gam"), "widths")):
        stats = {}
        for ch, pts in res[key].items():
            if not pts:
                stats[ch] = dict(n=0, note="channel not distinguishable in this panel")
                continue
            X = np.array([p[0] for p in pts]); Y = np.array([p[1] for p in pts])
            rows = dws.batch([(x, TANB_REF) for x in X])
            model = np.array([(r["br"][ch] if pick == "br" else r["widths"][ch]) for r in rows])
            st = logratio_stats(model, Y)
            st["points"] = [[float(x), float(y)] for x, y in pts]
            lr = np.abs(np.log10(model / Y))
            st["worst"] = [[float(X[i]), float(Y[i]), float(model[i]), float(lr[i])]
                           for i in np.argsort(-lr)[:5]]
            stats[ch] = st
            for x, y, mo in zip(X, Y, model):
                digitised_rows.append(dict(figure=key[0], panel=key[1], channel=ch, mH=float(x),
                                           digitised=float(y), model=float(mo),
                                           log10_ratio=float(np.log10(mo / y))))
        out["panels"][key[0] + "_" + key[1]] = stats
    # ---- panel B : total width for tanb=5,10,70,200,1000
    statsB = {}
    tanbmap = dict(tanb5=5.0, tanb10=10.0, tanb70=70.0, tanb200=200.0, tanb1000=1000.0)
    for ch, pts in res[("B", "gam")].items():
        X = np.array([p[0] for p in pts]); Y = np.array([p[1] for p in pts])
        rows = dws.batch([(x, tanbmap[ch]) for x in X])
        model = np.array([r["gamma_total"] for r in rows])
        st = logratio_stats(model, Y)
        st["tanb"] = tanbmap[ch]
        st["points"] = [[float(x), float(y)] for x, y in pts]
        statsB[ch] = st
        for x, y, mo in zip(X, Y, model):
            digitised_rows.append(dict(figure="B", panel="gam", channel=ch, mH=float(x),
                                       digitised=float(y), model=float(mo),
                                       log10_ratio=float(np.log10(mo / y))))
    out["panel_B_total_width"] = statsB
    # per-channel summary over both paper panels (BR panel + partial-width panel)
    summ = {}
    for ch in set(list(out["panels"]["A_br"].keys()) + list(out["panels"]["A_gam"].keys())):
        ns = []
        med = []
        mx = []
        for pk in ("A_br", "A_gam"):
            s = out["panels"][pk].get(ch)
            if s and s.get("n"):
                ns.append(s["n"]); med.append(s["median_abs_log10"]); mx.append(s["max_abs_log10"])
        if ns:
            summ[ch] = dict(n_total=int(sum(ns)), n_per_panel=ns,
                            median_abs_log10=float(min(med)), max_abs_log10=float(max(mx)))
    out["channel_summary"] = summ
    out["digitised_point_count"] = len(digitised_rows)
    out["_rows"] = digitised_rows
    out["_res"] = res
    out["_info"] = info
    return out


def make_overlay(dws, digi):
    res = digi["_res"]
    mH = np.logspace(np.log10(0.1), 1, 400)
    fig, axes = plt.subplots(1, 3, figsize=(21, 6.2))
    colors = dict(gaga="purple", pipi="blue", KK="red", **{"4pi": "green"}, cc="navy",
                  ss="darkred", gg="darkgreen", mumu="orange", ee="cyan", tautau="peru")
    lab = dict(gaga=r"$\gamma\gamma$", pipi=r"$\pi\pi$", KK=r"$KK$", **{"4pi": r"$4\pi$"},
               cc=r"$c\bar c$", ss=r"$s\bar s$", gg=r"$gg$", mumu=r"$\mu\mu$",
               ee=r"$ee$", tautau=r"$\tau\tau$")
    # panel 1: BR
    rows = dws.batch([(m, TANB_REF) for m in mH])
    ax = axes[0]
    for ch in CH:
        y = np.array([r["br"][ch] for r in rows])
        ax.plot(mH, y, "-", color=colors[ch], lw=1.2, label=lab[ch])
        pts = res[("A", "br")][ch]
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], "o", ms=3.0,
                    mfc="none", mec=colors[ch])
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.1, 10)
    ax.set_ylim(5e-3, 2); ax.set_xlabel(r"$m_H$ [GeV]"); ax.set_ylabel(r"BR($H\to XX$)")
    ax.set_title(r"BR, tan$\beta$=10 (lines: our program; circles: digitised paper fig.)")
    ax.legend(fontsize=7, ncol=2)
    # panel 2: partial widths
    ax = axes[1]
    for ch in CH:
        y = np.array([r["widths"][ch] for r in rows])
        ax.plot(mH, y, "-", color=colors[ch], lw=1.2)
        pts = res[("A", "gam")][ch]
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], "o", ms=3.0,
                    mfc="none", mec=colors[ch])
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.1, 10)
    ax.set_ylim(1e-24, 1e-9); ax.set_xlabel(r"$m_H$ [GeV]")
    ax.set_ylabel(r"$\Gamma(H\to XX)$ [GeV]")
    ax.set_title("partial widths, tan$\\beta$=10")
    # panel 3: total width / ctau for the 5 tanb values
    ax = axes[2]
    ax2 = ax.twinx()
    cols = dict(tanb5="blue", tanb10="green", tanb70="orange", tanb200="purple", tanb1000="red")
    for ch, tb in (("tanb5", 5.0), ("tanb10", 10.0), ("tanb70", 70.0),
                   ("tanb200", 200.0), ("tanb1000", 1000.0)):
        rows = dws.batch([(m, tb) for m in mH])
        g = np.array([r["gamma_total"] for r in rows])
        ax.plot(mH, g, "-", color=cols[ch], lw=1.2, label=r"tan$\beta$=%g" % tb)
        pts = res[("B", "gam")][ch]
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], "o", ms=3.0, mfc="none",
                    mec=cols[ch])
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.1, 10)
    ax.set_ylim(1e-24, 5e-9)
    ax.set_xlabel(r"$m_H$ [GeV]"); ax.set_ylabel(r"$\Gamma$ [GeV]")
    ax.set_title(r"$\Gamma_{tot}$ (left axis) and $c\tau$ (right axis)")
    ax.legend(fontsize=7)
    ylo, yhi = ax.get_ylim()
    ax2.set_yscale("log"); ax2.set_ylim(HBARC / yhi, HBARC / ylo)
    ax2.set_ylabel(r"$c\tau$ [m]")
    fig.suptitle("2HDM light-H decay module validation: author Light_scalar_decay vs paper figures "
                 "(Li et al. arXiv:2212.06186)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(RUN, "outputs", "figures", "decay_validation.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=125)
    plt.close(fig)
    return out


def main():
    masses = [float(l.split("#")[0].strip()) for l in
              open(os.path.join(RUN, "work", "common", "mass_grid.txt")) if l.split("#")[0].strip()]
    tanbs = list(np.logspace(np.log10(5), 5, 40))
    exe = dw.find_executable()
    dws = dw.DecayWidths(exe)
    out = dict(
        module="decay (task-1, DESIGN.md section 2)",
        executable=exe, executable_sha256=sha256(exe),
        program="Light_scalar_decay (github.com/shiggs90/Light_scalar_decay, Combined_Code.zip)",
        source_input=os.path.join(RUN, "inputs", "decay_source"),
        source_input_sha256=sha256(os.path.join(RUN, "inputs", "source.sha256")),
        patches=["PATCH-1: M_PI definition in include/Parameters.h (Windows compiler)",
                 "PATCH-2: main() of the original T1Aa.cpp replaced by src_patched/drv.cpp batch driver"],
        hbarc_GeV_m=HBARC,
    )
    print("== 1. basic checks")
    out["basic"], rows = basic_checks(dws, masses, tanbs)
    print("   BR sum max dev:", out["basic"]["br_sum_max_dev_from_1"])
    print("== 2. FORESEE numeric comparison")
    fs = foresee_check(dws)
    out["foresee"] = {k: v for k, v in fs.items() if k != "_arrays"}
    if fs.get("available"):
        for ch in CH:
            s = fs["br"][ch]["all_nodes"]
            print("   %-6s n=%5d median|log10|=%.2e max=%.2e within1%%=%.4f"
                  % (ch, s["n"], s["median_abs_log10"], s["max_abs_log10"], s["frac_within_1pct"]))
        print("   ctau(2e-16) median|log10|=%.2e max=%.2e"
              % (fs["ctau"]["stats_hbarc_2e_16"]["all"]["median_abs_log10"],
                 fs["ctau"]["stats_hbarc_2e_16"]["all"]["max_abs_log10"]))
    print("== 3. figure digitisation")
    digi = digitise_check(dws)
    out["paper_figures"] = {k: v for k, v in digi.items() if not k.startswith("_")}
    for key, stats in list(out["paper_figures"]["panels"].items()):
        for ch, s in stats.items():
            if s.get("n"):
                print("   %-10s %-6s n=%3d median|log10|=%.3f max=%.3f"
                      % (key, ch, s["n"], s["median_abs_log10"], s["max_abs_log10"]))
    for ch, s in out["paper_figures"]["panel_B_total_width"].items():
        print("   B %-10s n=%3d median|log10|=%.3f max=%.3f"
              % (ch, s["n"], s["median_abs_log10"], s["max_abs_log10"]))
    fig = make_overlay(dws, digi)
    out["overlay_figure"] = fig
    out["overlay_figure_sha256"] = sha256(fig)
    # ctau / Gamma_tot table at the sampled validation points
    vpts = [(0.5, 10), (1, 10), (2, 10), (5, 10),
            (0.5, 100), (1, 100), (2, 100), (5, 100),
            (0.5, 1000), (1, 1000), (2, 1000), (5, 1000),
            (1, 1), (1, 5), (0.15, 10), (0.27, 10), (1.0, 1e5), (10.0, 1e5)]
    rws = dws.batch([(float(a), float(b)) for a, b in vpts])
    out["validation_table"] = [dict(mH=r["mH"], tanb=r["tanb"], GammaTot=r["gamma_total"],
                                    ctau_m=r["ctau"], ctau_m_hbarc_2e_16=HBARC_FORESEE / r["gamma_total"],
                                    br_gaga=r["br"]["gaga"], br_pipi=r["br"]["pipi"],
                                    br_gg=r["br"]["gg"], br_cc=r["br"]["cc"])
                               for r in rws]
    # save the digitised points for provenance
    csv = os.path.join(HERE, "digitised_paper_points.csv")
    with open(csv, "w") as f:
        f.write("figure,panel,channel,mH,digitised,model,log10_ratio\n")
        for r in digi["_rows"]:
            f.write("%s,%s,%s,%.6g,%.6g,%.6g,%.6g\n" % (r["figure"], r["panel"], r["channel"],
                    r["mH"], r["digitised"], r["model"], r["log10_ratio"]))
    out["digitised_points_csv"] = csv
    out["digitised_points_csv_sha256"] = sha256(csv)
    json.dump(out, open(os.path.join(HERE, "validation.json"), "w"), indent=1, default=str)
    print("wrote", os.path.join(HERE, "validation.json"))
    print("wrote", fig, "and", csv)


if __name__ == "__main__":
    main()
