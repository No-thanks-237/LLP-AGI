#!/usr/bin/env python
"""make_grids.py -- produce the (mH, tanb) width/BR tables of the decay module.

Outputs
  <run>/outputs/data/decay_widths_grid.csv   full width+BR table, 67 mass nodes x 40 tanb nodes
  <run>/work/analysis/grids/scan_points_tanb40.txt   "<mH> <tanb>" one line per point (scan input)
  <run>/work/analysis/grids/widths_img_13cols.txt    13 numbers per point, same order as the batch
                                             executable stdout (mH tanb Gtot Ggaga Gee Gmumu
                                             Gtautau Gpipi GKK G4pi Gcc Gss Ggg)
  <run>/work/decay/scan_points_tanb40.txt    identical copy (owner's area)
  <run>/work/decay/widths_img_13cols.txt     identical copy (owner's area)
Provenance (exe sha256, grids, commands) is written to work/decay/grids_provenance.json.
"""
import hashlib, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
import decay_widths as dw

HBARC = 1.973269804e-16
TANB = np.logspace(np.log10(5.0), np.log10(1e5), 40)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_mass_grid():
    out = []
    for ln in open(os.path.join(RUN, "work", "analysis", "grids", "mass_grid.txt")):
        s = ln.split("#")[0].strip()
        if s:
            out.append(float(s.split()[0]))
    return out


def main():
    masses = read_mass_grid()
    exe = dw.find_executable()
    dws = dw.DecayWidths(exe)
    pts = [(m, float(t)) for m in masses for t in TANB]
    rows = dws.batch(pts)

    datadir = os.path.join(RUN, "outputs", "data")
    common = os.path.join(RUN, "work", "analysis", "grids")
    os.makedirs(datadir, exist_ok=True)
    os.makedirs(common, exist_ok=True)

    exe_sha = sha256(exe)
    hdr = ("# 2HDM light-H decay widths and branching fractions\n"
           "# Type-I 2HDM, m_h=125, v=246, m_A=m_Hpm=600 GeV, lambda v^2=0, cos(beta-alpha)=1/tanb\n"
           "# AUTHOR Light_scalar_decay program (github.com/shiggs90/Light_scalar_decay), batch driver,\n"
           "# only adaptations: M_PI definition + main() -> drv.cpp batch driver. No physics edits.\n"
           "# executable sha256: %s\n"
           "# widths in GeV; ctau = hbar*c/Gamma_tot with hbar*c = %.10e GeV m; BR = Gamma_i/Gamma_tot\n"
           "# columns: mH,tanb,GammaTot_GeV,ctau_m,ctau_m_if_hbarc2.0e-16,Gamma_{gaga,ee,mumu,tautau,pipi,KK,4pi,cc,ss,gg},BR_{same order}\n"
           "# last column mH_exact_2GeV_double_counted: 1 only at mH=2.0 where the author program\n"
           "# evaluates BOTH the hadronic and the partonic channels (strict-inequality matching),\n"
           "# so Gamma_tot is too large there (+4.9%% at tanb=10, +43%% at tanb=5, +85%% at tanb=1).\n"
           "# Interpolate across that node (grid neighbours 1.81693 / 2.07163 GeV) instead.\n"
           % (exe_sha, HBARC))
    csv_path = os.path.join(datadir, "decay_widths_grid.csv")
    with open(csv_path, "w") as f:
        f.write(hdr)
        f.write("mH,tanb,GammaTot_GeV,ctau_m,ctau_m_hbarc2e-16," +
                ",".join("Gamma_%s_GeV" % c for c in dw.CHANNELS) + "," +
                ",".join("BR_%s" % c for c in dw.CHANNELS) + ",mH_exact_2GeV_double_counted\n")
        for (m_in, t_in), r in zip(pts, rows):
            f.write("%.17g,%.17g,%.10e,%.10e,%.10e," % (m_in, t_in, r["gamma_total"],
                                                          r["ctau"], 2.0e-16 / r["gamma_total"]))
            f.write(",".join("%.10e" % r["widths"][c] for c in dw.CHANNELS) + ",")
            f.write(",".join("%.10e" % r["br"][c] for c in dw.CHANNELS) +
                    ",%d\n" % (1 if abs(m_in - 2.0) < 1e-9 else 0))

    # raw 13-column block for the scan (Lead's interface)
    ptxt = "".join("%.17g %.17g\n" % (m, t) for m, t in pts)   # exact round-trip
    wtxt = ("# scan input block for the 2HDM light-H decay stage\n"
            "# executable sha256: %s\n"
            "# mass grid: work/analysis/grids/mass_grid.txt (%d nodes); tanb grid: numpy logspace(log10(5),log10(1e5),40)\n"
            "# one line per point, same order as scan_points_tanb40.txt, 13 whitespace-separated numbers:\n"
            "#   mH tanb GammaTot Gamma_gaga Gamma_ee Gamma_mumu Gamma_tautau Gamma_pipi Gamma_KK\n"
            "#           Gamma_4pi Gamma_cc Gamma_ss Gamma_gg      (all widths in GeV, %%.10e)\n"
            % (exe_sha, len(masses)))
    for (m_in, t_in), r in zip(pts, rows):
        wtxt += "%.17g %.17g" % (m_in, t_in) + " %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n" % (
            r["gamma_total"],
            *[r["widths"][c] for c in dw.CHANNELS])
    for path in (os.path.join(common, "scan_points_tanb40.txt"),
                 os.path.join(HERE, "scan_points_tanb40.txt")):
        open(path, "w").write(ptxt)
    for path in (os.path.join(common, "widths_img_13cols.txt"),
                 os.path.join(HERE, "widths_img_13cols.txt")):
        open(path, "w").write(wtxt)

    # ---- variant without the mH = 2.0 node (matches work/analysis/grids/scan_mass_grid.txt,
    # where the Lead dropped the double-counted node).  Same 40-point tanb grid.
    scan_grid_file = os.path.join(common, "scan_mass_grid.txt")
    variant = None
    if os.path.isfile(scan_grid_file):
        masses2 = [float(l.split("#")[0].strip()) for l in open(scan_grid_file)
                   if l.split("#")[0].strip()]
        pts2 = [(m, float(t)) for m in masses2 for t in TANB]
        rows2 = dws.batch(pts2)
        ptxt2 = "".join("%.17g %.17g\n" % (m, t) for m, t in pts2)
        wtxt2 = ("# scan input block for the 2HDM light-H decay stage: NO mH=2.0 node\n"
                 "# (the author program double-counts hadron+parton at exactly 2.0 GeV;\n"
                 "#  this file matches work/analysis/grids/scan_mass_grid.txt, 66 mass nodes)\n"
                 "# executable sha256: %s\n"
                 "# tanb grid: numpy logspace(log10(5), log10(1e5), 40); %d points\n"
                 "# 13 numbers per line: mH tanb GammaTot Gamma_gaga Gamma_ee Gamma_mumu\n"
                 "#   Gamma_tautau Gamma_pipi Gamma_KK Gamma_4pi Gamma_cc Gamma_ss Gamma_gg [GeV]\n"
                 % (exe_sha, len(pts2)))
        for (m_in, t_in), r in zip(pts2, rows2):
            wtxt2 += ("%.17g %.17g" % (m_in, t_in) +
                      " %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n" % (
                          r["gamma_total"], *[r["widths"][c] for c in dw.CHANNELS]))
        paths2 = [os.path.join(common, "scan_points_tanb40_no2GeV.txt"),
                  os.path.join(HERE, "scan_points_tanb40_no2GeV.txt")]
        for p in paths2:
            open(p, "w").write(ptxt2)
        paths2b = [os.path.join(common, "widths_img_13cols_no2GeV.txt"),
                   os.path.join(HERE, "widths_img_13cols_no2GeV.txt")]
        for p in paths2b:
            open(p, "w").write(wtxt2)
        variant = dict(mass_grid_file=scan_grid_file,
                       mass_grid_sha256=sha256(scan_grid_file),
                       mass_grid_nodes=len(masses2), n_points=len(rows2),
                       points=paths2[0], points_sha256=sha256(paths2[0]),
                       widths13=paths2b[0], widths13_sha256=sha256(paths2b[0]))
        print("wrote", paths2b[0], len(rows2), "rows (no mH=2.0 variant)")

    prov = dict(
        executable=exe, executable_sha256=exe_sha,
        wrapper=os.path.join(HERE, "decay_widths.py"),
        wrapper_sha256=sha256(os.path.join(HERE, "decay_widths.py")),
        mass_grid_file=os.path.join(RUN, "work", "analysis", "grids", "mass_grid.txt"),
        mass_grid_sha256=sha256(os.path.join(RUN, "work", "analysis", "grids", "mass_grid.txt")),
        mass_grid_nodes=len(masses),
        tanb_grid="logspace(log10(5), log10(1e5), 40)",
        tanb_grid_sha256=hashlib.sha256("".join("%.17g\n" % t for t in TANB).encode()).hexdigest(),
        n_points=len(rows), hbarc_GeV_m=HBARC,
        outputs=dict(csv=csv_path, csv_sha256=sha256(csv_path),
                     points=os.path.join(common, "scan_points_tanb40.txt"),
                     points_sha256=sha256(os.path.join(common, "scan_points_tanb40.txt")),
                     widths13=os.path.join(common, "widths_img_13cols.txt"),
                     widths13_sha256=sha256(os.path.join(common, "widths_img_13cols.txt"))),
        variant_no_2GeV=variant,
        timestamp_utc=__import__("datetime").datetime.utcnow().isoformat() + "Z",
    )
    json.dump(prov, open(os.path.join(HERE, "grids_provenance.json"), "w"), indent=1)
    print("wrote", csv_path, len(rows), "rows")
    print("exe sha256", exe_sha)
    print(json.dumps(prov["outputs"], indent=1))


if __name__ == "__main__":
    main()
