#!/usr/bin/env python
"""decay_widths.py -- python wrapper around the AUTHOR's Light_scalar_decay CP-even-H
batch executable (DESIGN.md section 2, task-1).

Physics model: Type-I 2HDM light CP-even H benchmark of Li et al. arXiv:2212.06186v1
  m_h=125, v=246, m_A=m_Hpm=600, lambda v^2 = 0, cos(beta-alpha)=1/tanb.
Every partial width is the author program's output; nothing is recomputed here.
Units: GeV (widths), metres (c*tau).

Public API (contract of DESIGN.md section 2):
    widths(mH, tanb)        -> dict(channel -> Gamma [GeV])  10 channels
    widths_and_br(mH,tanb)  -> (widths dict, BR dict)        BR sums to 1
    gamma_total(mH,tanb)    -> Gamma_tot [GeV]
    ctau(mH,tanb)           -> c*tau [m] = hbar*c / Gamma_tot
    batch(points)           -> list of row dicts, ONE process invocation
Channel names (fixed by the contract):
    gaga, ee, mumu, tautau, pipi, KK, 4pi, cc, ss, gg

NOTE on mH = 2 GeV exactly: the author program switches from hadronic
(pi pi, K K, 4pi) to partonic (s sbar, g g) descriptions with *strict* cuts
(hadronic channels are floored for mH > 2, partonic ones are zero for mH < 2),
so at exactly mH = 2.0 BOTH sets are evaluated and Gamma_tot is too large
(+4.9 % at tanb=10, +43 % at tanb=5, +85 % at tanb=1; negligible for tanb>=100).
This wrapper is faithful to the program and returns that value; consumers should
not use the mH=2.0 node as-is but interpolate across it (the grid nodes 1.81693
and 2.07163 bracket it) or use the mH -> 2 limit from either side, which agree
to <0.5 % for tanb >= 5.
Executable contract: stdin "<mH> <tanb>" per line; stdout 13 whitespace-separated
numbers per line (%.10e):
    mH tanb GammaTot Gamma_gaga Gamma_ee Gamma_mumu Gamma_tautau
            Gamma_pipi Gamma_KK Gamma_4pi Gamma_cc Gamma_ss Gamma_gg

CLI:
    python decay_widths.py --mass-grid <file> --tanb-grid "5:1e5:40" --out <csv>
    python decay_widths.py --points <file with two columns> [--json out.json]
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys
import numpy as np

HBARC_GEV_M = 1.973269804e-16          # DESIGN.md section 0
CHANNELS = ("gaga", "ee", "mumu", "tautau", "pipi", "KK", "4pi", "cc", "ss", "gg")
_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_EXE_CANDIDATES = (
    os.path.join(_HERE, "build", "light_scalar_decay.exe"),
    os.path.join(_HERE, "build", "light_scalar_decay"),
)


def find_executable():
    env = os.environ.get("LIGHT_SCALAR_DECAY_EXE")
    cands = ([env] if env else []) + list(DEFAULT_EXE_CANDIDATES)
    for c in cands:
        if c and os.path.isfile(c):
            return c
    raise FileNotFoundError(
        "batch executable not found; looked for: " + ", ".join(str(c) for c in cands))


def _key(mH, tanb):
    return (float(mH), float(tanb))


class DecayWidths:
    """batched, cached access to the author executable."""

    def __init__(self, exe=None):
        self.exe = exe or find_executable()
        self._cache = {}          # (mH,tanb) -> tuple(13 floats)
        self.ncalls = 0
        self.npoints = 0

    # ---------------------------------------------------------------- raw
    def _run(self, points):
        """one subprocess invocation for a list of (mH,tanb)."""
        if not points:
            return
        inp = "".join("%.17g %.17g\n" % (m, t) for m, t in points)
        self.ncalls += 1
        self.npoints += len(points)
        p = subprocess.run([self.exe], input=inp, capture_output=True, text=True)
        if p.returncode != 0:
            raise RuntimeError("decay executable failed rc=%d stderr=%s" % (p.returncode, p.stderr[-2000:]))
        rows = [ln.split() for ln in p.stdout.strip().splitlines() if ln.strip()]
        if len(rows) != len(points):
            raise RuntimeError("decay executable returned %d rows for %d points" % (len(rows), len(points)))
        for (m, t), r in zip(points, rows):
            if len(r) != 13:
                raise RuntimeError("bad row (%d cols): %s" % (len(r), r))
            self._cache[(m, t)] = tuple(float(x) for x in r)

    # ------------------------------------------------------------- public
    def batch(self, points):
        """evaluate a list of (mH,tanb); one process call for all uncached points."""
        pts = [_key(m, t) for m, t in points]
        need = [p for p in dict.fromkeys(pts) if p not in self._cache]
        self._run(need)
        return [self._row(p) for p in pts]

    def _row(self, p):
        r = self._cache[p]
        mH, tanb, Gtot = r[0], r[1], r[2]
        g = dict(zip(CHANNELS, r[3:13]))
        return dict(mH=mH, tanb=tanb, gamma_total=Gtot, ctau=HBARC_GEV_M / Gtot,
                    br={k: (v / Gtot) for k, v in g.items()}, widths=g)

    def row(self, mH, tanb):
        p = _key(mH, tanb)
        if p not in self._cache:
            self._run([p])
        return self._row(p)

    def widths(self, mH, tanb):
        return dict(self.row(mH, tanb)["widths"])

    def widths_and_br(self, mH, tanb):
        r = self.row(mH, tanb)
        return dict(r["widths"]), dict(r["br"])

    def gamma_total(self, mH, tanb):
        return self.row(mH, tanb)["gamma_total"]

    def ctau(self, mH, tanb):
        return self.row(mH, tanb)["ctau"]


_DEFAULT = None
def default():
    global _DEFAULT
    if _DEFAULT is None:
        _DEFAULT = DecayWidths()
    return _DEFAULT

def widths(mH, tanb):        return default().widths(mH, tanb)
def widths_and_br(mH, tanb): return default().widths_and_br(mH, tanb)
def gamma_total(mH, tanb):   return default().gamma_total(mH, tanb)
def ctau(mH, tanb):          return default().ctau(mH, tanb)
def batch(points, exe=None): return (DecayWidths(exe) if exe else default()).batch(points)


# ------------------------------------------------------------------- CLI
def _read_grid(path):
    out = []
    for ln in open(path):
        ln = ln.split("#")[0].strip()
        if ln:
            out.append(float(ln.split()[0]))
    return out


def _grid_spec(spec):
    parts = spec.split(":")
    if len(parts) == 3:
        a, b, n = float(parts[0]), float(parts[1]), int(parts[2])
        return list(np.logspace(np.log10(a), np.log10(b), n))
    return [float(x) for x in spec.split(",")]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mass-grid", help="file with one mH per line (comments allowed)")
    ap.add_argument("--tanb-grid", default="5:1e5:40", help="'lo:hi:n' (log) or comma list")
    ap.add_argument("--points", help="file with '<mH> <tanb>' per line instead of grids")
    ap.add_argument("--out", help="output CSV")
    ap.add_argument("--json", help="output JSON")
    ap.add_argument("--exe")
    a = ap.parse_args(argv)
    dw = DecayWidths(a.exe)
    if a.points:
        pts = []
        for ln in open(a.points):
            ln = ln.split("#")[0].strip()
            if not ln:
                continue
            m, t = ln.split()[:2]
            pts.append((float(m), float(t)))
    else:
        masses = _read_grid(a.mass_grid)
        tanbs = _grid_spec(a.tanb_grid)
        pts = [(m, t) for m in masses for t in tanbs]
    rows = dw.batch(pts)
    if a.out:
        with open(a.out, "w") as f:
            f.write("# 2HDM light-H decay widths from the AUTHOR Light_scalar_decay program "
                    "(Type-I, cos(ba)=1/tanb, mA=mHC=600, lambdav=0)\n")
            f.write("# channels order: " + ",".join(CHANNELS) + "\n")
            f.write("mH,tanb,GammaTot_GeV,ctau_m,ctau_m_hbarc2e-16," +
                    ",".join("Gamma_%s_GeV" % c for c in CHANNELS) + "," +
                    ",".join("BR_%s" % c for c in CHANNELS) + "\n")
            for r in rows:
                f.write("%.10g,%.10g,%.10e,%.10e,%.10e," % (r["mH"], r["tanb"], r["gamma_total"],
                        r["ctau"], 2.0e-16 / r["gamma_total"]))
                f.write(",".join("%.10e" % r["widths"][c] for c in CHANNELS) + ",")
                f.write(",".join("%.10e" % r["br"][c] for c in CHANNELS) + "\n")
    if a.json:
        json.dump([dict(mH=r["mH"], tanb=r["tanb"], gamma_total=r["gamma_total"], ctau=r["ctau"],
                        widths=r["widths"], br=r["br"]) for r in rows],
                  open(a.json, "w"), indent=1)
    if not a.out and not a.json:
        for r in rows:
            print("%.6g %.6g %.6e %.6e" % (r["mH"], r["tanb"], r["gamma_total"], r["ctau"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
