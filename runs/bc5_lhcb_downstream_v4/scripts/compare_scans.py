#!/usr/bin/env python
"""Compare two signal_scan.json files value by value (regression / provenance checks).

    python scripts/compare_scans.py outputs/data/signal_scan.json <other>/signal_scan.json

Exit code 0 when the physics payload is identical (metadata such as model_tag and input
fingerprints is allowed to differ), 1 otherwise.
"""
import json
import sys


def _payload(path):
    with open(path) as fh:
        d = json.load(fh)
    return d


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 2:
        raise SystemExit(__doc__)
    a, b = _payload(argv[0]), _payload(argv[1])
    problems = []
    for key in ("masses", "theta2"):
        if a.get(key) != b.get(key):
            problems.append("%s differs" % key)
    if a.get("nsamp_by_mass") != b.get("nsamp_by_mass"):
        print("note: nsamp_by_mass differs (the two scans used different sample sizes)")
    sa, sb = a.get("signal_per_fb_per_channel", {}), b.get("signal_per_fb_per_channel", {})
    if set(sa) != set(sb):
        problems.append("mass keys differ")
    compared = different = 0
    max_rel = 0.0
    for mass in sorted(set(sa) & set(sb)):
        ra, rb = sa[mass], sb[mass]
        if len(ra) != len(rb):
            problems.append("row length differs at m=%s" % mass)
            continue
        for x, y in zip(ra, rb):
            for key in set(x) | set(y):
                compared += 1
                va, vb = x.get(key), y.get(key)
                if isinstance(va, float) and isinstance(vb, float):
                    if va != vb:
                        different += 1
                        max_rel = max(max_rel, abs(va - vb) / max(abs(va), 1e-300))
                elif va != vb:
                    different += 1
    print("masses: %d  compared values: %d  differing: %d  max relative difference: %.3e"
          % (len(sa), compared, different, max_rel))
    for p in problems:
        print("PROBLEM:", p)
    return 0 if (not problems and different == 0) else 1


if __name__ == "__main__":
    sys.exit(main())
