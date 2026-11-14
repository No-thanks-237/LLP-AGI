# -*- coding: utf-8 -*-
"""Compare two signal_scan.json files point by point (restructure reproducibility check)."""
from __future__ import annotations
import json, math, sys
import numpy as np


def load(path):
    with open(path) as f:
        return json.load(f)


def index(scan):
    out = {}
    for p in scan['points']:
        out[(round(float(p['mH']), 9), round(float(p['tanb']), 9))] = p
    return out


def s_of(p, lumi):
    for k, v in p.get('s', {}).items():
        if abs(float(k) - float(lumi)) < 1e-9:
            return float(v)
    return float('nan')


def ds_of(p, lumi):
    for k, v in p.get('ds', {}).items():
        if abs(float(k) - float(lumi)) < 1e-9:
            return float(v)
    return float('nan')


def main():
    a = load(sys.argv[1]); b = load(sys.argv[2])
    ia, ib = index(a), index(b)
    keys = sorted(set(ia) & set(ib))
    only_a = sorted(set(ia) - set(ib)); only_b = sorted(set(ib) - set(ia))
    print('points: A=%d B=%d common=%d onlyA=%d onlyB=%d' % (len(ia), len(ib), len(keys), len(only_a), len(only_b)))
    for L in (25.0, 300.0):
        diffs = []
        for k in keys:
            sa, sb = s_of(ia[k], L), s_of(ib[k], L)
            if sa > 0 and sb > 0 and np.isfinite(sa) and np.isfinite(sb):
                diffs.append(abs(sb / sa - 1.0))
        d = np.array(diffs)
        print(' L=%g  n=%d  max|ds/s|=%.3e  median=%.3e' % (L, len(d), d.max(), np.median(d)))
    # boundaries
    import csv
    pa, pb = sys.argv[3], sys.argv[4]
    ra = {(round(float(r['mass_GeV']), 9), float(r['lumi_fb'])): r for r in csv.DictReader(open(pa))}
    rb = {(round(float(r['mass_GeV']), 9), float(r['lumi_fb'])): r for r in csv.DictReader(open(pb))}
    bad = 0
    for k in sorted(set(ra) & set(rb)):
        for tag in ('short', 'long'):
            va, vb = ra[k]['tanb_' + tag], rb[k]['tanb_' + tag]
            if (va == '') != (vb == ''):
                bad += 1
            elif va != '':
                rel = abs(float(vb) / float(va) - 1.0)
                if rel > 1e-6:
                    bad += 1
                    print(' boundary differs', k, tag, va, vb)
    print('boundary entries compared: %d, mismatches(>1e-6): %d' % (len(set(ra) & set(rb)) * 2, bad))


if __name__ == '__main__':
    main()
