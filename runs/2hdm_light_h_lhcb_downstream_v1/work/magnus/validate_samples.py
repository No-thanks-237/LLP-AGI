
"""Validate the T3 samples (host-side, no magnus needed).

    python validate_samples.py            # writes work/magnus/validation.json

Checks (acceptance criteria of task-3):
  * minbias_pik.txt.gz: n_gen from max event id, 10 columns, status 91/92 only (no status column in
    the 10-column format -> checked structurally), all rows are one of {211,-211,321,-321,311,-311},
    rows/event, fraction of events with a pi/K, eta range (no cut), on-shell mass per species,
    charge balance, Pythia sigmaGen vs the 79.5 mb normalisation, n_gen vs normalisation.json;
  * bbar_heavy.txt.gz and higgs.txt.gz: row counts, n_gen, B hadrons/event, status-62 count;
  * hadron_lib_*.txt.gz: distinct mass nodes EXACTLY equal to work/analysis/grids/mass_list_*.txt,
    sample index 0..1499, on-shell masses, charge conservation per (mass,sample) decay,
    rows and empty-sample fraction per node.
"""
import gzip, json, math, os, sys
from collections import Counter, defaultdict

RUN = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # work/magnus -> run
MG = os.path.join(RUN, "work", "magnus")
COMMON = os.path.join(RUN, "work", "common")
SAMP = os.path.join(MG, "samples")
HL = os.path.join(MG, "hadron_lib")
OUT = {}

PID_NAME = {211: "pi+", -211: "pi-", 321: "K+", -321: "K-", 311: "K0", -311: "K0bar"}


def eta(pz, px, py):
    p = math.sqrt(px * px + py * py + pz * pz)
    if p <= abs(pz):
        return math.copysign(1e9, pz)
    return 0.5 * math.log((p + pz) / (p - pz))


def check_pik():
    p = os.path.join(SAMP, "minbias_pik.txt.gz")
    if not os.path.isfile(p):
        return {"present": False}
    r = {"present": True, "file": "samples/minbias_pik.txt.gz"}
    n = 0; cols = Counter(); pdg = Counter(); maxev = -1; evrows = Counter()
    mass_ok = Counter(); mass_bad = Counter(); charge = 0
    eta_abs_max = 0.0; n_eta6 = 0; ev_eta6 = set(); badcol = []
    for line in gzip.open(p, "rt"):
        if line[0] == "#":
            continue
        v = line.split()
        if len(v) != 10:
            if len(badcol) < 5:
                badcol.append(line[:80])
            continue
        n += 1
        ev = int(v[0]); idd = int(v[1]); px = float(v[2]); py = float(v[3]); pz = float(v[4]); e = float(v[5])
        pdg[idd] += 1; evrows[ev] += 1; maxev = max(maxev, ev); charge += 1 if idd in (211, 321) else (-1 if idd in (-211, -321) else 0)
        m = math.sqrt(max(e * e - px * px - py * py - pz * pz, 0.0))
        ref = {211: 0.13957, 321: 0.49368, 311: 0.49761}.get(abs(idd), 0.0)
        if abs(m - ref) < 0.01 * ref + 0.002:
            mass_ok[idd] += 1
        else:
            mass_bad[idd] += 1
        et = eta(pz, px, py)
        eta_abs_max = max(eta_abs_max, abs(et))
        if abs(et) < 6.0:
            n_eta6 += 1; ev_eta6.add(ev)
    r.update({"rows": n, "n_gen_maxev_plus1": maxev + 1, "columns_bad": badcol[:5],
              "pdg": {PID_NAME.get(k, str(k)): v for k, v in sorted(pdg.items())},
              "charge_sum_piK": charge,
              "mass_on_shell_fraction": {PID_NAME.get(k, str(k)): mass_ok[k] / float(mass_ok[k] + mass_bad[k])
                                          for k in pdg},
              "rows_per_event_mean": n / float(maxev + 1),
              "events_with_piK": len(evrows), "events_with_piK_fraction": len(evrows) / float(maxev + 1),
              "rows_abs_eta_lt6": n_eta6, "events_with_piK_eta6": len(ev_eta6),
              "events_with_piK_eta6_fraction": len(ev_eta6) / float(maxev + 1),
              "abs_eta_max": eta_abs_max,
              "multiplicity_per_event": {"min": min(evrows.values()), "max": max(evrows.values()),
                                          "mean": sum(evrows.values()) / float(len(evrows))}})
    norm = json.load(open(os.path.join(SAMP, "normalisation.json")))
    r["n_gen_normalisation"] = norm["inel"]["n_gen"]
    r["n_gen_matches"] = (maxev + 1) == norm["inel"]["n_gen"]
    r["sigma_pb_normalisation"] = norm["inel"]["sigma_pb"]
    return r


def check_reused():
    out = {}
    for name, key in (("bbar_heavy.txt.gz", "B"), ("higgs.txt.gz", "H")):
        p = os.path.join(SAMP, name)
        if not os.path.isfile(p):
            out[key] = {"present": False}
            continue
        n = 0; maxev = -1; pdg = Counter(); n62 = 0; ncol = Counter()
        for line in gzip.open(p, "rt"):
            if line[0] == "#":
                continue
            v = line.split(); n += 1; ncol[len(v)] += 1
            maxev = max(maxev, int(v[0])); pdg[int(v[1])] += 1
            if len(v) >= 11 and abs(int(v[10])) == 62:
                n62 += 1
        norm = json.load(open(os.path.join(SAMP, "normalisation.json")))[key]
        out[key] = {"present": True, "file": name, "rows": n, "n_gen": maxev + 1,
                    "n_gen_matches_normalisation": (maxev + 1) == norm["n_gen"],
                    "sigma_pb": norm["sigma_pb"], "columns": dict(ncol),
                    "rows_per_event": n / float(maxev + 1),
                    "pdg_top": {str(k): v for k, v in pdg.most_common(8)}}
        if key == "H":
            out[key]["status62_rows"] = n62
            out[key]["one_status62_per_event"] = n62 == maxev + 1
    return out


LIB_MASS = {211: 0.13957, 321: 0.49368, 311: 0.49761, 2212: 0.93827, 11: 0.000511,
            13: 0.10566, 2112: 0.93957, 3222: 1.18937, 3112: 1.19745, 3312: 1.32171,
            3334: 1.67245, 15: 1.77686, 3212: 1.19264, 3122: 1.11568, 221: 0.54786,
            2224: 1.231, 1114: 1.234, 223: 0.78266, 113: 0.77526, 213: 0.77526,
            22: 0.0, 111: 0.13498, 130: 0.49761, 310: 0.49761}


def read_mass_list(name):
    p = os.path.join(COMMON, name)
    vals = []
    for line in open(p):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        vals.append(float(s))
    return vals


def check_hadron_lib():
    out = {}
    chans = {"gg": "mass_list_gg.txt", "ss": "mass_list_ss.txt",
             "cc": "mass_list_cc.txt", "tautau": "mass_list_tau.txt"}
    for ch, mfile in chans.items():
        p = os.path.join(HL, "hadron_lib_%s.txt.gz" % ch)
        if not os.path.isfile(p):
            out[ch] = {"present": False}
            continue
        r = {"present": True, "file": "hadron_lib/hadron_lib_%s.txt.gz" % ch}
        n = 0; nodes = Counter(); samples_per_node = defaultdict(set)
        mass_ok = 0; mass_bad = 0; bad_pid = Counter()
        charge_by_group = defaultdict(int); rows_by_group = Counter()
        eta_norm = [0.0, 0.0]
        hdr = None
        for line in gzip.open(p, "rt"):
            if line[0] == "#":
                if hdr is None:
                    hdr = line.strip()
                continue
            v = line.split()
            n += 1
            m0 = float(v[0]); isamp = int(v[1]); pid = int(v[2]); chg = int(v[3])
            e = float(v[4]); px = float(v[5]); py = float(v[6]); pz = float(v[7])
            nodes[round(m0, 6)] += 1
            samples_per_node[round(m0, 6)].add(isamp)
            m = math.sqrt(max(e * e - px * px - py * py - pz * pz, 0.0))
            ref = LIB_MASS.get(abs(pid))
            if ref is None:
                bad_pid[pid] += 1
            elif abs(m - ref) < max(0.02, 0.02 * ref):
                mass_ok += 1
            else:
                mass_bad += 1
            g = (round(m0, 6), isamp)
            charge_by_group[g] += chg
            rows_by_group[g] += 1
        want = [round(x, 6) for x in read_mass_list(mfile)]
        got = sorted(nodes)
        r.update({"header": hdr, "rows": n,
                  "n_mass_nodes": len(got), "mass_nodes_expected": len(want),
                  "mass_nodes_exact_match": got == sorted(want),
                  "mass_nodes_missing": [x for x in want if x not in nodes][:10],
                  "mass_nodes_extra": [x for x in got if x not in want][:10],
                  "sample_index_max": max((max(s) for s in samples_per_node.values()), default=-1),
                  "n_attempt_expected": 1500,
                  "rows_per_node_min": min(nodes.values()), "rows_per_node_max": max(nodes.values()),
                  "rows_per_node_mean": sum(nodes.values()) / float(len(nodes)),
                  "nonzero_charge_groups": sum(1 for g, c in charge_by_group.items() if c != 0),
                  "groups": len(charge_by_group),
                  "empty_attempt_fraction": 1.0 - len(charge_by_group) / float(len(nodes) * 1500),
                  "mass_on_shell_fraction": mass_ok / float(mass_ok + mass_bad),
                  "unknown_pids": dict(bad_pid)})
        out[ch] = r
    return out


if __name__ == "__main__":
    OUT["created"] = "validate_samples.py"
    OUT["pik"] = check_pik()
    OUT["reused"] = check_reused()
    OUT["hadron_lib"] = check_hadron_lib()
    json.dump(OUT, open(os.path.join(MG, "validation.json"), "w"), indent=1)
    print(json.dumps(OUT, indent=1))
