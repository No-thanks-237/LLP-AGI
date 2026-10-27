
"""Interface checks that the response module will rely on.

1. work/analysis/grids/mass_grid.txt x channel: for every scan mass where the channel is OPEN, the hadron
   library must have a node within the consumer tolerance (bc5_lib NODE_MTOL, default 0.6 % in
   mass).  Reports the nearest-node distance for every (mass, channel) pair.
2. numpy.loadtxt (the exact call in bc5_lib.HadronLib) parses each library file, and the
   (mass*1e6 + sample) key/group bookkeeping gives n_attempt = 1500 for every channel.
3. minbias_pik.txt.gz: on-shell mass reconstruction quality as a function of energy (the 10-column
   %8g text format limits e^2-p^2 for high-energy rows; quantify how many rows are affected).
"""
import gzip, json, math, os
import numpy as np

RUN = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
COMMON = os.path.join(RUN, "work", "common")
MG = os.path.join(RUN, "work", "magnus")
HL = os.path.join(MG, "hadron_lib")
out = {}

# channel opening thresholds used by the decay module (DESIGN 2/5)
THRESH = {"gg": 2.0, "ss": 2.0, "cc": 2.0 * 3.729686 / 2.0, "tautau": 2.0 * 1.7768}
THRESH["cc"] = 3.729686
LIBMAP = {"gg": "mass_list_gg.txt", "ss": "mass_list_ss.txt", "cc": "mass_list_cc.txt",
          "tautau": "mass_list_tau.txt"}
grid = [float(l) for l in open(os.path.join(COMMON, "mass_grid.txt")) if l.strip() and not l.startswith("#")]
NODE_MTOL = 0.006  # bc5_lib tolerance (relative mass difference accepted without error)

nodes = {}
for ch, f in LIBMAP.items():
    nodes[ch] = np.array([float(l) for l in open(os.path.join(COMMON, f)) if l.strip() and not l.startswith("#")])

prob = []
for ms in grid:
    for ch in nodes:
        if ms < THRESH[ch]:
            continue
        n = nodes[ch]
        rel = float(np.min(np.abs(n - ms)) / max(ms, 1e-30))
        if rel > NODE_MTOL:
            prob.append({"mass": ms, "channel": ch, "nearest_rel": rel,
                         "nearest": float(n[np.argmin(np.abs(n - ms))])})
out["open_channel_nearest_node_problems"] = prob
out["open_channel_pairs_checked"] = sum(1 for ms in grid for ch in nodes if ms >= THRESH[ch])

# 2) parse with the exact bc5_lib.HadronLib code path
parse = {}
for ch in nodes:
    p = os.path.join(HL, "hadron_lib_%s.txt.gz" % ch)
    d = np.loadtxt(p, comments='#')
    if d.ndim == 1:
        d = d.reshape(1, -1)
    mass = d[:, 0]; samp = d[:, 1].astype(np.int64)
    key = mass * 1e6 + samp
    _, start, counts = np.unique(key, return_index=True, return_counts=True)
    order = np.argsort(start)
    parse[ch] = {"rows": int(d.shape[0]), "cols": int(d.shape[1]),
                 "n_attempt": int(samp.max()) + 1, "n_groups": int(len(start)),
                 "mass_min": float(mass.min()), "mass_max": float(mass.max())}
out["loadtxt_parse"] = parse

# 3) minbias mass reconstruction vs energy
p = os.path.join(MG, "samples", "minbias_pik.txt.gz")
bins = [(0, 10), (10, 50), (50, 100), (100, 200), (200, 500), (500, 1e9)]
stat = {("%g-%g" % b): [0, 0] for b in bins}   # [n, n_ok]
with gzip.open(p, "rt") as fh:
    for line in fh:
        if line[0] == "#":
            continue
        v = line.split()
        idd = int(v[1]); px = float(v[2]); py = float(v[3]); pz = float(v[4]); e = float(v[5])
        m = math.sqrt(max(e * e - px * px - py * py - pz * pz, 0.0))
        ref = {211: 0.13957, 321: 0.49368, 311: 0.49761}[abs(idd)]
        for b in bins:
            if b[0] <= e < b[1]:
                k = "%g-%g" % b
                stat[k][0] += 1
                if abs(m - ref) < 0.01 * ref + 0.002:
                    stat[k][1] += 1
                break
out["minbias_mass_reconstruction_by_energy"] = {k: {"rows": v[0], "on_shell_ok_fraction": v[1] / float(v[0]) if v[0] else None}
                                                 for k, v in stat.items()}
json.dump(out, open(os.path.join(MG, "validation_interface.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
