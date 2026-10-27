
import json
d = json.load(open(r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\review_check\kaon\f04_results.json"))
for m, r in d.items():
    print("m=%s BR_KK=%.3f" % (m, r['br_KK']))
    for edge in ('long','short'):
        tot = 0.0
        line = []
        for ch in ('B','K','H'):
            if ch in r.get(edge, {}):
                v = r[edge][ch]
                line.append("%s: n=%d pc=%.5f pa=%.5f ratio=%.4f dss=%.4f" % (
                    ch, v['n'], v['pc'], v['pa'], v['pa']/v['pc'] if v['pc'] else float('nan'), v['dss']))
        print("   %-5s %s" % (edge, " | ".join(line)))
