
import json, os, shutil, tarfile, gzip, sys
import magnus

BASE = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\magnus"
SAMP = os.path.join(BASE, "samples")

def fetch(tag):
    rec = json.load(open(os.path.join(BASE, tag + "_result.json")))
    res = rec["raw"]
    if isinstance(res, str):
        res = json.loads(res)
    if not res.get("success"):
        raise SystemExit("%s failed: %r" % (tag, res))
    sec = res["output"]
    tgt = os.path.join(BASE, tag + "_out.tar.gz")
    p = magnus.download_file(sec, target_path=tgt)
    ex = os.path.join(BASE, tag + "_extract")
    with tarfile.open(p) as t:
        names = [m.name for m in t.getmembers()]
        t.extractall(ex)
    print(tag, "->", names, os.path.getsize(p))
    out = os.path.join(ex, "out")
    print(open(os.path.join(out, "summary.jsonl")).read().strip()[-300:])
    return out

which = sys.argv[1]
out = fetch(which)
if which == "bc5_higgs4":
    src = os.path.join(out, "higgs.txt.gz")
    dst = os.path.join(SAMP, "higgs.txt.gz")
    if os.path.exists(dst) and not os.path.exists(os.path.join(SAMP, "higgs_v34_backup.txt.gz")):
        shutil.copy2(dst, os.path.join(SAMP, "higgs_v34_backup.txt.gz"))
    shutil.copy2(src, dst)
    n = sum(1 for _ in gzip.open(dst, 'rt'))
    print("installed higgs.txt.gz :", n, "rows")
    n62 = sum(1 for l in gzip.open(dst, 'rt') if not l.startswith('#') and l.rstrip().endswith('62'))
    print("rows with status 62:", n62)
elif which == "bc5_minbias4":
    for name in ("minbias_light.txt.gz", "minbias_heavy.txt.gz"):
        src = os.path.join(out, name)
        if not os.path.exists(src):
            print("missing", name); continue
        dst = os.path.join(SAMP, name)
        old = os.path.join(SAMP, name.replace(".txt.gz", ".txt"))
        if os.path.exists(old) and not os.path.exists(old + ".v34"):
            os.rename(old, old + ".v34")
        shutil.copy2(src, dst)
        n = sum(1 for _ in gzip.open(dst, 'rt'))
        print("installed", name, ":", n, "rows")
