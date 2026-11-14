
"""Fetch a finished llp-job output tarball and show its log.

usage (inside llp-magnus container):  python fetch_job.py <tag>
Reads  work/magnus/jobs/<tag>_result.json  -> downloads the file secret ->
extracts into work/magnus/jobs/<tag>_extract/ and prints job.log + file listing.
"""
import json, os, sys, tarfile
import magnus

MG = "/workspace/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream/work/magnus"
tag = sys.argv[1]
rec = json.load(open(os.path.join(MG, "jobs", tag + "_result.json")))
res = rec["raw"]
if isinstance(res, str):
    res = json.loads(res)
print("RESULT:", res)
if not isinstance(res, dict) or "output" not in res:
    raise SystemExit("no output secret")
sec = res["output"]
tgt = os.path.join(MG, "jobs", tag + "_out.tar.gz")
p = magnus.download_file(sec, target_path=tgt)
ex = os.path.join(MG, "jobs", tag + "_extract")
with tarfile.open(p) as t:
    names = [m.name for m in t.getmembers()]
    t.extractall(ex)
print("tarball:", p, os.path.getsize(p), "members:", names[:20])
out = os.path.join(ex, "out")
if not os.path.isdir(out):
    print("NO out/ DIR"); raise SystemExit(1)
for f in sorted(os.listdir(out)):
    fp = os.path.join(out, f)
    print("  %-28s %10d" % (f, os.path.getsize(fp)))
logp = os.path.join(out, "job.log")
if os.path.exists(logp):
    print("===== job.log (tail) =====")
    print(open(logp).read()[-8000:])
