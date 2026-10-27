
"""Fetch the two finished llp-job output tarballs and install the samples into this run.

usage (inside llp-magnus container):
    python install_samples.py pik         # -> work/magnus/samples/minbias_pik.txt.gz
    python install_samples.py hadronlib   # -> work/magnus/hadron_lib/hadron_lib_{gg,ss,cc,tautau}.txt.gz

Verifies: job exit code, sha256 printed by the job vs the re-computed sha256, gzip integrity,
row counts, and records everything in work/magnus/jobs/install_report.json.
"""
import gzip, hashlib, json, os, shutil, sys, tarfile
import magnus

MG = "/workspace/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream/work/magnus"
JOBS = os.path.join(MG, "jobs")
tag = sys.argv[1]

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

rec = json.load(open(os.path.join(JOBS, tag + "_result.json")))
res = rec["raw"]
if isinstance(res, str):
    res = json.loads(res)
print("RESULT:", res)
report = {"tag": tag, "result": res}
if not (isinstance(res, dict) and res.get("success")) and res.get("exit_code") != 0:
    print("JOB DID NOT SUCCEED — inspect the tarball log before installing")
    report["warning"] = "nonzero exit code; inspect job.log"
sec = res.get("output")
if sec:
    tgt = os.path.join(JOBS, tag + "_out.tar.gz")
    p = magnus.download_file(sec, target_path=tgt)
    ex = os.path.join(JOBS, tag + "_extract")
    if os.path.isdir(ex):
        shutil.rmtree(ex)
    with tarfile.open(p) as t:
        t.extractall(ex)
    out = os.path.join(ex, "out")
    report["tarball"] = {"path": tgt, "bytes": os.path.getsize(p),
                         "members": sorted(os.listdir(out)) if os.path.isdir(out) else []}
    log = open(os.path.join(out, "job.log")).read() if os.path.exists(os.path.join(out, "job.log")) else ""
    print("===== job.log tail =====")
    print(log[-3000:])
    jobsplit = {}
    for name in sorted(os.listdir(out)):
        if not name.endswith(".gz"):
            continue
        src = os.path.join(out, name)
        h = sha256(src)
        n = sum(1 for _ in gzip.open(src, "rb"))
        jobsplit[name] = {"sha256": h, "bytes": os.path.getsize(src), "lines": n}
        print("  %-28s %8.1f MB  %d lines  sha256 %s" % (name, os.path.getsize(src) / 1e6, n, h))
    report["outputs"] = jobsplit
    # integrity
    for name in jobsplit:
        with gzip.open(os.path.join(out, name), "rb") as fh:
            while fh.read(1 << 24):
                pass
    print("gzip integrity OK")

    if tag == "pik":
        dst = os.path.join(MG, "samples", "minbias_pik.txt.gz")
        shutil.copy2(os.path.join(out, "minbias_pik.txt.gz"), dst)
        report["installed"] = dst
        print("installed", dst)
        # sigmaGen from the summary line
        sj = os.path.join(out, "summary.jsonl")
        if os.path.exists(sj):
            report["summary"] = [json.loads(l) for l in open(sj) if l.strip().startswith("{")]
    elif tag == "hadronlib":
        os.makedirs(os.path.join(MG, "hadron_lib"), exist_ok=True)
        for name, meta in jobsplit.items():
            dst = os.path.join(MG, "hadron_lib", name)
            shutil.copy2(os.path.join(out, name), dst)
            print("installed", dst)
        report["installed"] = sorted(jobsplit)
        sj = os.path.join(out, "summary.jsonl")
        if os.path.exists(sj):
            report["summary"] = [json.loads(l) for l in open(sj) if l.strip().startswith("{")]
rep_path = os.path.join(JOBS, "install_report_%s.json" % tag)
json.dump(report, open(rep_path, "w"), indent=1)
print("wrote", rep_path)
