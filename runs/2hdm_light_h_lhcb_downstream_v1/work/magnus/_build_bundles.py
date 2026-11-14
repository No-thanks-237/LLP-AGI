
import os, json, hashlib, subprocess, tarfile

RUN = r"E:\AAAwork\LLP-AGI\runs\20261025_1325_UTC_2HDM_lightH_lhcb_downstream"
MG = os.path.join(RUN, "work", "magnus")
PIK = os.path.join(MG, "bundles_stage_pik")
HLB = os.path.join(MG, "bundles_stage_hadronlib")
BUND = os.path.join(MG, "bundles")

script_pik = open(os.path.join(PIK, "script.sh.tpl")).read()
script_hlb = open(os.path.join(HLB, "script.sh.tpl")).read()
with open(os.path.join(PIK, "script.sh"), "w", newline="\n") as fh:
    fh.write(script_pik)
with open(os.path.join(HLB, "script.sh"), "w", newline="\n") as fh:
    fh.write(script_hlb)

def tar_bundle(stage, name):
    out = os.path.join(BUND, name)
    if os.path.exists(out):
        os.remove(out)
    r = subprocess.run(["tar", "-czf", out, "-C", stage, "."], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("tar failed: " + r.stderr)
    return out

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

res = {}
for stage, name in ((PIK, "pik.tar.gz"), (HLB, "hadronlib.tar.gz")):
    p = tar_bundle(stage, name)
    with tarfile.open(p) as t:
        members = sorted(m.name for m in t.getmembers())
    res[name] = {"path": p, "sha256": sha256(p), "bytes": os.path.getsize(p),
                 "stage": os.path.basename(stage), "members": members,
                 "staged_files": {f: sha256(os.path.join(stage, f)) for f in sorted(os.listdir(stage))
                                  if not f.endswith(".tpl")}}
print(json.dumps(res, indent=1))
jobs = os.path.join(MG, "jobs")
prev = json.load(open(os.path.join(jobs, "bundles.json")))
prev.update(res)
json.dump(prev, open(os.path.join(jobs, "bundles.json"), "w"), indent=1)
