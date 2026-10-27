
"""Submit the two T3 Pythia samples to the local Magnus service (serial queue).

Launch order matters: the long minbias job is launched first, the short hadron-library job second,
so the queue never idles.  Both are launched immediately (launch_blueprint returns a job id), then
polled to completion.  Results are written to work/magnus/jobs/<tag>_result.json in the same
{raw: <result string>} format used by the BC5 run, so fetch_job.py works unchanged.
"""
import json, os, pathlib, time, datetime
import magnus

MG = "/workspace/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream/work/magnus"
JOBS = os.path.join(MG, "jobs")
BUNDLES = os.path.join(MG, "bundles")

def utc():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

def launch(tag, bundle_name, job_name, expire, maxdl):
    b = os.path.join(BUNDLES, bundle_name)
    import hashlib
    h = hashlib.sha256()
    with open(b, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    sec = magnus.custody_file(b, expire_minutes=expire, max_downloads=maxdl)
    jid = magnus.launch_blueprint('llp-job', args={'bundle': sec, 'job_name': job_name},
                                  expire_minutes=expire, max_downloads=maxdl, timeout=120.0)
    rec = {"tag": tag, "job_id": jid, "job_name": job_name, "bundle": bundle_name,
           "bundle_sha256": h.hexdigest(), "bundle_secret": sec,
           "expire_minutes": expire, "max_downloads": maxdl, "launched_utc": utc()}
    print("launched", tag, jid, flush=True)
    return rec

state = {"started_utc": utc(), "jobs": []}
state["jobs"].append(launch("pik", "pik.tar.gz", "2hdm_pik", 1440, 3))
pathlib.Path(JOBS, "submitted.json").write_text(json.dumps(state, indent=1))
state["jobs"].append(launch("hadronlib", "hadronlib.tar.gz", "2hdm_hadronlib", 600, 3))
pathlib.Path(JOBS, "submitted.json").write_text(json.dumps(state, indent=1))

# poll: report status changes; each job's clock starts at its own launch record
for j in state["jobs"]:
    jid = j["job_id"]
    t0 = time.time()
    last = None
    deadline = t0 + 7200
    while True:
        st = magnus.get_job(jid)["status"]
        if st != last:
            print("%s %s -> %s (t+%.0fs)" % (utc(), jid, st, time.time() - t0), flush=True)
            last = st
        if st in ("Success", "Failed", "Terminated", "Cancelled"):
            break
        if time.time() > deadline:
            print("POLL TIMEOUT", jid, flush=True)
            break
        time.sleep(10.0)
    j["final_status"] = st
    j["observed_seconds"] = round(time.time() - t0, 1)
    j["finished_utc"] = utc()
    try:
        res = magnus.get_job_result(jid)
    except Exception as e:  # noqa: BLE001
        res = json.dumps({"error": str(e)})
    print("%s RESULT %s: %s" % (utc(), jid, res), flush=True)
    pathlib.Path(JOBS, j["tag"] + "_result.json").write_text(json.dumps({"raw": res}, indent=1))
    pathlib.Path(JOBS, "submitted.json").write_text(json.dumps(state, indent=1))

print("ALL DONE", json.dumps(state, indent=1), flush=True)
