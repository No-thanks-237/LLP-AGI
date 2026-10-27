
import os, json, pathlib, time
MG = "/workspace/runs/20261025_1325_UTC_2HDM_lightH_lhcb_downstream/work/magnus"
import magnus

B = pathlib.Path(MG) / "bundles" / "smoke.tar.gz"
sec = magnus.custody_file(str(B), expire_minutes=240, max_downloads=3)
print("bundle secret:", sec, flush=True)
res = magnus.run_blueprint('llp-job', args={'bundle': sec, 'job_name': '2hdm_smoke'},
                           expire_minutes=240, max_downloads=3, timeout=1200,
                           poll_interval=5.0, execute_action=False)
print("RESULT:", res, flush=True)
pathlib.Path(MG, "jobs", "smoke_result.json").write_text(json.dumps({'raw': res}, indent=1))
