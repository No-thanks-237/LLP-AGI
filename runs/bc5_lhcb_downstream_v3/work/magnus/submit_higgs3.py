
import json, pathlib
import magnus
B = pathlib.Path(r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\magnus\bundles\higgs3.tar.gz")
sec = magnus.custody_file(str(B), expire_minutes=360, max_downloads=3)
print("bundle secret:", sec, flush=True)
res = magnus.run_blueprint('llp-job', args={'bundle': sec, 'job_name': 'bc5_higgs3'},
                           expire_minutes=360, max_downloads=3, timeout=1800,
                           poll_interval=5.0, execute_action=False)
print("RESULT:", res, flush=True)
open(r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\magnus\higgs3_result.json", 'w').write(json.dumps({'raw': res}, indent=1))
