
import json, pathlib, sys, time
import magnus
B = pathlib.Path(r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\magnus\bundles\hadronlib3.tar.gz")
sec = magnus.custody_file(str(B), expire_minutes=360, max_downloads=3)
print("bundle secret:", sec)
args = {'bundle': sec, 'job_name': 'bc5_hadronlib3'}
res = magnus.run_blueprint('llp-job', args=args, expire_minutes=360, max_downloads=3,
                           timeout=2400, poll_interval=5.0, execute_action=False)
print("RESULT:", res)
open(r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\magnus\hadronlib3_result.json", 'w').write(json.dumps({'raw': res}, indent=1))
