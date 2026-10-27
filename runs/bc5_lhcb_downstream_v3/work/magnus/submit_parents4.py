
import json, pathlib, sys
import magnus
which = sys.argv[1]
tag = {'higgs4': 'bc5_higgs4', 'minbias4': 'bc5_minbias4'}[which]
B = pathlib.Path(r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\magnus\bundles\%s.tar.gz" % which)
sec = magnus.custody_file(str(B), expire_minutes=600, max_downloads=3)
print("bundle secret:", sec, flush=True)
res = magnus.run_blueprint('llp-job', args={'bundle': sec, 'job_name': tag},
                           expire_minutes=600, max_downloads=3, timeout=3600,
                           poll_interval=8.0, execute_action=False)
print("RESULT:", res, flush=True)
open(r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\magnus\%s_result.json" % tag, 'w').write(json.dumps({'raw': res}, indent=1))
