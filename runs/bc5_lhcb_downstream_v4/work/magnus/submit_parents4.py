"""Submit a parent-sample bundle (bc5_higgs4 / bc5_minbias4) to Magnus."""
import json
import pathlib
import sys

import magnus

HERE = pathlib.Path(__file__).resolve().parent
which = sys.argv[1]
B = HERE / "bundles" / ("%s.tar.gz" % which)
sec = magnus.custody_file(str(B), expire_minutes=360, max_downloads=3)
print("bundle secret:", sec)
res = magnus.run_blueprint("llp-job", args={"bundle": sec, "job_name": which},
                           expire_minutes=360, max_downloads=3,
                           timeout=2400, poll_interval=5.0, execute_action=False)
print("RESULT:", res)
(HERE / ("%s_result.json" % which)).write_text(json.dumps({"raw": res}, indent=1))
