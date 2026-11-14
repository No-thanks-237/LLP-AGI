"""Submit the hadronlib3 bundle to Magnus and record the raw result."""
import json
import pathlib
import sys
import time

import magnus

HERE = pathlib.Path(__file__).resolve().parent
B = HERE / "bundles" / "hadronlib3.tar.gz"
sec = magnus.custody_file(str(B), expire_minutes=360, max_downloads=3)
print("bundle secret:", sec)
args = {"bundle": sec, "job_name": "bc5_hadronlib3"}
res = magnus.run_blueprint("llp-job", args=args, expire_minutes=360, max_downloads=3,
                           timeout=2400, poll_interval=5.0, execute_action=False)
print("RESULT:", res)
(HERE / "hadronlib3_result.json").write_text(json.dumps({"raw": res}, indent=1))
