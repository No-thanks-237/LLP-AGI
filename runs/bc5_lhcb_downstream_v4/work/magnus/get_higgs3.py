"""Download and unpack the higgs3 sample job output (one-off helper)."""
import os
import pathlib
import tarfile

import magnus

HERE = pathlib.Path(__file__).resolve().parent
sec = "magnus-secret:9491-peach-tend-kill"
tgt = HERE / "higgs3_out.tar.gz"
p = magnus.download_file(sec, target_path=str(tgt))
print("downloaded:", p, os.path.getsize(p))
ex = HERE / "higgs3_extract"
with tarfile.open(p) as t:
    print([m.name for m in t.getmembers()][:8])
    t.extractall(str(ex))
d = ex / "out"
print(sorted(os.listdir(d)))
print((d / "job_status.json").read_text())
print((d / "job.log").read_text()[-500:])
