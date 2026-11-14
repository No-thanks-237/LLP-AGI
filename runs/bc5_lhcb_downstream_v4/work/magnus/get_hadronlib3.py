"""Download and unpack the hadronisation-library job output (one-off helper)."""
import os
import pathlib
import tarfile

import magnus

HERE = pathlib.Path(__file__).resolve().parent
sec = "magnus-secret:41641-touch-odds-mood"
tgt = HERE / "hadronlib3_out.tar.gz"
p = magnus.download_file(sec, target_path=str(tgt))
print("downloaded:", p, os.path.getsize(p))
ex = HERE / "hadronlib3_extract"
with tarfile.open(p) as t:
    print([m.name for m in t.getmembers()][:12])
    t.extractall(str(ex))
d = ex / "out"
print(sorted(os.listdir(d)))
print((d / "job_status.json").read_text())
print((d / "job.log").read_text()[-1200:])
