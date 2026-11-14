"""Download and unpack the higgs4 sample job output (one-off helper)."""
import os
import pathlib
import tarfile

import magnus

HERE = pathlib.Path(__file__).resolve().parent
sec = "magnus-secret:49057-west-lure-ozone"
tgt = HERE / "higgs4_out.tar.gz"
ex = HERE / "higgs4_extract"
p = magnus.download_file(sec, target_path=str(tgt))
with tarfile.open(p) as t:
    t.extractall(str(ex))
    print([m.name for m in t.getmembers()])
print((ex / "out" / "job.log").read_text()[-2500:])
