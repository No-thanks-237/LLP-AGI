
import json, os, tarfile
import magnus
sec = "magnus-secret:49057-west-lure-ozone"
tgt = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\magnus\higgs4_out.tar.gz"
p = magnus.download_file(sec, target_path=tgt)
ex = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\magnus\higgs4_extract"
with tarfile.open(p) as t:
    t.extractall(ex)
    print([m.name for m in t.getmembers()])
print(open(os.path.join(ex, "out", "job.log")).read()[-2500:])
