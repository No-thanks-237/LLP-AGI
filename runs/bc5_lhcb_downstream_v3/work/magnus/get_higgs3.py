
import json, os, shutil, tarfile, glob, gzip
import magnus
sec = "magnus-secret:9491-peach-tend-kill"
tgt = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\magnus\higgs3_out.tar.gz"
p = magnus.download_file(sec, target_path=tgt)
print("downloaded:", p, os.path.getsize(p))
ex = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\magnus\higgs3_extract"
with tarfile.open(p) as t:
    print([m.name for m in t.getmembers()][:8])
    t.extractall(ex)
d = os.path.join(ex, "out")
print(sorted(os.listdir(d)))
print(open(os.path.join(d, "job_status.json")).read())
print(open(os.path.join(d, "job.log")).read()[-500:])