
import json, os, shutil, tarfile, pathlib
import magnus
sec = "magnus-secret:41641-touch-odds-mood"
tgt = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\magnus\hadronlib3_out.tar.gz"
p = magnus.download_file(sec, target_path=tgt)
print("downloaded:", p, os.path.getsize(p))
with tarfile.open(p) as t:
    print([m.name for m in t.getmembers()][:12])
    t.extractall(r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\magnus\hadronlib3_extract")
d = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3\work\magnus\hadronlib3_extract\out"
print(sorted(os.listdir(d)))
print(open(os.path.join(d, "job_status.json")).read())
log = open(os.path.join(d, "job.log")).read()
print(log[-1200:])
