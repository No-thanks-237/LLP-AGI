"""Two tiny no-physics jobs validate the local single-job adapter."""
import json
import pathlib
import time
import magnus
run=pathlib.Path(json.loads(pathlib.Path('/workspace/deployment/container/smoke-run.json').read_text())['run'])
record=run/'logs/concurrency.json'
if record.exists():raise SystemExit('Concurrency attempt already exists; inspect its jobs, do not resubmit')
entry='python3 -c \'import os,pathlib,time; time.sleep(12); pathlib.Path(os.environ["MAGNUS_RESULT"]).write_text("nonphysical-queue-ok")\''
data={'job_ids':[], 'snapshots':[], 'observed_serial_queue':False}
record.write_text(json.dumps(data))
pins=json.loads(pathlib.Path('/opt/llp/blueprint-pins.json').read_text())
for i in range(2):
    job=magnus.submit_job(task_name=f'Deployment serial queue {i+1}',entry_command=entry,
        namespace='HET-AGI',repo_name='ColliderAgent',branch='main',commit_sha=pins['collideragent_commit'],
        container_image=pins['images']['docker://git.pku.edu.cn/het-agi/collider:latest'],cpu_count=4,memory_demand='8G')
    data['job_ids'].append(job)
    record.write_text(json.dumps(data,indent=2))
for _ in range(80):
    states=[magnus.get_job(j)['status'] for j in data['job_ids']]
    data['snapshots'].append({'time':time.time(),'states':states})
    if states==['Running','Pending']:data['observed_serial_queue']=True
    assert states.count('Running')<=1,states
    record.write_text(json.dumps(data,indent=2))
    if all(s=='Success' for s in states):break
    if any(s=='Failed' for s in states):raise RuntimeError(states)
    time.sleep(0.5)
assert data['observed_serial_queue'] and states==['Success','Success'],data
print(json.dumps({'job_ids':data['job_ids'],'serial_queue_passed':True}))
