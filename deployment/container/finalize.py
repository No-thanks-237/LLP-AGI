"""Record current deployment evidence after validation; run inside service."""
import datetime
import hashlib
import json
import pathlib
import shutil
import subprocess
import sys
import urllib.request
import magnus

root=pathlib.Path('/workspace')
here=root/'deployment/container'
run=pathlib.Path(json.loads((here/'smoke-run.json').read_text())['run'])
def docker(*args):return subprocess.check_output(['/usr/local/bin/docker',*args],text=True)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
now=datetime.datetime.now(datetime.timezone.utc)
evidence={'utc':now.isoformat(),'beijing':now.astimezone(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
    'physics_events':0,'physical_end_to_end_validated':False}
evidence['backend']=json.loads(urllib.request.urlopen('http://127.0.0.1:8017/health').read())
evidence['frontend_status']=urllib.request.urlopen('http://127.0.0.1:3011/jobs').status
assert evidence['backend']=={'status':'ok'} and evidence['frontend_status']==200
registration=json.loads((here/'registered/registration.json').read_text())
schemas={}
for name in registration['blueprints']:
    code=magnus.get_blueprint(name)['code']
    assert hashlib.sha256(code.encode()).hexdigest()==registration['sha256'][name]
    schemas[name]=len(magnus.get_blueprint_schema(name))
evidence['blueprint_schemas']=schemas
evidence['tools']={}
for tool in ['collider','mma','micromegas']:
    output=json.loads((run/'outputs'/f'{tool}.json').read_text())
    assert output['success'] and output['physics_events']==0
    evidence['tools'][tool]=output
queue=json.loads((run/'logs/concurrency.json').read_text())
assert queue['observed_serial_queue'] and queue['snapshots'][-1]['states']==['Success','Success']
evidence['serial_queue']=queue['job_ids']
evidence['successful_tool_jobs']={'collider':'95484066ce9fa228','mma':'e4f8420cceeadef9','micromegas':'4de590e0cfdcd166'}
job_inputs=run/'inputs/jobs'
job_inputs.mkdir(exist_ok=True)
for job in magnus.list_jobs(limit=50)['items']:
    detail=magnus.get_job(job['id'])
    snapshot=job_inputs/(job['id']+'.json')
    if not snapshot.exists():snapshot.write_text(json.dumps(detail,indent=2))
    log_path=run/'logs'/(job['id']+'_stdout.json')
    if not log_path.exists():log_path.write_text(json.dumps(magnus.get_job_logs(job['id']),indent=2))
(run/'logs/collider_job.json').write_text(json.dumps({'state':'rejected_before_job_creation',
    'job_id':None,'reason':'HTTP 400: unescaped dictionary braces in deployment probe f-string; no job created'},indent=2))
unit=subprocess.check_output([sys.executable,str(here/'check_guard.py')],text=True)
evidence['guard_unit']=json.loads(unit)
for job in [*evidence['successful_tool_jobs'].values(),*queue['job_ids']]:
    assert magnus.get_job(job)['status']=='Success'
evidence['resource_limits']={}
audit=pathlib.Path('/var/lib/docker/volumes/llp-agi-magnus-20261016/_data/runtime/resource-audit')
for job in ['95484066ce9fa228',*queue['job_ids']]:
    caps=json.loads((audit/f'magnus-job-{job}.json').read_text())
    assert caps=={'Memory':8589934592,'MemorySwap':8589934592,'NanoCpus':4000000000,'PidsLimit':1024}
    evidence['resource_limits'][job]=caps
sources={}
for name,commit in [('ColliderAgent','1140f39e8730889422a64a141fbd3ca10529e13b'),('Magnus','6baa667853412397aab4631e2ac184d0a79a4c91')]:
    path=root/'deployment/sources'/name
    assert subprocess.check_output(['git','-C',str(path),'rev-parse','HEAD'],text=True).strip()==commit
    assert subprocess.check_output(['git','-C',str(path),'status','--porcelain'],text=True)==''
    sources[name]=commit
images={}
for name in ['local/llp-magnus:20261016','git.pku.edu.cn/het-agi/collider:latest','git.pku.edu.cn/het-agi/mma-het:latest','git.pku.edu.cn/rise-agi/micromegas:latest','local/wolfram-engine:13.3-clean']:
    info=json.loads(docker('image','inspect',name))[0]
    images[name]={'id':info['Id'],'digests':info.get('RepoDigests',[]),'size_bytes':info['Size']}
service=json.loads(docker('inspect','llp-magnus'))[0]
assert service['Image']==images['local/llp-magnus:20261016']['id']
assert service['HostConfig']['Memory']==4294967296
assert service['HostConfig']['NanoCpus']==2000000000
assert service['NetworkSettings']['Networks']['llp-agi-services']['IPAddress']=='172.30.11.2'
shutil.copyfile('/opt/llp/python-environment.txt',here/'python-environment.txt')
(run/'logs/deployment_validation.json').write_text(json.dumps(evidence,indent=2,ensure_ascii=False))
files=[p for p in here.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='smoke-run.json']
files+=[root/'README.md',root/'docs/container_deployment.md',root/'configs/magnus_container.private.json',root/'deployment/.dockerignore']
lock={'schema':1,'utc':evidence['utc'],'beijing':evidence['beijing'],'sources':sources,'images':images,
    'service_image':service['Image'],'validation_run':str(run.relative_to(root)),
    'official_blueprints':registration,'physical_end_to_end_validated':False,
    'files':{str(p.relative_to(root)):sha(p) for p in files}}
lock['runtime']={'memory':service['HostConfig']['Memory'],'nano_cpus':service['HostConfig']['NanoCpus'],
    'ip':'172.30.11.2','mounts':{m['Destination']:{k:m[k] for k in ['Type','Source','RW']} for m in service['Mounts']},
    'port_bindings':service['HostConfig']['PortBindings'],
    'path':next(e for e in service['Config']['Env'] if e.startswith('PATH='))}
(root/'configs/container_deployment.lock.json').write_text(json.dumps(lock,indent=2,ensure_ascii=False))
print(json.dumps({'recorded':True,'run':str(run),'official_blueprints':len(schemas),'physics_events':0}))
