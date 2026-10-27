"""Nonphysical deployment validation. Execute in service container, one submit per stage."""
import datetime
import hashlib
import json
import pathlib
import shlex
import sys
import magnus

base=pathlib.Path('/workspace/runs')
pointer=pathlib.Path('/workspace/deployment/container/smoke-run.json')
if not pointer.exists():
    name=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d_%H%M%S_UTC_deployment-smoke')
    run=base/name
    for d in ['inputs','logs','outputs']:(run/d).mkdir(parents=True,exist_ok=True)
    pointer.write_text(json.dumps({'run':str(run)}))
run=pathlib.Path(json.loads(pointer.read_text())['run'])
pins=json.loads(pathlib.Path('/opt/llp/blueprint-pins.json').read_text())
mode=sys.argv[1]
if mode=='register':
    payload=run/'inputs/payload.txt'
    payload.write_text('LLP-AGI deployment roundtrip; no physics events.\n')
    for tool,image in [('collider','git.pku.edu.cn/het-agi/collider'),('mma','git.pku.edu.cn/het-agi/mma-het'),('micromegas','git.pku.edu.cn/rise-agi/micromegas')]:
        python='python3'
        install=('uv pip install --python /usr/local/bin/python3 --break-system-packages --quiet magnus-sdk==0.8.2' if tool=='mma' else 'python3 -m pip install --quiet magnus-sdk==0.8.2')
        if tool=='collider': install='source /opt/root/bin/thisroot.sh\n'+install
        probe='''import hashlib,json,os,pathlib,subprocess,shutil,urllib.request
import magnus
p=magnus.download_file(SECRET,target_path="/magnus/workspace/payload.txt")
assert pathlib.Path(p).read_text()=="LLP-AGI deployment roundtrip; no physics events.\\n"
r={"success":True,"tool":TOOL,"sdk":magnus.__version__,"physics_events":0}
assert json.loads(urllib.request.urlopen(os.environ["MAGNUS_ADDRESS"]+"/health").read())["status"]=="ok"
'''
        if tool=='collider':
            probe+='''import ROOT
r["ROOT"]=ROOT.gROOT.GetVersion()
r["MG5"]=str(pathlib.Path(os.environ["MG5_DIR"])/"bin/mg5_aMC")
assert pathlib.Path(r["MG5"]).is_file()
assert pathlib.Path(os.environ["MG5_DIR"],"Delphes/DelphesHepMC2").is_file()
r["DelphesHepMC2"]=True
r["Pythia8"]=subprocess.check_output(["pythia8-config","--version"],text=True).strip()
'''
        elif tool=='mma':
            probe+='''w='Get["FeynRules`"]; Print["LLP_DEPLOY_OK:",2+2,":",$Version,":",$FeynRulesVersion];'
c=subprocess.run(["wolframscript","-code",w],capture_output=True,text=True,timeout=120)
assert c.returncode==0 and "LLP_DEPLOY_OK:4:" in c.stdout,(c.returncode,c.stdout,c.stderr)
r["wolfram"]=c.stdout
'''
        else:
            probe+='''assert pathlib.Path("/opt/micromegas_6.3.0").is_dir()
r["micromegas"]="6.3.0"
r["compiler"]=subprocess.check_output(["gcc","--version"],text=True).splitlines()[0]
'''
        probe+='''out=pathlib.Path("/magnus/workspace/probe.json")
out.write_text(json.dumps(r,indent=2))
secret=magnus.custody_file(str(out),max_downloads=3)
pathlib.Path(os.environ["MAGNUS_RESULT"]).write_text(json.dumps({"success":True,"output":secret}))
'''
        probe=probe.replace('TOOL',repr(tool))
        # Runtime FileSecret is passed through the blueprint argument.
        encoded=probe.replace('SECRET','{secret!r}')
        entry=install+'\npython3 - <<\'LLP_PROBE\'\n'+encoded+'\nLLP_PROBE'
        entry=entry.replace('{','{{').replace('}','}}').replace('{{secret!r}}','{secret!r}')
        # repr creates safe Python literal; only the secret placeholder is dynamic.
        code='''Payload = Annotated[FileSecret, {"description":"Synthetic deployment input"}]
def blueprint(payload: Payload):
    secret = payload
    entry = ENTRY
    submit_job(task_name="Deployment probe TOOL", namespace="HET-AGI", repo_name="ColliderAgent", branch="main", commit_sha="COMMIT", container_image="IMAGE", cpu_count=4, memory_demand="8G", system_entry_command=SYSTEM, entry_command=entry)
'''
        code=code.replace('ENTRY', 'f'+repr(entry)).replace('TOOL',tool).replace('COMMIT',pins['collideragent_commit']).replace('IMAGE',pins['images']['docker://'+image+':latest'])
        system=('mounts=("'+pins['license_bind']+'")') if tool=='mma' else ''
        code=code.replace('SYSTEM',repr(system))
        (run/'inputs'/f'{tool}_blueprint.py').write_text(code)
        magnus.save_blueprint('deployment-probe-'+tool,'Deployment probe '+tool,'Nonphysical tool and file roundtrip validation',code)
    print(json.dumps({'run':str(run),'registered':True}))
elif mode=='launch':
    tool=sys.argv[2]
    stage=sys.argv[3] if len(sys.argv)>3 else tool
    record=run/'logs'/f'{stage}_job.json'
    if record.exists():raise SystemExit('Already submitted: collect existing job')
    record.write_text(json.dumps({'state':'submission_started','tool':tool}))
    job=magnus.launch_blueprint('deployment-probe-'+tool,{'payload':str(run/'inputs/payload.txt')},save_preference=False,max_downloads=2,timeout=30)
    record.write_text(json.dumps({'job_id':job,'tool':tool},indent=2))
    print(record.read_text())
elif mode=='collect':
    tool=sys.argv[2]
    stage=sys.argv[3] if len(sys.argv)>3 else tool
    record=json.loads((run/'logs'/f'{stage}_job.json').read_text())
    job=magnus.get_job(record['job_id'])
    (run/'logs'/f'{tool}_status.json').write_text(json.dumps(job,indent=2))
    print(json.dumps({'job_id':record['job_id'],'status':job['status']}))
    if job['status']=='Success':
        result=json.loads(magnus.get_job_result(record['job_id']))
        assert result['success']
        out=run/'outputs'/f'{tool}.json'
        magnus.download_file(result['output'],str(out),overwrite=False)
        assert json.loads(out.read_text())['success']
        print(out.read_text())
