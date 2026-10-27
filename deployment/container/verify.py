"""Read-only current deployment check. No registration, submissions or repairs."""
import hashlib
import json
import pathlib
import subprocess
import urllib.request
import magnus
root=pathlib.Path('/workspace')
lock=json.loads((root/'configs/container_deployment.lock.json').read_text())
fail=[]
for name,expected in lock['files'].items():
    p=root/name
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=expected:fail.append(name)
for name,commit in lock['sources'].items():
    p=root/'deployment/sources'/name
    actual=subprocess.check_output(['git','-C',str(p),'rev-parse','HEAD'],text=True).strip()
    dirty=subprocess.check_output(['git','-C',str(p),'status','--porcelain'],text=True)
    if actual!=commit or dirty:fail.append('source:'+name)
service=json.loads(subprocess.check_output(['/usr/local/bin/docker','inspect','llp-magnus']))[0]
if service['Image']!=lock['service_image']:fail.append('service image')
runtime=lock['runtime']
if service['HostConfig']['Memory']!=runtime['memory'] or service['HostConfig']['NanoCpus']!=runtime['nano_cpus']:fail.append('service resources')
mounts={m['Destination']:{k:m[k] for k in ['Type','Source','RW']} for m in service['Mounts']}
if mounts!=runtime['mounts']:fail.append('service mounts')
if service['NetworkSettings']['Networks'].get('llp-agi-services',{}).get('IPAddress')!=runtime['ip']:fail.append('service network')
if service['HostConfig']['PortBindings']!=runtime['port_bindings']:fail.append('service ports')
if runtime['path'] not in service['Config']['Env']:fail.append('Docker adapter PATH')
for name,expected in lock['official_blueprints']['sha256'].items():
    actual=hashlib.sha256(magnus.get_blueprint(name)['code'].encode()).hexdigest()
    if actual!=expected:fail.append('blueprint:'+name)
assert json.loads(urllib.request.urlopen('http://127.0.0.1:8017/health').read())=={'status':'ok'}
assert urllib.request.urlopen('http://127.0.0.1:3011/jobs').status==200
print(json.dumps({'read_only':True,'ok':not fail,'mismatches':fail,'physical_validation':False},ensure_ascii=False))
raise SystemExit(bool(fail))
