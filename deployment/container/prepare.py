"""Prepare current container manifests; never reads or prints license content."""
import json
import pathlib
import subprocess

root = pathlib.Path(__file__).resolve().parents[2]
here = root / 'deployment/container'
def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')
def digest(image):
    return json.loads(subprocess.check_output(['docker','image','inspect',image]))[0]['RepoDigests'][0]
images = {f'docker://{name}:latest': 'docker://' + digest(name+':latest') for name in [
    'git.pku.edu.cn/het-agi/collider', 'git.pku.edu.cn/het-agi/mma-het', 'git.pku.edu.cn/rise-agi/micromegas']}
pins = dict(collideragent_commit='1140f39e8730889422a64a141fbd3ca10529e13b',
    images=images, license_bind='/var/lib/docker/volumes/llp-agi-wolfram-license-20261016/_data:/root/.WolframEngine/Licensing')
dump(here/'blueprint-pins.json', pins)
cfgpath = root/'configs/magnus_container.private.json'
cfg = json.loads(cfgpath.read_text())
cfg['cluster'].update(max_cpu_count=4, max_memory_demand='8G', default_cpu_count=4,
    default_memory_demand='8G', default_container_image='docker://'+digest('python:3.14-slim-bookworm'))
dump(cfgpath,cfg)
dump(here/'docker-client.json', {'proxies': {'default': {
    'httpProxy':'http://http.docker.internal:3128','httpsProxy':'http://http.docker.internal:3128',
    'noProxy':'localhost,127.0.0.1,172.30.11.2,host.docker.internal'}}})
data='/var/lib/docker/volumes/llp-agi-magnus-20261016/_data'
service_image=json.loads(subprocess.check_output(['docker','image','inspect','local/llp-magnus:20261016']))[0]['Id']
services={'services':{'magnus': {
    'image':service_image,'pull_policy':'never','container_name':'llp-magnus',
    'restart':'unless-stopped','mem_limit':'4g','cpus':2,
    'environment': {'PATH':'/opt/llp/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin'},
    'ports':['127.0.0.1:8017:8017','127.0.0.1:3011:3011'],
    'networks':{'services':{'ipv4_address':'172.30.11.2'}},
    'volumes': [
        '/var/run/docker.sock:/var/run/docker.sock',
        f'llp-data:{data}', str(root).replace('\\','/')+':/workspace',
        '../../configs/magnus_container.private.json:/opt/magnus/configs/magnus_config.yaml:ro',
        './blueprint-pins.json:/opt/llp/blueprint-pins.json:ro',
        './register.py:/opt/llp/register.py:ro',
        './docker-guard:/opt/llp/bin/docker:ro',
        './docker-client.json:/root/.docker/config.json:ro'
    ]}},
    'networks':{'services':{'external':True,'name':'llp-agi-services'}},
    'volumes':{'llp-data':{'external':True,'name':'llp-agi-magnus-20261016'}}}
dump(here/'compose.json',services)
print('Prepared pinned manifests (private config not printed).')
