"""Install reviewed official blueprints with pinned sources and image digests."""
import json
import pathlib
import re
import hashlib
from ruamel.yaml import YAML
import magnus

root = pathlib.Path('/opt/collideragent/src/blueprints')
pins = json.loads(pathlib.Path('/opt/llp/blueprint-pins.json').read_text())
yaml = YAML(typ='safe')
ids = []
prepared = []
for path in sorted(root.glob('*.yaml')):
    doc = yaml.load(path.read_text())
    code = doc['code']
    code = code.replace('commit_sha="HEAD"', 'commit_sha="' + pins['collideragent_commit'] + '"')
    code = code.replace('commit_sha = "HEAD"', 'commit_sha = "' + pins['collideragent_commit'] + '"')
    for old, new in pins['images'].items():
        code = code.replace(old, new)
    code = code.replace('$HOME/.wolfram-container-license:$HOME/.WolframEngine/Licensing', pins['license_bind'])
    # Pin SDK requirements in registration copies, without editing upstream files.
    code = re.sub(r'magnus-sdk>=\d+\.\d+\.\d+', 'magnus-sdk==0.8.2', code)
    code = re.sub(r'cpu_count=\d+', 'cpu_count=4', code)
    code = re.sub(r'memory_demand="\d+G"', 'memory_demand="8G"', code)
    code = code.replace('repo_name="ColliderAgent",', 'repo_name="ColliderAgent",\n        branch="main",')
    if path.stem.startswith('micromegas-'):
        code = code.replace('export PYTHONPATH=$PWD', 'python3 -m pip install --quiet "magnus-sdk==0.8.2"\n  export PYTHONPATH=$PWD')
    if 'collider@sha256:' in code:
        code = code.replace('entry_command = f"""', 'entry_command = f"""\n  source /opt/root/bin/thisroot.sh')
    assert 'HEAD' not in code and ':latest' not in code and 'magnus-sdk>=' not in code, path
    assert pins['collideragent_commit'] in code and '@sha256:' in code, path
    assert '$HOME/.wolfram-container-license' not in code, path
    prepared.append((path, doc, code))
out = pathlib.Path('/workspace/deployment/container/registered')
out.mkdir(exist_ok=True)
for path, doc, code in prepared:
    magnus.save_blueprint(path.stem, doc.get('title', path.stem), doc.get('description', ''), code)
    actual = magnus.get_blueprint(path.stem)
    assert actual['code'] == code, (path, actual.keys())
    (out/(path.stem+'.py')).write_text(code)
    ids.append(path.stem)
(out/'registration.json').write_text(json.dumps({'blueprints': ids, 'sha256': {
    path.stem: hashlib.sha256(code.encode()).hexdigest() for path, doc, code in prepared}}, indent=2))
print(json.dumps({'registered': ids}))
