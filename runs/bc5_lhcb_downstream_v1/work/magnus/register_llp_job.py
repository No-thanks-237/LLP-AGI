"""Register the project-local llp-job blueprint inside the service container."""
import hashlib, json, pathlib, sys
import magnus

src = pathlib.Path('/workspace/runs/20261016_144605_UTC_bc5_lhcb_downstream/work/magnus/llp_job_blueprint.py')
code = src.read_text()
record = {
    'blueprint_id': 'llp-job',
    'title': 'LLP generic job (pinned collider image)',
    'description': 'Run a shell bundle in the pinned collider image; returns an output tarball file secret.',
    'source_file': str(src.relative_to('/workspace')),
    'sha256': hashlib.sha256(code.encode()).hexdigest(),
}
magnus.save_blueprint('llp-job', record['title'], record['description'], code)
actual = magnus.get_blueprint('llp-job')
assert actual['code'] == code, 'roundtrip mismatch'
schema = magnus.get_blueprint_schema('llp-job')
record['parameters'] = list(schema.keys()) if isinstance(schema, dict) else schema
out = pathlib.Path('/workspace/runs/20261016_144605_UTC_bc5_lhcb_downstream/work/magnus/llp_job_registration.json')
out.write_text(json.dumps(record, indent=2))
print(json.dumps(record, indent=2))
