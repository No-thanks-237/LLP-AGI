"""Magnus job entry: synthetic fixtures only, with outputs retained on failure."""
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shlex
import subprocess
import tarfile
import magnus

inp=Path('test_input').resolve()
out=Path('test_output').resolve()
out.mkdir()
report={'scope':'synthetic_boundary_events_only','model_events':0,'sensitivity_computed':False,
        'tests':[], 'runtime_packages': {d.metadata['Name']:d.version for d in importlib.metadata.distributions()}}
env=os.environ.copy()
def run(command, log):
    p=subprocess.run(command,env=env,capture_output=True,text=True)
    (out/log).write_text(p.stdout+'\n'+p.stderr)
    if p.returncode:
        raise RuntimeError('Command failed: '+command[0]+' (see '+log+')')
    return p.stdout
try:
    archive=inp/'downstream-delphes.tar.gz'
    expected=(inp/'package.sha256').read_text().strip()
    assert hashlib.sha256(archive.read_bytes()).hexdigest()==expected, 'Detector package SHA mismatch'
    package=Path('detector_package').resolve()
    package.mkdir()
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            dest=(package/member.name).resolve()
            if not dest.is_relative_to(package) or member.issym() or member.islnk():
                raise RuntimeError('Unsafe archive member')
        tar.extractall(package)
    env['LD_LIBRARY_PATH']=str(package/'lib')+':'+env.get('LD_LIBRARY_PATH','')
    source=Path('/opt/MG5_aMC_v3_7_0/Delphes')
    # ROOT dictionaries embed source includes; the relocated package must provide
    # both the unchanged base headers and new module headers to Cling.
    env['ROOT_INCLUDE_PATH']=':'.join([str(package/'source'),str(source),str(source/'external')])
    flags=shlex.split(subprocess.check_output(['root-config','--cflags'],text=True))
    libs=shlex.split(subprocess.check_output(['root-config','--libs'],text=True))
    checker=Path('check_boundary_output').resolve()
    run(['g++',*flags,'-I'+str(source),'-I'+str(source/'external'),
         str(inp/'check_boundary_output.cpp'),'-L'+str(package/'lib'),'-lDelphes',
         *libs,'-lEG','-o',str(checker)],'compile_checker.log')
    fixture=inp/'fixture'
    spec=json.loads((fixture/'expected.json').read_text())
    report.update({'seed':spec['seed'],'events_per_file':spec['event_count'],
                   'fixture_sha256':spec['hashes'],'package_sha256':expected})
    geometry_card=package/'cards/delphes_card_downstream_geometry_only.tcl'
    if not geometry_card.exists():
        # Compatibility with historical geometry-only packages; alter only a
        # fresh job-local card, never the archived package or source inputs.
        geometry_card=package/'cards/geometry_only_override.tcl'
        # ExRootConfReader's source is relative to the top card directory.
        geometry_card.write_text('source delphes_card_downstream.tcl\n'
                                 'set Downstream::MinBirthEnergy 0.0\n'
                                 'set Downstream::MinBirthPT 0.0\n')
    report['kinematic_mode']='geometry_only_birth_thresholds_disabled'
    for eventfile in sorted(fixture.glob('*.hepmc')):
        name=eventfile.stem
        version=3 if 'hepmc3' in name else 2
        rootfile=out/(name+'.root')
        record={'input':eventfile.name,'reader_version':version}
        try:
            run([str(package/'bin'/f'DelphesHepMC{version}'),
                 str(geometry_card),str(rootfile),str(eventfile)],name+'_reader.log')
            record['checker_stdout']=run([str(checker),str(rootfile),str(fixture/'expected.tsv')],name+'_check.log')
            record['status']='passed'
        except Exception as e:
            record.update(status='failed',error=str(e))
        report['tests'].append(record)
    report['success']=bool(report['tests']) and all(t['status']=='passed' for t in report['tests'])
except Exception as e:
    report.update(success=False,error=str(e))
finally:
    (out/'test_report.json').write_text(json.dumps(report,indent=2)+'\n')
    secret=magnus.custody_file(str(out),expire_minutes=1440,max_downloads=2)
    Path(os.environ['MAGNUS_ACTION']).write_text('magnus receive '+secret+' --output boundary-results\n')
    Path(os.environ['MAGNUS_RESULT']).write_text(json.dumps({'success':report.get('success',False),
        'passed_files':sum(t['status']=='passed' for t in report['tests']),
        'total_files':len(report['tests']),'scope':report['scope']}))
# Preserve diagnostics through Magnus even when validation fails. Job completion
# means the test runner returned artifacts; test_report.success is the verdict.
print(json.dumps({'validation_success':report.get('success',False),'files':len(report['tests'])}),flush=True)
