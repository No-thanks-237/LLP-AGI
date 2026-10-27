"""No Docker calls: validate resource arguments and audit-failure semantics."""
import contextlib
import io
import json
import pathlib
import runpy
import subprocess
import sys
from unittest.mock import patch

guard=pathlib.Path(__file__).with_name('docker-guard')
captured=[]
def launch(args,**kwargs):
    captured.append(args)
    return subprocess.CompletedProcess(args,0,'synthetic-container-id\n','')
output=io.StringIO()
with patch.object(sys,'argv',['docker','run','--name','magnus-job-unit','example']), patch('subprocess.run',launch), patch('subprocess.check_output',side_effect=OSError('synthetic audit failure')), contextlib.redirect_stdout(output), contextlib.redirect_stderr(io.StringIO()):
    try:runpy.run_path(str(guard),run_name='__main__')
    except SystemExit as e:assert e.code==0
assert output.getvalue()=='synthetic-container-id\n'
for flag in ['--memory=8g','--memory-swap=8g','--cpus=4','--pids-limit=1024']:assert flag in captured[0]
with patch.object(sys,'argv',['docker','run','--name','magnus-job-unit','--privileged','example']), patch('subprocess.run') as forbidden:
    try:runpy.run_path(str(guard),run_name='__main__')
    except SystemExit as e:assert 'cannot override' in str(e)
    else:raise AssertionError('Privileged request accepted')
    forbidden.assert_not_called()
print(json.dumps({'resource_argument_test':True,'audit_failure_preserves_started_job':True,'privileged_request_rejected':True,'actual_jobs_submitted':0}))
