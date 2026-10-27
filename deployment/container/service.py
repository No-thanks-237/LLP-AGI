"""Container listener adapter; the official Magnus source is left unchanged."""
import sys
sys.argv = ["magnus-container", "--deliver", "--config", "/opt/magnus/configs/magnus_config.yaml"]
sys.path.insert(0, "/opt/magnus/back_end")
import uvicorn
from server.main import app
from server._scheduler._submit import _SubmitMixin
import subprocess
import threading
# Official local scheduling does not reserve RAM/CPU across sibling jobs.
# Serialize managed jobs on this memory-limited workstation; retain Pending state.
_submission_lock = threading.Lock()
_original_submit = _SubmitMixin._submit_to_docker
def _serial_submit(self, db, job):
    with _submission_lock:
        active = subprocess.check_output(['/usr/local/bin/docker','ps','--filter','name=magnus-job-','--format','{{.Names}}'],text=True)
        if any(name.startswith('magnus-job-') for name in active.splitlines()):
            return False
        return _original_submit(self, db, job)
_SubmitMixin._submit_to_docker = _serial_submit
uvicorn.run(app, host="0.0.0.0", port=8017)
