"""Project-local generic execution blueprint: run a shell bundle in the pinned collider image.

This is the controlled compile/response entry point required by docs/downstream_import.md.
The upstream ColliderAgent blueprints are left untouched; this registers an additional
project-local blueprint under the same pinned image digest.

Input : a .tar.gz ("bundle") containing script.sh (entry point) plus any other inputs.
Output: a .tar.gz of <work>/out (log, status json and whatever the script produced),
        returned as a Magnus file secret.
"""
Bundle = Annotated[FileSecret, {
    "description": "tar.gz containing script.sh and any input files for the job",
}]

JobName = Annotated[str, {
    "placeholder": "bc5-pythia-minbias",
    "description": "short label; also the output tarball name",
}]


def blueprint(bundle: Bundle, job_name: JobName):
    safe = job_name.replace("'", "'\\''")
    entry_command = (
        "python3 -m pip install --quiet magnus-sdk==0.8.2\n"
        "python3 - <<'LLP_JOB'\n"
        "import json, os, pathlib, shutil, subprocess, tarfile\n"
        "import magnus\n"
        "work = pathlib.Path('/magnus/workspace/job')\n"
        "if work.exists(): shutil.rmtree(work)\n"
        "(work / 'out').mkdir(parents=True)\n"
        "b = magnus.download_file(" + repr(bundle) + ", target_path='/magnus/workspace/bundle.tar.gz')\n"
        "with tarfile.open(b) as t: t.extractall(work)\n"
        "log = work / 'out' / 'job.log'\n"
        "rc = 1\n"
        "with open(log, 'w') as fh:\n"
        "    p = subprocess.run(['bash', str(work / 'script.sh')], cwd=str(work), stdout=fh, stderr=subprocess.STDOUT)\n"
        "    rc = p.returncode\n"
        "(work / 'out' / 'job_status.json').write_text(json.dumps({'exit_code': rc, 'job_name': " + repr(job_name) + "}, indent=2))\n"
        "tar_path = '/magnus/workspace/' + " + repr(safe) + " + '-out.tar.gz'\n"
        "with tarfile.open(tar_path, 'w:gz') as t: t.add(str(work / 'out'), arcname='out')\n"
        "secret = magnus.custody_file(tar_path, max_downloads=3)\n"
        "pathlib.Path(os.environ['MAGNUS_RESULT']).write_text(json.dumps({'success': rc == 0, 'exit_code': rc, 'output': secret}))\n"
        "LLP_JOB\n"
    )
    submit_job(
        task_name="[Blueprint] LLP generic job",
        description="Pinned-collider execution entry point for the LLP-AGI BC5 LHCb downstream analysis",
        namespace="HET-AGI",
        repo_name="ColliderAgent",
        branch="main",
        commit_sha="1140f39e8730889422a64a141fbd3ca10529e13b",
        system_entry_command="",
        entry_command=entry_command,
        container_image="docker://git.pku.edu.cn/het-agi/collider@sha256:1bf6d8db2cb9348f965c3c6650f1c6af8234c1d8be0b49040ddd55e79a33c1a5",
        job_type=JobType.A2,
        runner="magnus",
        cpu_count=4,
        memory_demand="8G",
    )
