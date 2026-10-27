Payload = Annotated[FileSecret, {"description": "Immutable synthetic boundary-test inputs and detector package"}]

def blueprint(payload: Payload):
    safe = payload.replace("'", "'\\''")
    entry = """set -e
exec 2>&1
echo 'BOOTSTRAP: install isolated SDK'
python3 -m pip install --disable-pip-version-check --target .test-sdk 'magnus-sdk==0.8.2' --quiet
export PYTHONPATH="$PWD/.test-sdk"
echo 'BOOTSTRAP: download immutable payload'
python3 -c "import magnus; magnus.download_file(file_secret='SECRET', target_path='test_input', overwrite=False)"
echo 'BOOTSTRAP: execute detector tests'
python3 test_input/run_boundary_job.py
""".replace('SECRET', safe)
    submit_job(
        task_name="Downstream synthetic boundary test",
        description="Synthetic HepMC reader/transport/ROOT tests; no scalar model or sensitivity calculation",
        namespace="HET-AGI", repo_name="ColliderAgent", branch="main",
        commit_sha="1140f39e8730889422a64a141fbd3ca10529e13b",
        entry_command=entry, runner="magnus", job_type=JobType.A2,
        container_image="docker://git.pku.edu.cn/het-agi/collider@sha256:1bf6d8db2cb9348f965c3c6650f1c6af8234c1d8be0b49040ddd55e79a33c1a5",
        cpu_count=2, memory_demand="4G", time_limit=10,
    )
