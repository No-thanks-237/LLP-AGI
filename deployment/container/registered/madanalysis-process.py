Events = Annotated[FileSecret, {
    "placeholder": "114514-apple-banana-cat",
    "description": "事件目录的传输密钥（madgraph-launch 的产出）",
}]

Script = Annotated[str, {
    "placeholder": "import ... as sample\nset sample.type = signal\nplot PT(mu+) 50 0 500",
    "description": "MA5 分析脚本（不含 submit）",
    "multi_line": True,
}]

Output = Annotated[str, {
    "placeholder": "e.g. tmp/ma5_output",
    "description": "分析结果的下载路径",
}]

Level = Annotated[Literal["parton", "hadron", "reco"], {
    "description": "分析级别",
    "options": {
        "parton": {"label": "parton", "description": "LHE"},
        "hadron": {"label": "hadron", "description": "HepMC"},
        "reco": {"label": "reco", "description": "LHCO/ROOT"},
    },
}]

def blueprint(
    events: Events,
    script: Script,
    output: Output,
    level: Level = "parton",
):
    description = "## MadAnalysis5 分析任务"

    safe_events = events.replace("'", "'\\''")
    safe_script = script.replace("'", "'\\''")
    safe_target = output.replace("'", "'\\''")
    safe_level = level.replace("'", "'\\''")

    system_entry_command = """
# No system entry command needed.
"""

    entry_command = f"""
  source /opt/root/bin/thisroot.sh
pip3 install "magnus-sdk==0.8.2" --quiet
python3 scripts/run_madanalysis_process.py --events_secret '{safe_events}' --script '{safe_script}' --level '{safe_level}' --target_path '{safe_target}'
"""

    submit_job(
        task_name="[Blueprint] MadAnalysis Process",
        description=description,
        namespace="HET-AGI",
        repo_name="ColliderAgent",
        branch="main",
        commit_sha="1140f39e8730889422a64a141fbd3ca10529e13b",
        system_entry_command=system_entry_command,
        entry_command=entry_command,
        container_image="docker://git.pku.edu.cn/het-agi/collider@sha256:1bf6d8db2cb9348f965c3c6650f1c6af8234c1d8be0b49040ddd55e79a33c1a5",
        job_type=JobType.A2,
        runner="magnus",
        cpu_count=4,
        memory_demand="8G",
    )
