Process = Annotated[FileSecret, {
    "placeholder": "114514-apple-banana-cat",
    "description": "已编译 process 目录的传输密钥（madgraph-compile 的产出）",
}]

Commands = Annotated[str, {
    "placeholder": "done\nset nevents 10000\ndone",
    "description": "MG5 launch 指令（set 命令、shower/detector 开关、done）",
    "multi_line": True,
}]

Output = Annotated[str, {
    "placeholder": "e.g. tmp/pp_ttbar",
    "description": "事件输出目录的下载路径",
}]

Pdf = Annotated[Optional[str], {
    "placeholder": "LUXlep-NNPDF31_nlo_as_0118_luxqed",
    "description": "LHAPDF PDF 集名称（可选，从 CERN 下载安装）",
}]

Interactive = Annotated[Optional[bool], {
    "description": "使用 MG5 launch -i 模式，对已有事例文件执行后续处理（如 Pythia/Delphes）",
}]

def blueprint(
    process: Process,
    commands: Commands,
    output: Output,
    pdf: Pdf = None,
    interactive: Interactive = None,
):
    description = "## MadGraph5 事件生成任务"

    safe_process = process.replace("'", "'\\''")
    safe_commands = commands.replace("'", "'\\''")
    safe_target = output.replace("'", "'\\''")

    pdf_arg = ""
    if pdf is not None:
        safe_pdf = pdf.replace("'", "'\\''")
        pdf_arg = f" --pdf '{safe_pdf}'"

    interactive_arg = ""
    if interactive:
        interactive_arg = " --interactive"

    system_entry_command = """
# No system entry command needed.
"""

    entry_command = f"""
  source /opt/root/bin/thisroot.sh
export PYTHONUSERBASE=/tmp/magnus-pip-$MAGNUS_JOB_ID
export PIP_CACHE_DIR=/tmp/magnus-pip-cache-madgraph-launch
pip3 install "magnus-sdk==0.8.2" --quiet
python3 scripts/run_madgraph_launch.py --process_secret '{safe_process}' --launch_commands '{safe_commands}' --target_path '{safe_target}'{pdf_arg}{interactive_arg}
"""

    submit_job(
        task_name="[Blueprint] MadGraph Launch",
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
