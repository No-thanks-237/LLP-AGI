Project = Annotated[FileSecret, {
    "placeholder": "114514-apple-banana-cat",
    "description": "已编译的 micrOmegas 项目目录（micromegas-compile 的产出）",
}]

TargetPath = Annotated[str, {
    "placeholder": "tmp/dm_run_output",
    "description": "运行输出目录的下载路径（含 results.json 与 main.c 写出的附加文件）",
}]

Slha = Annotated[Optional[FileSecret], {
    "placeholder": "114514-apple-banana-cat",
    "description": "SLHA 参数卡传输密钥（若 main.c 用 slhaRead 读参数时填写）",
}]

ExtraArgs = Annotated[Optional[str], {
    "placeholder": "e.g. 700 0.1",
    "description": "传给 ./main 的额外命令行参数（空格分隔，需 main.c 支持）",
}]

def blueprint(
    project: Project,
    output: TargetPath,
    slha: Slha = None,
    extra_args: ExtraArgs = None,
):
    description = "## micrOmegas 计算任务"

    safe_project = project.replace("'", "'\\''")
    safe_target = output.replace("'", "'\\''")

    slha_arg = ""
    if slha is not None:
        safe_slha = slha.replace("'", "'\\''")
        slha_arg = f" --slha_secret '{safe_slha}'"

    extra_arg = ""
    if extra_args is not None and extra_args.strip():
        safe_extra = extra_args.replace("'", "'\\''")
        extra_arg = f" --extra_args '{safe_extra}'"

    system_entry_command = """
# No system entry command needed.
"""

    entry_command = f"""
python3 -m pip install --quiet "magnus-sdk==0.8.2"
  export PYTHONPATH=$PWD
python3 scripts/run_micromegas_calc.py --project_secret '{safe_project}' --target_path '{safe_target}'{slha_arg}{extra_arg}
"""

    submit_job(
        task_name="[Blueprint] micrOmegas Calc",
        description=description,
        namespace="HET-AGI",
        repo_name="ColliderAgent",
        branch="main",
        commit_sha="1140f39e8730889422a64a141fbd3ca10529e13b",
        system_entry_command=system_entry_command,
        entry_command=entry_command,
        container_image="docker://git.pku.edu.cn/rise-agi/micromegas@sha256:070498baf78f4efc29a48ba67d0eaec2730af006b386a2d7c9047f65a123e62a",
        job_type=JobType.A2,
        runner="magnus",
        cpu_count=4,
        memory_demand="8G",
    )
