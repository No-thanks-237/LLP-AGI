CalcHEP = Annotated[FileSecret, {
    "placeholder": "114514-apple-banana-cat",
    "description": "CalcHEP 模型目录的传输密钥（generate-calchep 的产出，含 Z2-odd 标记的 prtcls1.mdl）",
}]

Main = Annotated[FileSecret, {
    "placeholder": "114514-apple-banana-cat",
    "description": "用户 main.c（调用 micrOmegas API，约定向 results.json 写出结构化结果）",
}]

TargetPath = Annotated[str, {
    "placeholder": "tmp/my_dm_project",
    "description": "编译完成的 micrOmegas 项目目录下载路径",
}]

ProjectName = Annotated[Optional[str], {
    "placeholder": "my_dm_project",
    "description": "micrOmegas 项目名（newProject 使用）；留空则用 'dm_project'",
}]

def blueprint(
    calchep: CalcHEP,
    main: Main,
    output: TargetPath,
    project: ProjectName = None,
):
    description = "## micrOmegas 项目编译任务"

    safe_calchep = calchep.replace("'", "'\\''")
    safe_main = main.replace("'", "'\\''")
    safe_target = output.replace("'", "'\\''")

    project_arg = ""
    if project is not None and project.strip():
        safe_project = project.replace("'", "'\\''")
        project_arg = f" --project '{safe_project}'"

    system_entry_command = """
# No system entry command needed.
"""

    entry_command = f"""
python3 -m pip install --quiet "magnus-sdk==0.8.2"
  export PYTHONPATH=$PWD
python3 scripts/run_micromegas_compile.py --calchep_secret '{safe_calchep}' --main_secret '{safe_main}' --target_path '{safe_target}'{project_arg}
"""

    submit_job(
        task_name="[Blueprint] micrOmegas Compile",
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
