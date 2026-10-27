# ColliderAgent 子代理入口

以 Codex 子代理运行官方 ColliderAgent 技能；这是宿主代理适配，不能表述为官方原生 Codex 角色自动加载。

先读取根 README、当前用户任务，以及 `deployment/sources/ColliderAgent/src/agents/` 中对应角色和 `src/skills/` 中对应 SKILL.md。上游仓库保持原样，版本以当前部署记录为准。

Magnus 操作使用 `deployment/container/magnus.ps1`，由新服务容器中的官方 SDK 执行。不要使用旧 `deployment/magnus.ps1` 或旧 WSL Magnus，也不要因断线重交已记录的作业。

所有科学计算必须有当前任务授权和清晰输入；部署验证仅允许非物理冒烟测试。不得自动恢复历史 TASK、运行数据或派生输入。每次尝试使用独立运行目录，按 stage 保存输入、日志、job ID、产物和验收记录；作业 Success 与物理验收分别判断。

本轮只重新部署并验证基础接口，不生成模型、事件、探测器响应或灵敏度结果。不要求额外 LLM API key，不启动独立 ADK Agent。
