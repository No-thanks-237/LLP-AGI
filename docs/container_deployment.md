# 当前容器部署

部署日期：2026-10-16（北京时间）。此文档及 `configs/container_deployment.lock.json` 描述重新部署的环境；旧部署文件和旧锁只作历史材料，不用于恢复本环境。

## 来源与结构

- ColliderAgent 官方仓库：<https://github.com/HET-AGI/ColliderAgent>，固定提交 `1140f39e8730889422a64a141fbd3ca10529e13b`。
- Magnus 官方仓库：<https://github.com/Rise-AGI/magnus>，固定提交 `6baa667853412397aab4631e2ac184d0a79a4c91`。
- 两份全新源码在 `deployment/sources/`；没有恢复首轮模型、UFO、事件或数据库。
- `llp-magnus` 是新服务容器，内含 Magnus 前后端、SDK 和 ColliderAgent 源码／技能。工具作业由它通过 Docker socket 启动独立的兄弟容器；不是把所有科学程序放进常驻服务进程。
- 已获取官方 `collider`、`mma-het`、`micromegas` 工具镜像，注册 8 个官方蓝图的本机适配副本。完整镜像摘要与实际环境见当前锁及 `deployment/container/python-environment.txt`。
- 保留的 `local/wolfram-engine:13.3-clean` 未修改。新工具镜像另外包含 FeynRules、MG5 等程序。

## 启停与使用

在项目根目录使用 PowerShell：

```powershell
./deployment/container/start.ps1
./deployment/container/magnus.ps1 config
./deployment/container/magnus.ps1 blueprint list
docker exec llp-magnus python /workspace/deployment/container/verify.py
```

界面：<http://localhost:3011/jobs>；健康检查：<http://127.0.0.1:8017/health>。Windows 仅绑定回环地址，不开放局域网访问。容器自动重启策略为 `unless-stopped`，需要 Docker Desktop 正在运行。

确认没有运行中的作业后，可用 `./deployment/container/stop.ps1` 停止服务。它不删除卷；不要使用 `down -v`、volume prune 或重跑初始化来解决一般连接问题。服务容器可重建，数据库与作业文件在独立数据卷中持久保存。

CLI 包装器在服务容器内运行新版 SDK。项目内绝对 Windows 路径会转成 `/workspace/...`；相对文件路径建议从项目根目录使用正斜杠，例如 `runs/<run-id>/inputs/model.fr`。不要调用旧 `deployment/magnus.ps1` 或旧 WSL 服务入口。

ColliderAgent 继续由宿主 Codex 的子代理执行。显式让子代理读取 `deployment/container/collideragent_role.md`、当前获准任务及 `deployment/sources/ColliderAgent/src/skills/` 中相关技能。该角色适配文件不自动加载；未额外启动 ADK 或配置 LLM API key。官方源码／技能保存在服务镜像和项目目录中。

## 本机适配与资源边界

官方源码保持原样，适配集中在 `deployment/container/`：

1. `service.py` 将后端监听地址设为 `0.0.0.0`，并在本机串行提交工具作业。官方 local 后端原本允许多容器并行，且不执行申请的 CPU／内存硬限制。
2. `docker-guard` 对 Magnus 作业容器实际添加 4 CPU、8 GiB 内存、8 GiB 内存加交换总量（即禁止额外 swap）、1024 PID 限制。服务容器另限 2 CPU／4 GiB。超过限制仍会 OOM 失败；这些限制不能代替生产样本的分块、流式处理与磁盘预算。
3. 注册副本固定源码提交、镜像 digest、SDK 0.8.2，资源申请统一为 4 CPU／8G。对 collider 入口显式加载 ROOT 环境；micrOmegas 入口将内置旧 SDK 升为固定版本。上游蓝图源码没有改写。已注册代码保存在 `deployment/container/registered/`。
4. 控制面使用 Python 3.14、CPU 版 Torch 2.9.1／Torchvision 0.24.1。保留上游 lock，但构建时移除了其 CUDA 依赖并安装相应 CPU 包，因此不声称是上游完整锁逐字复现。最终镜像 ID 和包清单是本机环境依据；重新构建不保证字节一致。
5. 容器内 Git 将两个指定 GitHub 仓库映射到镜像中的固定源码。依赖安装仍可能联网；不宣称完全离线。

Docker socket 允许服务控制本机 Docker，仅用于可信本地任务；不要把本地免登录服务公开到外网。工具容器使用官方 Linux host 网络，回调固定为 `http://172.30.11.2:8017`。

## 持久化与文件管理

- 网络：`llp-agi-services`，`172.30.11.0/24`；服务固定 `172.30.11.2`。
- 数据卷：`llp-agi-magnus-20261016`。必须按 `/var/lib/docker/volumes/llp-agi-magnus-20261016/_data` **同一绝对路径**挂入服务；这是兄弟作业的宿主绑定源，不能随意改为 `/data`。
- 许可卷：`llp-agi-wolfram-license-20261016`。仅从原许可复制 mathpass；原文件不修改。官方 mma-het 首次运行会在新卷内缓存 FeynArts 补丁。不要输出许可内容。
- 私有配置：`configs/magnus_container.private.json`，包含本地服务密钥，不展示内容或上传公共仓库。恢复备份时一并保留；公开锁只记录其摘要。
- `compose.json` 定义挂载、端口、固定网络和资源。`docker-client.json` 使用 Docker Desktop 代理；代理失效时先定位联网故障，不擅自升级依赖。
- 部署测试在独立 `runs/20261016_051119_UTC_deployment-smoke/`；失败尝试与作业 ID 保留。临时下载／构建记录在 `temp/redeploy_20261016/`。不要把测试当作物理验收。

以后正式计算仍须独立 `runs/<时间戳_任务>/`，输入、代码快照、分阶段作业 ID、日志与结果归属该 run。断线只收集已有作业，不自动重交。旧自定义物理接入及首轮数据未恢复，本次也没有创建模型、生成事件或扫描。

## 验证和恢复

当前锁与 `runs/.../logs/deployment_validation.json` 记录健康检查、注册代码比对、非物理工具探针、文件往返、资源限制和串行队列测试。作业 Success 需要同时核对返回结果与实际输出，不能单独作为科学程序成功的判据。

普通恢复：启动 Docker Desktop，运行 `start.ps1`，再检查健康、蓝图及已有作业。不要重新注册蓝图或重新提交测试来代替检查。

迁机恢复需要保存服务镜像和工具镜像、两个命名卷、项目配置和源码。先导入镜像，再用同名空卷恢复数据和许可，创建同名网络与子网，核对锁后启动 compose。卷内数据库、服务密钥和镜像／注册代码必须对应；不要混用旧 WSL 数据库。跨机器 Wolfram 许可可能需要重新激活，复制文件不保证许可可迁移。


