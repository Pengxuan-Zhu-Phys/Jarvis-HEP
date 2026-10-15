# Jarvis-HEP 2.0.17

1. **Roadmap phase**：Redis 环境引导和标准 CLI 命令组完成后的 HEP 2.0.17 提交、候选验收与正式发布准备。
2. **Affected modules**：Redis 环境/缺失程序诊断、CLI 路由/help、agent catalog/search、对应测试及维护文档；HEP 版本声明。
3. **Public interfaces**：新增正式 `Jarvis redis info/check/install/start` 子命令、独立 help、JSON 和 agent 发现；裸命令继续默认 info。安装/启动指引只打印，YAML 结构和扫描算法不变。
4. **Main risks**：误纳入工作区 Sampler/runtime/Diver 改动；以本地 SDK 安装冒充公开索引可安装；所需 `Jarvis-HEP-Sampler>=0.1.0,<0.2` 尚未公开，不能先发布依赖不可安装的 HEP。
5. **Verification performed**：待本阶段测试结果补录；此前 Redis 202 项定向、25 项安装 CLI 和 fast gate 结果属于开发阶段历史证据，不能替代本次发行源码/安装包验收。

## 授权与发行范围

2026-10-15 用户明确授权提交 Redis 改进并发布 HEP 新版本，同时要求忽略当前 Sampler 改动。本版本为 2.0.17；既有 2.0.15/2.0.16 候选和历史记录保留。

本次提交/候选只纳入 Redis 环境引导、正式子命令及帮助/agent 发现、对应测试与文档、HEP 两处版本声明。已有 committed baseline 中的 sampler SDK 集成仍是依赖关系；本轮未提交的 sampler/runtime 修补、Diver 源码/测试、Diver package-data 与 sampler 审阅资料均排除并保留。`pyproject.toml` 只纳入本次版本变更，不带入其它未验收工作区差异。

## 用户可见变化

- 首次 CLI 使用后缺少 Redis-compatible executable 时，stderr 给出系统适用的安装/启动指引；它是 PATH 检查，不是 pip 安装钩子。
- `Jarvis redis` / `info` 展示本机系统、executable/version 与说明；`check` 每次执行 PATH 和有界版本检查，缺失、不可用/过旧或 native Windows 退出 1。
- `install` / `start` 只打印指导，不执行系统安装、服务操作、网络/端口探测或扫描。普通扫描仍由 Jarvis 自动管理本次 Redis；Windows 引导在 WSL 内安装/运行 HEP 与 Redis。
- Redis 是正式命令组：父级 help 列出四个 Commands，子级有独立说明和 `--json`；JSON 可在子命令前后使用，裸命令保留 info 默认，agent 公开 exact child commands。
- 安装建议按发行版/版本处理 Redis 或 Valkey；旧 Ubuntu/Debian 使用已核验的官方 signed APT recipe，RHEL 8 说明 `redis:6` 与 updated AppStream 要求。

## 本次验收与源码状态

本阶段应先完成 Redis-only 提交，再从该干净 committed source 构建 wheel/sdist；核对内容没有混入工作区未授权改动，并在独立环境安装全部声明依赖。执行 `pip check`、installed Redis help/JSON/agent gate、真实 Redis/spawn smoke 和发行回归，保存源码 commit、命令、结果、产物 checksum 与文件范围。

验收状态：待本阶段测试结果补录。提交 SHA、push 状态、产物路径及 draft Release URL 均待实际操作结果补录；此处不声称已经提交、推送或创建 Release。

开发阶段历史证据为 [Redis 环境引导验收](redis-environment-verification.json) 和 [标准命令组验收](redis-cli-help-verification.json)。后者 202 项定向、25 项安装检查通过；工作区 fast gate 为 1509 passed / 478 subtests，保留 1 项既有外部 iDM 描述编码失败。这些结果不证明新的版本源码或发行包已经通过。

## 公开依赖与发布状态

2026-10-15 公共 PyPI 查询快照：HEP 最新公开版本为 2.0.14，Operas 为 1.5.0，`Jarvis-HEP-Sampler` 项目查询返回 HTTP 404。查询入口分别为 [HEP JSON](https://pypi.org/pypi/Jarvis-HEP/json)、[Operas JSON](https://pypi.org/pypi/Jarvis-Operas/json)、[Sampler JSON](https://pypi.org/pypi/Jarvis-HEP-Sampler/json)；这些是当日观测，发布前须重新核验。

既有 baseline 声明 `Jarvis-HEP-Sampler>=0.1.0,<0.2`，本次不改变该依赖。允许使用已验收且保持原样的本地 SDK 0.1.0 wheel 验证候选，其余依赖按正常完整安装解析；这属于候选验收，不是仅依赖公开索引即可安装的证明。

正式 HEP 发布仍等待公开依赖可用。可以准备 draft Release 及完整验收证据；本次 Redis-only 授权不包括修改、提交或发布 SDK。依赖问题未解决前，不将 draft/本地候选标记为公开发布完成，不上传不可正常安装的 HEP。后续公开发行沿用 GitHub Release/Trusted Publishing。
