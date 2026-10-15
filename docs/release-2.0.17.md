# Jarvis-HEP 2.0.17

1. **Roadmap phase**：Redis 环境引导和标准 CLI 命令组完成后的 HEP 2.0.17 提交、候选验收与正式发布准备。
2. **Affected modules**：Redis 环境/缺失程序诊断、CLI 路由/help、agent catalog/search、对应测试及维护文档；HEP 版本声明。
3. **Public interfaces**：新增正式 `Jarvis redis info/check/install/start` 子命令、独立 help、JSON 和 agent 发现；裸命令继续默认 info。安装/启动指引只打印，YAML 结构和扫描算法不变。
4. **Main risks**：误纳入工作区 Sampler/runtime/Diver 改动；以本地 SDK 安装冒充公开索引可安装；所需 `Jarvis-HEP-Sampler>=0.1.0,<0.2` 尚未公开，不能先发布依赖不可安装的 HEP。
5. **Verification performed**：干净 detached 源码 `1c8c82fc` 全量 fast gate：1489 passed、3 skipped（缺少 sibling Jarvis-Examples）、166 deselected、478 subtests；完整隔离安装与 pip check、真实 Redis/spawn/project reload、24 项 installed agent、39 项 installed Redis 和 22 项发行回归全部通过。完整命令和历史结果见 [发行验收证据](release-2.0.17-verification.json)。

## 授权与发行范围

2026-10-15 用户明确授权提交 Redis 改进并发布 HEP 新版本，同时要求忽略当前 Sampler 改动。本版本为 2.0.17；既有 2.0.15/2.0.16 候选和历史记录保留。

同一会话中用户随后选择“暂只完成 HEP 候选，保留 SDK 发布待办”。本次以提交、候选验收与 draft 交付收尾，正式 HEP/SDK 公开发行作为后续工作。

本次提交/候选只纳入 Redis 环境引导、正式子命令及帮助/agent 发现、对应测试与文档、HEP 两处版本声明。已有 committed baseline 中的 sampler SDK 集成仍是依赖关系；本轮未提交的 sampler/runtime 修补、Diver 源码/测试、Diver package-data 与 sampler 审阅资料均排除并保留。`pyproject.toml` 只纳入本次版本变更，不带入其它未验收工作区差异。

## 用户可见变化

- 首次 CLI 使用后缺少 Redis-compatible executable 时，stderr 给出系统适用的安装/启动指引；它是 PATH 检查，不是 pip 安装钩子。
- `Jarvis redis` / `info` 展示本机系统、executable/version 与说明；`check` 每次执行 PATH 和有界版本检查，缺失、不可用/过旧或 native Windows 退出 1。
- `install` / `start` 只打印指导，不执行系统安装、服务操作、网络/端口探测或扫描。普通扫描仍由 Jarvis 自动管理本次 Redis；Windows 引导在 WSL 内安装/运行 HEP 与 Redis。
- Redis 是正式命令组：父级 help 列出四个 Commands，子级有独立说明和 `--json`；JSON 可在子命令前后使用，裸命令保留 info 默认，agent 公开 exact child commands。
- 安装建议按发行版/版本处理 Redis 或 Valkey；旧 Ubuntu/Debian 使用已核验的官方 signed APT recipe，RHEL 8 说明 `redis:6` 与 updated AppStream 要求。

## 本次验收与源码状态

本阶段应先完成 Redis-only 提交，再从该干净 committed source 构建 wheel/sdist；核对内容没有混入工作区未授权改动，并在独立环境安装全部声明依赖。执行 `pip check`、installed Redis help/JSON/agent gate、真实 Redis/spawn smoke 和发行回归，保存源码 commit、命令、结果、产物 checksum 与文件范围。

验收状态：候选验收完成。Redis 与版本提交 `569d3865f60486da048b69c74c4ef3730713265d`、发行测试补正 `1c8c82fc41efbc732f533cf78f74034da611033d` 已推送至 `Jarvis/feat/sampler-plugins-s1`。最终发行源码是干净 detached `1c8c82fc`；本次仅在后续提交中补录不打包的维护文档与验收证据。

产物位于 `dist/2.0.17/`：wheel SHA256 `2c76bb9658604b5bb764223176c9594a981e739f3300cca0b26681cd1d1d810a`，sdist SHA256 `40bba3ff02bf9ebb47530f237091a802e3e7a4fa6c29a1b2d95d7cbeff10ad89`。wheel 的全部 289 个文件（含 dist-info）与先前完整安装验收版本逐字节一致；283 个运行文件和 sdist 412 项 tracked 文件与最终源码相符。

[GitHub draft Release](https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP/releases/tag/untagged-c8534defe6d5e2c72938) 保存 HEP 产物和验收证明，保持未公开状态。完整验收：1489 passed、3 skipped（依赖不在干净检出中的 sibling Jarvis-Examples）、166 deselected、478 subtests；隔离完整安装、pip check、真实 Redis/spawn/project reload、installed agent 24 项、installed Redis 39 项、发行回归 22 项均通过。没有使用 `--no-deps` 或 system-site-packages；唯一非公共依赖来源是保持原样的已验收 SDK 0.1.0 wheel。

开发阶段历史证据为 [Redis 环境引导验收](redis-environment-verification.json) 和 [标准命令组验收](redis-cli-help-verification.json)。后者 202 项定向、25 项安装检查通过；工作区 fast gate 为 1509 passed / 478 subtests，保留 1 项既有外部 iDM 描述编码失败。这些结果不证明新的版本源码或发行包已经通过。

### 发行门禁的检出环境补正

首次干净源码全量运行暴露了既有 `test_splash_lists_ps_rows_and_quit` 的环境依赖：测试无条件要求分支标志，但 detached HEAD 没有分支名，生产 topbar 按设计只显示路径。GitHub Release 标签检出同样是 detached HEAD；仅给本地验收创建临时分支不能解决发行 CI。

本阶段将该 UI 测试自己的 `git_branch` 输入固定为测试分支，并断言显示该分支。此补正只让发行测试自包含，不改变 monitor 生产代码，也不纳入 Sampler 修补。最初失败与修正后的验收分别保留，不放宽或跳过断言。

## 公开依赖与发布状态

2026-10-15 公共 PyPI 查询快照：HEP 最新公开版本为 2.0.14，Operas 为 1.5.0，`Jarvis-HEP-Sampler` 项目查询返回 HTTP 404。查询入口分别为 [HEP JSON](https://pypi.org/pypi/Jarvis-HEP/json)、[Operas JSON](https://pypi.org/pypi/Jarvis-Operas/json)、[Sampler JSON](https://pypi.org/pypi/Jarvis-HEP-Sampler/json)；这些是当日观测，发布前须重新核验。

既有 baseline 声明 `Jarvis-HEP-Sampler>=0.1.0,<0.2`，本次不改变该依赖。允许使用已验收且保持原样的本地 SDK 0.1.0 wheel 验证候选，其余依赖按正常完整安装解析；这属于候选验收，不是仅依赖公开索引即可安装的证明。

正式 HEP 发布仍等待公开依赖可用。HEP 候选、draft Release 及完整验收证据已完成；按用户后续选择保留 SDK 发布待办，本次不修改、提交或发布 SDK。依赖问题未解决前，不将 draft/本地候选标记为公开发布完成，不上传不可正常安装的 HEP。后续公开发行沿用 GitHub Release/Trusted Publishing。
