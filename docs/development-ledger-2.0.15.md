# Jarvis-HEP 开发台账（2.0.15—2.0.17）

本发行从公开 master `0bbcf675` 移植独立的可靠性、CLI agent/YAML 与 Redis 功能。原集成分支 `30f2467` 保留早期阶段完整记录；旧SDK候选`1c8c82fc`已被替代，不复用其测试数或产物。以下为当前生效阶段。


## 2026-10-15 2.0.17 SDK-free 公开发行（最新授权）

1. **Roadmap phase**：按用户最新明确要求公开发布 HEP，所有 Sampler 工作延后；取代此前仅完成候选的选择。
2. **Affected modules**：从公开 master `0bbcf675` 移植不依赖 SDK 的可靠性、CLI agent/YAML、Redis 及测试文档、发行版本/工作流。
3. **Public interfaces**：发布 agent/YAML 与 Redis 增量能力，保留内置采样；不声明外部 SDK/plugin 依赖或 CLI。
4. **Main risks**：误带入已提交 SDK integration `0d859502` / `52876060` 或 dirty Sampler/runtime/Diver；复用旧 SDK 候选验收；未核验公共 PyPI 全依赖安装即发布。
5. **Verification performed**：见下方已完成的源码/安装验收与发布状态；必须对当前隔离源码/wheel/sdist重新执行门禁与公开发布后安装验证。

旧 `1c8c82fc` 候选曾在本地 SDK 安装安排下完成测试与 draft，但因其 SDK 依赖不适合本次公开发行，现已被 SDK-free master 派生源替代。旧证据/产物作为历史保留；不能沿用其 checksum、1489 计数或当时的 SDK 阻断结论判断新版本状态。

当前隔离源保留 HEP 内置采样和公共依赖，移植可靠性/CLI agent/Redis 独立功能，不合并 Sampler 分支、不改动或代发布 SDK。完整验收应包括 clean source fast gate、公共 PyPI 全依赖安装/pip check、installed agent/Redis/help/JSON、真实 Redis/spawn/project reload 与相关发行回归；公开 wheel 已核验 PyPI 2.0.17 和全新安装；本次 OIDC 配置失败后的上传例外与 sdist 待办见下方。当前事实见 [2.0.17 发行记录](release-2.0.17.md)。

## 本次源码、CI 与公开安装状态

源码功能提交 `3b9d41d0`、验证接口兼容修正 `e8af5221` 已进入 master；公开 tag `v2.0.17` 指向 `8e0add06`。macOS 与 Linux CI 均为 1426 passed、3 skipped（缺少 sibling Jarvis-Examples）、164 deselected、477 subtests；分别用时 205.99s / 265.40s。Linux CI 的完整安装/pip check、真实 Redis/spawn/双 Worker/monitor/项目重跑、agent 24、Redis 39、发行回归 22 均通过。

[GitHub Release](https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP/releases/tag/v2.0.17) 已公开。OIDC publisher 因 PyPI 未匹配 Trusted Publisher 配置而报 `invalid-publisher`；本次使用既有本地 PyPI 凭证上传成功 CI 的原样产物，没有重建。PyPI wheel 已公开，全新环境仅从公共索引安装完整依赖，pip check、无 SDK、276 个运行文件身份比对、Redis 39 / agent 24 均通过。PyPI 源码包上传未获服务器响应；最后一次有界请求在 36.913 秒后 ReadTimeout，全部上传进程已停止。PyPI wheel 与公开安装已完成，已验收源码包可从 GitHub Release 下载。

发布 wheel SHA256 为 `25257b646c1a7ce828dae136e1e111de7c70d749d7b5b6857b7c687f227fb3ae`；CI sdist 的预期 SHA256 为 `1843e006921161dd6a46507b7b7c12b52b48990b9e3af1c76541788dafbe9915`，未在 PyPI 公开；GitHub Release 另提供通过本地完整验收的源码包。发布前本地产物的 a823… / 9171… hash 仅作本地验收记录，不能当成公开 CI 产物 hash。完整状态见 [发行记录](release-2.0.17.md) 与 [验收证据](release-2.0.17-verification.json)。
