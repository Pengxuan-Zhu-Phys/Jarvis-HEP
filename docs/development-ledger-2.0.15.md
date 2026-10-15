# Jarvis-HEP 开发台账（2.0.15—2.0.17）

本发行从公开 master `0bbcf675` 移植独立的可靠性、CLI agent/YAML 与 Redis 功能。原集成分支 `30f2467` 保留早期阶段完整记录；旧SDK候选`1c8c82fc`已被替代，不复用其测试数或产物。以下为当前生效阶段。


## 2026-10-15 2.0.17 SDK-free 公开发行（最新授权）

1. **Roadmap phase**：按用户最新明确要求公开发布 HEP，所有 Sampler 工作延后；取代此前仅完成候选的选择。
2. **Affected modules**：从公开 master `0bbcf675` 移植不依赖 SDK 的可靠性、CLI agent/YAML、Redis 及测试文档、发行版本/工作流。
3. **Public interfaces**：发布 agent/YAML 与 Redis 增量能力，保留内置采样；不声明外部 SDK/plugin 依赖或 CLI。
4. **Main risks**：误带入已提交 SDK integration `0d859502` / `52876060` 或 dirty Sampler/runtime/Diver；复用旧 SDK 候选验收；未核验公共 PyPI 全依赖安装即发布。
5. **Verification performed**：待本阶段测试结果补录；必须对当前隔离源码/wheel/sdist重新执行门禁与公开发布后安装验证。

旧 `1c8c82fc` 候选曾在本地 SDK 安装安排下完成测试与 draft，但因其 SDK 依赖不适合本次公开发行，现已被 SDK-free master 派生源替代。旧证据/产物作为历史保留；不能沿用其 checksum、1489 计数或当时的 SDK 阻断结论判断新版本状态。

当前隔离源保留 HEP 内置采样和公共依赖，移植可靠性/CLI agent/Redis 独立功能，不合并 Sampler 分支、不改动或代发布 SDK。完整验收应包括 clean source fast gate、公共 PyPI 全依赖安装/pip check、installed agent/Redis/help/JSON、真实 Redis/spawn/project reload 与相关发行回归；通过正常 OIDC 工作流公开后，还需核验 PyPI 2.0.17 和全新安装。当前事实见 [2.0.17 发行记录](release-2.0.17.md)。
