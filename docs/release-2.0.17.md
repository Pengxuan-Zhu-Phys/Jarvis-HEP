# Jarvis-HEP 2.0.17

1. **Roadmap phase**：按最新授权公开发布不依赖外部 Sampler SDK 的 HEP 2.0.17；从公开 master 隔离移植可靠性、CLI agent/YAML 与 Redis 改进。
2. **Affected modules**：不依赖 SDK 的 Worker/表达式/monitor/CLI 可靠性、YAML 诊断、agent 辅助、Redis 引导、对应测试文档、版本与发行流程。
3. **Public interfaces**：CLI agent 合同/查询/schema/模板/inspect/草稿与 Redis 正式子命令/help/JSON；保留内置采样、YAML 结构和正常扫描恢复，不包含外部 SDK/plugin CLI。
4. **Main risks**：误带入已提交或未提交的 Sampler 集成；沿用旧 SDK 候选证据；发行包依赖或 metadata 残留 SDK；公开安装未经核验。
5. **Verification performed**：干净源码 `e8af5221`：1426 passed、3 skipped（缺少 sibling Jarvis-Examples）、164 deselected、477 subtests；335 项代码/测试文件与提交逐字节一致。全新隔离环境从公共 PyPI 解析全部依赖，pip check 与 Twine 通过，SDK 模块/发行包均不存在；真实 Redis/spawn/双 Worker/monitor/项目打包重跑、installed agent 24 项、installed Redis 39 项、发行回归 22 项均通过。

## 当前授权与源码范围

用户最新明确要求公开 HEP，并将所有 Sampler 工作后移，取代此前“暂只完成候选”的选择。本版从公开 `master` 基线 `0bbcf675` 的独立 worktree 准备，只移植不依赖外部 SDK 的 2.0.15 可靠性、2.0.16 CLI agent/YAML 和 2.0.17 Redis 改进。已提交 SDK integration `0d859502` / `52876060` 与原工作区未提交的 Sampler/runtime/Diver/审阅资料都延后，不改动原工作区，不合并或发布 SDK。

采样继续使用 HEP 内置实现，发行包不声明 `Jarvis-HEP-Sampler` 依赖，不提供其插件注册/加载 CLI。旧 `1c8c82fc` SDK-dependent 候选和 draft 被这份公开发行源取代；旧产物/checksum/测试保留为历史，不能复用为当前结果。

## 用户可见变化

- Worker 初始化故障有界反馈、表达式预检诊断、monitor JSON/缺失值与 CPU snapshot 等不依赖 SDK 的可靠性修复；正常运行期恢复保留。
- Agents 可查询 CLI 合同、命令效果和参数，导出本地 schema/准确模板，解释卡片生效值及来源，并沿 validate JSON 导航修正 YAML；显式重复键拒绝，合法 merge override 保留。
- 功能反馈生成可填写邮件草稿，不自动发送；查询不隐含授权执行 check/run 或用户代码。
- `Jarvis redis info/check/install/start` 为正式命令组，含独立 help、JSON 与 exact agent 发现；缺失程序时提供按发行版/版本选择的安装指导。安装/启动命令只打印，正常扫描继续自动管理本次 Redis。

## 当前验收与发布状态

干净源码 `e8af5221`：1426 passed、3 skipped（缺少 sibling Jarvis-Examples）、164 deselected、477 subtests；335 项代码/测试文件与提交逐字节一致。全新隔离环境从公共 PyPI 解析全部依赖，pip check 与 Twine 通过，SDK 模块/发行包均不存在；真实 Redis/spawn/双 Worker/monitor/项目打包重跑、installed agent 24 项、installed Redis 39 项、发行回归 22 项均通过。

源码：分离提交 `3b9d41d0dc2171bdaf79990a4e23faddd9ff6ce9`，验证接口兼容修正/最终打包源码 `e8af5221eb84b4c1feccf492b2676df764505b3f`。独立分支为 `codex/release-2.0.17`，Sampler 原分支和未提交修改保留。

发行产物：`dist/2.0.17-sdk-free/`。wheel SHA256 `a823b9c71a7959fdb2ca5a45b4d03e9190058d0a47ed14f97fe4f76393c53c1b`；sdist SHA256 `91711348da271fe6356791c0c945033b1c8bae2f03de466599242418ff10a3a2`。完整命令、源码范围、初次移植发现的 notes 接口问题与最终验收见 [结构化证据](release-2.0.17-verification.json)。

当前本地源码和安装包验收完成；正式 GitHub/PyPI 发布执行中，尚未在此预先标记成功。后续仅补录不参与打包的维护文档，不改变已验收的产品源码和产物。

通过 GitHub Release/Trusted Publishing 公开后，继续核验工作流成功、PyPI 2.0.17 wheel/sdist 与发行元数据，并在全新环境从公共索引安装全部依赖和运行 smoke。当前没有 SDK 发布前置条件，也不能用本地 SDK、`--no-deps` 或 system-site-packages 替代此验收。

## 历史候选说明

此前 `1c8c82fc` 候选使用本地SDK完成测试并创建draft，但其外部SDK依赖不适合本次公开发行，已被公开 master `0bbcf675` 派生的无SDK依赖版本取代。旧源码、产物、checksum和验收记录保留在原集成分支，不构成本次发行证据；本次不再等待或发布SDK。
