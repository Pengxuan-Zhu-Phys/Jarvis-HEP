# Jarvis-HEP 项目上下文

维护基线：2026-10-11。本文件在原约定路径及 Git 历史中均无旧版本，现依据当前代码、README、安装与发布文档重建；不代表恢复了遗失的历史原文。

## 产品与执行模型

发行包为 `Jarvis-HEP`，Python 包为 `jarvishep2`，用户命令为 `Jarvis`。源码版本为 2.0.17。2.0.17 已公开不依赖外部 Sampler SDK 的 GitHub Release 与 PyPI wheel；完整公共索引安装已核验，PyPI 源码包上传超时，源码包由 GitHub Release 提供。当前 master 在该发行之后要求 `Jarvis-HEP-Sampler>1`。Python 最低版本为 3.10。

每次扫描由控制进程管理独立 Redis、spawn Worker 和 Archiver。Worker 执行 Operas 或外部 calculator；Archiver 持久化样本结果，控制进程负责采样、等待、恢复和清理。监控按已知扫描与进程读取运行数据。

Jarvis-Operas 提供函数与常量注册，Jarvis-HEP-Portal 提供 calculator I/O，JarvisPLOT 提供绘图，Jarvis-HEP-Sampler 提供插件采样器，均为完整发行包依赖。`Jarvis sampler` 与 `jsampler` 是同一条命令。2.0.17 那次公开发行不依赖 Jarvis-HEP-Sampler。

## 本轮授权与范围

用户最新明确要求公开发布 HEP 2.0.17，并将所有 Sampler 工作后移，取代此前仅交付候选的选择。发行从 master 0bbcf675 隔离构建，仅移植不依赖 SDK 的运行可靠性、CLI agent/YAML 和 Redis 改进。已提交的 SDK 集成 0d859502 / 52876060 与未提交的 Sampler/runtime/Diver 工作都不纳入；不修改或代为发布 SDK。必须完成本发行源码、公共索引全依赖安装与公开发布核验，历史候选结果不代替当前验收。

工作区开始时已有 Diver 源码、测试、package-data 改动及 sampler 审阅资料。保留这些改动；发行候选明确记录其是否纳入，并以干净源码与完整发行包验证，不能误把未验收内容包含进本轮发行。

## 维护规则

代码修改前同时阅读本文件、PRD、架构审阅、代码地图和实施路线。阶段完成后同步这五份维护文件与台账。

保留现有命令、YAML 结构和正常 Worker 恢复行为；诊断应给出路径、原因和纠正方法。修复后运行相关定向测试、真实运行对照和完整 fast gate；发行包还需在干净环境安装全部依赖并核验。

## 关联维护文件

- [产品要求](PRD_JARVIS_HEP_MASTER.md)
- [架构审阅](ARCHITECTURE_REVIEW_JARVIS_HEP.md)
- [代码地图](CODE_MAP_JARVIS_HEP.md)
- [实施路线](IMPLEMENTATION_ROADMAP.md)
- [发布流程](../RELEASING.md)



## 2026-10-15 2.0.17 SDK-free 公开发行（最新授权）

用户明确要求这次公开发布 HEP，把所有 Sampler 工作后移，取代仅完成候选的选择。发行从公开 master `0bbcf675` 隔离构建，移植不依赖外部 SDK 的 2.0.15 可靠性、2.0.16 CLI agent/YAML 与 2.0.17 Redis 改进；保留内置采样。已提交 SDK 集成 `0d859502` / `52876060` 和原工作区所有未提交 Sampler/runtime/Diver 编辑都延后，HEP 不声明 Sampler SDK 依赖，不需要 SDK 先公开。

旧 `1c8c82fc` 候选及其 draft/本地 SDK 安装证据被此发行源替代；它们保留为历史，不能复用其 checksum 或测试计数。当前必须对隔离源码和完整 wheel/sdist 重新运行 fast gate、公共 PyPI 全依赖安装/pip check、installed agent/Redis、真实 Redis/spawn，发布并核验公共索引安装；本次 OIDC 配置失败后的原样 CI 产物上传例外见下方记录。

当前结果：见下方已完成的源码/安装验收与发布状态。提交、产物和公开发布状态见 [2.0.17 发行记录](release-2.0.17.md)。原 Sampler 分支与用户工作区保持原样。

## 本次源码、CI 与公开安装状态

源码功能提交 `3b9d41d0`、验证接口兼容修正 `e8af5221` 已进入 master；公开 tag `v2.0.17` 指向 `8e0add06`。macOS 与 Linux CI 均为 1426 passed、3 skipped（缺少 sibling Jarvis-Examples）、164 deselected、477 subtests；分别用时 205.99s / 265.40s。Linux CI 的完整安装/pip check、真实 Redis/spawn/双 Worker/monitor/项目重跑、agent 24、Redis 39、发行回归 22 均通过。

[GitHub Release](https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP/releases/tag/v2.0.17) 已公开。OIDC publisher 因 PyPI 未匹配 Trusted Publisher 配置而报 `invalid-publisher`；本次使用既有本地 PyPI 凭证上传成功 CI 的原样产物，没有重建。PyPI wheel 已公开，全新环境仅从公共索引安装完整依赖，pip check、无 SDK、276 个运行文件身份比对、Redis 39 / agent 24 均通过。PyPI 源码包上传未获服务器响应；最后一次有界请求在 36.913 秒后 ReadTimeout，全部上传进程已停止。PyPI wheel 与公开安装已完成，已验收源码包可从 GitHub Release 下载。

发布 wheel SHA256 为 `25257b646c1a7ce828dae136e1e111de7c70d749d7b5b6857b7c687f227fb3ae`；CI sdist 的预期 SHA256 为 `1843e006921161dd6a46507b7b7c12b52b48990b9e3af1c76541788dafbe9915`，未在 PyPI 公开；GitHub Release 另提供通过本地完整验收的源码包。发布前本地产物的 a823… / 9171… hash 仅作本地验收记录，不能当成公开 CI 产物 hash。完整状态见 [发行记录](release-2.0.17.md) 与 [验收证据](release-2.0.17-verification.json)。
