# Jarvis-HEP 项目上下文

维护基线：2026-10-11。本文件在原约定路径及 Git 历史中均无旧版本，现依据当前代码、README、安装与发布文档重建；不代表恢复了遗失的历史原文。

## 产品与执行模型

发行包为 `Jarvis-HEP`，Python 包为 `jarvishep2`，用户命令为 `Jarvis`。当前源码版本为 2.0.17，正在从公开 master 基线准备不依赖外部 Sampler SDK 的正式发行。公开发布结果待本阶段测试结果补录。Python 最低版本为 3.10。

每次扫描由控制进程管理独立 Redis、spawn Worker 和 Archiver。Worker 执行 Operas 或外部 calculator；Archiver 持久化样本结果，控制进程负责采样、等待、恢复和清理。监控按已知扫描与进程读取运行数据。

Jarvis-Operas 提供函数与常量注册，Jarvis-HEP-Portal 提供 calculator I/O，JarvisPLOT 提供绘图，均为完整发行包依赖。采样使用 HEP 内置实现；本次公开发行不依赖 Jarvis-HEP-Sampler，不提供外部 SDK 插件注册/加载 CLI。

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

旧 `1c8c82fc` 候选及其 draft/本地 SDK 安装证据被此发行源替代；它们保留为历史，不能复用其 checksum 或测试计数。当前必须对隔离源码和完整 wheel/sdist 重新运行 fast gate、公共 PyPI 全依赖安装/pip check、installed agent/Redis、真实 Redis/spawn，再经正常 Trusted Publishing 发布并核验公共索引安装。

当前结果：待本阶段测试结果补录。提交、产物和公开发布状态见 [2.0.17 发行记录](release-2.0.17.md)。原 Sampler 分支与用户工作区保持原样。
