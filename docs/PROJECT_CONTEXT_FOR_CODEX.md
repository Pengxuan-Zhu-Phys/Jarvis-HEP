# Jarvis-HEP 项目上下文

维护基线：2026-10-11。本文件在原约定路径及 Git 历史中均无旧版本，现依据当前代码、README、安装与发布文档重建；不代表恢复了遗失的历史原文。

## 产品与执行模型

发行包为 `Jarvis-HEP`，Python 包为 `jarvishep2`，用户命令为 `Jarvis`。当前源码版本为 2.0.15，发行候选验证中，尚未公开发布。Python 最低版本为 3.10。

每次扫描由控制进程管理独立 Redis、spawn Worker 和 Archiver。Worker 执行 Operas 或外部 calculator；Archiver 持久化样本结果，控制进程负责采样、等待、恢复和清理。监控按已知扫描与进程读取运行数据。

Jarvis-Operas 提供函数与常量注册，Jarvis-HEP-Portal 提供 calculator I/O，JarvisPLOT 提供绘图。Jarvis-HEP-Sampler 提供当前代码使用的 sampler SDK。这些均为完整发行包依赖。

## 本轮授权与范围

执行 [2.0.15 开发台账](development-ledger-2.0.15.md)的七项修复与验收，准备完整 wheel、sdist 和发布说明。**最终发布由用户执行**：本轮不上传 PyPI、不发布 GitHub Release、不推送发行 tag。

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

## 2026-10-11 第一阶段实施

已实现表达式编译诊断 JV2-EXP-001（不执行数值函数）；Worker 使用 spawn pipe 报告 ready/error，Factory 在 120 秒内等待启动，初始化异常直接终止且不重启；adaptive feedback 与结果等待检查后续启动错误。monitor 按 OS Worker 清单读取已知 boards，JSON ref 输出单一对象，缺失数量为 null、采集时间为当前时间。CLI 补上 intent 冲突与 task-card plot 提示。Operas 外部源码改用 Styled 保留 Rich renderable。

相关 fast 回归 148 passed、30 subtests passed；完整 fast gate 与真实 spawn/恢复/分布式回归进行中。当前 HEP 版本维持 2.0.14，全部验收后才准备 2.0.15。最终公开发布由用户执行。

## 2026-10-12 验收与发行范围核对

更正第一阶段的启动顺序：Worker 与 Archiver 并行启动，supervisor 在两者 spawn 后等待 Worker ready/error（默认有界等待 120 秒）；Factory.start_workers 仍异步返回。启动错误通过 spawn pipe 返回原始类型和消息；后续 adaptive feedback 与结果等待也检查失败。创建部分 Worker 后 spawn 失败时停止已启动进程。正常样本执行失败仍进入原有样本计数，正常运行期崩溃仍由 watchdog 恢复。

新增发行阻断项 QA-02（macOS/Windows CPU 活进程估算，PID+create_time 增量，不重复过期读数）、SDK-01（nested 正常返回只表示 engine returned，预算/custom/convergence 未被引擎明确报告时不猜测收敛）及 QA-03（吞吐 gate 预热管线、16 个测量点，保留 1.25 倍门槛）。相关控制面 140 passed、30 subtests；运行/恢复/分布式主体 97 passed，吞吐与最新启动失败独立验收 2 passed。SDK 完整包安装测试 146 passed，Operas CLI 43 passed。第一版独立发行源码 fast gate 1265 passed、3 skipped、141 subtests；最终代码快照 gate 进行中。

发行包含当前已跟踪 sampler plug-in 功能及本轮修复；未跟踪 Diver/F1/F2/F5 与其 package-data 改动保留在工作区且排除发行。独立源码不依赖 sibling Examples 的未发布非 ASCII 卡片；该卡片在工作区语料测试仍有诊断，不修改外部案例。依赖候选为 SDK 0.1.0 和 Operas 1.4.2，用户先发布依赖，再发布 HEP 2.0.15。

## 2026-10-12 发行候选完成

所有台账修复与 QA 项完成。源码版本 2.0.15，完整 wheel/sdist 在干净 venv 安装所有依赖后通过 Twine、pip check、真实双 Worker/ref JSON/CLI/项目重载 smoke，以及 22 项安装包回归。最终源码 fast gate 1272 passed、3 skipped、141 subtests。依赖候选 SDK 0.1.0（146 tests）、Operas 1.4.2（43 CLI tests、安装包帮助验证）。

交付位置 `dist/2.0.15/`，持久证据与发布顺序见开发台账、bugfix-2.0.15/verification.json 和 release-notes.md。仅候选准备完成，用户先发布依赖，再自行发布 HEP；代理未上传、推送或创建 tag。保留未验收 Diver 与其它用户改动，实际 wheel/sdist 检查确认未包含 Diver。
