# Jarvis-HEP 项目上下文

维护基线：2026-10-11。本文件在原约定路径及 Git 历史中均无旧版本，现依据当前代码、README、安装与发布文档重建；不代表恢复了遗失的历史原文。

## 产品与执行模型

发行包为 `Jarvis-HEP`，Python 包为 `jarvishep2`，用户命令为 `Jarvis`。当前源码版本为 2.0.16，发行候选准备中，尚未公开发布。Python 最低版本为 3.10。

每次扫描由控制进程管理独立 Redis、spawn Worker 和 Archiver。Worker 执行 Operas 或外部 calculator；Archiver 持久化样本结果，控制进程负责采样、等待、恢复和清理。监控按已知扫描与进程读取运行数据。

Jarvis-Operas 提供函数与常量注册，Jarvis-HEP-Portal 提供 calculator I/O，JarvisPLOT 提供绘图。Jarvis-HEP-Sampler 提供当前代码使用的 sampler SDK。这些均为完整发行包依赖。

## 本轮授权与范围

2.0.15 阶段只授权准备候选，由用户最终发布。2026-10-13 用户新增授权：提交两轮 CLI agent 改进并发布新版本。当前 2.0.16 可执行 commit、push 和正式发行；必须先通过发行验收、确认依赖已公开。历史候选准备记录保留。

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

## 2026-10-13 CLI agent 集成阶段

用户新授权：先审计形成 report，再补齐 Codex/Claude 等通过 CLI 编写 YAML、查询能力和问答所需的显式行为约定，并仿照 Jarvis-PLOT 增加邮件需求反馈。实施前报告见 [CLI agent 集成报告](cli-agent-integration-report.md)。

新增 `Jarvis agent` 自举合同（与 `man agent-readme` 同源）、离线 schema、模板/前置条件导出、邮件草稿；validate JSON 增加机器诊断导航和验证范围，JSON 参数/I/O/内部失败如实输出。接口为增量扩展，不改 YAML 结构或扫描算法。查询仍可能初始化依赖缓存和首次提示 marker；check/run 是真实执行，现有用户授权持续有效。PLOT 邮件功能实际为反馈模板，无自动发送，本次沿用该边界。

纠正旧示例元数据：四个仅含 Sampling 的示例是 fragment，需合入有 Scan/EnvReqs 及 Calculators/Operas 的完整卡。loader 仍对重复 YAML 键使用 last-value 语义，合同要求唯一键并说明验证限制。保留任务开始前的 Diver、pyproject 与 sampler 审阅改动；版本和发行候选不变。本阶段验收以 report 和 `cli-agent-verification.json` 为准。

本阶段实施及文档同步完成：最终定向 194 tests/316 subtests、安装测试 wheel 13 项 CLI 检查通过；工作区 fast gate 1345 tests/437 subtests 通过，唯一失败是已知 sibling iDM 卡片非 ASCII 描述，未修改该外部案例。无真实邮件发送或公开发布。

## 2026-10-13 CLI agent 第二轮实际体验改进

按用户新请求通过 CLI 实测后，先将问题写入原报告，再修正。新增 brief/section/command 定向合同、关键词 search、用户卡片 inspect 生效配置/来源/路径、schema 小索引与本地 JSONPointer；修正 fragment 验证目标、import-path 算子查询、未知手册路径假成功及 token/阈值/heartbeat 说明。邮件草稿可填写证据段并报告缺项，仍不发信。

**更正第一轮边界**：task/default YAML 现在拒绝重复显式键（JV2-YAML-003，路径/两处位置），不再 last-wins；合法 YAML merge overrides/aliases 仍允许。外部旧 MultiNest 卡的重复 seed 保留但被准确拒绝，语料测试只为该已知文件断言新诊断，其它断言不放宽。inspect 仅反映已加载/规范化设置；实际 workers/port 未选择，不导入用户 operator。

最终定向 332 tests/373 subtests、22 项安装包 CLI 检查通过；fast gate 1416/478 通过，仍 1 个已知外部 iDM 编码失败。原默认完整合同兼容，新短合同约 10 KB、runtime 查询约 3 KB。证据见 `cli-agent-review2-verification.json`；版本、原候选及 Diver 改动不变，未发布。


## 2026-10-13 2.0.16 提交与发行

用户已授权 commit 与新版本发布，本阶段取代此前仅准备候选的权限边界。版本 2.0.16 包含两轮 CLI agent 改进及已提交的 2.0.15 可靠性/sampler integration 基线；保持原有 YAML 结构，显式重复键改为 JV2-YAML-003 错误。Diver 草稿、其 package-data 及未验收的 sampler 审阅文件保留且排除。

从干净 commit 副本构建完整 wheel/sdist，执行源码 fast gate、全依赖安装、pip check、真实 Redis/spawn smoke、旧可靠性回归和新增 installed agent CLI smoke。发布正常使用 GitHub Release/Trusted Publishing，不使用本地 token 上传。当前依赖检查：Operas 1.4.2 已公开；Jarvis-HEP-Sampler 0.1.0 在 PyPI 缺失，须先解决才能正式发布 HEP。具体状态与证据见 [2.0.16 发行记录](release-2.0.16.md)。
