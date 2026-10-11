# Jarvis-HEP 架构审阅

2026-10-11 根据当前实现重建，未恢复历史审阅原文。

## 当前职责

`cli/client.py` 路由命令并转发依赖 CLI；`cardload` 加载和验证 task card；`runtime/core.py` 通过 supervisor 与 scan driver 管理运行。

Worker 在子进程中建立表达式上下文、mapper、Operas、calculator 和 likelihood，再拉取 Redis 任务。样本级异常处理位于 `process_task` 内；初始化失败在此边界之前，不能依赖样本失败统计发现。

Factory 管理 Worker 启停、watchdog 和进程内 monitor snapshot。新建 Factory 不拥有其它控制进程的内存快照。外部 monitor 应从 OS inventory 获得 Worker owner IDs，再读取对应 Redis boards。

Archiver 负责结果持久化；scan driver 等待持久化与运行计数。初始化失败必须另有控制进程可观察的信号，不能一直走 archive 等待或普通 respawn。

## 2.0.15 修复约束

- 表达式验证使用与运行时一致的函数与常量命名空间，只编译，不做数值调用。
- Worker 启动错误与样本级失败分别表达，不虚构未执行样本的终态计数。
- Redis 快照读取只读，只读取明确的 Worker IDs；不扫描 Worker namespace 或队列 payload。
- snapshot 的采集时间与 Worker heartbeat 时间是不同字段。
- 依赖 help 修复必须进入可安装的发行包，不能仅修改本机 site-packages。
- 发行构建只使用明确纳入的源码，保留开始任务前已有的无关改动。

## 阶段状态

确认了验证漏检、初始化失败后持续等待、带 ref 的 JSON 丢失、monitor 的默认零值、plot 诊断不足、intent 参数位置差异与 Operas help 对象字符串化。修复与验收状态见 [台账](development-ledger-2.0.15.md)。

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
