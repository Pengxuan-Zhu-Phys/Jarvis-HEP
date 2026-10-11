# Jarvis-HEP 代码地图

2026-10-11 按当前 checkout 重建。顶层若干模块为兼容导入 shim，编辑时定位实际实现。

| 职责 | 实现 | 主要验证 |
| --- | --- | --- |
| CLI 路由与参数 | `jarvishep2/cli/client.py`、`cli/help.py`、`cli/run.py` | `tests/test_cli.py`、`test_diagnostic_guidance.py` |
| 卡片加载与验证 | `jarvishep2/cardload/task_config.py`、`cardload/task_validation.py`、`cardload/contracts/` | `tests/test_task_validation.py`、`test_task_schema.py` |
| 表达式与 Operas 注册 | `jarvishep2/expression.py`、`operas_functions.py`、`operas.py`、`likelihood.py` | `tests/test_expression.py`、`test_operas_functions.py`、`test_operas_bridge.py` |
| Worker 初始化与任务 | `jarvishep2/runtime/worker.py`、`runtime/sample_executor.py` | `tests/test_worker_mvp.py`、`test_worker_failure.py` |
| 控制进程与恢复 | `jarvishep2/runtime/core.py`、`runtime/factory.py`、`runtime/_runtime_supervisor.py`、`runtime/_scan_driver.py` | `tests/test_core_collaborators.py`、`test_graceful_stop.py`、`test_distributed_acceptance.py` |
| Redis 协议 | `jarvishep2/queue/redis_queue.py`、`queue/_redis_proc_board.py`、`queue/_redis_control.py` | `tests/test_redis_proc_board.py`、`test_redis_queue.py` |
| 快照与进程清单 | `jarvishep2/dashboard.py`、`process_cleanup.py`、`monitor/collector.py` | `tests/test_dashboard_reader.py`、`test_monitor_snapshot.py`、`test_process_cleanup.py` |
| plot 桥接 | `jarvishep2/cli/client.py`、`plot_bridge.py`、`plot_scene.py` | `tests/test_cli.py`、`test_plot_bridge.py` |
| 发行准备 | `pyproject.toml`、`jarvishep2/__init__.py`、`scripts/full-release.sh`、`scripts/smoke_install.py`、`jhrel` | `tests/test_distribution_version.py`、`test_package_layout.py` |

Jarvis-Operas 的 help renderer 位于其仓库 `jarvis_operas/cli.py`，对应 `tests/test_cli.py`。不得把安装路径中的临时修补当成源码修复。

相关文档：验证见 `docs/validation-diagnostics.md`；监控见 `docs/monitor-tui.md` 与 `docs/monitor-inflight.md`；本轮任务与发行门槛见 [开发台账](development-ledger-2.0.15.md)。

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
