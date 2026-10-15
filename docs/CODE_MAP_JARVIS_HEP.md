# Jarvis-HEP 代码地图

2026-10-11 重建；本次 2.0.17 当前范围为公开 master 0bbcf675 + 不依赖 SDK 的可靠性/CLI agent/Redis 移植。顶层若干模块为兼容导入 shim，编辑时定位实际实现。外部 SDK/plugin 模块不在本次发行范围；旧 SDK 候选记录在原集成分支归档。

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

本轮行为约束见 [产品要求](PRD_JARVIS_HEP_MASTER.md) 和 [架构审阅](ARCHITECTURE_REVIEW_JARVIS_HEP.md)；任务与发行门槛见 [开发台账](development-ledger-2.0.15.md)。

## 新增 CLI 与 YAML 职责

| 职责 | 实现 | 主要验证 |
| --- | --- | --- |
| agent 合同/发现/离线 schema/模板 | `cli/agent.py`、`agent_search.py`、`agent_artifacts.py` | agent CLI/artifacts/experience 回归 |
| 卡片解释与草稿 | `cli/agent_inspect.py`、`agent_request.py` | agent inspect/request 回归 |
| 显式重复键和诊断导航 | `cardload/yaml_loader.py`、`task_config.py`、`cli/run.py`、`man_codes.py` | YAML/validate 回归 |
| Redis 系统/版本/安装启动说明 | `queue/redis_notice.py`、兼容 `redis_notice.py`、`cli/redis.py` | `test_redis_notice.py`、`test_redis_environment.py` |
| Redis 父子 help/JSON/agent 发现 | `cli/help.py`、`client.py`、`agent.py`、`agent_search.py` | parser/help/agent 与 installed Redis gate |



## 2026-10-15 SDK-free 公开发行代码地图（最新）

| 职责 | 当前实现范围 | 验收 |
| --- | --- | --- |
| 内置采样与 runtime 基线 | 公开 master `0bbcf675`；不移植 SDK integration `0d859502` / `52876060` | 内置采样、真实 Redis/spawn、checkpoint/恢复回归 |
| 可靠性与 YAML | 不依赖 SDK 的 worker 启动/表达式/monitor/CLI 修复、duplicate YAML 诊断 | 定向可靠性/YAML 与完整 fast gate |
| CLI agent 与 Redis | agent 查询/离线 schema/模板/inspect/邮件草稿；Redis 正式子命令/help/JSON | installed agent/Redis gate，保持仅打印安装/启动指导 |
| 公开发行 | `pyproject.toml`、`jarvishep2/__init__.py`、构建/安装 smoke、Trusted Publishing | 2.0.17；无 SDK Requires-Dist；公共 PyPI 全依赖安装/pip check |

当前验收：见下方已完成的源码/安装验收与发布状态。旧 `1c8c82fc` 候选不再是发行源码，其测试/checksum只作历史；最新产物、源码与公开安装结果见 [发行记录](release-2.0.17.md)。

## 本次源码与安装验收状态

干净源码 `e8af5221`：1426 passed、3 skipped（缺少 sibling Jarvis-Examples）、164 deselected、477 subtests；335 项代码/测试文件与提交逐字节一致。全新隔离环境从公共 PyPI 解析全部依赖，pip check 与 Twine 通过，SDK 模块/发行包均不存在；真实 Redis/spawn/双 Worker/monitor/项目打包重跑、installed agent 24 项、installed Redis 39 项、发行回归 22 项均通过。

源码分离提交 `3b9d41d0`，基础 ValidationReport 兼容修正 `e8af5221`；修正只调用公开基线的 warnings 接口，新增正常文本验证和警告回归。wheel 的 276 个运行文件、sdist 的 401 项 tracked 文件与提交一致，无 SDK/plugin bridge、对应测试或 Diver。当前本地验收完成，正式发布执行中；实际 PyPI 成功与公开安装结果随后记录。见 [发行记录](release-2.0.17.md) 与 [验收证据](release-2.0.17-verification.json)。
