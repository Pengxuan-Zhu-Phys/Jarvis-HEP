# Jarvis-HEP 代码地图

2026-10-11 重建。2.0.17 公开发行范围是公开 master 0bbcf675 加上不依赖 SDK 的可靠性、CLI agent 和 Redis 移植。当前 master 另含 Jarvis-HEP-Sampler 插件运行时：依赖 `Jarvis-HEP-Sampler>=1.0.3`，不设上限，命令入口是 `cli/sampler.py`。顶层若干模块为兼容导入 shim，编辑时定位实际实现。

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
| 发行准备 | `pyproject.toml`、`jarvishep2/__init__.py`、`jhrel`（独立脚本，唯一选项 `--dry`）、`scripts/smoke_install.py` | `tests/test_distribution_version.py`、`test_package_layout.py` |

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



## 2026-10-15 SDK-free 公开发行代码地图（历史发布阶段）

| 职责 | 当前实现范围 | 验收 |
| --- | --- | --- |
| 内置采样与 runtime 基线 | 公开 master `0bbcf675`；不移植 SDK integration `0d859502` / `52876060` | 内置采样、真实 Redis/spawn、checkpoint/恢复回归 |
| 可靠性与 YAML | 不依赖 SDK 的 worker 启动/表达式/monitor/CLI 修复、duplicate YAML 诊断 | 定向可靠性/YAML 与完整 fast gate |
| CLI agent 与 Redis | agent 查询/离线 schema/模板/inspect/邮件草稿；Redis 正式子命令/help/JSON | installed agent/Redis gate，保持仅打印安装/启动指导 |
| 公开发行 | `pyproject.toml`、`jarvishep2/__init__.py`、构建/安装 smoke、Trusted Publishing | 2.0.17；无 SDK Requires-Dist；公共 PyPI 全依赖安装/pip check |

当前验收：见下方已完成的源码/安装验收与发布状态。旧 `1c8c82fc` 候选不再是发行源码，其测试/checksum只作历史；最新产物、源码与公开安装结果见 [发行记录](release-2.0.17.md)。

## 本次源码、CI 与公开安装状态

源码功能提交 `3b9d41d0`、验证接口兼容修正 `e8af5221` 已进入 master；公开 tag `v2.0.17` 指向 `8e0add06`。macOS 与 Linux CI 均为 1426 passed、3 skipped（缺少 sibling Jarvis-Examples）、164 deselected、477 subtests；分别用时 205.99s / 265.40s。Linux CI 的完整安装/pip check、真实 Redis/spawn/双 Worker/monitor/项目重跑、agent 24、Redis 39、发行回归 22 均通过。

[GitHub Release](https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP/releases/tag/v2.0.17) 已公开。OIDC publisher 因 PyPI 未匹配 Trusted Publisher 配置而报 `invalid-publisher`；本次使用既有本地 PyPI 凭证上传成功 CI 的原样产物，没有重建。PyPI wheel 已公开，全新环境仅从公共索引安装完整依赖，pip check、无 SDK、276 个运行文件身份比对、Redis 39 / agent 24 均通过。PyPI 源码包上传未获服务器响应；最后一次有界请求在 36.913 秒后 ReadTimeout，全部上传进程已停止。PyPI wheel 与公开安装已完成，已验收源码包可从 GitHub Release 下载。

发布 wheel SHA256 为 `25257b646c1a7ce828dae136e1e111de7c70d749d7b5b6857b7c687f227fb3ae`；CI sdist 的预期 SHA256 为 `1843e006921161dd6a46507b7b7c12b52b48990b9e3af1c76541788dafbe9915`，未在 PyPI 公开；GitHub Release 另提供通过本地完整验收的源码包。发布前本地产物的 a823… / 9171… hash 仅作本地验收记录，不能当成公开 CI 产物 hash。完整状态见 [发行记录](release-2.0.17.md) 与 [验收证据](release-2.0.17-verification.json)。

## 2026-10-16 HEP 修复代码地图（历史阶段，CLI 边界见 2026-10-17）

| 责任 | 实现 | 回归 |
| --- | --- | --- |
| 当前运行策略与恢复回退 | `sampling/checkpointed_sampler.py`、`runtime/_runtime_supervisor.py`、`runtime/_resume_service.py` | fresh/resume/拒绝 checkpoint/check × Grid/Random/Bridson/CSV |
| 仅 resume 过滤旧结果 | `sampling/fixed_set_sampler.py`、`sampling/csv_sampler.py` | `test_hep_review_regressions.py`；真实 Redis/spawn/HDF5 旧数据重算 |
| 当时 HEP 的保留名和旧注册冲突处理（现已替换） | `sampler_catalog.py`、`distributor.py`、`runtime/core.py` | 现行入口见下文 SDK 默认保留名与 HEP 注册表报告 |
| 当时 HEP 的 CLI 调用边界（现已替换） | `cli/sampler.py` | 现行 dimension、file-only check、登记诊断由 SDK 处理 |
| 既有 resume 故障矩阵 | `tests/test_plugin_resume_fault_matrix.py` | 独立 sampler 单测显式开启 resume 模式，保留稀疏 Failed/prefix/背压断言 |

发布脚本、SDK 仓库和 Diver 移植均不属于本次修改。完整命令与结果见 [本轮验证证据](hep-review-2026-10-16-verification.json)。

## 2026-10-16 续跑协议代码地图（四种内置方法已实现）

| 责任 | 实现位置 | 回归 |
| --- | --- | --- |
| 没有 checkpoint 就不续跑 | `runtime/_resume_service.py` `_preload_resume_checkpoint` | `enabled: false` / 无 state.pkl 的 `--resume` |
| check 模式关闭 checkpoint | `runtime/_scan_driver.py` `_apply_check_modules_runtime_policy` | check 不写 state.pkl / points.npy |
| fresh 初始 checkpoint、`replacing_existing`、resume 进度来源 | `runtime/_runtime_supervisor.py`（Archiver 启动后的 DATABASE 载入与 `advance_to_persisted_prefix`）、`runtime/_resume_service.py` `save_runtime_checkpoint` | 心跳前/后中断的真实 Redis 续跑 |
| Archiver 初始水位与清尾 | `io/archiver.py`（按 fresh/marked resume 设置初始前缀和 `replace_duplicates`）、`io/database.py` `trim_scan_records`（跨 shard 清尾与去重） | fresh 前缀从 0；短扫描覆盖长扫描无尾部 |
| 卡片与点集指纹 | `sampling/runtime_checkpoint.py`（`build_payload` 与核对函数）、`sampling/fixed_set_sampler.py`、`sampling/grid.py`、`sampling/csv_sampler.py` | 改 Bounds/selection/CSV 拒绝；改 workers 通过 |
| Bridson 点集文件 | `sampling/bridson.py`（`initialize`、export/import）、`sampling/runtime_checkpoint.py` `checkpoint_path` 同目录 | points.npy 只写一次；resume 不生成；损坏/缺失从头扫描 |
| 不累积提交 UUID | `sampling/stateless_batch.py`、`sampling/fixed_set_sampler.py`、`runtime/_scan_driver.py` 结束统计 | state.pkl 与内存不随点数增长 |

新增验收位于 `tests/test_resume_protocol.py`：4 × 2 真实中断续跑、真实长短扫描、check 禁写文件、卡片/点集拒绝、百万点 mmap 和 checkpoint 大小。归档 writer 的 `drop_matching_records` 仅处理 live shard，故清尾定位已更正为跨 shard 的 `trim_scan_records`。

## 2026-10-17 注册表与命名代码地图（N1–N6 已完成）

| 责任 | 实现位置 | 回归 |
| --- | --- | --- |
| SDK 默认保留名单（N1） | SDK `metadata.py`（`RESERVED_METHOD_NAMES`、`check_method_name`）；HEP `sampler_catalog.py` `reserved_method_names` | 两个 load 拒绝重名；HEP builtin 与 builtin-plugin ⊆ SDK 名单；普通 init/load/check 无豁免 |
| 官方 builtin-plugin 的精确低层豁免（N1） | SDK `metadata.py`、`loader.py` 的 `trusted_builtin_method`；HEP `sampler_catalog.py` 官方来源加载、`sampling/plugin_identity.py` 官方旧绑定判定、`sampling/plugin_sampler.py` 扫描副本恢复 | `LatinHypercube` 的旧扫描副本及保存 hash 在包路径/字节变化后仍可续跑；项目本地绑定可编辑，不能防有意伪造 |
| check 只接受文件、名称第一项、CLI 诊断（N2） | SDK `cli.py`（`_run_check`、`--dimension` 解析、`_run_info` 与统一渲染）、`testing.check_plugin` | 方法名参数 exit 2；重名文件后续项未执行；`--dimension 0` 与负 `--seed` exit 2；损坏登记的 list/info/remove exit 2；`--version --json` |
| remove 只删外部登记及托管快照 | SDK `registry.py` `remove`、`entries`（不做名字检查）；HEP `cli/sampler.py` 原样转发，不动内置方法或用户源文件 | 重名条目 `grid` 可读取并 `remove grid` |
| HEP CLI 转发（N4） | `cli/sampler.py`：`_load_sampler`、`_guard_check`、info 中的重名拦截已删除；load/check 转发 SDK | `Jarvis sampler load/check` 与 `jsampler` 输出一致 |
| HEP 帮助解析器（N4） | `cli/help.py` 的 SDK 解析器拷贝 `build_sampler_help_parser` 已删除；帮助与 `cli/agent.py` 命令目录使用 SDK `build_parser("Jarvis sampler")`；错误用 SDK 面板，SDK CLI 不可导入时由 HEP 输出等价面板或 JSON。规则见 [CLI 规范](JARVIS_CLI_DESIGN.md) §10 | 无扩展子命令替换程序名后帮助一致；根帮助及 list/info 通过 SDK 扩展点增加 HEP 内置说明 |
| 注册表报告（N5 已实现） | `sampler_catalog.py`：`RegistryReport` / `registry_report()`；`runtime/core.py` `preflight_sampler_plugin` 在 Redis 前核查；`cli/run.py` 将注册表 preflight 错误映射 exit 2 | 重名任意卡片退出；无关快照问题警告；真实 run/check 入口回归 |
| 删除分散拦截（N5 已实现） | 已删除 `assert_no_reserved_name_collision` 及 `distributor.py` `set_method`、`runtime/core.py`、`cli/sampler.py` 的调用；直接 Distributor 内置优先 | 新报告回归替代旧冲突拦截回归 |
| 报告展示（N6 已实现） | `cli/sampler.py` 的 list/info 用 SDK 登记 JSON 载荷与 HEP `registry_report()` 合并；文本主表和 JSON 共用报告状态 | `Grid`/`grid` 旧登记的 name-conflict、精确 remove 命令、移除后内置方法恢复 |
| 发布与依赖下限（N3 已完成） | 维护者已公开 SDK 1.0.3；HEP `pyproject.toml`、`sampler_sdk.py`、两份 `environment_default.yaml` 均设 `Jarvis-HEP-Sampler>=1.0.3`，不加上限；SDK README 同步 | 干净 Python 3.11 环境从本地 HEP 源码实际安装 66 个包，依赖从公共 PyPI 解析并选中 SDK 1.0.3；`pip check`、已安装运行时版本门禁及双入口 JSON 版本通过；HEP 未公开新版本 |
| Diver 移植（暂留） | `sampling/Source/Diver/`、`tests/test_diver_engine.py`、`pyproject.toml` package-data | 未注册为 `Sampling.Method`；未做数值验证 |

本地验证：SDK 源码 234 passed、Ruff 通过；HEP 最终定向集成 169 passed / 4 deselected，快速全套门禁 1597 passed / 196 deselected（先于最后 CLI 边界修复，随后定向复跑）。HEP N5 非 CLI 专测 5 passed，旧冲突的真实 run/check 与 list/info/remove 入口 4 passed。N4 的原生帮助、错误与 JSON 转发已有定向回归。

## 2026-10-17 复查遗留清理代码地图（C1 已实现）

| 责任 | 实现位置 | 回归 |
| --- | --- | --- |
| sampler 持久化绑定（resume 模式、persisted UUID/prefix/pending、按前缀推进） | `runtime/_runtime_supervisor.py` 原两处调用同一 `Jarvis2Core._bind_sampler_persistence()`；`runtime/_resume_service.py` `apply_resume_checkpoint` 及 Core 包装已删除。为兼容现有 MagicMock 回归，两处显式调用 Core 类方法体。队列清理仍在 supervisor 的 `prepare_resume` 调用处，重新提议标记仍在 checkpoint 导入处 | 指定三组非 slow 回归 58 passed；Core collaborator 17 passed；真实 Redis 慢测 20 passed / 54 deselected（278.16 秒） |
