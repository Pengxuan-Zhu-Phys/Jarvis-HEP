# Jarvis-HEP 2.0.15 开发台账

1. **Roadmap phase**：2.0.15 运行可靠性与 CLI 缺陷修复。
2. **Affected modules**：task validation、Worker 初始化与控制进程、monitor snapshot、CLI 参数和 plot 诊断；关联 Jarvis-Operas help 与 Jarvis-PLOT 卡片诊断。
3. **Public interfaces**：保留现有命令与 YAML 结构；修正验证结果、退出行为、JSON 输出和错误提示。具体 JSON 字段契约在 MON-01 实施时补充。
4. **Main risks**：误拒绝合法 Operas 表达式、误终止正常 Worker 启动、把缺失监控数据显示为零、修复只存在于本机依赖而未进入发行包。
5. **Verification performed**：严格验证、Worker 初始化、真实双 Worker 对照运行、受控 Redis monitor、plot 输入诊断、intent flag 顺序和 Operas help 的最小复现；回归与完整发行包安装验收通过；证据见 verification.json。

更新日期：2026-10-11，Australia/Adelaide。当前源码版本：**2.0.15（候选，未发布）**。目标发行版：**Jarvis-HEP 2.0.15**。

本轮确认的五类问题拆为七个修复任务，另设回归验收与发行准备任务。**全部修复及验收完成后，准备 2.0.15 发行候选，最终发布由用户执行。** 用户在 2026-10-11 授权开始实施，并明确保留最终发布责任；代理不上传 PyPI、不发布 GitHub Release、不推送发行 tag。

## 实施前提

上级 `AGENTS.md` 要求代码修改前阅读五份 canonical documents。下列文件在约定的 `docs/` 路径均缺失，代码实施前应定位或恢复权威原件，并按其要求同步后续阶段结论：

- `PROJECT_CONTEXT_FOR_CODEX.md`
- `PRD_JARVIS_HEP_MASTER.md`
- `ARCHITECTURE_REVIEW_JARVIS_HEP.md`
- `CODE_MAP_JARVIS_HEP.md`
- `IMPLEMENTATION_ROADMAP.md`

开始实施后再次检索工作区与本仓 Git 历史，未找到旧原件。上述五份文件已依据当前代码、现有文档和用户授权重建为当前维护基线，逐份标注重建来源，不声称还原历史原文。代码修改前阅读这些重建文件。工作区已有的 sampler、Diver 和 `pyproject.toml` 改动，应在发行范围核对时明确归属；若新增任务进入 2.0.15，先补入台账，再执行同样的完成门槛。

## 任务总表

“待开发”表示问题已确认但未修复。“待验收”表示实现完成，验收尚未全部通过。“完成”必须同时具备实现、验证证据和文档同步。“待发布”只有在全部前置任务完成后才能进入。

| ID | 优先级 | 任务 | 责任仓库与主要模块 | 前置任务 | 状态 |
| --- | --- | --- | --- | --- | --- |
| VAL-01 | P1 | 拒绝不存在的 Operas 限定名与不可编译表达式 | Jarvis-HEP：`cardload/task_validation.py`、`expression.py`、`operas_functions.py` | 实施前提 | 完成 |
| RUN-01 | P1 | Worker 初始化失败应反馈给控制进程并有界退出 | Jarvis-HEP：`runtime/worker.py`、`runtime/factory.py`、`runtime/_runtime_supervisor.py`、`runtime/_scan_driver.py` | 实施前提；与 VAL-01 配套 | 完成 |
| MON-01 | P2 | 带 scan ref 时兑现 `monitor --once --json` | Jarvis-HEP：`cli/client.py`、`dashboard.py` | 实施前提 | 完成 |
| MON-02 | P2 | Redis attach 正确读取 Worker 信息与快照时间 | Jarvis-HEP：`cli/client.py`、`dashboard.py`、`queue/redis_queue.py` | 实施前提；与 MON-01 配套 | 完成 |
| PLOT-01 | P2 | task card 传给 plot 时给出明确用途与纠正命令 | Jarvis-PLOT：卡片加载诊断；Jarvis-HEP：`cli/client.py` 转发诊断 | 实施前提 | 完成 |
| CLI-01 | P3 | 统一不同参数位置下的 intent 冲突诊断 | Jarvis-HEP：`cli/client.py`、`cli/help.py` | 实施前提 | 完成 |
| OPR-01 | P3 | 修复 Operas help 中 Rich 对象被字符串化 | Jarvis-Operas：CLI `_render_help_options`；Jarvis-HEP：依赖集成验收 | 实施前提 | 完成 |
| QA-02 | P2 | 修复 macOS/Windows scan CPU 旧值冻结并标明活进程估算 | Jarvis-HEP：monitor collector、resources | fast gate 暴露的已记录问题 | 完成 |
| SDK-01 | P2 | nested 停止日志不能将预算或自定义停止误报为收敛 | Jarvis-HEP：dynesty_sampler | 发行范围核对发现已合并 F3 | 完成 |
| QA-03 | P2 | 吞吐验收预热完整管线并增加测量点数，保持 1.25 倍门槛 | Jarvis-HEP：distributed acceptance tests | 短批次冷启动重复失败 | 完成 |
| QA-01 | 发布门槛 | 完成定向回归、真实运行和完整发行包验收 | Jarvis-HEP 与相关依赖 | VAL-01 至 OPR-01 全部完成 | 完成 |
| REL-01 | 发布门槛 | 准备 Jarvis-HEP 2.0.15 发行候选与用户发布步骤 | Jarvis-HEP：版本、`jhrel`、PyPI workflow | QA-01 完成；外部依赖修复已可打包 | 发行候选就绪，待用户发布 |

## 复现基线

可持久复现结果见 [reproduction-evidence.json](bugfix-2.0.15/reproduction-evidence.json)。环境为 Python 3.13.13、Jarvis-HEP 2.0.14、Jarvis-Operas 1.4.1、JarvisPLOT 2.1.9、Typer 0.25.1、Click 8.4.2、Rich 15.0.0。

| 现象 | 已观察结果 | 证据边界 |
| --- | --- | --- |
| LogLikelihood 写入 `nosuch.f(z)` | `validate --strict --json` 退出 0，`ok: true`，`issues: []`；Worker 初始化报 `AttributeError: 'Symbol' object has no attribute 'f'` | 同类错误放在 `Sampling.Mapper` 中已被 `JV2-MAP-007` 拒绝，不能概括为所有表达式位置都漏检 |
| 真实双 Worker 运行错误卡片 | 35.15 秒仍未退出；completed、failed、running 均为 0，task queue 为 2；Worker 已重启 4 次，scan mode 为 degraded | 测试主动发送 Ctrl+C，退出 130，并完成控制锁释放与 managed Redis 关闭；用户报告的三分钟时长未重复跑满 |
| 同环境正常卡片 | 完成 2 个样本、failed 为 0、自然退出 0 | 作为真实运行对照，排除整个测试环境都无法运行 |
| monitor 带 ref 与 JSON | 返回文本、`timestamp: None`、`workers: 0/0 alive` | 使用受控进程清单与内存 Redis，其中已有两份 Worker heartbeat board；不等同于原始现场的进程盘点 |
| monitor 不带 ref 与 JSON | 返回可解析的 JSON 扫描列表 | 之前“不带 ref 也坏”的判断已撤回，不列为缺陷 |
| task card 传给 plot | 退出 1，`Jarvis plot failed: 'DataSet'` | 需要专用绘图 YAML，缺陷在输入类型与纠正提示不足 |
| 冲突 flag 的位置 | 子命令前报 `conflicting CLI intents`，子命令后报 argparse `unrecognized arguments` | 两种方式均拒绝执行，差异在诊断一致性 |
| `jopera init -h` | 退出 0，但 `--force` 说明包含 `<rich.columns.Columns object at ...>` | 属于已报告 Operas help 问题的同一修复项，不重复计算任务 |

## 修复要求与验收

### Operas 表达式验证

VAL-01 应覆盖 Worker 使用的 LogLikelihood、Operas 输入等表达式字段。对未知命名空间、已有命名空间下未知成员以及编译失败，返回带 YAML 字段路径和纠正建议的验证错误；验证过程不能执行用户的数值函数。

验收应同时证明 `nosuch.f(z)`、`helper.nosuch(z)` 被拒绝，合法注册函数与常量仍能通过，原有 Mapper 诊断保持有效。`validate --json` 的 stdout 必须可解析，验证失败退出非零，且不启动 Redis 或 Worker。

### Worker 初始化失败退出

RUN-01 应让控制进程识别 Worker 启动失败的原因。确定性的配置或初始化失败应结束本次运行，不能只进入进程恢复与 archive 等待。正常运行中的可恢复 Worker 崩溃仍应保留恢复能力。

验收需独立注入初始化失败，覆盖“验证后配置变化”或直接 Python 调用等绕过预检查的路径：运行在明确的启动失败期限内退出非零，保留原始异常与相关模块信息，终止已启动的 Worker 和 FileOperation，释放控制锁并关闭本次 managed Redis。启动失败应体现为运行失败，不能虚构未执行样本的完成或失败计数。另验证正常双 Worker 扫描与既有崩溃恢复仍通过。

### monitor JSON 输出

MON-01 应让带 ref 与不带 ref 的 JSON 模式都有稳定输出。带 ref 时输出一个完整 JSON 文档，不能混入 `R1: ...` 或文本快照标题；不带 ref 时保留扫描列表语义。具体字段结构及无效 ref、元数据不匹配、Redis 连接失败时的错误输出契约，应写入监控文档并由测试固定。

验收至少包含 `Jarvis monitor --once --json R1`、按扫描名称选择、不带 ref 列表以及无效 ref。命令结束后 stdout 能按约定解析，错误情况退出非零。纯文本 `--once` 与 TUI 入口继续可用。

### monitor Worker 数据与时间

MON-02 应从选中扫描的 OS 进程清单获取 Worker owner IDs，再读取明确对应的 Redis boards；不要依赖新建 TaskFactory 的进程内快照。快照时间应表示本次采集时间，缺失或不确定的信息应明确表达，不能默认伪装成 `0/0 alive`。

验收需用真实双 Worker 扫描确认数量、PID 和状态一致，并覆盖进程退出、缺失 board、过期 board、零 Worker 与 Redis 故障。快照读取保持只读，不额外扫描所有 Worker Redis keys，也不启用 TUI inflight telemetry。

### plot 卡片诊断

PLOT-01 应识别扫描 task card 与绘图 YAML 的差异。缺少 `DataSet` 时应说明所需文件类型；识别出 task card 时，提示先执行 `Jarvis gen-plot-yaml TASK.yaml`，再将生成的绘图 YAML 交给 `Jarvis plot`。保留缺失文件、损坏 YAML 和真正绘图错误各自的原因。

验收覆盖 task card、合法绘图 YAML、损坏 YAML 及缺失文件，确认诊断不再是孤立的 `'DataSet'`，已有 `jplot` 参数转发继续有效。若修复位于 Jarvis-PLOT，须验证 HEP 使用可分发的依赖版本也得到相同结果。

### CLI intent 冲突

CLI-01 应在 argparse 退出之前或适当的参数归一化阶段，识别 Jarvis 自己的冲突 intent。已识别的冲突在子命令前后位置均返回 `conflicting CLI intents` 与 usage 退出码；真正未知的参数仍按未知参数诊断处理。

验收覆盖 `Jarvis --monitor validate TASK.yaml`、`Jarvis validate TASK.yaml --monitor`、legacy bare YAML 与冲突 flag 的组合，以及正常单 intent。plot、portal、operas 等依赖 CLI 的自由参数转发需要专项回归，不能把它们合法的参数错误当成 HEP intent。

### Operas help 渲染

OPR-01 应保留 Rich renderable 的正常渲染，修正 `_render_help_options` 对危险提示执行 `Text(str(help_text))` 的路径。该项归属 Jarvis-Operas；同一 renderer 中其它触发该路径的选项应一起验证。

验收包括 `jopera init -h`、已有受影响 help 页面、`Jarvis operas init -h`，确认说明可读且不包含 `<rich.* object at ...>`。修复后的 Operas 版本必须能够从发行包安装，再确定 HEP 所需的最低依赖版本。

## 回归验收与发布顺序

QA-01 需保存每个修复项的提交、测试命令、结果和证据路径。最小验收路径包括相关 CLI、validation、expression、dashboard、Worker failure 与 distributed acceptance 测试，以及本台账的真实运行对照；全部 fast gate 也必须通过。历史审阅中的失败记录应在当前发行候选上重新核实，不能拿历史复现完成或测试脚本退出 0 代替回归通过。

REL-01 按现有 [完整发行脚本](../scripts/full-release.sh)与 [PyPI workflow](../.github/workflows/publish-pypi.yml)执行：

1. VAL-01 至 OPR-01、QA-01 全部完成；核对这次发行包含的所有变更，并排除尚未验收的改动。相关外部依赖修复须已进入可分发版本。
2. 更新 `pyproject.toml`、`jarvishep2/__init__.py` 与 README 的版本为 2.0.15，补充发行说明，确认 tag、运行时版本和 distribution metadata 一致。
3. 使用 `jhrel` 构建完整 wheel 和 sdist，在干净环境安装全部依赖，运行 `pip check` 与安装 smoke，保留 `release-verification.json` 和 `SHA256SUMS`。完整依赖安装方式按现有脚本执行。
4. 从干净安装结果复验本轮 CLI、双 Worker 运行和外部依赖 help 修复。构建源码来自已纳入发行的 tracked files，不能依靠当前源码 checkout、editable sibling 或未跟踪文件才能通过。
5. 代理交付 wheel、sdist、验证报告、SHA256SUMS、依赖候选及发布步骤，将 REL-01 标为“发行候选就绪，待用户发布”。如果依赖修复尚未公开发行，列明先发布依赖的顺序与包版本。
6. 用户自行创建 `v2.0.15` tag 和 GitHub Release，执行现有 PyPI 发布流程。公开发布与安装核验成功后才将正式发布记为完成。

在步骤 1 满足前，版本维持 2.0.14。候选准备与正式发布分别记录；代理本轮工作完成条件是交付通过验证的发行候选，正式发布由用户执行。

## 完成记录

| 日期 | 任务 | 实现提交 | 验证命令与证据 | 结果 |
| --- | --- | --- | --- | --- |
| 2026-10-11 | 台账建立与缺陷确认 | 文档阶段，无代码修复 | [复现证据](bugfix-2.0.15/reproduction-evidence.json) | 七个修复项待开发；QA-01、REL-01 待执行 |

实施过程中持续更新任务状态与完成记录，阶段结束时同步权威维护文档。任何新增的 2.0.15 任务都要具备明确归属、验收条件和状态。

### 第一阶段（2026-10-11）

七个修复已实现，相关 fast 回归 `148 passed, 30 subtests passed`；真实 Worker/恢复/分布式回归及完整 fast gate 进行中。MON JSON ref 成功契约为 `{ok:true, reference, scan_name, sampler, factory, workers, calculators, occupancy, samples, resources, queues, op_counts, proc_core, proc_archiver, proc_redis, timestamp}`；错误为 `{ok:false, reference, error:{message}}`，退出 1 或无效 ref 退出 2。无 ref 保留扫描数组。OS 可见进程决定 alive/PID；未知数量为 null；timestamp 为 Unix 采集时间，board age 单独表示。

QA-02 与 SDK-01 为当前 tracked 发行范围内已有审阅问题，纳入本次发布门槛。F1/F2/F5 属未跟踪 Diver，不纳入候选包，保留原始工作区内容。CPU 验收区分 Linux reaped-child accounting 与 macOS/Windows live-tree estimate，不声称捕获两次刷新之间完全结束的子进程。

2026-10-12：最终发行源码 fast gate 1272 passed、3 skipped、141 subtests；吞吐及启动失败独立验收 2 passed。开始 2.0.15 完整安装包验收。

## 最终交付（2026-10-12）

VAL-01 至 OPR-01、QA-02、SDK-01、QA-03、QA-01 完成。REL-01 的候选准备完成，正式发布留给用户，未创建 tag/Release、未推送、未上传。

最终发行源码 fast gate：**1272 passed、3 skipped、141 subtests**。运行/恢复/分布式主体 97 passed；最终吞吐与启动失败 2 passed；控制面 140 passed、30 subtests。完整 HEP wheel 在干净环境安装全部依赖，Twine、pip check、双 Worker、ref/name JSON、坏表达式、plot/intent/Operas help、项目重载 smoke 通过，另有 **22 项安装包回归通过**。SDK 完整安装测试 146 passed；Operas CLI 43 passed。

发行文件位于 `dist/2.0.15/`；依赖候选位于 `dependencies/sampler` 和 `dependencies/operas`，每个目录都有 wheel/sdist、校验和和验证报告。见 [验证结果](bugfix-2.0.15/verification.json)、[发行源码清单](bugfix-2.0.15/source-manifest.json)、[用户发布步骤](bugfix-2.0.15/release-notes.md)。实现统一进入本次 2.0.15 本地修复提交（以本地 Git log 为准）。

限制：本机 macOS/Python 3.13 验证不代表跨平台 CI 已运行；CPU 活进程估算可能遗漏短命进程；nested 未报告终止原因时只报告 engine returned；启动 ready 上限默认 120 秒。Diver 原始改动及其 F1/F2/F5 未进入发行包。
