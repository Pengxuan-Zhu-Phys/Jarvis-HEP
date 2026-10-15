# Jarvis-HEP 产品要求

2026-10-11 重建的维护基线，来源为当前 README、安装文档、运行接口与 [2.0.15 台账](development-ledger-2.0.15.md)，不声称还原旧 PRD。

## 核心行为

- `validate` 只验证卡片，不启动 Redis 或 Worker；配置错误退出 2，并给出 YAML 字段路径与纠正建议。
- `run` 返回真实运行结果：成功为 0，运行错误或样本失败为 1，使用或配置错误为 2，中断为 130。
- Worker 启动失败必须反馈给控制进程，结束不可运行的扫描并清理本次资源；正常扫描中的可恢复进程故障保留恢复能力。
- 监控读取选中扫描的数据。缺失或未知信息不伪装成零；JSON 模式输出可解析的 JSON。
- plot 使用绘图卡片。收到扫描卡片时说明用途并指向 `gen-plot-yaml`。
- help 和冲突参数的错误提示必须可读、稳定；依赖 CLI 的自由参数保持转发。

## 2.0.15 验收

VAL-01、RUN-01、MON-01、MON-02、PLOT-01、CLI-01、OPR-01 的具体验收以台账为准。合法 Operas 注册函数与常量、Mapper 校验、正常双 Worker 扫描和既有 Worker 恢复必须保持可用。

发行候选需通过定向回归、真实运行、fast gate、完整 wheel 与 sdist 检查、干净安装与依赖检查。最终发布由用户执行；候选准备完成与公开发布成功是两个独立状态。

## 边界

本轮保持现有采样算法和 YAML 结构。额外发现的问题先记录范围、证据和验收，再决定是否作为发行阻断项修复。不能把历史测试结果或源码环境成功等同于当前发行包通过。

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

## 2026-10-13 AI agent CLI 要求

- 顶层 help、man 索引及 README 必须可发现 `Jarvis agent --json`/`Jarvis man agent-readme --json`；合同包含安装版本、contract version、argv、查询路由、命令效果、授权边界和退出码。
- agent 应能导出完整本地 bundled schema、读取准确示例正文/前置条件，区分完整 task card 和 Sampling fragment；插件字段 schema 由 sampler man 暴露，不声称 bundled schema 覆盖任意扩展运行语义。
- 编写闭环为查询→编辑→validate JSON→manual_argv→修正→复验。不得因编写任务自动运行 check/run、安装/加载代码或改变物理假设；已获授权不反复确认。
- JSON-capable 查询的参数、加载和内部失败给出单个可解析结果与非零退出码，保留既有成功字段；help 输出遵循各命令约定。
- 邮件反馈收件人 `Jarvis.HEP.AI@mail.grokbot.com`、精确 marker `Resquest@Jarvis-HEP`。输出含版本、环境、最小 YAML、命令诊断和实际/期望/约束的草稿；不自动发送、不推断维护者身份。

验收及已知重复键/缓存限制见 [报告](cli-agent-integration-report.md)。本阶段不改变版本、既有候选或公开发布状态。

验收完成并记录外部例外：194 定向 tests/316 subtests、13 安装包 CLI 检查通过；full fast gate 1345 tests/437 subtests 通过，已知 sibling Examples 非 ASCII 描述导致 1 项失败。外部案例保持原文，不声称整个工作区 gate 全绿。

## 第二轮新增要求与完成情况

- 短合同保留全部规则与工作流，旧默认完整 JSON 继续可用；section/command/search/schema-index 让 agent 按需获取证据，search 只能返回查询/帮助 argv。
- inspect 解释具体卡片的实际默认文件、任务覆盖、规范化值和来源；明确未做 runtime 资源选择或用户算子可用性检查，大输入截断显式可见。
- fragment 的组合/验证目标、算子注册/导入边界及手册子路径必须准确。未知后缀和冲突 selector 返回失败及可用父页，不默默丢参数。
- 重复显式 YAML 键必须在任务与默认环境加载时拒绝，并给字段和两处位置；合法 merge override 保持。此要求取代第一轮仅约定唯一键的限制。
- request 可填写 YAML/命令/诊断/实际/期望/约束/安装方式；字面文本不执行/不隐含读文件，draft_ready 仅代表内容齐全。

验收 332/373 定向、22 安装包检查通过；fast gate 1416/478 通过及 1 个已知外部编码例外。报告与第二轮 verification 已记录严格重复键对旧案例的兼容影响；无发布动作。


## 2026-10-13 2.0.16 提交与发行

用户已授权 commit 与新版本发布，本阶段取代此前仅准备候选的权限边界。版本 2.0.16 包含两轮 CLI agent 改进及已提交的 2.0.15 可靠性/sampler integration 基线；保持原有 YAML 结构，显式重复键改为 JV2-YAML-003 错误。Diver 草稿、其 package-data 及未验收的 sampler 审阅文件保留且排除。

从干净 commit 副本构建完整 wheel/sdist，执行源码 fast gate、全依赖安装、pip check、真实 Redis/spawn smoke、旧可靠性回归和新增 installed agent CLI smoke。发布正常使用 GitHub Release/Trusted Publishing，不使用本地 token 上传。当前依赖检查：Operas 1.4.2 已公开；Jarvis-HEP-Sampler 0.1.0 在 PyPI 缺失，须先解决才能正式发布 HEP。具体状态与证据见 [2.0.16 发行记录](release-2.0.16.md)。

2.0.16 发行候选验收完成（源码 commit `50631bae` 已推送到 `feat/sampler-plugins-s1`）：带 Git 元数据的干净 checkout fast gate 1396 tests/478 subtests、3 skipped、166 slow deselected；完整包全依赖干净安装、pip check、真实 Redis/spawn/项目重载 smoke、24 项 agent CLI 和 22 项安装包回归均通过。wheel/sdist 与已提交源码逐文件一致并排除 Diver。SDK 使用本地已验收候选，其余依赖来自官方 PyPI，因此仍不能声称仅用公开索引安装成功。SDK 0.1.0 尚未公开，正式 HEP 发布等待依赖；准备 draft Release。详见 `release-2.0.16.md` 与 `release-2.0.16-verification.json`。


## 2026-10-15 Redis 环境引导要求

- 安装后的首次 CLI 调用完成时，复用一次 PATH-only 提示；不声称 `pip install` 执行 post-install 检测。缺失时将系统适用的安装/启动说明写入 stderr，保持 JSON stdout 可解析，marker 兼容既有路径。
- `Jarvis redis` 默认 `info`；`info` 展示系统/发行版、兼容 executable 和版本及安装/启动说明。`check` 每次重新查 PATH，以有超时的 `--version` 查询验证本地 executable；缺失、执行失败、未知/过旧版本（Redis < 6.2）或原生 Windows 退出 1。
- `install` / `start` 仅打印可复制指引，不运行包管理器、不启动服务、不执行版本查询；所有 topic 支持 `--json`。环境查询不连接 Redis、不检查服务或端口、不启动 Worker/扫描，成功不代表真实扫描已验收。
- 平台说明按 `/etc/os-release` 的实际 ID 优先于 ID_LIKE 选择，再按发行版版本和可用包管理器处理；原生 Windows 使用 WSL，HEP 与 Redis 都在 WSL 中运行。Redis/Valkey 选择必须匹配发行版维护的软件包。
- 默认 Redis 过旧的 Ubuntu 20.04/22.04、Debian 11 应打印官方签名 APT 仓库步骤，未核验的更早版本只给文档/升级建议；RHEL 8 选择 `redis:6` 并解释 updated AppStream / 旧 stream 切换要求，不能把仍为 Redis 6.0 的固定旧仓库当作合格。
- 正常扫描由 Jarvis 自动启动/清理本次 Redis，不要求全局服务。可选手动调试命令仅绑定 127.0.0.1、示例端口 6380、关闭持久化，并说明 Ctrl+C 停止。

验收状态：最终定向 177 项通过（含真实 Redis 生命周期）；临时 wheel 安装 9 项 CLI 检查通过。完整 fast gate 为 1464 passed / 478 subtests、1 项已知 sibling iDM 非 ASCII 描述失败；最终增补平台用例由定向复验覆盖。详见 [验收记录](redis-environment-verification.json)。此阶段仅增量 CLI/引导与文档修正，YAML 结构、扫描生命周期及发行状态不变。


## 2026-10-15 Redis 命令组一致性要求

- 更正第一阶段 positional `topic` 仅能选择行为、未形成真正子命令的接口理解。Redis 必须使用与现有 CLI 一致的命令组结构，`Jarvis redis --help` 显示 info/check/install/start 四个 Commands。
- `Jarvis redis COMMAND --help` 必须展示该子命令自己的说明和参数（包括 `--json`），不再返回同一父级帮助页；帮助发现不执行环境检查。
- 保留裸命令与 `--json` 的 info 默认；`--json` 在子命令前后均有效。agent catalog/search/`--command` 应提供四个 exact child command，而非只有父级。
- 本次只修正命令组织、帮助和发现；安装平台逻辑、报告含义、退出码、仅打印 install/start 指导及扫描边界不变。

第一阶段的 177 项定向、9 项安装 smoke 与 fast gate 仅为历史验收，记录见 [Redis 环境验收记录](redis-environment-verification.json)。命令组补正验收：202 项定向通过，临时 wheel 安装 25 项检查通过；完整 fast gate 为 1509 passed / 478 subtests、1 项既有外部 iDM 描述编码失败。详见 [本次验收](redis-cli-help-verification.json)；版本和公开发行状态不变。


## 2026-10-15 HEP 2.0.17 Redis-only 发行要求

用户授权只提交 Redis 改进并发布新的 HEP，明确忽略本轮 Sampler 改动。发行包含 Redis 环境引导、正式子命令/help/JSON/agent 发现及对应测试/文档，版本更新为 2.0.17；工作区 Sampler/runtime/Diver 改动不得混入，baseline 已有的 SDK 集成不因本次授权移除。

验收必须针对干净 committed source 的完整 wheel/sdist，覆盖安装包 Redis 帮助/JSON/agent、pip check、真实 Redis/spawn 运行和相关回归。允许以保持原样的本地已验收 SDK 0.1.0 wheel 验证候选，明确区别于公开索引完整安装；Sampler 仍返回 PyPI 404 时，只能准备候选/draft，正式发布等待依赖可用。本次不改动或发布 SDK。

本阶段验收：干净 detached 源码 `1c8c82fc` 全量 fast gate：1489 passed、3 skipped（缺少 sibling Jarvis-Examples）、166 deselected、478 subtests；完整隔离安装与 pip check、真实 Redis/spawn/project reload、24 项 installed agent、39 项 installed Redis 和 22 项发行回归全部通过。代码提交 `569d3865` 与发行测试环境补正 `1c8c82fc` 均已 push；用户随后选择“暂只完成 HEP 候选，保留 SDK 发布待办”。本次完成候选和 draft，不公开发布 HEP/SDK。详见 [2.0.17 发行记录](release-2.0.17.md) 与 [发行验收证据](release-2.0.17-verification.json)。
