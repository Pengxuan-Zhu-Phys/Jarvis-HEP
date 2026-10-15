# Jarvis-HEP 实施路线

2026-10-11 重建的维护基线。旧路线文件在约定路径与 Git 历史均未找到，本文件记录已获用户授权的当前阶段。

## 2.0.15 阶段顺序

| 阶段 | 内容 | 状态 | 结束门槛 |
| --- | --- | --- | --- |
| 维护基线 | 重建五份维护文件，调整发布责任 | 完成 | 文件保存并阅读；明确最终发布由用户执行 |
| 配置与启动 | VAL-01、RUN-01 | 完成 | 无效表达式提前拒绝；绕过校验的初始化故障有界退出；正常运行与恢复通过 |
| monitor | MON-01、MON-02 | 完成 | ref JSON、Worker 信息、采集时间和错误输出通过只读验收 |
| CLI 与依赖 | PLOT-01、CLI-01、OPR-01 | 完成 | 诊断和参数转发通过，Operas 源码修复可打包 |
| 回归 | QA-01、QA-02、SDK-01、QA-03 | 完成 | 定向测试、真实运行、fast gate 与发行包验收通过 |
| 发行候选 | REL-01 准备部分 | 完成，待用户发布 | 2.0.15 wheel、sdist、校验和、验证报告及用户发布步骤齐全 |
| 正式发布 | 用户执行 | 待用户发布 | 用户发布后再核验公开发行状态 |

所有任务的细项、依赖与验收以 [开发台账](development-ledger-2.0.15.md)为准。任务状态只能随实际证据推进；发行候选准备完毕不等同于已发布。

## 范围控制

保留已有 sampler 与 Diver 工作；本轮未直接授权的功能变更不自动并入发行候选。若完整回归发现已有失败，记录当前证据、原因与发行影响，再处理必要的发行阻断问题。

每个完成阶段同步上下文、PRD、架构审阅、代码地图、路线和台账。2.0.15 阶段最终发布权限保留给用户；2026-10-13 的新增发行授权见下方 2.0.16 阶段。

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

## 2026-10-13 CLI agent 集成（独立新增阶段）

| 阶段 | 范围 | 状态/结束门槛 |
| --- | --- | --- |
| 先行审计 | CLI/man/schema/诊断与 PLOT 邮件模式 | 已完成，报告先于代码写入 |
| 显式合同 | 自举入口、问题路由、命令效果、授权/证据边界 | 已实现，同源人类/JSON入口可发现 |
| 编写与反馈辅助 | 本地 schema、准确模板/fragment/前置条件、邮件草稿 | 已实现，无扫描执行或自动发信 |
| 诊断与集成 | manual_argv、范围metadata、JSON参数/I/O/内部失败 | 已实现，真实查询→错误→导航→修正→复验闭环通过 |
| 最终回归 | 定向测试、完整fast gate、canonical 文档同步 | 完成并记录例外：194 定向测试/316 子测试、13 安装包 CLI 检查通过；fast gate 1345/437 通过，1 个已知外部 Examples 编码失败保留 |

此阶段由用户新请求授权，与旧 2.0.15 公开发布状态独立。保留现有 Diver/pyproject 改动，不重建/覆盖旧发行候选，不发布。重复键 loader 行为及依赖缓存边界明确记录，不静默改变兼容行为。

## 第二轮 CLI agent 实际体验 review

阶段已完成：先实测并追加报告，再实现短合同/定向查询、关键词搜索、卡片解释、精确 schema/模板导航、重复键拒绝及可填写草稿。未知路径和冲突 selector 精确拒绝；重复键从第一轮的行为约定改为加载约束，合法 merge 保留，旧外部 seed 重复的测试预期已明确更新。

结束门槛：332 tests/373 subtests 定向通过；22 安装包 CLI 检查通过；完整 fast gate 1416/478 通过、166 slow deselected、1 个已知 sibling iDM 编码失败保留。五份 canonical、报告、开发台账与第二轮 verification 同步。版本和旧公开发布状态不变。


## 2026-10-13 2.0.16 提交与发行

用户已授权 commit 与新版本发布，本阶段取代此前仅准备候选的权限边界。版本 2.0.16 包含两轮 CLI agent 改进及已提交的 2.0.15 可靠性/sampler integration 基线；保持原有 YAML 结构，显式重复键改为 JV2-YAML-003 错误。Diver 草稿、其 package-data 及未验收的 sampler 审阅文件保留且排除。

从干净 commit 副本构建完整 wheel/sdist，执行源码 fast gate、全依赖安装、pip check、真实 Redis/spawn smoke、旧可靠性回归和新增 installed agent CLI smoke。发布正常使用 GitHub Release/Trusted Publishing，不使用本地 token 上传。当前依赖检查：Operas 1.4.2 已公开；Jarvis-HEP-Sampler 0.1.0 在 PyPI 缺失，须先解决才能正式发布 HEP。具体状态与证据见 [2.0.16 发行记录](release-2.0.16.md)。

2.0.16 发行候选验收完成（源码 commit `50631bae` 已推送到 `feat/sampler-plugins-s1`）：带 Git 元数据的干净 checkout fast gate 1396 tests/478 subtests、3 skipped、166 slow deselected；完整包全依赖干净安装、pip check、真实 Redis/spawn/项目重载 smoke、24 项 agent CLI 和 22 项安装包回归均通过。wheel/sdist 与已提交源码逐文件一致并排除 Diver。SDK 使用本地已验收候选，其余依赖来自官方 PyPI，因此仍不能声称仅用公开索引安装成功。SDK 0.1.0 尚未公开，正式 HEP 发布等待依赖；准备 draft Release。详见 `release-2.0.16.md` 与 `release-2.0.16-verification.json`。


## 2026-10-15 Redis 环境引导（独立新增阶段）

| 阶段 | 范围 | 状态/结束门槛 |
| --- | --- | --- |
| 现状核对 | 首次命令 PATH 提示、managed Redis 生命周期、YAML/test 文档误述 | 已核对；区分 CLI 首次使用与 pip 安装钩子 |
| Redis CLI | `redis [info/check/install/start]`、JSON、重新检测与退出码 | 完成；仅输出指导；查询不执行安装/启动，不探测服务/网络/端口 |
| 平台说明 | 具体发行版优先、版本对应 Redis/Valkey、旧 Ubuntu/Debian 官方签名 APT recipe、RHEL 8 redis:6/updated AppStream、Windows WSL、可选前台调试 | 实现与平台矩阵同步；未核验的旧版不输出不合格默认包命令，安装/启动建议必须与本地检测事实区分 |
| 回归与维护同步 | 定向/既有 CLI 回归，README/INSTALL/五份 canonical/台账 | 最终定向 177 项通过（含真实 Redis 生命周期）；临时 wheel 安装 9 项 CLI 检查通过。完整 fast gate 为 1464 passed / 478 subtests、1 项已知 sibling iDM 非 ASCII 描述失败；最终增补平台用例由定向复验覆盖。详见 [验收记录](redis-environment-verification.json)；不得用旧发行验收代替 |

本阶段来自用户 Redis 环境体验请求，与 2.0.16 发行准备独立。不调整版本或既有发行候选，不公开发布；仅在临时目录构建安装测试 wheel，排除已有 Diver/pyproject 和 sampler 审阅改动。运行期自动启动/清理 Redis 的边界不变。


## 2026-10-15 Redis CLI 帮助布局补正（后续阶段）

| 阶段 | 范围 | 状态/结束门槛 |
| --- | --- | --- |
| 接口理解更正 | 第一阶段 positional topic 未建立真正子命令，child help 仍显示父页 | 已确认并显式记录，第一阶段验收保留为历史 |
| 标准命令组 | info/check/install/start subparsers、父级 Commands、独立 child help | 完成；裸 info 默认及 JSON 子命令前后兼容 |
| agent 发现同步 | exact child catalog / search / 定向查询 | 完成；必须发现可直接查询的子命令 |
| 回归与文档 | parser/help/JSON/agent 定向、安装包 help smoke、README/INSTALL 与维护层 | 202 项定向通过，临时 wheel 安装 25 项检查通过；完整 fast gate 为 1509 passed / 478 subtests、1 项既有外部 iDM 描述编码失败。详见 [本次验收](redis-cli-help-verification.json)；不得沿用第一阶段计数充当通过 |

[前一阶段验收](redis-environment-verification.json) 的 177 项定向、9 项安装 smoke 与 fast gate 作为历史证据保留。本次不变更平台安装逻辑、服务/扫描执行边界、版本、既有发行候选或公开发布状态。


## 2026-10-15 HEP 2.0.17 Redis-only 发行（新增授权阶段）

| 阶段 | 范围 | 状态/结束门槛 |
| --- | --- | --- |
| 范围与版本 | 只提交 Redis 两阶段改进、对应测试文档、HEP 2.0.17 声明 | 已获授权；排除本轮 Sampler/runtime/Diver，实际 commit/push 结果待记录 |
| 干净源码候选 | 从 committed source 构建完整 wheel/sdist 并核对文件 | 待本阶段测试结果补录；不沿用 dirty checkout 或开发阶段计数 |
| 安装包验收 | 完整依赖、pip check、installed Redis gate、真实 Redis/spawn、发行回归 | 可用保持原样的本地已验收 SDK 0.1.0 wheel，仅证明候选；结果待补录 |
| draft/公开发布 | 正常 GitHub Release/Trusted Publishing | 2026-10-15 Sampler PyPI 404 阻断正式 HEP 发布；可先准备 draft，不修改或发布 SDK |

本次用户授权取代此前 Redis 阶段“不 bump/不发布”范围，但没有扩大到当前 Sampler 工作。以前的 Redis 202/25 与 fast gate 为历史；完整状态和证据随实际操作更新于 [2.0.17 发行记录](release-2.0.17.md)。不将准备完成误记为公开发布成功。
