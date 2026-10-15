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

## 2026-10-13 agent 辅助层

`cli/agent.py` 提供版本化自举合同和命令效果；命令参数来自当前 help parser，转发 CLI/项目兼容参数明确标记不完整。`cli/man.py` 的 agent-readme 返回同一合同，原有手册和 schema 继续作为字段事实源。插件 OPTIONS_SCHEMA 原样公开，metadata 查询不导入插件实现。

`cli/agent_artifacts.py` 只输出本地 schema bundle 和模板/前置条件；$id URI 为本地 registry 标识，不触发网络。模板保持原文，fragment 必须组合为完整卡。`cli/agent_request.py` 只生成草稿和反馈协议，不实现 SMTP 或其它发送通道。

`cli/run.py` 在原 validate 结构上增加导航和证据范围；loader/validator 内部错误与配置错误区分，JSON stdout 不混入诊断。validate 无扫描执行，但依赖导入缓存/CLI marker 并非绝对无写入。重复键兼容行为本阶段未改。边界和验证记录见 [报告](cli-agent-integration-report.md)。

最终验收确认 authoring smoke 不实例化 Core/Redis；安装包从独立 venv 的 site-packages 导入，13 项 CLI 检查通过。最终定向 194/316 通过；运行回归通过，fast gate 唯一失败来自已知 sibling Examples 描述编码。完整例子引用按实际正文筛选，fragment 不再被 complete-card 导航误标。

## 第二轮架构补充

`cli/agent_search.py` 从安装手册/nav 和 help metadata 建立有界词汇索引，返回查询 argv，不将命中或无命中作为能力支持证明。`cli/agent_inspect.py` 复用 loader/validator/runtime policy 读模型，输出生效值、来源和计划路径；CLI section 投影和 preview 限制展示，不改变验证数据，不构造 Core/Redis 或导入用户算子。

`cardload/yaml_loader.py` 在 SafeLoader flatten merges 前检查原始 mapping 显式键，捕获 silent overwrite，并保留合法继承覆盖、alias 和 quoted literal merge-like keys。task_config 将重复诊断转为 TaskCardLoadError 的 path/location/related_locations；run/inspect 输出实际位置。

手册 exact lookup 不再丢弃后缀或覆盖冲突 context selector；schema fragment 查询只解析 manifest 的本地文档，不访问 URI。模板识别已知 import-path 边界，草稿输入保持数据。原完整 agent JSON 不改为短格式，新增投影避免破坏消费者。最终 332/373 定向、22 安装 CLI 检查通过；fast gate 1416/478 通过及已知外部编码失败。


## 2026-10-13 2.0.16 提交与发行

用户已授权 commit 与新版本发布，本阶段取代此前仅准备候选的权限边界。版本 2.0.16 包含两轮 CLI agent 改进及已提交的 2.0.15 可靠性/sampler integration 基线；保持原有 YAML 结构，显式重复键改为 JV2-YAML-003 错误。Diver 草稿、其 package-data 及未验收的 sampler 审阅文件保留且排除。

从干净 commit 副本构建完整 wheel/sdist，执行源码 fast gate、全依赖安装、pip check、真实 Redis/spawn smoke、旧可靠性回归和新增 installed agent CLI smoke。发布正常使用 GitHub Release/Trusted Publishing，不使用本地 token 上传。当前依赖检查：Operas 1.4.2 已公开；Jarvis-HEP-Sampler 0.1.0 在 PyPI 缺失，须先解决才能正式发布 HEP。具体状态与证据见 [2.0.16 发行记录](release-2.0.16.md)。

2.0.16 发行候选验收完成（源码 commit `50631bae` 已推送到 `feat/sampler-plugins-s1`）：带 Git 元数据的干净 checkout fast gate 1396 tests/478 subtests、3 skipped、166 slow deselected；完整包全依赖干净安装、pip check、真实 Redis/spawn/项目重载 smoke、24 项 agent CLI 和 22 项安装包回归均通过。wheel/sdist 与已提交源码逐文件一致并排除 Diver。SDK 使用本地已验收候选，其余依赖来自官方 PyPI，因此仍不能声称仅用公开索引安装成功。SDK 0.1.0 尚未公开，正式 HEP 发布等待依赖；准备 draft Release。详见 `release-2.0.16.md` 与 `release-2.0.16-verification.json`。


## 2026-10-15 Redis 环境查询与引导

`queue/redis_notice.py` 是平台识别、安装建议和环境报告的共同事实源；顶层 `redis_notice.py` 保持兼容导出。首次 CLI 退出提示沿用 marker 和 stderr，只查 PATH，不执行 server。独立的 `cli/redis.py` 渲染文本或一个 JSON 报告，`cli/client.py` / `cli/help.py` 提供 topic 路由，agent 合同公开查询范围。

`info` / `check` 以有界 `--version` subprocess 查询本地 executable；`check` 将不存在、不可用、版本无法确认/低于要求和原生 Windows 判为退出 1。`install` / `start` 不执行该 subprocess，只形成安装或启动指引。所有查询与 runtime/Redis 网络连接、端口探测、服务管理和扫描创建分离；显式查询每次重新检查，不使用首次提示 marker 作为环境事实。

平台识别先用具体发行版 ID，再用 ID_LIKE，再退回可用包管理器；按已知发行版版本选 Redis/Valkey，避免衍生系统被其它 family 提前吞掉。原生 Windows 给 WSL 安装说明，不能虚构原生 Windows Redis 包。正常扫描继续复用现有 managed Redis 生命周期；可选手动调试只打印 loopback/6380/无持久化的前台命令。

`_legacy_apt_release` 只为明确的 Ubuntu 20.04/22.04 与 Debian 11 选择已核验的官方签名 APT recipe，不拿衍生系统 VERSION_ID 推断上游仓库；更老未核验版本不输出不合格的默认包命令。RHEL 8 选择 `redis:6` module 并说明 AppStream 更新/已有旧 stream 切换，版本仍由显式 check 复验。

验收状态：最终定向 177 项通过（含真实 Redis 生命周期）；临时 wheel 安装 9 项 CLI 检查通过。完整 fast gate 为 1464 passed / 478 subtests、1 项已知 sibling iDM 非 ASCII 描述失败；最终增补平台用例由定向复验覆盖。详见 [验收记录](redis-environment-verification.json)。未更改队列协议、运行资源选择或公开发行状态。


## 2026-10-15 Redis 正式子命令架构补正

更正第一阶段 `topic` 路由的理解：当时 help parser 使用 positional choices，并非 argparse 的正式命令组，因此不同 `redis COMMAND --help` 无独立子级帮助。现由 `cli/help.py:build_parser` 为 Redis 创建 info/check/install/start subparsers；parent 和 child 共同接受 `--json`，子级缺省不得覆盖父级已设的 JSON 状态，裸 Redis 继续落到 info。

`cli/client.py` 从正式子命令选择分派给既有报告/render 层；`cli/agent.py` 命令 catalog 与 `cli/agent_search.py` 搜索发现公开 exact child commands，并复用子 parser metadata，使帮助、参数和 agent 查询一致。帮助解析不调用 server；既有平台安装建议、版本检查和 managed Redis 生命周期不受此结构调整影响。

第一阶段 [验收记录](redis-environment-verification.json) 中 177 项定向、9 项安装 smoke 及 fast gate 是历史证据。本次正式命令组验收：202 项定向通过，临时 wheel 安装 25 项检查通过；完整 fast gate 为 1509 passed / 478 subtests、1 项既有外部 iDM 描述编码失败。详见 [本次验收](redis-cli-help-verification.json)；不以历史结果代替。


## 2026-10-15 2.0.17 发行源码隔离

Redis-only 授权包含环境/CLI 层及对应测试文档、HEP 版本声明，不包括当前工作区 Sampler/runtime/Diver 修补。发行由明确提交的干净源码构建，`pyproject.toml` 只纳入版本变更；wheel/sdist 逐项核对范围，不能由 dirty checkout 隐式收集其它源码或 package-data。

现有 committed baseline 仍需要 SDK >=0.1.0,<0.2；使用不变的本地已验收 SDK 0.1.0 wheel 完整安装候选是测试安排，不改变依赖声明也不证明公共索引解析可用。2026-10-15 Sampler PyPI 404 阻断正式 HEP 发布，本次不修改/代发布 SDK；依赖解决前保留候选/draft 状态。

本阶段干净源码/安装包 Redis gate、pip check、真实 Redis/spawn 与发行回归结果：待本阶段测试结果补录。此前 Redis 202 项定向/25 项安装检查为历史，不替代当前源码验收。状态详见 [2.0.17 发行记录](release-2.0.17.md)。
