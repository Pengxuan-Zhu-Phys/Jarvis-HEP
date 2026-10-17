# Jarvis-HEP 产品要求

2026-10-11 重建的维护基线，来源为当前 README、安装文档、运行接口与 [2.0.15 台账](development-ledger-2.0.15.md)，不声称还原旧 PRD。

> 当前源码状态（2026-10-17 提交复查后）：本轮复查修复使用 SDK 源码的 `HOST_CONTRACT_VERSION >= 2`，包含入口保留名策略和 `CatalogExtension` 数据接口；已发布的 SDK 1.0.3 尚无该契约。两仓版本号和 `>=1.0.3` 的历史最低声明保持不变，不加上限。HEP 运行时额外核查契约，旧 SDK 给出匹配源码构建的安装诊断；维护者后续发布 SDK 后再将 HEP 版本下限对齐真实发行号。本轮不发布。

## 当前范围

当前 master 要求 `Jarvis-HEP-Sampler>=1.0.3`，不设上限。`Jarvis sampler` 复用 `jsampler` 命令面，并增加 HEP 的内置方法和边界检查。插件由 SDK 定义，HEP 负责运行、检查点和项目包。

## 2.0.17 公开发行范围

2.0.17 从公开 master 0bbcf675 基线移植不依赖外部 SDK 的可靠性、CLI agent/YAML 和 Redis 功能，使用内置采样器；该次发行不声明外部 Sampler SDK 依赖、插件加载或注册命令。用户已授权正式公开发布，验收与发布结果见下方已完成的源码/安装验收与发布状态。后续历史阶段中的 SDK 要求和候选限定不适用于该次发行。

## 核心行为

- `validate` 只验证卡片，不启动 Redis 或 Worker；配置错误退出 2，并给出 YAML 字段路径与纠正建议。
- `run` 返回真实运行结果：成功为 0，运行错误或样本失败为 1，使用或配置错误为 2，中断为 130。
- Worker 启动失败必须反馈给控制进程，结束不可运行的扫描并清理本次资源；正常扫描中的可恢复进程故障保留恢复能力。
- 监控读取选中扫描的数据。缺失或未知信息不伪装成零；JSON 模式输出可解析的 JSON。
- plot 使用绘图卡片。收到扫描卡片时说明用途并指向 `gen-plot-yaml`。
- help 和冲突参数的错误提示必须可读、稳定；依赖 CLI 的自由参数保持转发。

## 当前验收边界

保持现有内置采样算法、YAML 结构与正常 Worker 恢复。CLI agent 的合同/schema/模板/inspect 是查询和编写辅助，不能把 validate 成功当作用户算子、外部计算器或物理结果正确。重复显式 YAML 键给出路径与两处位置；合法 merge overrides 保留。邮件功能仅生成草稿。

Redis 父级列出 info/check/install/start，各子级有独立 help/JSON；裸命令默认 info，JSON 支持子命令前后。info/check 仅做本地 PATH/有界版本检测，check 对缺失/不可用/Redis<6.2/原生 Windows 退出1。install/start只打印按系统选择的安装/启动指导，不连接服务、不探测端口、不启动扫描；普通扫描自动管理本次 Redis。



## 2026-10-15 2.0.17 SDK-free 公开发行要求（历史发布要求）

- 从公开 master `0bbcf675` 保留内置采样，移植不依赖 SDK 的可靠性、CLI agent/YAML 与 Redis 改进；本次没有外部 SDK/plugin 能力或依赖。
- 用户授权公开 HEP，推迟已提交 SDK 集成 `0d859502` / `52876060` 及未提交 Sampler/runtime/Diver 工作；不为发布 HEP 修改或发布 SDK。
- 公开候选应在干净环境仅从公共 PyPI 安装全部依赖，并通过 pip check、内置采样真实 Redis/spawn、installed agent/Redis/help/JSON、相关可靠性和完整 fast gate。
- `1c8c82fc` 旧候选的本地 SDK 依赖和候选限定已被替代；旧验收只作历史。当前公开成功须核验 GitHub 工作流、PyPI 版本/产物和全新公共索引安装。

当前验收与发布结果：见下方已完成的源码/安装验收与发布状态。以 [2.0.17 发行记录](release-2.0.17.md) 为当前事实源。

## 本次源码、CI 与公开安装状态

源码功能提交 `3b9d41d0`、验证接口兼容修正 `e8af5221` 已进入 master；公开 tag `v2.0.17` 指向 `8e0add06`。macOS 与 Linux CI 均为 1426 passed、3 skipped（缺少 sibling Jarvis-Examples）、164 deselected、477 subtests；分别用时 205.99s / 265.40s。Linux CI 的完整安装/pip check、真实 Redis/spawn/双 Worker/monitor/项目重跑、agent 24、Redis 39、发行回归 22 均通过。

[GitHub Release](https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP/releases/tag/v2.0.17) 已公开。OIDC publisher 因 PyPI 未匹配 Trusted Publisher 配置而报 `invalid-publisher`；本次使用既有本地 PyPI 凭证上传成功 CI 的原样产物，没有重建。PyPI wheel 已公开，全新环境仅从公共索引安装完整依赖，pip check、无 SDK、276 个运行文件身份比对、Redis 39 / agent 24 均通过。PyPI 源码包上传未获服务器响应；最后一次有界请求在 36.913 秒后 ReadTimeout，全部上传进程已停止。PyPI wheel 与公开安装已完成，已验收源码包可从 GitHub Release 下载。

发布 wheel SHA256 为 `25257b646c1a7ce828dae136e1e111de7c70d749d7b5b6857b7c687f227fb3ae`；CI sdist 的预期 SHA256 为 `1843e006921161dd6a46507b7b7c12b52b48990b9e3af1c76541788dafbe9915`，未在 PyPI 公开；GitHub Release 另提供通过本地完整验收的源码包。发布前本地产物的 a823… / 9171… hash 仅作本地验收记录，不能当成公开 CI 产物 hash。完整状态见 [发行记录](release-2.0.17.md) 与 [验收证据](release-2.0.17-verification.json)。

## 2026-10-16 HEP 修复验收约束（历史阶段；CLI 边界已由 N1–N6 更新）

- 正式 fresh 不得因旧 DATABASE 的 index/UUID 而跳过重新求值，包括 checkpoint 拒绝后的 fresh 回退；resume 保留对连续 prefix、稀疏终态行和 Failed 行的去重。check 保持独立 DATABASE/test 的追加烟测，不开启 resume 批次过滤。fresh 中断后的续跑规则见下节。
- 当时要求 `Jarvis sampler load` 拒绝内置方法及内置插件名，大小写不敏感；旧冲突在运行前报错。现由 SDK 默认保留名单与统一注册表报告实现，见下文。
- 当时由 HEP 校验 `check --dimension`、损坏注册表以及按名字 check 的快照 SHA-256。现 check 只接受文件，维度和 CLI 错误均由 SDK 处理；按名字 check 已删除，见下文。
- `check --card` 保留卡片预算、seed、Variables 与选项；显式 `--points` 改变总候选预算，文本模式说明这一点，JSON 模式保持可解析。
- 不改 SDK 项目和发布脚本，不提升版本或发布。原有 Diver 移植暂留，见 [架构审阅](ARCHITECTURE_REVIEW_JARVIS_HEP.md)"Diver 移植"。本轮门禁见 [验证证据](hep-review-2026-10-16-verification.json)。

## 2026-10-16 续跑协议验收要求（已实现）

设计见 [架构审阅](ARCHITECTURE_REVIEW_JARVIS_HEP.md)"续跑以 checkpoint 为准"。对 Grid/Random/Bridson/CSV，没有 checkpoint 就没有 resume；进度只以 checkpoint 为准。Dynesty 独立引擎文件与插件续跑边界保持原状。

- **fresh 中断后续跑。** Grid/Random/Bridson/CSV 在旧 DATABASE 上重新开始、中途中断、再 `--resume`，最终库中只有当前卡片的结果，旧行一条不剩（真实 Redis/spawn/HDF5）。第一次心跳前中断也一样。
- **清尾。** 新扫描点数少于旧扫描时，结束后不留 sample_index ≥ 新总数的旧行或重复 index；保留原分片边界与名称，不合成大文件。无变化分片及 CSV 不重写。仅确认 Archiver 正常退出后清尾；超时/退出失败跳过清尾。磁盘或停止错误仍尝试 summary、锁释放、Worker 与 Redis 关闭，最后返回失败。部分清尾失败不丢最后有效副本，可重试，重试应修复已替换 HDF5 对应的旧 CSV；清空 live 时不留旧 CSV；覆盖扫描出图使用清尾后的库。
- **没有 checkpoint。** `checkpoint.enabled: false` 加 `--resume`、或删除 state.pkl 后 `--resume`，都给出警告并从头扫描，不读取 DATABASE 推进进度。
- **check 模式。** 不写 state.pkl 和 points.npy；DATABASE/test 追加行为不变。
- **核对卡片。** 改动 Bounds（seed、point_number、radius、max_attempt）或 selection 后续跑，拒绝 checkpoint 并从头扫描，日志说明原因；只改 workers 时正常续跑。
- **核对点集。** Grid 点集或 CSV 文件内容与 checkpoint 指纹不符时拒绝续跑。
- **Bridson。**
  - fresh 写一次 points.npy，扫描期间不重写；
  - resume 读取 points.npy，不调用点集生成；
  - points.npy 缺失、损坏或哈希不符时从头扫描；
  - 用大点集（≥10⁶ 点）验证 resume 不重新生成，并记录读取与校验耗时。
- **checkpoint 大小。** 同一 sampler 在 10 点与 10⁵ 点下 state.pkl 大小基本相同。indexed sampler 不再累积提交 UUID。
- **插件。** 不变：在已有数据的 scan 上重新开始仍被 preflight 拒绝。

续跑协议验收结果：4 种内置方法各在首次心跳前、后中断，再经真实 Redis/spawn/HDF5 续跑，8 项通过；实际长 Grid 扫描由短 Grid 扫描覆盖并清尾通过；Bridson check 不生成 checkpoint/points.npy 通过；无 state.pkl 与关闭 checkpoint 的真实续跑各通过 1 项。无 checkpoint、关闭 checkpoint、Sampling/CSV/Grid 指纹、Bridson 点集缺失/损坏/哈希不符及 10/10⁵ 点 state.pkl 大小均有定向回归。完整数值见 [验证证据](hep-review-2026-10-16-verification.json)。

## 2026-10-17 注册表与命名验收要求（N1–N6 已完成）

设计见 [架构审阅](ARCHITECTURE_REVIEW_JARVIS_HEP.md)"采样器注册表与命名"。外部 sampler 不能与内置方法重名，这是硬约束。

- **保留名单。** SDK `RESERVED_METHOD_NAMES` 包含 HEP 全部内置方法；HEP 测试断言 builtin 与 builtin-plugin 都在名单内。
- **官方 builtin-plugin。** 元数据解析只查格式，官方 `LatinHypercube` 直接经通用 loader 加载，无低层保留名豁免。外部 init/load/check 入口执行统一保留名单。已绑定扫描从原扫描副本和保存的 SHA-256 续跑，即使当前安装包路径或源码字节变化。项目内绑定及扫描副本可由用户编辑，不能将这项校验当作防止有意伪造的安全边界。
- **load。** `jsampler load` 与 `Jarvis sampler load` 对 `Grid`、`grid`、`LATINHYPERCUBE` 都拒绝，退出码和信息相同，不写登记和快照。
- **check。**
  - 只接受文件；传方法名时 exit 2，提示 `Jarvis sampler info <Method>`。
  - 方法名是报告第一项；重名文件在 options schema/prepare/draw 之前失败，exit 2，后续项标为未执行；算法检查失败 exit 1。
  - `jsampler check` 与 `Jarvis sampler check` 结果相同。
  - `--dimension` 为 0 或负数时 exit 2，不再回落为 2。
  - `--seed` 为负数时在载入插件前 exit 2，不进入检查报告。
- **info。** 查询内置 `Grid` / `GRID` 同时显示旧 `grid` 等大小写无关的冲突及精确修复命令；文本与 JSON 信息一致。登记损坏时 `jsampler info` 和 `Jarvis sampler info`（包括查询内置 Grid）给出可读错误，JSON 模式输出单个 JSON 错误，没有 traceback。
- **run/check 启动。**
  - 用旧版 SDK 或手写 `sources.json` 造出 `Grid` 登记后，`Jarvis run` 和 `Jarvis check` 对任何卡片（包括用 Random 的卡片）都在启动 Redis 之前退出，列出重名条目与修复命令（登记中的确切名字）；按提示 `Jarvis sampler remove` 只删除该外部登记及其托管快照，不动内置方法或用户源文件，之后恢复正常。登记中的 `grid` 用 `remove grid` 可删除。
  - 卡片所用插件的快照被改：启动前退出，不导入快照。
  - 无关条目快照缺失：警告并继续。
  - 登记无法读取：用内置方法的卡片警告并继续；用插件的卡片退出。
- **展示。** `Jarvis sampler list` 列出重名条目并标明状态，不再隐藏。
- **CLI 风格。** `jsampler` 与 `Jarvis sampler` 都用 Jarvis CLI 风格，load/check 等无扩展命令在替换程序名后帮助与错误一致，退出码一致；list/info 经 SDK `CatalogExtension` 只扩充数据，文本和 JSON 均由 SDK 渲染，不在 HEP 再复制展示。登记损坏时 list/info/remove 都 exit 2；SDK 不可用时 HEP 仍输出统一错误面板或单个 JSON 错误对象。全生态规则见 [CLI 规范](JARVIS_CLI_DESIGN.md)。
- **不变。** 没有重名且登记正常时，run/check 的行为与输出不变；生成报告不导入插件代码。

08b2dab5 / d9d2853 阶段的历史验收：SDK N1/N2 源码测试 234 passed、Ruff 通过；HEP 最终定向集成 169 passed、4 deselected，快速全套门禁 1597 passed、196 deselected，续跑真实 Redis 慢测 20 passed、54 deselected。HEP N5 非 CLI 专测 5 passed，旧 `Grid`/`grid` 登记的真实 `Jarvis run`/`Jarvis check`、list/info/remove 入口回归 4 passed，均覆盖精确修复命令与移除后内置方法恢复。HEP 官方 `LatinHypercube` 扫描副本在包路径和字节变化后仍可读取绑定并按旧 hash 续跑。N4 的转发与 CLI 一致性已在本地实现。快速全套门禁先于最后三项 CLI 边界修复完成，定向集成在修复后复跑。N3 的 SDK 1.0.3 已公开并包含 N1/N2；HEP 最低依赖及运行时门禁为 `>=1.0.3`，不设上限。全新 Python 3.11 环境从本地 HEP 源码完成全依赖实际安装（66 个包），依赖从公共 PyPI 解析并选中 SDK 1.0.3；`pip check`、已安装模块与两个 sampler 入口的 JSON 版本输出通过。HEP 包来自本地源码，未公开发布新版本。完整结果见 [验证证据](hep-review-2026-10-16-verification.json)。

2026-10-17 C1 持久化绑定收拢已实现，行为保持不变；指定三组续跑测试的非 slow 部分 58 passed，Core collaborator 17 passed；真实 Redis 慢测 20 passed、54 deselected，用时 278.16 秒。

## 2026-10-17 提交复查验收

必须覆盖：Archiver 超时/关闭异常不清尾，磁盘异常不跳过其他收尾；多分片边界和无变化文件身份保持；每分片原子替换及中途失败重试；真实 Redis/spawn/HDF5 的长短扫描强制轮转；官方插件加载不需豁免；外部保留名优先于 schema/prepare/draw，check exit 2；算法检查失败仍 exit 1；内置 info 显示大小写不同的旧冲突；两个入口调用同一 SDK 渲染器；旧 SDK 契约给可读诊断。数值见现有 [验证记录](hep-review-2026-10-16-verification.json)。本轮不以未发布源码替代公共索引发行验收。

提交复查验证完成：SDK 全套 253 passed、Ruff 通过；HEP 最终快速门禁 1627 passed / 197 deselected / 478 subtests，真实 Redis/spawn/HDF5 慢测 21 passed，最后边界补丁后长短扫描两种分片模式复跑 2 passed，收尾定向 40 passed。匹配暂存源码的两个 wheel 在干净 Python 3.11 环境完整安装 66 个包，pip check、官方插件、分片清尾与双入口 CLI 通过；公开 SDK 1.0.3 的旧契约诊断也经真实换装核验。wheel 内 HEP 227 / SDK 16 个 Python 文件逐字匹配暂存内容。详细命令、耗时、构件哈希和范围保护见现有 [验证 JSON](hep-review-2026-10-16-verification.json) 的 `review_followup_2026_10_17`；本轮未改发行版本、未发布。
