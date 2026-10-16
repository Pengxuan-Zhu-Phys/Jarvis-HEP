# Jarvis-HEP 架构审阅

2026-10-11 根据当前实现重建，未恢复历史审阅原文。

## 当前发行架构

2.0.17 的公开发行不带外部 Sampler SDK，记录见 [2.0.17 发行记录](release-2.0.17.md)。当前 master 在该发行之后接入 Jarvis-HEP-Sampler：依赖为 `Jarvis-HEP-Sampler>=1.0.3`，不设上限；`Jarvis sampler` 复用 `jsampler` 命令面，HEP 适配层另负责内置方法、保留名和调用前校验。插件格式在 SDK，运行时、检查点和项目包在 HEP。

## 2.0.17 发行架构

2.0.17 使用 master 0bbcf675 的内置采样器与依赖边界，移植不依赖 SDK 的可靠性、CLI agent/YAML 和 Redis 功能。该次发行没有外部 Sampler SDK 依赖、插件注册表或插件 CLI。0d859502 / 52876060 的 SDK 集成与其它工作区改动当时留在原分支；下方保留的旧候选说明仅作历史记录。

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

## CLI 辅助层

`cli/agent.py`、`agent_search.py`、`agent_artifacts.py`、`agent_inspect.py`、`agent_request.py` 提供自举合同、发现、离线资源、生效卡片解释与反馈草稿。查询不构造 Core/Redis 或导入用户算子；bundled schema只描述内置结构词汇。`cardload/yaml_loader.py` 拒绝显式重复键并保留合法 merge。`cli/run.py` 输出验证导航和范围。

`queue/redis_notice.py` 是平台/安装建议/报告的事实源；首次CLI结束后的提示仅查PATH并写stderr，显式info/check有界执行--version。`cli/redis.py` 渲染报告，`cli/help.py` 的正式subparsers提供独立help/JSON；agent以相同parser metadata发现四个exact子命令。install/start不调用server；查询不连接服务或探测端口。



## 2026-10-15 SDK-free 发行隔离架构（历史发布阶段）

本次源自公开 master `0bbcf675`，只移植不依赖外部 SDK 的 Worker 启动/表达式诊断/monitor/CLI 可靠性、CLI agent/YAML 与 Redis 引导。内置 sampler 的 factory、checkpoint/恢复和调用边界沿用公开基线，不带入 SDK 集成 `0d859502` / `52876060` 或后续未提交 Sampler/runtime/Diver 编辑。

HEP 的 wheel/sdist 不应包含 SDK/plugin registry/CLI，不应有 `Jarvis-HEP-Sampler` Requires-Dist。构建后的完整依赖安装使用公共 PyPI，不接受本地 SDK 或 `--no-deps` 作为发行验收。实际 fast gate、installed CLI、pip check 和真实 Redis/spawn 重新执行，不能引用已被替代的 `1c8c82fc` SDK 候选作为当前证明。

当前验证与公开发布状态：见下方已完成的源码/安装验收与发布状态；证据见 [2.0.17 发行记录](release-2.0.17.md)。原工作区/分支的 Sampler 开发继续独立保留。

## 本次源码、CI 与公开安装状态

源码功能提交 `3b9d41d0`、验证接口兼容修正 `e8af5221` 已进入 master；公开 tag `v2.0.17` 指向 `8e0add06`。macOS 与 Linux CI 均为 1426 passed、3 skipped（缺少 sibling Jarvis-Examples）、164 deselected、477 subtests；分别用时 205.99s / 265.40s。Linux CI 的完整安装/pip check、真实 Redis/spawn/双 Worker/monitor/项目重跑、agent 24、Redis 39、发行回归 22 均通过。

[GitHub Release](https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP/releases/tag/v2.0.17) 已公开。OIDC publisher 因 PyPI 未匹配 Trusted Publisher 配置而报 `invalid-publisher`；本次使用既有本地 PyPI 凭证上传成功 CI 的原样产物，没有重建。PyPI wheel 已公开，全新环境仅从公共索引安装完整依赖，pip check、无 SDK、276 个运行文件身份比对、Redis 39 / agent 24 均通过。PyPI 源码包上传未获服务器响应；最后一次有界请求在 36.913 秒后 ReadTimeout，全部上传进程已停止。PyPI wheel 与公开安装已完成，已验收源码包可从 GitHub Release 下载。

发布 wheel SHA256 为 `25257b646c1a7ce828dae136e1e111de7c70d749d7b5b6857b7c687f227fb3ae`；CI sdist 的预期 SHA256 为 `1843e006921161dd6a46507b7b7c12b52b48990b9e3af1c76541788dafbe9915`，未在 PyPI 公开；GitHub Release 另提供通过本地完整验收的源码包。发布前本地产物的 a823… / 9171… hash 仅作本地验收记录，不能当成公开 CI 产物 hash。完整状态见 [发行记录](release-2.0.17.md) 与 [验收证据](release-2.0.17-verification.json)。

## 2026-10-16 fresh/resume 与 CLI 边界修正（历史阶段）

旧理解“DATABASE 中已持久化的 index 就可以在 FixedSetSampler.flush_batch 中跳过”不成立：bootstrap 在 fresh 运行也会读取旧库，save_runtime_checkpoint 还会不断更新持久化水位。改为显式 `set_resume_mode` 绑定最终策略，不从 checkpoint 恢复此标志；只有 resume 过滤 prefix/pending。CSV 的 legacy UUID 过滤使用同一策略。supervisor 在 sampler 初始化后绑定模式，resume import/拒绝回退再次同步，避免策略切换后继续跳过。

当时保留名由 HEP catalog 提供，通过 SDK `SamplerRegistry.load(reserved_names=...)` 公共 API 校验；旧冲突由 Core preflight、Distributor 和 info 分散拦截，HEP 还在 SDK check 前校验 dimension 与按名字 check 的快照哈希。这些实现现已删除：SDK 自带唯一保留名单与 file-only check，HEP 用注册表报告和 SDK CLI 转发，见下文。

发布脚本按用户指示排除；SDK 库内 CLI 行为保持原样。该轮修复只保证 fresh 完成时旧行被重新计算。fresh 中断后再 `--resume` 仍会把未覆盖的旧行算作已完成；当时把这一点写成"不引入 DATABASE 世代或跨卡片中断恢复协议"的可接受限制，这个判断不成立，改由下节协议解决。验证证据见 [本轮记录](hep-review-2026-10-16-verification.json)。

## 2026-10-16 续跑以 checkpoint 为准（R1–R6 已实现）

### 缺陷

同名 scan 在已有 DATABASE 上重新开始（fresh）时，旧行要等新批次到达才被覆盖。中断后有三条路径把旧行当成本次进度：

- `--resume` 导入 checkpoint 游标后调用 `advance_to_persisted_prefix`，跳到 `io/database.py` `read_persisted_sample_index_state` 按全库 sample_index 算出的连续前缀，前缀包含旧行。
- 修复前 fresh 期间 Archiver 启动时以旧库前缀作为持久化水位（`io/archiver.py`），`CheckpointedSampler` "只保存已落盘快照"的判断因此总是成立，checkpoint 可能超前于真正落盘的进度。
- fresh 开始时删除旧 state.pkl，第一次心跳前中断时没有 checkpoint；`--resume` 走 DATABASE-only 恢复（`runtime/_resume_service.py` `_preload_resume_checkpoint`），同样跳到旧前缀。

受影响：Grid、Random、Bridson、CSV。插件不受影响：`plugin_sampler.preflight_plugin` 拒绝在已有数据、state.pkl 或身份文件的 scan 上重新开始，该行为保留。

### 原则

本节的协议实现范围为 Grid、Random、Bridson、CSV；插件路径及 Dynesty 独立引擎续跑保持原状。

- 进度只以 checkpoint 为准；DATABASE 只证明本次运行写入的行。没有 checkpoint 就没有 resume。
- checkpoint 只保存游标与固定大小的少量状态，不随扫描规模增长。大数组不进 state.pkl（`assert_checkpoint_attribute_contract` 已强制每个属性声明存或不存）。

### 规则

1. **没有 checkpoint 就从头扫描。** `checkpoint.enabled: false` 或 state.pkl 不存在时，`--resume` 给出警告后从头扫描。删除无 checkpoint 时的 DATABASE-only 恢复分支；`advance_to_persisted_prefix` 只保留给无 `replacing_existing` 标记且有有效 checkpoint 的扫描。
2. **check 模式不走本协议。** `_scan_driver._apply_check_modules_runtime_policy` 强制关闭 checkpoint，与已强制的单 Worker、不打包一起；DATABASE/test 追加烟测流程不变。
3. **fresh 开始即写 checkpoint。** 正式 fresh 在启动 Archiver 后立刻写游标为 0 的 state.pkl；库中已有行时写入 `replacing_existing`。该标记随后续每次 checkpoint 保留，直到本次扫描正常结束。
4. **fresh 不用旧库作水位。** Archiver 的初始持久化前缀从 0 开始；supervisor 只在 resume 时载入 DATABASE 事实（check 维持自身对 DATABASE/test 的读取）。
5. **resume 的进度来源。**
   - 带 `replacing_existing`：只信 checkpoint 游标以下的行。Archiver 初始前缀取游标；不调用 `advance_to_persisted_prefix`，不按 DATABASE 跳过；继续 `replace_duplicates`。代价是重算最近一次 checkpoint 之后的点（默认心跳 30 秒）。
   - 不带标记：DATABASE 只含本扫描自身的行，保留按 DATABASE 前缀推进，只重放 checkpoint 之后的空档。
6. **fresh 结束时清尾。** 带标记的扫描正常结束、Archiver 关闭后，按本次提交总数跨 live/封存 HDF5 shards 清理旧尾部与重复 index。原设计指定 writer `drop_matching_records`，但该 writer 只能处理 live shard；实现改为临时 SQLite 按 index 保留最后落盘行，再原子替换 live shard 和刷新 CSV，以覆盖封存 shard 的旧行，内存不随行数增长。
7. **核对卡片。** checkpoint 记录 Sampling 部分（Bounds 与 selection）的哈希；续跑时不一致，按 Mapper 指纹不符处理：拒绝 checkpoint，从头扫描。workers 等运行参数不参与比较。现有写入但无人读取的 `config_hash` 已由 `sampling_hash` 取代。现在 checkpoint 中的 seed/selection 会在 `FixedSetSampler._import_common_fields` 静默覆盖卡片值，核对后这种覆盖不再可能。
8. **核对点集。** checkpoint 记录点集指纹：Grid 为 `_P` 的 SHA-256（含 shape、dtype），CSV 为数据文件的 SHA-256 与大小，Bridson 见下节。续跑时重新生成或读取后比对，不符则拒绝 checkpoint。Random 不生成点集，NumPy 随机状态本已在 checkpoint 内（固定约 2.5 KB）。插件沿用现有身份绑定。
9. **不累积提交 UUID。** `stateless_batch.flush_batch` 每批调用 `CheckpointedSampler.record_submitted_batch`，Grid/Random/Bridson 原先把整次扫描的 UUID 留在内存；CSV 已改为空操作，现已提升到 `FixedSetSampler`（插件保持原行为）；`_scan_driver` 结束时只为非 indexed sampler 构造提交集合。

### Bridson 点集（重点）

Bridson 的 `_P` 不能当作可以随时重算的数据：

- **重算不保证一致。** 生成依赖全局 `np.random.seed(seed)` 和逐点接受/拒绝，接受判断比较浮点距离与 r²（`bridson.py` `in_neighborhood`）。NumPy 版本或 CPU/SIMD 路径带来的末位差异，只要翻转一次临界判断，之后的所有点都不同。
- **卡片变化不被发现。** `radius`、`max_attempt` 位于 `Sampling.Bounds`，现有续跑核对不覆盖；改了再续跑会在同一台机器上生成另一套点，游标却照旧前进。
- **空间可以很大。** 点数随 (length/radius)^d 增长，4 维小半径可达 10⁶–10⁷ 点。`_P` 为 N×d float64（10⁷×4 约 320 MB），纯 Python 逐点生成的耗时也随 N 增长。

设计：

- **`_P` 与游标分开保存。** fresh 生成 `_P` 后，以临时文件加 rename 原子写入 `checkpoints/<scan>/Bridson/points.npy`（`np.save`，不用 pickle），整个扫描期间不再重写。随后写入初始 state.pkl，记录 `points_sha256`（覆盖 shape、dtype 与字节）、点数和维度。
- **state.pkl 大小固定。** 只存游标、seed 与进度信息。
- **resume 不重新生成。** 用 `np.load(path, mmap_mode="r", allow_pickle=False)` 打开 points.npy，校验 SHA-256 与 checkpoint 记录一致后按游标继续。NumPy 或机器变化不影响续跑，大空间也省去重新生成的时间。内存映射避免把整份点集复制进进程内存；校验需要完整读一遍文件。
- **点集文件属于 checkpoint。** points.npy 缺失、损坏或哈希不符时 checkpoint 不完整，拒绝续跑并从头扫描。
- **卡片核对在读取点集之前。** radius、max_attempt、length、seed 任一变化，按第 7 条直接拒绝。
- **磁盘。** 每个 scan 一份，约 N×d×8 字节；新的 fresh 原子覆盖旧文件。check 模式没有 checkpoint，不写 points.npy。

Grid 不保存 `_P`：网格只由各变量的 `num` 决定，没有随机性，重算后比对 SHA-256 即可。

### 待查

Dynesty 的 `nested_engine.pkl` 看起来不受 `checkpoint.enabled` 控制（`dynesty_sampler.py` 无条件设置引擎 checkpoint 路径，`--resume` 时 `arm_engine_resume` 会读取）。不属于本协议，尚未核实。

实现核对：正式 fresh 在检测到既有 HDF5 shard 时设置 `replacing_existing`；新 scan 无旧库时不设置，仍可按自身 DATABASE 前缀优化续跑。初始 checkpoint 在 Archiver 启动及 Worker 就绪后、提交任何样本前写入。带标记续跑核实 DATABASE 连续前缀至少覆盖 checkpoint 游标，否则拒绝 checkpoint 并 fresh。百万点验证使用合成 4 维 points.npy，检验 mmap/哈希/不调用生成，未将合成坐标当作 Bridson 算法输出。

## 2026-10-17 采样器注册表与命名（N1–N6 已完成）

### 修复前缺陷

- `jsampler` 与 `Jarvis sampler` 共用 `~/.jarvis/samplers`（或 `$JARVIS_SAMPLERS_HOME`），但只有 HEP 重写的 `load`（`cli/sampler.py` `_load_sampler`）传入保留名；`jsampler load` 和两边的 `check` 都传空名单。直接用 `jsampler load` 就能登记 `Grid`、`grid`。
- 重名靠三处分散的拦截：`runtime/core.py` `preflight_sampler_plugin`、`distributor.py` `set_method`、`cli/sampler.py` info。只有卡片恰好用到重名方法时才报错；`_user_layer` 把重名条目从 catalog 隐藏，`list` 看不到。
- `check` 同时接受文件和已登记名。按名字 check 读的是快照，SDK `_check_target` 不校验快照 hash，HEP 为此另加 `_guard_check`。快照就是 load 时的源文件字节，按名字 check 不提供新信息。

### 原则（用户决定）

- 外部 sampler 不能与内置方法重名，大小写不敏感。这是硬约束，没有绕过选项。
- 保留名单只有一份，在 SDK。
- `Jarvis run` 与 `Jarvis check` 每次启动都生成注册表报告，内置与外部在内存里合并，内置方法不写进登记文件。
- 发现重名，运行时退出。

### 规则

1. **保留名单。** SDK 已定义 `RESERVED_METHOD_NAMES`，`check_method_name`（`metadata.py`）默认使用，不再需要调用方传入名单。`init`、`load`、`check` 在两个入口下拒绝同样的名字。读取登记和 `remove` 不做名字检查，否则重名条目既不能读也删不掉。HEP `reserved_method_names()` = SDK 名单 ∪ HEP catalog；HEP 测试断言 builtin 与 builtin-plugin 都在 SDK 名单内。新增内置方法先改 SDK 名单并发布，HEP 再提高 SDK 版本下限，不加上限。
2. **load 等价。** `Jarvis sampler load` 已直接转发 `jsampler load`；HEP `_load_sampler` 已删除。
3. **check 只接受文件路径。** 参数不是文件时 exit 2，提示用 `Jarvis sampler info <Method>` 找到源文件。导入文件后第一项检查方法名，单独列在报告首位；重名即失败，options/prepare/draw 不执行。SDK `_check_target` 的登记分支和 HEP `_guard_check` 已删除。
4. **注册表报告。** `sampler_catalog.registry_report()` 合并 builtin、builtin-plugin、user 三层。user 条目状态为 ok、与内置重名、快照缺失或被改（复用 SDK `SamplerRegistry.problems`）；登记整体无法读取时记一个全局错误。只读文件，不导入插件代码。登记部分沿用 `_user_layer` 的签名缓存，快照校验每次都做。
5. **启动时处理。** 报告在 `core.preflight_sampler_plugin` 开头生成，先于 Redis、Worker、Archiver 和 DATABASE 整理。
   - 任一条目与内置重名：退出（exit 2），不论卡片用哪个方法；列出全部重名条目和 `Jarvis sampler remove <X> --yes` 修复命令，`<X>` 用登记中的确切名字。`remove` 只删除外部登记及其托管快照，不动内置方法或用户源文件。
   - 卡片所用插件的快照缺失或被改：退出。
   - 其他条目的快照问题：警告并继续。
   - 登记无法读取：卡片用内置方法时警告并继续，用插件时退出。
6. **删除分散拦截。** `assert_no_reserved_name_collision` 及其在 `distributor.py`、`runtime/core.py`、`cli/sampler.py` 的调用已删除。不经 Core、直接在 Python 中使用 Distributor 时内置方法优先；当前 SDK 源码拒绝新重名登记，旧 SDK 或手写登记留下的条目由 run/check 启动时拦截。
7. **展示。** `Jarvis sampler list` / `info` 读同一份报告，重名条目照常列出并标明状态。
8. **CLI 风格。** `jsampler` 已改用 Jarvis CLI 风格，HEP 已删除 `cli/help.py` 自带的 SDK 解析器拷贝；`Jarvis sampler` 的帮助与错误由 SDK 按程序名 `Jarvis sampler` 渲染。损坏登记的 list/info/remove 统一 exit 2，`--version --json` 输出规范 JSON。全生态的帮助、错误、退出码与 JSON 规则见 [CLI 规范](JARVIS_CLI_DESIGN.md)。

SDK README 中"check 也接受已登记 sampler"的旧说明已修改。完整命令语义见 [技术设计](jarvis-hep-sampler/technical-design.md) §5.1、§5.2、§5.4。

官方打包的 builtin-plugin `LatinHypercube` 必须使用其 SDK 保留名。HEP 加载已核对的官方插件时，向 SDK 低层 metadata/loader 传入精确 `trusted_builtin_method="LatinHypercube"`；仅该规范大小写名称获得豁免，普通 SDK init/load/check 没有豁免入口。已有扫描的 `sampler_identity.json` 保存原来源与 SHA-256，preflight 以原扫描副本的字节及保存的 hash 验证并续跑；识别官方旧绑定使用 `jarvishep2/sampling/plugins/latin_hypercube.py` 的包内相对路径，不要求后来安装的包仍位于同一绝对路径或具有相同字节。项目中的绑定文件和扫描副本是用户可编辑状态：该协议能发现偶然损坏，不能防御有意伪造。

本地实现核对：注册表报告将旧重名登记保留为 `name-conflict` 行，每次读取重新核验快照；Core 在初始化 Redis 前针对全部重名、所选插件坏快照及不可读登记执行上述分流。`Jarvis sampler list/info` 的文本主表与 JSON 使用同一报告状态，显示登记中的精确名字及移除命令。N1/N2 SDK 本地源码 234 passed 且 Ruff 通过；HEP 最终定向集成 169 passed / 4 deselected，快速全套门禁 1597 passed / 196 deselected（最后三项 CLI 修复经后续定向复跑）；N5 非 CLI 专测 5 passed，N6 及真实 run/check 入口 4 passed。N4 已转发 SDK 并移除 HEP 自有 `_load_sampler`、`_guard_check` 和解析器拷贝。N3 使用已公开且包含 N1/N2 的 SDK 1.0.3；HEP 打包依赖、运行时版本门禁及默认环境卡同步采用 `>=1.0.3`，以免安装后由运行时接受低于该协议的 SDK。2.x 继续受支持。全新 Python 3.11 环境从本地 HEP 源码完成全依赖安装（66 个包），依赖来自公共 PyPI，已安装 SDK 1.0.3 与 HEP 运行时版本门禁一致，`pip check` 及双入口 JSON 版本输出通过。此安装验证不代表 HEP 已公开新版本。详情见 [验证证据](hep-review-2026-10-16-verification.json)。

已公开的 SDK 1.0.3 包的 PyPI 长描述内嵌了旧 README，仍写 1.0.2 和 HEP `>1`。这是发布描述的历史快照，不能据此推断 wheel 缺少 N1/N2；本地 SDK README 已改正，公开页须待后续版本才能同步。当前 HEP 依赖约束以其 `pyproject.toml` 为准。

### 复查遗留清理 C1（已实现）

续跑 sampler 的运行策略、持久化 UUID/prefix/pending 和按前缀推进，在 sampler 初始化后及 Archiver 启动后都通过同一个 `Jarvis2Core._bind_sampler_persistence()` 方法体绑定。原无调用的 `apply_resume_checkpoint` 与 Core 包装已删除。现有回归以 MagicMock 模拟 Core，故 supervisor 两处采用显式类方法调用，不改变实际 Core 行为。队列清理与重新提议标记仍保留各自原位置。指定三组非 slow 回归 58 passed、Core collaborator 17 passed；真实 Redis 慢测 20 passed、54 deselected，用时 278.16 秒。

## Diver 移植（暂留，未验证）

`jarvishep2/sampling/Source/Diver/` 是 Diver 1.3.0 的 NumPy 移植，没有注册为 `Sampling.Method`，目前只有 `tests/test_diver_engine.py` 使用。移植后没有做过数值验证。目录内 `LICENSE` 为仅限学术使用的条款，与 HEP 的 MIT 发行不一致；pyproject package-data 列入了它的 LICENSE 与 VENDOR.md，按当前工作区构建的 wheel 会带上这份源码。已发布的 2.0.16、2.0.17 不含它。

审阅指出了这一点。用户决定暂时保留现状：完成验证后联系原作者，说明所做的工作并确认许可，同时请对方审阅源码或建议测试。在此之前不把 Diver 作为已验证功能对外描述。
