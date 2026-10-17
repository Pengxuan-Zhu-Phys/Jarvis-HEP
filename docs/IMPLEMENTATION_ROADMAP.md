# Jarvis-HEP 实施路线

当前 master 要求 `Jarvis-HEP-Sampler>=1.0.3`，不设上限。`Jarvis sampler` 复用 `jsampler` 命令面，HEP 侧增加内置方法与安全校验。下方 2.0.17 记录是已发布的 SDK-free 发行，不改写。


> 当前源码状态（2026-10-17 提交复查后）：本轮复查修复使用 SDK 源码的 `HOST_CONTRACT_VERSION >= 2`，包含入口保留名策略和 `CatalogExtension` 数据接口；已发布的 SDK 1.0.3 尚无该契约。两仓版本号和 `>=1.0.3` 的历史最低声明保持不变，不加上限。HEP 运行时额外核查契约，旧 SDK 给出匹配源码构建的安装诊断；维护者后续发布 SDK 后再将 HEP 版本下限对齐真实发行号。本轮不发布。

## 2026-10-15 2.0.17 SDK-free 公开发行（已发布）

| 阶段 | 范围 | 当前状态/结束门槛 |
| --- | --- | --- |
| 源码隔离 | 从公开 master `0bbcf675` 移植不依赖 SDK 的可靠性/CLI agent/Redis | 隔离源码已提交 `3b9d41d0` / `e8af5221`；不带入已提交 SDK integration `0d859502` / `52876060` 或 dirty Sampler |
| 依赖与打包 | 内置采样；无外部 SDK Requires-Dist / plugin CLI | 已通过：wheel 276 个运行文件 / sdist 401 项 tracked 文件；公共依赖完整安装与 pip check |
| 当前源码验收 | fast gate、installed agent/Redis、pip check、真实 Redis/spawn、发行回归 | 见下方已完成的源码/安装验收与发布状态；旧 `1c8c82fc` 候选不作替代 |
| 公开发布 | GitHub Release/Trusted Publishing、PyPI、全新安装 | GitHub Release / PyPI wheel / 公开完整安装已通过；OIDC 配置例外已记录，PyPI 源码包上传超时，源码包见 GitHub Release |

此阶段取代旧 SDK 依赖阻断与候选限定方案，所有 Sampler 工作延后且保留原分支/工作区。当前状态以 [2.0.17 发行记录](release-2.0.17.md) 为准。

## 本次源码、CI 与公开安装状态

源码功能提交 `3b9d41d0`、验证接口兼容修正 `e8af5221` 已进入 master；公开 tag `v2.0.17` 指向 `8e0add06`。macOS 与 Linux CI 均为 1426 passed、3 skipped（缺少 sibling Jarvis-Examples）、164 deselected、477 subtests；分别用时 205.99s / 265.40s。Linux CI 的完整安装/pip check、真实 Redis/spawn/双 Worker/monitor/项目重跑、agent 24、Redis 39、发行回归 22 均通过。

[GitHub Release](https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP/releases/tag/v2.0.17) 已公开。OIDC publisher 因 PyPI 未匹配 Trusted Publisher 配置而报 `invalid-publisher`；本次使用既有本地 PyPI 凭证上传成功 CI 的原样产物，没有重建。PyPI wheel 已公开，全新环境仅从公共索引安装完整依赖，pip check、无 SDK、276 个运行文件身份比对、Redis 39 / agent 24 均通过。PyPI 源码包上传未获服务器响应；最后一次有界请求在 36.913 秒后 ReadTimeout，全部上传进程已停止。PyPI wheel 与公开安装已完成，已验收源码包可从 GitHub Release 下载。

发布 wheel SHA256 为 `25257b646c1a7ce828dae136e1e111de7c70d749d7b5b6857b7c687f227fb3ae`；CI sdist 的预期 SHA256 为 `1843e006921161dd6a46507b7b7c12b52b48990b9e3af1c76541788dafbe9915`，未在 PyPI 公开；GitHub Release 另提供通过本地完整验收的源码包。发布前本地产物的 a823… / 9171… hash 仅作本地验收记录，不能当成公开 CI 产物 hash。完整状态见 [发行记录](release-2.0.17.md) 与 [验收证据](release-2.0.17-verification.json)。

## 2026-10-16 HEP bug 复查与修复

- 已确认 `c3544160` 删除无效 catalog 参数，原子命令 TypeError 不再存在。
- 已实现：fresh/check 不过滤旧持久化结果，resume 继续去重；checkpoint 拒绝回退同步 fresh 策略。
- 已实现：HEP 保留名注册拒绝、旧冲突运行前拒绝；HEP CLI 的 dimension/注册表/快照哈希与 card 预算诊断。
- 验收：新增矩阵与 CLI 回归，以及四种内置方法的真实 Redis/spawn/HDF5 fresh 重算。最终命令和门禁结果记入 [证据](hep-review-2026-10-16-verification.json)。
- 范围外：SDK 项目内问题、发布脚本、原有 Diver 移植（暂留，见架构审阅"Diver 移植"）。本轮不提交、不提升版本、不发布。

本轮验收完成：定向 78 passed；真实 fresh 重算 4 passed；真实 resume 故障矩阵 4 passed；最终快速门禁 1566 passed / 184 deselected / 478 subtests passed。外部 Examples 非 ASCII 用例与原有 Diver 草稿测试的排除原因已记录在验证证据中。

## 2026-10-16 续跑以 checkpoint 为准（R1–R6 已完成）

来源：复查 issue 1（fresh 中断后续跑保留旧行）。设计见 [架构审阅](ARCHITECTURE_REVIEW_JARVIS_HEP.md)，验收见 [产品要求](PRD_JARVIS_HEP_MASTER.md)。

| 步骤 | 范围 | 结束门槛 |
| --- | --- | --- |
| R1 已完成：没有 checkpoint 就从头扫描 | 删除 DATABASE-only 恢复；check 模式强制关闭 checkpoint | `enabled: false` / 无 state.pkl 的 `--resume` 从头扫描；check 不写 state.pkl |
| R2 已完成：fresh 水位与初始 checkpoint | fresh 立刻写游标 0 与 `replacing_existing`；Archiver 前缀从 0 开始；supervisor 仅 resume 载入 DATABASE | 第一次心跳前中断后续跑不跳到旧前缀 |
| R3 已完成：带标记的 resume 与清尾 | 只信游标以下的行；不按 DATABASE 推进；确认 Archiver 正常退出后跨 shard 删除超出总数的旧行并去重，保留分片，异常完成其他收尾 | 四种内置方法真实 Redis 中断续跑后库中无旧行；短扫描覆盖长扫描后无尾部 |
| R4 已完成：卡片与点集指纹 | Sampling（Bounds、selection）哈希；Grid `_P` 哈希；CSV 文件哈希 | 改 Bounds/selection/CSV 后拒绝续跑；只改 workers 照常续跑 |
| R5 已完成：Bridson 点集文件 | points.npy 原子写一次；resume 内存映射读取并校验，不重新生成 | ≥10⁶ 点 resume 不调用生成；缺失/损坏/哈希不符时从头扫描 |
| R6 已完成：内存 | indexed sampler 不累积 `_submitted_uuids` | 10 点与 10⁵ 点的 state.pkl 大小基本相同 |

R1–R6 阶段未改插件路径或 SDK。Dynesty `nested_engine.pkl` 与 `checkpoint.enabled` 的关系另列待查，不在本阶段；已有无 state.pkl 的 Dynesty resume 行为保留。

完成证据：扩展定向 162 passed；真实 Redis/spawn/HDF5 12 passed（四种方法的心跳前后中断续跑、长短扫描清尾、check）；Bridson 10⁶×4 点集 mmap 读取 0.00495 秒，SHA-256 校验 0.01155 秒，生成调用 0。详见 [现有验证记录](hep-review-2026-10-16-verification.json)。

## 2026-10-17 注册表与命名（N1–N6 已完成）

来源：复查 issue 3、5、6、8 及用户的硬约束（外部 sampler 不能与内置方法重名）。设计见 [架构审阅](ARCHITECTURE_REVIEW_JARVIS_HEP.md)"采样器注册表与命名"，命令语义见 [技术设计](jarvis-hep-sampler/technical-design.md) §5.1、§5.2、§5.4，验收见 [产品要求](PRD_JARVIS_HEP_MASTER.md)。

| 步骤 | 范围 | 结束门槛 |
| --- | --- | --- |
| N1 已完成：SDK 保留名单 | SDK `RESERVED_METHOD_NAMES` 由入口策略 `check_external_method_name` 默认使用；通用 metadata/loader 仅校验格式，init/load/check 无需调用方传名单，读取登记与 remove 不查保留名。remove 只删除外部登记及其托管快照，不动内置方法或用户源文件。官方打包的 `LatinHypercube` 通过通用 loader；已删除全部 trusted_builtin 豁免与路径判断 | `jsampler`/`Jarvis sampler` load 拒绝 `Grid`/`grid`/`LATINHYPERCUBE`；旧 `grid` 登记可读可删；HEP builtin 与 builtin-plugin ⊆ SDK 名单 |
| N2 已完成：SDK file-only check 与 CLI 诊断 | `_check_target` 登记分支已删除，报告首项为方法名；`--dimension` 非正数及 `--seed` 负数在导入前 exit 2；`_run_info` 捕获 RegistryError；SDK README 已改为 file-only。帮助、错误面板、退出码、JSON 字段按 [CLI 规范](JARVIS_CLI_DESIGN.md) §3–§8；帮助中旧版本上限文字已删除，README 保持无版本上限 | 方法名参数 exit 2；重名文件在 options/prepare/draw 前失败；损坏登记的 list/info/remove exit 2 且 JSON 含 `error.type`；`--version --json` 输出规范对象；SDK 本地源码 234 passed、Ruff 通过 |
| N3 已完成：SDK 发版与 HEP 下限 | 维护者已发布含 N1/N2 的 SDK 1.0.3；HEP `pyproject.toml`、运行时版本门禁及默认环境卡提高至 `>=1.0.3`，不加上限；SDK README 同步新下限 | 全新 Python 3.11 环境从本地 HEP 源码完成全依赖实际安装（66 个包），依赖从公共 PyPI 解析并选中 SDK 1.0.3；`pip check`、运行时版本门禁及双入口 JSON 版本通过；HEP 未公开新版本 |
| N4 已完成：HEP load/check 转发与 CLI 风格 | `cli/sampler.py` 的 `_load_sampler`、`_guard_check` 与 info 旧重名拦截已删除；load/check 直接转发 SDK；`cli/help.py` 的 SDK 解析器拷贝已删除，帮助和 agent 目录用 SDK `build_parser("Jarvis sampler")`；错误用 SDK 面板 | load/check 等无扩展命令替换程序名后帮助与错误一致、退出码一致；list/info 通过 SDK 扩展点增加 HEP 内置说明与报告；损坏登记含内置 info 均 exit 2，SDK 缺失时仍输出统一错误面板或单个 JSON 错误；按名字 check 的旧用例已换成 file-only 断言 |
| N5 已完成：注册表报告 | `sampler_catalog.registry_report()` 在内存合并三层且每次核验快照；`core.preflight_sampler_plugin` 开头调用；`assert_no_reserved_name_collision` 及三处调用已删除；`cli/run.py` 将错误映射 exit 2 | 真实 `Jarvis run`/`Jarvis check` 对 Random 卡片与旧 Grid 冲突均在 Redis 前 exit 2；N5 非 CLI 专测 5 passed |
| N6 已完成：展示 | `Jarvis sampler list/info` 从同一报告呈现外部登记和内置方法，文本与 JSON 都显示冲突及修复命令 | 旧 `Grid`/`grid` 登记的 list/info/remove 和移除后内置方法恢复，4 个入口测试通过 |

N1/N2/N4/N5/N6 已在本地源码整合，维护者公开了含 N1/N2 的 SDK 1.0.3，N3 的 HEP 依赖变更与公共索引验收已完成。最终 HEP 定向集成 169 passed / 4 deselected，快速全套门禁 1597 passed / 196 deselected；全套门禁先于最后三项 CLI 边界修复，定向集成在修复后复跑。官方 `LatinHypercube` 已绑定扫描继续使用原扫描副本与保存的 hash，容忍安装包绝对路径或字节变化；本地绑定和扫描副本可编辑，不能抵御有意伪造。CLI 风格的全生态要求见 [CLI 规范](JARVIS_CLI_DESIGN.md)，本节只覆盖 SDK 与 `Jarvis sampler`；`Jarvis` 自身及 Operas/PLOT/Portal 的差距按规范 §13–§14 另立步骤。Diver 不在本阶段：验证完成后由用户联系原作者。

## 2026-10-17 复查遗留清理（C1 已实现）

来源：复查 issue 9。行为不变，只收拢重复代码。

| 步骤 | 范围 | 结束门槛 |
| --- | --- | --- |
| C1 已完成：持久化绑定收拢 | 无调用的 `_resume_service.py` `apply_resume_checkpoint`、Core 包装与两处 docstring 引用已删除；`_runtime_supervisor.py` 原两处 sampler 持久化绑定共用 `Jarvis2Core._bind_sampler_persistence()`，各调用一次 | 指定三组非 slow 回归 58 passed；Core collaborator 17 passed；真实 Redis 慢测 20 passed / 54 deselected，用时 278.16 秒 |

## 2026-10-17 提交复查修复

| 范围 | 源码结果 | 验收 |
| --- | --- | --- |
| 清尾与关闭 | 按原分片清尾；无变化不改文件；CSV 导出失败重试修复、空 live 清除旧 CSV；确认 Archiver 正常退出，异常逐步收尾并报错；清尾后出图 | 故障注入、分片与重试回归、真实 Redis/spawn/HDF5 强制轮转 |
| 命名策略 | 通用 metadata/loader 只查格式；init/load/check 执行默认保留名单；删除所有 trusted_builtin 与路径判定；重名 check exit 2 | SDK 全仓、官方插件与旧绑定恢复 |
| list/info | SDK `CatalogExtension` 提供数据扩展与共享渲染；HEP 只提供状态与内置数据；大小写不同的旧冲突也显示 | SDK 与 HEP 文本/JSON、共享渲染器回归 |
| SDK 接口 | 两仓不改版本、不发布；host contract 2 在源码可用；发布 1.0.3 缺少时明确拒绝新宿主接口 | 本地匹配构件安装；维护者下一次 SDK 发布后再提高 HEP 依赖下限 |

这轮更正了“合并清尾仍满足分片约束”和“metadata 层必须拒绝保留名”的旧理解。先前 N3 的公共安装证据仍属于 08b2dab5 / d9d2853 阶段，不代表这轮源码已发布。实际测试结果记入现有验证 JSON。

提交复查验证完成：SDK 全套 253 passed、Ruff 通过；HEP 最终快速门禁 1627 passed / 197 deselected / 478 subtests，真实 Redis/spawn/HDF5 慢测 21 passed，最后边界补丁后长短扫描两种分片模式复跑 2 passed，收尾定向 40 passed。匹配暂存源码的两个 wheel 在干净 Python 3.11 环境完整安装 66 个包，pip check、官方插件、分片清尾与双入口 CLI 通过；公开 SDK 1.0.3 的旧契约诊断也经真实换装核验。wheel 内 HEP 227 / SDK 16 个 Python 文件逐字匹配暂存内容。详细命令、耗时、构件哈希和范围保护见现有 [验证 JSON](hep-review-2026-10-16-verification.json) 的 `review_followup_2026_10_17`；本轮未改发行版本、未发布。
