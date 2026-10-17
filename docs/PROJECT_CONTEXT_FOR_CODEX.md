# Jarvis-HEP 项目上下文

维护基线：2026-10-11。本文件在原约定路径及 Git 历史中均无旧版本，现依据当前代码、README、安装与发布文档重建；不代表恢复了遗失的历史原文。

> 当前源码状态（2026-10-17 提交复查后）：本轮复查修复使用 SDK 源码的 `HOST_CONTRACT_VERSION >= 2`，包含入口保留名策略和 `CatalogExtension` 数据接口；已发布的 SDK 1.0.3 尚无该契约。两仓版本号和 `>=1.0.3` 的历史最低声明保持不变，不加上限。HEP 运行时额外核查契约，旧 SDK 给出匹配源码构建的安装诊断；维护者后续发布 SDK 后再将 HEP 版本下限对齐真实发行号。本轮不发布。

## 产品与执行模型

发行包为 `Jarvis-HEP`，Python 包为 `jarvishep2`，用户命令为 `Jarvis`。源码版本为 2.0.17。2.0.17 已公开不依赖外部 Sampler SDK 的 GitHub Release 与 PyPI wheel；完整公共索引安装已核验，PyPI 源码包上传超时，源码包由 GitHub Release 提供。当前 master 在该发行之后要求 `Jarvis-HEP-Sampler>=1.0.3`，不设上限。Python 最低版本为 3.10。

每次扫描由控制进程管理独立 Redis、spawn Worker 和 Archiver。Worker 执行 Operas 或外部 calculator；Archiver 持久化样本结果，控制进程负责采样、等待、恢复和清理。监控按已知扫描与进程读取运行数据。

Jarvis-Operas 提供函数与常量注册，Jarvis-HEP-Portal 提供 calculator I/O，JarvisPLOT 提供绘图，Jarvis-HEP-Sampler 提供插件采样器，均为完整发行包依赖。`Jarvis sampler` 转发 SDK 命令并通过 HEP 注册表报告补充内置方法的 list/info 展示与运行前校验；方法保留名单由 SDK 维护。2.0.17 那次公开发行不依赖 Jarvis-HEP-Sampler。

## 本轮授权与范围

当前 2026-10-17 阶段落实注册表与命名 N1–N6、续跑代码清理 C1：SDK 兄弟仓库与 HEP 本地源码均已修改。维护者已公开 SDK 1.0.3，HEP 已将最低依赖版本提高至该版且不设上限；HEP 本阶段不改自身版本、不发布。此前 2026-10-16 HEP 单仓修复基于 master `cdb0a116`，当时 SDK 不在范围内。发布脚本维持每个项目一个独立脚本：`jhrel` 为 `61781536` 的独立版本，唯一选项 `--dry`（不加即发布），`scripts/full-release.sh` 已删除；本阶段不修改发布脚本。

工作区原有 Diver 源码、测试和 pyproject package-data 改动保持原样，不纳入本次修复。Diver 移植未做验证，其 LICENSE 仅限学术使用；用户决定暂留，验证后联系原作者确认许可，见 [架构审阅](ARCHITECTURE_REVIEW_JARVIS_HEP.md)"Diver 移植"。2.0.17 SDK-free 发布记录保留为历史，不能将当时“推迟 Sampler”的范围误当成当前 master 的依赖状态。

## 维护规则

代码修改前同时阅读本文件、PRD、架构审阅、代码地图和实施路线。阶段完成后同步这五份维护文件与台账。

保留现有命令、YAML 结构和正常 Worker 恢复行为；诊断应给出路径、原因和纠正方法。修复后运行相关定向测试、真实运行对照和完整 fast gate；发行包还需在干净环境安装全部依赖并核验。

## 关联维护文件

- [产品要求](PRD_JARVIS_HEP_MASTER.md)
- [架构审阅](ARCHITECTURE_REVIEW_JARVIS_HEP.md)
- [代码地图](CODE_MAP_JARVIS_HEP.md)
- [实施路线](IMPLEMENTATION_ROADMAP.md)
- [Jarvis 生态 CLI 设计规范](JARVIS_CLI_DESIGN.md)（跨项目）
- [发布流程](../RELEASING.md)



## 2026-10-15 2.0.17 SDK-free 公开发行（历史发布授权）

用户明确要求这次公开发布 HEP，把所有 Sampler 工作后移，取代仅完成候选的选择。发行从公开 master `0bbcf675` 隔离构建，移植不依赖外部 SDK 的 2.0.15 可靠性、2.0.16 CLI agent/YAML 与 2.0.17 Redis 改进；保留内置采样。已提交 SDK 集成 `0d859502` / `52876060` 和原工作区所有未提交 Sampler/runtime/Diver 编辑都延后，HEP 不声明 Sampler SDK 依赖，不需要 SDK 先公开。

旧 `1c8c82fc` 候选及其 draft/本地 SDK 安装证据被此发行源替代；它们保留为历史，不能复用其 checksum 或测试计数。当前必须对隔离源码和完整 wheel/sdist 重新运行 fast gate、公共 PyPI 全依赖安装/pip check、installed agent/Redis、真实 Redis/spawn，发布并核验公共索引安装；本次 OIDC 配置失败后的原样 CI 产物上传例外见下方记录。

当前结果：见下方已完成的源码/安装验收与发布状态。提交、产物和公开发布状态见 [2.0.17 发行记录](release-2.0.17.md)。原 Sampler 分支与用户工作区保持原样。

## 本次源码、CI 与公开安装状态

源码功能提交 `3b9d41d0`、验证接口兼容修正 `e8af5221` 已进入 master；公开 tag `v2.0.17` 指向 `8e0add06`。macOS 与 Linux CI 均为 1426 passed、3 skipped（缺少 sibling Jarvis-Examples）、164 deselected、477 subtests；分别用时 205.99s / 265.40s。Linux CI 的完整安装/pip check、真实 Redis/spawn/双 Worker/monitor/项目重跑、agent 24、Redis 39、发行回归 22 均通过。

[GitHub Release](https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP/releases/tag/v2.0.17) 已公开。OIDC publisher 因 PyPI 未匹配 Trusted Publisher 配置而报 `invalid-publisher`；本次使用既有本地 PyPI 凭证上传成功 CI 的原样产物，没有重建。PyPI wheel 已公开，全新环境仅从公共索引安装完整依赖，pip check、无 SDK、276 个运行文件身份比对、Redis 39 / agent 24 均通过。PyPI 源码包上传未获服务器响应；最后一次有界请求在 36.913 秒后 ReadTimeout，全部上传进程已停止。PyPI wheel 与公开安装已完成，已验收源码包可从 GitHub Release 下载。

发布 wheel SHA256 为 `25257b646c1a7ce828dae136e1e111de7c70d749d7b5b6857b7c687f227fb3ae`；CI sdist 的预期 SHA256 为 `1843e006921161dd6a46507b7b7c12b52b48990b9e3af1c76541788dafbe9915`，未在 PyPI 公开；GitHub Release 另提供通过本地完整验收的源码包。发布前本地产物的 a823… / 9171… hash 仅作本地验收记录，不能当成公开 CI 产物 hash。完整状态见 [发行记录](release-2.0.17.md) 与 [验收证据](release-2.0.17-verification.json)。

## 2026-10-16 HEP 侧审阅修复

正式 fresh 即使读取到旧 DATABASE，也必须重新提交当前卡片的所有样本；只有最终运行策略为 resume 时，才允许在批次提交中按持久化 prefix/pending 或 CSV legacy UUID 跳过工作。check 维持独立 DATABASE/test 的追加烟测流程，不使用 resume 批次过滤。运行策略单独保存在 sampler 的进程状态中，不进入 checkpoint；checkpoint 被拒转 fresh 时同步关闭跳过逻辑。

该轮曾由 HEP load 传入保留名，并在 preflight/Distributor/info 分散处理旧冲突；HEP CLI 曾单独校验 dimension、注册表与按名字 check 的快照 SHA-256。这些边界已由 2026-10-17 的 SDK 默认保留名、file-only check、统一注册表报告及 SDK CLI 转发取代。历史验收见 [本轮证据](hep-review-2026-10-16-verification.json)。

## 2026-10-16 续跑协议（R1–R6 已实现）

上一段记录的是前一轮修复；本轮已修正 fresh 中断后 `--resume` 把旧行算作已完成的缺陷。现行规则：

- Grid/Random/Bridson/CSV 的进度只以 checkpoint 为准，DATABASE 只证明本次运行写入的行；这四种方法没有 checkpoint 就没有 resume，`checkpoint.enabled: false` 时每次都从头扫描。Dynesty 的独立引擎文件仍属待查，未纳入此协议。
- check 模式只用于测试，不走本协议，checkpoint 强制关闭。
- 插件 SDK 没有自己的 checkpoint：independent 插件的点只由 (seed, 候选编号) 决定，续跑完全由 HEP 的 indexed resume 负责；SDK 只保证可复现（`check_plugin` 的 resume 项）。
- checkpoint 只存游标与固定大小的指纹，不随扫描规模增长。Bridson 的 `_P` 单独存为 points.npy，fresh 原子写一次，续跑内存映射读取并校验哈希，不重新生成。带 `replacing_existing` 的运行完成后会跨 HDF5 shards 清理旧尾部和重复 index。

完整规则与 Bridson 细节见 [架构审阅](ARCHITECTURE_REVIEW_JARVIS_HEP.md)，验收要求见 [产品要求](PRD_JARVIS_HEP_MASTER.md)。

本轮验收：扩展定向 162 passed；真实 Redis/spawn/HDF5 中断续跑（4 方法 × 心跳前后）、长短扫描清尾、check 与无 checkpoint 回退共 12 passed；百万点（4 维、32,000,128 字节）点集读取与哈希耗时记录在 [验证证据](hep-review-2026-10-16-verification.json)。百万点文件为合成坐标，检验的是 Bridson 续跑读取路径，不声称执行了百万点生成。

2026-10-17 C1 已删除无调用的续跑包装，sampler 初始化后与 Archiver 启动后的持久化绑定现共用 Core 的同一方法体；指定三组非 slow 续跑回归 58 passed，Core collaborator 17 passed。两处以显式类方法形式调用，以保持现有 MagicMock 回归不变。

## 2026-10-17 注册表与命名（N1–N6 已完成）

现行规则：

- 外部 sampler 不能与内置方法重名（大小写不敏感）；SDK 的 `RESERVED_METHOD_NAMES` 是唯一保留名单。`check_method_name` 只查语法；`check_external_method_name` 只由 init/load/check 入口调用。通用 metadata/loader 不拒绝内置插件名字，全部 `trusted_builtin_method` 参数及路径尾缀判定已删除。
- `Jarvis sampler load/check` 直接转发 SDK。check 只接受文件路径，报告首项检查方法名，重名 exit 2；按名字 check、HEP `_load_sampler`/`_guard_check` 和自带 SDK 解析器均已删除。SDK 统一处理非正 dimension、损坏登记的 list/info/remove（exit 2）及 `--version --json`。
- `Jarvis run` 与 `Jarvis check` 每次启动都生成只读注册表报告，内置与外部在内存中合并；发现任一重名即在 Redis 前退出。`Jarvis sampler list/info` 使用同一报告并显示旧冲突与修复命令。
- `Jarvis sampler remove` 只删除外部登记及其托管快照，不动内置方法或用户源文件；旧重名登记仍可按登记中的确切名字移除。
- 官方 `LatinHypercube` 已绑定扫描用原扫描副本和保存的 SHA-256 续跑，容忍后来安装的包路径或源码字节变化。项目中的 `sampler_identity.json` 和扫描副本属于用户可编辑状态；校验能发现偶然损坏，不能防御有意伪造。
- SDK 与 `Jarvis sampler` 已按 [CLI 规范](JARVIS_CLI_DESIGN.md) 对齐；生态中其他 CLI 的差距另列步骤。

08b2dab5 / d9d2853 阶段的历史源码验收：SDK 234 passed，Ruff 通过；HEP 最终定向集成 169 passed / 4 deselected，快速全套门禁 1597 passed / 196 deselected（之后对最后三项 CLI 修复复跑定向）。HEP N5 非 CLI 专测 5 passed、旧冲突真实 run/check 与 list/info/remove 入口 4 passed。C1 的真实 Redis 慢测 20 passed、54 deselected，用时 278.16 秒。N3 的 SDK 1.0.3 已公开，HEP 的 `pyproject.toml`、运行时版本门禁和环境默认卡均要求 `>=1.0.3`。全新 Python 3.11 环境从本地 HEP 源码构建并完整安装 66 个包，依赖均从公共 PyPI 解析，选中 SDK 1.0.3；`pip check` 通过，已安装的 HEP 与 SDK 模块来自该环境，`Jarvis sampler` / `jsampler --version --json` 都返回 1.0.3，HEP 运行时版本门禁返回 1.0.3。公开 SDK 的 N1/N2 行为也已在干净环境核验。HEP 本身仍未发布新版本。详细结果见 [验证证据](hep-review-2026-10-16-verification.json)。完整规则见 [架构审阅](ARCHITECTURE_REVIEW_JARVIS_HEP.md) 与 [实施路线](IMPLEMENTATION_ROADMAP.md) N1–N6。

## 2026-10-17 提交复查收尾

前次提交的清尾将所有分片合并，且未确认 Archiver 已退出；此前验收遗漏这两个条件。现按原分片清尾，未变化文件不重写；停止失败跳过清尾，磁盘/关闭错误不会跳过其他资源收尾，最终报告失败并保留原扫描异常。覆盖扫描在清尾后出图。

保留名政策现仅在 init/load/check 入口执行，元数据和运行时 loader 不再查保留名；已删除 trusted_builtin 豁免及路径尾缀判断。SDK 负责 list/info 全部展示，HEP 的 CatalogExtension 只提供数据；内置 info 也列出大小写不同的旧冲突。重名 check 统一 exit 2，算法检查失败 exit 1。验证以现有 JSON 的提交复查记录为准。

提交复查验证完成：SDK 全套 253 passed、Ruff 通过；HEP 最终快速门禁 1627 passed / 197 deselected / 478 subtests，真实 Redis/spawn/HDF5 慢测 21 passed，最后边界补丁后长短扫描两种分片模式复跑 2 passed，收尾定向 40 passed。匹配暂存源码的两个 wheel 在干净 Python 3.11 环境完整安装 66 个包，pip check、官方插件、分片清尾与双入口 CLI 通过；公开 SDK 1.0.3 的旧契约诊断也经真实换装核验。wheel 内 HEP 227 / SDK 16 个 Python 文件逐字匹配暂存内容。详细命令、耗时、构件哈希和范围保护见现有 [验证 JSON](hep-review-2026-10-16-verification.json) 的 `review_followup_2026_10_17`；本轮未改发行版本、未发布。
