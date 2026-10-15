# Jarvis-HEP 架构审阅

2026-10-11 根据当前实现重建，未恢复历史审阅原文。

## 当前发行架构

2.0.17 使用 master 0bbcf675 的内置采样器与依赖边界，移植不依赖 SDK 的可靠性、CLI agent/YAML 和 Redis 功能。当前没有外部 Sampler SDK 依赖、插件注册表或插件 CLI。0d859502 / 52876060 的 SDK 集成与其它工作区改动留在原分支；下方保留的旧候选说明仅作历史记录。

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



## 2026-10-15 SDK-free 发行隔离架构（最新）

本次源自公开 master `0bbcf675`，只移植不依赖外部 SDK 的 Worker 启动/表达式诊断/monitor/CLI 可靠性、CLI agent/YAML 与 Redis 引导。内置 sampler 的 factory、checkpoint/恢复和调用边界沿用公开基线，不带入 SDK 集成 `0d859502` / `52876060` 或后续未提交 Sampler/runtime/Diver 编辑。

HEP 的 wheel/sdist 不应包含 SDK/plugin registry/CLI，不应有 `Jarvis-HEP-Sampler` Requires-Dist。构建后的完整依赖安装使用公共 PyPI，不接受本地 SDK 或 `--no-deps` 作为发行验收。实际 fast gate、installed CLI、pip check 和真实 Redis/spawn 重新执行，不能引用已被替代的 `1c8c82fc` SDK 候选作为当前证明。

当前验证与公开发布状态：见下方已完成的源码/安装验收与发布状态；证据见 [2.0.17 发行记录](release-2.0.17.md)。原工作区/分支的 Sampler 开发继续独立保留。

## 本次源码与安装验收状态

干净源码 `e8af5221`：1426 passed、3 skipped（缺少 sibling Jarvis-Examples）、164 deselected、477 subtests；335 项代码/测试文件与提交逐字节一致。全新隔离环境从公共 PyPI 解析全部依赖，pip check 与 Twine 通过，SDK 模块/发行包均不存在；真实 Redis/spawn/双 Worker/monitor/项目打包重跑、installed agent 24 项、installed Redis 39 项、发行回归 22 项均通过。

源码分离提交 `3b9d41d0`，基础 ValidationReport 兼容修正 `e8af5221`；修正只调用公开基线的 warnings 接口，新增正常文本验证和警告回归。wheel 的 276 个运行文件、sdist 的 401 项 tracked 文件与提交一致，无 SDK/plugin bridge、对应测试或 Diver。当前本地验收完成，正式发布执行中；实际 PyPI 成功与公开安装结果随后记录。见 [发行记录](release-2.0.17.md) 与 [验收证据](release-2.0.17-verification.json)。
