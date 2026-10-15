# Jarvis-HEP 产品要求

2026-10-11 重建的维护基线，来源为当前 README、安装文档、运行接口与 [2.0.15 台账](development-ledger-2.0.15.md)，不声称还原旧 PRD。

## 当前公开发行范围

2.0.17 从公开 master 0bbcf675 基线移植不依赖外部 SDK 的可靠性、CLI agent/YAML 和 Redis 功能，使用内置采样器；不声明外部 Sampler SDK 依赖、插件加载或注册命令。用户已授权正式公开发布，验收与发布结果见下方已完成的源码/安装验收与发布状态。后续历史阶段中的 SDK 要求和候选限定不适用于本次发行。

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



## 2026-10-15 2.0.17 SDK-free 公开发行要求（最新）

- 从公开 master `0bbcf675` 保留内置采样，移植不依赖 SDK 的可靠性、CLI agent/YAML 与 Redis 改进；本次没有外部 SDK/plugin 能力或依赖。
- 用户授权公开 HEP，推迟已提交 SDK 集成 `0d859502` / `52876060` 及未提交 Sampler/runtime/Diver 工作；不为发布 HEP 修改或发布 SDK。
- 公开候选应在干净环境仅从公共 PyPI 安装全部依赖，并通过 pip check、内置采样真实 Redis/spawn、installed agent/Redis/help/JSON、相关可靠性和完整 fast gate。
- `1c8c82fc` 旧候选的本地 SDK 依赖和候选限定已被替代；旧验收只作历史。当前公开成功须核验 GitHub 工作流、PyPI 版本/产物和全新公共索引安装。

当前验收与发布结果：见下方已完成的源码/安装验收与发布状态。以 [2.0.17 发行记录](release-2.0.17.md) 为当前事实源。

## 本次源码与安装验收状态

干净源码 `e8af5221`：1426 passed、3 skipped（缺少 sibling Jarvis-Examples）、164 deselected、477 subtests；335 项代码/测试文件与提交逐字节一致。全新隔离环境从公共 PyPI 解析全部依赖，pip check 与 Twine 通过，SDK 模块/发行包均不存在；真实 Redis/spawn/双 Worker/monitor/项目打包重跑、installed agent 24 项、installed Redis 39 项、发行回归 22 项均通过。

源码分离提交 `3b9d41d0`，基础 ValidationReport 兼容修正 `e8af5221`；修正只调用公开基线的 warnings 接口，新增正常文本验证和警告回归。wheel 的 276 个运行文件、sdist 的 401 项 tracked 文件与提交一致，无 SDK/plugin bridge、对应测试或 Diver。当前本地验收完成，正式发布执行中；实际 PyPI 成功与公开安装结果随后记录。见 [发行记录](release-2.0.17.md) 与 [验收证据](release-2.0.17-verification.json)。
