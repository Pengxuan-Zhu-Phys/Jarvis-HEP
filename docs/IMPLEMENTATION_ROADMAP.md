# Jarvis-HEP 实施路线

当前路线按用户2026-10-15最新授权：公开发布 HEP 2.0.17，所有外部 Sampler 工作延后。旧候选记录保留在原集成分支，不作为当前源码能力或验收结论。


## 2026-10-15 2.0.17 SDK-free 公开发行（最新授权）

| 阶段 | 范围 | 当前状态/结束门槛 |
| --- | --- | --- |
| 源码隔离 | 从公开 master `0bbcf675` 移植不依赖 SDK 的可靠性/CLI agent/Redis | 隔离源码已提交 `3b9d41d0` / `e8af5221`；不带入已提交 SDK integration `0d859502` / `52876060` 或 dirty Sampler |
| 依赖与打包 | 内置采样；无外部 SDK Requires-Dist / plugin CLI | 已通过：wheel 276 个运行文件 / sdist 401 项 tracked 文件；公共依赖完整安装与 pip check |
| 当前源码验收 | fast gate、installed agent/Redis、pip check、真实 Redis/spawn、发行回归 | 见下方已完成的源码/安装验收与发布状态；旧 `1c8c82fc` 候选不作替代 |
| 公开发布 | GitHub Release/Trusted Publishing、PyPI、全新安装 | 用户已授权；本地全部门禁通过，正式发布执行中，成功后记录公开安装结果 |

此阶段取代旧 SDK 依赖阻断与候选限定方案，所有 Sampler 工作延后且保留原分支/工作区。当前状态以 [2.0.17 发行记录](release-2.0.17.md) 为准。

## 本次源码与安装验收状态

干净源码 `e8af5221`：1426 passed、3 skipped（缺少 sibling Jarvis-Examples）、164 deselected、477 subtests；335 项代码/测试文件与提交逐字节一致。全新隔离环境从公共 PyPI 解析全部依赖，pip check 与 Twine 通过，SDK 模块/发行包均不存在；真实 Redis/spawn/双 Worker/monitor/项目打包重跑、installed agent 24 项、installed Redis 39 项、发行回归 22 项均通过。

源码分离提交 `3b9d41d0`，基础 ValidationReport 兼容修正 `e8af5221`；修正只调用公开基线的 warnings 接口，新增正常文本验证和警告回归。wheel 的 276 个运行文件、sdist 的 401 项 tracked 文件与提交一致，无 SDK/plugin bridge、对应测试或 Diver。当前本地验收完成，正式发布执行中；实际 PyPI 成功与公开安装结果随后记录。见 [发行记录](release-2.0.17.md) 与 [验收证据](release-2.0.17-verification.json)。
