# Jarvis-HEP 实施路线

当前路线按用户2026-10-15最新授权：公开发布 HEP 2.0.17，所有外部 Sampler 工作延后。旧候选记录保留在原集成分支，不作为当前源码能力或验收结论。


## 2026-10-15 2.0.17 SDK-free 公开发行（最新授权）

| 阶段 | 范围 | 当前状态/结束门槛 |
| --- | --- | --- |
| 源码隔离 | 从公开 master `0bbcf675` 移植不依赖 SDK 的可靠性/CLI agent/Redis | 隔离 worktree 准备中；不带入已提交 SDK integration `0d859502` / `52876060` 或 dirty Sampler |
| 依赖与打包 | 内置采样；无外部 SDK Requires-Dist / plugin CLI | 干净 wheel/sdist 内容和公共索引全依赖安装核验 |
| 当前源码验收 | fast gate、installed agent/Redis、pip check、真实 Redis/spawn、发行回归 | 待本阶段测试结果补录；旧 `1c8c82fc` 候选不作替代 |
| 公开发布 | GitHub Release/Trusted Publishing、PyPI、全新安装 | 用户已授权；实际发布结果待记录，完成后再标记公开成功 |

此阶段取代旧 SDK 依赖阻断与候选限定方案，所有 Sampler 工作延后且保留原分支/工作区。当前状态以 [2.0.17 发行记录](release-2.0.17.md) 为准。
