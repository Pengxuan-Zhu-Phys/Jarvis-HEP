# Jarvis-HEP 实施路线

当前路线按用户2026-10-15最新授权：公开发布 HEP 2.0.17，所有外部 Sampler 工作延后。旧候选记录保留在原集成分支，不作为当前源码能力或验收结论。


## 2026-10-15 2.0.17 SDK-free 公开发行（最新授权）

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
