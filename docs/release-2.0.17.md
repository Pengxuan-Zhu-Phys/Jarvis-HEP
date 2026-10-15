# Jarvis-HEP 2.0.17

1. **Roadmap phase**：按最新授权公开发布不依赖外部 Sampler SDK 的 HEP 2.0.17；从公开 master 隔离移植可靠性、CLI agent/YAML 与 Redis 改进。
2. **Affected modules**：不依赖 SDK 的 Worker/表达式/monitor/CLI 可靠性、YAML 诊断、agent 辅助、Redis 引导、对应测试文档、版本与发行流程。
3. **Public interfaces**：CLI agent 合同/查询/schema/模板/inspect/草稿与 Redis 正式子命令/help/JSON；保留内置采样、YAML 结构和正常扫描恢复，不包含外部 SDK/plugin CLI。
4. **Main risks**：误带入已提交或未提交的 Sampler 集成；沿用旧 SDK 候选证据；发行包依赖或 metadata 残留 SDK；公开安装未经核验。
5. **Verification performed**：macOS 与 Linux CI 1426 passed / 477 subtests、3 skipped、164 deselected；完整安装/pip check/真实运行、agent 24 / Redis 39 / 发行回归 22 通过。GitHub Release 与 PyPI wheel 已公开，全新公共索引全依赖安装、SDK 缺失和运行文件身份比对通过；OIDC 配置失败后的原样 CI 产物上传例外已记录，PyPI 源码包上传超时，GitHub Release 提供已验收源码包。

## 当前授权与源码范围

用户最新明确要求公开 HEP，并将所有 Sampler 工作后移，取代此前“暂只完成候选”的选择。本版从公开 `master` 基线 `0bbcf675` 的独立 worktree 准备，只移植不依赖外部 SDK 的 2.0.15 可靠性、2.0.16 CLI agent/YAML 和 2.0.17 Redis 改进。已提交 SDK integration `0d859502` / `52876060` 与原工作区未提交的 Sampler/runtime/Diver/审阅资料都延后，不改动原工作区，不合并或发布 SDK。

采样继续使用 HEP 内置实现，发行包不声明 `Jarvis-HEP-Sampler` 依赖，不提供其插件注册/加载 CLI。旧 `1c8c82fc` SDK-dependent 候选和 draft 被这份公开发行源取代；旧产物/checksum/测试保留为历史，不能复用为当前结果。

## 用户可见变化

- Worker 初始化故障有界反馈、表达式预检诊断、monitor JSON/缺失值与 CPU snapshot 等不依赖 SDK 的可靠性修复；正常运行期恢复保留。
- Agents 可查询 CLI 合同、命令效果和参数，导出本地 schema/准确模板，解释卡片生效值及来源，并沿 validate JSON 导航修正 YAML；显式重复键拒绝，合法 merge override 保留。
- 功能反馈生成可填写邮件草稿，不自动发送；查询不隐含授权执行 check/run 或用户代码。
- `Jarvis redis info/check/install/start` 为正式命令组，含独立 help、JSON 与 exact agent 发现；缺失程序时提供按发行版/版本选择的安装指导。安装/启动命令只打印，正常扫描继续自动管理本次 Redis。

## 当前验收与发布状态

功能提交 `3b9d41d0dc2171bdaf79990a4e23faddd9ff6ce9` 和 ValidationReport 兼容修正 `e8af5221eb84b4c1feccf492b2676df764505b3f` 已进入 master；公开 tag `v2.0.17` 为 `8e0add06dbcfb6020bd17f8af1a281105906c44a`。兼容修正只调用公开基线的 warnings 接口，并增加正常文本验证与警告回归。内置采样与原 Sampler 分支/工作区保持原边界。

macOS 本地和 Linux Python 3.10 CI full gate 均为 **1426 passed、3 skipped、164 deselected、477 subtests**；用时分别为 205.99s 和 265.40s。三个 skip 是干净检出没有 sibling Jarvis-Examples。Linux CI 构建/完整依赖安装/pip check/Twine、真实 Redis/spawn/双 Worker/monitor/项目打包重跑、installed agent 24、installed Redis 39 与发行回归 22 均通过。公开 wheel 的 276 个运行文件与发行 tag 相符，无 SDK 模块/发行包/依赖。

[GitHub Release v2.0.17](https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP/releases/tag/v2.0.17) 已公开。[发行工作流](https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP/actions/runs/38694564478) 的 build/验收 job 成功，但 OIDC publisher 因 PyPI 没有匹配 Trusted Publisher 配置而失败，错误为 `invalid-publisher`。本次使用既有本地 PyPI 凭证上传成功 CI 的原样产物，没有重建；不将该 job 误记为成功。后续正常 OIDC 发行前需补齐对应 PyPI 配置。

PyPI wheel 已公开。全新隔离环境以 `pip install --index-url https://pypi.org/simple Jarvis-HEP==2.0.17` 安装全部依赖，通过 pip check、2.0.17 版本、SDK 缺失、276 个运行文件身份比对、installed Redis 39 与 agent 24。没有使用本地 SDK、`--no-deps` 或 system-site-packages；公开安装结果不是由本地 candidate 安装替代。

| 产物 | SHA256 | 当前状态 |
| --- | --- | --- |
| CI / PyPI wheel | `25257b646c1a7ce828dae136e1e111de7c70d749d7b5b6857b7c687f227fb3ae` | PyPI 已公开，fresh public-index install 通过 |
| CI sdist | `1843e006921161dd6a46507b7b7c12b52b48990b9e3af1c76541788dafbe9915` | PyPI 请求超时，未公开；GitHub Release 提供已验收源码包 |
| 发布前本地 wheel | `a823b9c71a7959fdb2ca5a45b4d03e9190058d0a47ed14f97fe4f76393c53c1b` | 仅为本地预发布验收产物，不是公开 wheel hash |
| 发布前本地 sdist | `91711348da271fe6356791c0c945033b1c8bae2f03de466599242418ff10a3a2` | 仅为本地预发布验收产物，不是公开 sdist hash |

公开 wheel 的发布与完整安装验收已完成。PyPI 源码包上传在多次请求及最后一次 35 秒 read-timeout 的有界重试后仍无响应，上传进程均已停止；不声称源码包已在 PyPI 上架。GitHub Release 已提供通过完整验收的源码包，不影响正常 wheel 安装。完整命令、源码/产物范围、CI 和公开安装摘要见 [结构化证据](release-2.0.17-verification.json)。原始记录位于 `dist/2.0.17-sdk-free/` 的 `ci-publication-verification.json` 与 `public-install-verification.json`。本次补录维护文档不改变已验收/已公开产品源码或产物。

## 历史候选说明

此前 `1c8c82fc` 候选使用本地SDK完成测试并创建draft，但其外部SDK依赖不适合本次公开发行，已被公开 master `0bbcf675` 派生的无SDK依赖版本取代。旧源码、产物、checksum和验收记录保留在原集成分支，不构成本次发行证据；本次不再等待或发布SDK。
