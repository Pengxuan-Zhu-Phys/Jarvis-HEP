# Jarvis-HEP 历史瘦身 To-Do

Last updated: 2026-03-02
目标: 清理已无运行价值的历史代码与打包负担，保持当前主流程行为不变。
阶段总结: [HISTORY_SLIMMING_SUMMARY.md](./HISTORY_SLIMMING_SUMMARY.md)

## 约束

- 只删除“已确认无主链路引用”的内容。
- 每一轮瘦身后必须跑最小回归测试。
- 不改动核心计算语义（Sampling/Workflow/Likelihood 的行为保持不变）。

## Phase 0: 盘点与分级（已完成）

- [DONE] 主运行链路复核: `client -> core -> distributor -> Sampling/*`
- [DONE] 找到无引用的历史模块与资产候选。

## Phase 1: 安全瘦身（本轮执行）

- [DONE] 删除无引用且会误导维护的历史日志模块:
  - `src/logger.py`
  - `src/logger_old.py`

- [DONE] 从打包配置移除历史目录:
  - `pyproject.toml` 中的 `Sampling/Old/**/*`

- [DONE] 清理失活配置项:
  - `src/card/preference.json` 中的 `PoissonSchema`
  - 删除 `src/card/possion.json`

- [DONE] 删除明显无引用的大体积历史目录/资产:
  - `src/Sampling/Old/`
  - `src/Sampling/Source/rjdnn.zip`

- [DONE] 最小回归:
  - `python -m unittest discover -s tests -p 'test_*.py'`
  - `python -m py_compile`（改动文件）

## Phase 2: 中风险瘦身（后续）

- [DONE] 清理 `dataconvert.py` 中未完成且无调用的旧 API 分支。
- [DONE] 清理 `Sampling/dnn.py` 中无用 import 与历史残留函数块。
- [DONE] 将第三方源码快照的“示例/文档冗余文件”迁移到外部归档（若需要保留历史）。

## Phase 3: 发布体积优化（后续）

- [DONE] 对 sdist/wheel 做体积基线统计。
- [DONE] 校验 PyPI 包内文件白名单，仅保留运行所需资源。
- [DONE] 修复 `project.license` 旧式 TOML table 警告，切换到 SPDX 字符串。

## 验收标准

- 主流程功能不回退。
- 单测通过。
- 包体积和仓库噪音显著下降。
- 删除项均可追溯到“无引用证据”。
