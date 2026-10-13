# Jarvis-HEP 2.0.16

1. **Roadmap phase**：两轮 CLI agent 集成完成；2.0.16 提交、完整发行包验证和正式发布。
2. **Affected modules**：CLI agent/man/validate、cardload YAML loader、installed smoke、维护文档与版本声明。
3. **Public interfaces**：新增 agent 合同、search/inspect/schema/template/request；validate JSON 增量导航；未知手册后缀与显式重复 YAML 键现在报错。原默认完整 agent 合同保留。
4. **Main risks**：重复键旧卡须修正；schema/inspect 不证明用户算子或外部资源可运行；SDK 0.1.0 尚未上 PyPI，正式发布前必须先公开依赖。Diver 草稿不在发行范围。
5. **Verification**：两轮报告已保存；第二轮 332 tests/373 subtests、22 installed CLI checks 通过。工作区 fast gate 1416/478 通过，1 项已知外部 Examples 编码失败；本版本干净源码 fast gate 1396 tests/478 subtests、3 skipped，完整包安装及 pip check、真实 Redis/spawn smoke、24 项 agent CLI 检查和 22 项安装包回归全部通过。

## Changes

- Agents can discover the installed CLI contract, command effects, authorization rules, question routes and YAML authoring workflow through `Jarvis agent --json` or `Jarvis man agent-readme --json`.
- Brief/section/command views, bounded keyword search, bundled schema references and exact template text reduce unnecessary output and provide executable documentation follow-ups.
- `Jarvis agent inspect TASK.yaml --json` explains effective defaults, overrides, runtime policies, planned paths and value sources without starting a scan or importing user operators.
- Validation diagnostics include manual argv and evidence scope. JSON failures stay parseable. Unknown manual paths and conflicting selectors fail with parent-page guidance.
- Explicit duplicate YAML keys fail with `JV2-YAML-003`, field path and both source locations; legal merge overrides remain supported.
- `Jarvis agent request --json` creates a fillable feedback draft for `Jarvis.HEP.AI@mail.grokbot.com`, with the exact `Resquest@Jarvis-HEP` marker. It never sends mail or executes supplied text.

This release also carries the previously committed 2.0.15 reliability and sampler integration changes, which have not yet been published as 2.0.15. See the historical ledger and release notes for that baseline.

## Publication status

2026-10-13: user authorized commit and publication. Source version is 2.0.16. Complete candidate verification passed. HEP source commit `50631bae3512c4c04021e1db36bcaec1d5d4d63d` is pushed on `feat/sampler-plugins-s1`. Public dependency lookup confirms Operas 1.4.2 is available and Jarvis-HEP-Sampler is absent (HTTP 404). HEP must not be published until the SDK is available from the public index.


## Verified candidate

Complete wheel/sdist and evidence are in `dist/2.0.16/`; a verified SDK candidate is bundled under `dependencies/sampler/`. Local verification installs all declared dependencies into a fresh venv, using that SDK wheel and the official PyPI index for every other dependency. This is candidate verification, not evidence that public-index-only installation already works.

Clean Git checkout fast gate: 1396 passed, 478 subtests passed, 3 skipped, 166 slow deselected. Complete installed runtime/project reload smoke, 24 agent CLI checks and 22 reliability regression tests passed. Wheel has 282 source files matching the release commit; sdist has 409 matching files, including the new installed smoke. Both exclude Diver and its package-data edits. Detailed evidence: [verification](release-2.0.16-verification.json).

The required SDK is still absent from PyPI. Prepare a draft GitHub release while retaining the feature branch; publish SDK 0.1.0, fast-forward the verified branch into master, then publish HEP through the existing OIDC workflow and verify public availability. The SDK already has `publish.yml`; creating its GitHub repository/remote and configuring its PyPI pending Trusted Publisher are the remaining setup steps.
