# Jarvis-HEP 2.0.16

1. **Roadmap phase**：两轮 CLI agent 集成完成；2.0.16 提交、完整发行包验证和正式发布。
2. **Affected modules**：CLI agent/man/validate、cardload YAML loader、installed smoke、维护文档与版本声明。
3. **Public interfaces**：新增 agent 合同、search/inspect/schema/template/request；validate JSON 增量导航；未知手册后缀与显式重复 YAML 键现在报错。原默认完整 agent 合同保留。
4. **Main risks**：重复键旧卡须修正；schema/inspect 不证明用户算子或外部资源可运行；SDK 0.1.0 尚未上 PyPI，正式发布前必须先公开依赖。Diver 草稿不在发行范围。
5. **Verification**：两轮报告已保存；第二轮 332 tests/373 subtests、22 installed CLI checks 通过。工作区 fast gate 1416/478 通过，1 项已知外部 Examples 编码失败；本版本的干净源码与完整包验收进行中。

## Changes

- Agents can discover the installed CLI contract, command effects, authorization rules, question routes and YAML authoring workflow through `Jarvis agent --json` or `Jarvis man agent-readme --json`.
- Brief/section/command views, bounded keyword search, bundled schema references and exact template text reduce unnecessary output and provide executable documentation follow-ups.
- `Jarvis agent inspect TASK.yaml --json` explains effective defaults, overrides, runtime policies, planned paths and value sources without starting a scan or importing user operators.
- Validation diagnostics include manual argv and evidence scope. JSON failures stay parseable. Unknown manual paths and conflicting selectors fail with parent-page guidance.
- Explicit duplicate YAML keys fail with `JV2-YAML-003`, field path and both source locations; legal merge overrides remain supported.
- `Jarvis agent request --json` creates a fillable feedback draft for `Jarvis.HEP.AI@mail.grokbot.com`, with the exact `Resquest@Jarvis-HEP` marker. It never sends mail or executes supplied text.

This release also carries the previously committed 2.0.15 reliability and sampler integration changes, which have not yet been published as 2.0.15. See the historical ledger and release notes for that baseline.

## Publication status

2026-10-13: user authorized commit and publication. Source version is 2.0.16. Candidate verification is in progress. Public dependency lookup confirms Operas 1.4.2 is available and Jarvis-HEP-Sampler is absent (HTTP 404). HEP must not be published until the SDK is available from the public index.
