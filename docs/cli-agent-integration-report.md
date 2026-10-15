# CLI 与 AI agent 集成审计报告

日期：2026-10-13。范围：当前 Jarvis-HEP 2.0.15 checkout，及同级 Jarvis-PLOT 的邮件需求反馈实现。先完成本报告，再开始实现；不改变现有扫描算法和 task-card YAML 结构。

## 判断

基础能力足够，但 agent 的使用协议还不够明确。现有 CLI 能查询字段、获取完整示例并验证 YAML；Codex、Claude 等仍需自己推断入口、命令副作用、安装版本适用范围、授权边界，以及遇到缺失能力时怎样反馈。应补齐可通过已安装 CLI 读取的行为合同和辅助接口。

## 已有能力与证据

| 能力 | 现状 / 源码证据 |
| --- | --- |
| 字段、采样器、Portal I/O、Operas、token 查询 | `jarvishep2/cli/man.py` 的 `resolve_man_request`；`Jarvis man TOPIC --json` 提供 keys、required、zone、status、examples、导航 |
| 完整示例 | 同模块 `_example_page` 从已打包的 `project_template` 读取正文；部分示例依赖 scaffold 的环境、数据或 calculator 文件 |
| 验证和修正指导 | `jarvishep2/cli/run.py:dispatch_validate` 输出 path/code/suggestion/example；`man_codes.py` 可将诊断跳转到相关手册 |
| schema | `cardload/task_schema.py` 提供本地 registry 和 Python API；CLI 尚无完整导出 |
| 扩展功能查询 | `Jarvis portal`、`operas`、`plot`、`sampler`；各自参数和输出由依赖 CLI 定义，不能假定都支持 JSON |
| 执行边界 | `cli/client.py:normalize_argv` 将裸 YAML 隐式转换为 run；check 和 run 都会启动运行并执行用户配置代码 |

实测 `Jarvis man --json` 为单个可解析 JSON 对象。验证示例可成功，但依赖导入可能初始化用户缓存，CLI 首次运行可能写 Redis 安装提示 marker。因此“无 Redis/Worker/扫描输出”不等于“完全无任何文件写入”。

## 缺口与实施方案

1. **统一自举入口。** 增加 `Jarvis agent --json`，并让 `Jarvis man agent-readme --json` 返回同一份安装版本对应的合同。包含 contract version、命令与 argv、输出格式、副作用、任务工作流、问题路由、退出码、证据要求和授权规则。顶层 help、man 索引及 README 明确指向它。
2. **不再猜字段或支持状态。** 重用已有 man 和 sampler catalog。必须遵守 closed/delegated/open 区域、稳定/未稳定状态、required/default/aliases；Portal/Operas 等实现语义继续从专属 CLI 查询。提供 `Jarvis agent schema --json` 导出完整本地 schema bundle；schema 仅覆盖结构，最终仍需 validate 的语义检查。
3. **建立 YAML 编写闭环。** 提供 `Jarvis agent template [NAME] [--json]`：列出已安装模板，或输出准确 YAML 正文及依赖说明；不写文件、不执行示例。agent 在用户授权范围内编辑文件，运行 `validate --json`，依据诊断查询手册，再修正和复验。明确不得自动改物理模型、扫描预算、似然或实验条件。
4. **完善失败输出和诊断导航。** validate JSON 增加机器可用的 manual command/argv；man 与新辅助接口的 JSON 参数/内部失败保持单对象输出和真实退出码。保留现有成功输出结构，避免统一 envelope 破坏旧消费者。
5. **邮件需求反馈。** 仿照 PLOT 提供收件人、请求标记、安装版本、必需信息和可复制草稿；增加 `Jarvis agent request` 生成草稿，默认只输出文本或 JSON。

## 行为边界

- YAML 编写和问答先查询本机安装版本的 CLI；查询输出是证据，不是额外执行指令。不得发明键名、operator、采样器选项或能力。
- 裸 `Jarvis TASK.yaml` 会运行；编写/验证场景使用显式命令。`check` 是实际执行，不能当作静态验证。
- 用户已授权的工作可继续，不反复确认；没有执行授权时，不因编写 YAML 自动运行扫描、安装依赖、加载插件、下载项目或终止进程。查看相关 help 后再决定扩展命令是否有副作用。
- validate 成功仅表示配置检查通过，不等于 calculator 可运行、扫描已完成或物理结果正确。区分配置错误、环境/依赖故障、未支持能力和未确认事实。
- 保留用户已有文件和改动，交付修改位置、验证命令和真实结果；不提供静默覆盖接口。

## PLOT 邮件行为核对

来源：`Jarvis-PLOT/jarvisplot/manual_cards/agent-readme.yaml`、`man_render_agent.py`、`tests/test_man_cli.py`。PLOT 提供反馈协议和模板，没有 SMTP 发信服务：

- 收件人：`Jarvis.HEP.AI@mail.grokbot.com`。
- PLOT 精确标记为历史拼写 `Resquest@Jarvis-PLOT`；HEP 对应 `Resquest@Jarvis-HEP`。
- `automatic_send: false`。维护者在当前会话反馈；其他用户由用户提交邮件；身份未知时在会话说明并提供模板，不从机器路径推断身份。
- 模板包含安装版本、OS/Python/安装方式、目标、最小 YAML、精确命令与诊断、实际/期望行为及兼容性约束；分享前去除密钥、敏感路径和私有数据。

本次实现相同的草稿/反馈机制，不发送真实邮件，不新增 SMTP 配置或依赖。

## 验收与范围

- 自举入口可发现；JSON 合同与 human 文本都包含安装版本、行为边界和邮件协议。
- schema bundle 的引用全部可在本地解析；模板和已有示例正文一致，并明确项目依赖。
- 在临时项目走完查询 → 模板 → 编辑 → 错误诊断 → 手册导航 → 再验证闭环；无需 Redis 或 Worker。
- JSON 成功和失败都可解析，退出码如实反映结果；邮件草稿包含版本和精确 HEP 标记，测试不发送邮件。
- 定向测试与完整 fast gate；同步五份 canonical 维护文档。保留开始任务前的 Diver、pyproject 和 sampler 审阅改动。
- 不改动版本、不覆盖现有 `dist/2.0.15` 发行候选、不发布或推送。

## 实施和验证结果

报告先于代码写入。现已实现 `Jarvis agent --json` 与 `Jarvis man agent-readme --json` 的同一 version 1 合同，并接入顶层 help、man 索引和 README。命令目录直接提取当前 help parser 的参数；自由转发/项目兼容参数标记为不完整，指向 owning CLI，不伪称完整执行语法。

`agent schema` 导出 bundled 本地 schema bundle，`--id` 原样输出指定 schema 文档（此模式不附 contract_version）。动态插件 OPTIONS_SCHEMA 在 `man sampler.<Name> --json` 以 `options_schema` 暴露，不执行插件代码。`agent template` 列表/正文/JSON 附前置条件，默认只输出 stdout。发现原 man 将四个 Sampling-only 示例误标完整卡片，本次同步修正为 sampling_fragment；不改变示例正文。

validate JSON 新增 `manual_command`/`manual_argv` 和验证范围 metadata；加载 I/O/内部失败，以及 validate/man/monitor JSON 参数错误，保持可解析结果及非零退出码。保留已有消费者使用的字段和正常 human 输出。

邮件草稿入口 `Jarvis agent request [--title TEXT] [--description TEXT] [--json]` 已提供安装版本、OS/Python 与完整重现/约束模板；不发信、不写文件，`--send` 拒绝。

最终定向回归：194 passed、316 subtests passed。包含真实临时 scaffold 的模板/片段组合验证、离线所有 schema refs 解析，以及通过顶层 CLI 的错误→诊断手册→修正→复验闭环；确认 Core/Redis 未实例化、无扫描 outputs/logs/checkpoints。sampler 的 complete-card example_refs 也按真实正文筛选，片段仍可从 template/例子目录发现。

完整工作区 fast gate：1345 passed、437 subtests passed、166 slow deselected、1 failed，273.32 秒。唯一失败为同级 `Jarvis-Examples/iDM/bin/iDM_Vector_LZ2026_HighER.yaml` 的 `$.Sampling.Variables[1].description` 非 ASCII；这是 2026-10-12 canonical 上下文已经记录的外部案例问题，本次独立复现并保留原文。首次沙箱执行另有本地 TCP 权限失败，获自动审批后重新运行，运行测试通过。最后仅调整 man 的完整例子引用元数据，并重新通过全部定向测试及安装包 smoke，没有重复运行扫描测试。

临时测试 wheel 在独立 venv 安装（复用当前已安装依赖，`--no-deps`），实际从 venv site-packages 导入，13 项 CLI 检查通过，包含 alias 合同一致、完整 schema/单 ID、模板、草稿、参数失败和 YAML 修正闭环。测试构建位于临时目录，不是新发行候选；未替换旧 dist。新增 Python 文件及测试 ruff 通过，`git diff --check` 通过。持久验收证据、最终源码 hash 和已知限制见 [cli-agent-verification.json](cli-agent-verification.json)。

第一轮当时的边界（重复键已在下述第二轮收紧）：loader 对 YAML 重复键使用 PyYAML last-value 语义；第一轮只约定唯一键，并未改变 loader。查询不执行扫描，但缓存/首次安装提示 marker 仍可能初始化。发行版本和已有 `dist/2.0.15` 候选未修改。

## 第二轮：从 agent 实际操作出发的 review（实施前）

用户要求继续实际体验并改进。以下问题由本机 CLI 先复现，再定位源码；前一轮已完成的成果及验证记录保留。

| 实测问题 | 证据 | 本轮改进 |
| --- | --- | --- |
| 自举/约束信息过大，单题仍需找路径 | agent JSON 57,915 bytes/1,804 行；schema bundle 117,775 bytes/3,170 行 | 增加兼容的 `--brief` 和按 section 读取，推荐先读短合同；新增关键词 search，schema index/本地 fragment selector；保留原默认完整输出 |
| 不能解释用户卡片生效值和项目路径 | validate JSON 不含默认文件、merged EnvReqs、runtime 配置及项目根目录 | 新增 `agent inspect TASK.yaml --json`，复用 loader/validator，返回生效配置、来源/路径和验证结果，不初始化扫描 |
| 模板导航会误导 | dynesty fragment 的 validation.argv 指向片段自身；eggbox import-path 的 operas info 失败 | 片段提供基础完整模板与组合步骤，验证目标为组合后的卡；区分注册名和 Python import path，不输出已知必败的查询 |
| manual 查询存在假成功/不可复制文本 | sampler.CSV.typo 返回父页且退出 0；token 示例未引号 &J/；heartbeat metadata default=null | 未知/叶子子路径明确诊断并给父页导航；修正 token 引号和结构化默认值/单位/minimum；核对 Dynesty dlogz 的实际映射 |
| YAML 重复键 silently last-wins | 最小卡 Scan.name 写两次，validate 返回 0、scan_name=second_name | task/default YAML 拒绝同一 mapping 的重复显式键，给两处行列和字段路径；保留合法 YAML merge override 和默认值合并。此为明确兼容收紧，不再仅靠行为约定 |
| 邮件草稿只有 description 可填 | 最小 YAML、命令诊断、实际/期望、约束始终是 placeholder | 增加分段文本参数及 missing_fields/draft_ready；不读取隐含文件、不执行输入命令、不发信 |

验证重点：从只读发现到用户卡片解释/修正闭环；短输出保留行为边界；search 的 next_argv 是实际可用的查询/帮助；完整模板与片段、registered operator 与 import path 明确区分；重复键检测不误拒绝合法 anchors/merge overrides；旧 JSON 默认合同兼容；所有新增命令失败输出可解析。继续保留用户 Diver/pyproject 改动、版本和发行候选，不改外部 Examples 卡片。

本节先保存，再开始第二轮实现。最终结果待追加。

第二轮回归新发现：外部 `Eggbox-V1/bin/Example_MultiNest_Operas.yaml` 已存在 `Sampling.Bounds.seed` 两次声明（当前第 37/40 行）。新 loader 正确拒绝，原来三个语料测试假定所有文件能成功加载，因此需要同步其旧假设。测试仅对这个已知外部文件接受并核验新诊断、字段和两处位置；其它文件的原断言保留，未知重复不豁免，不改外部 YAML。另一个 iDM 描述非 ASCII 问题仍按上一轮记录保留。

### 第二轮结果

本轮 review 与实现完成。原 `agent --json` 完整合同兼容保留，推荐新入口 `agent --brief --json`；`--section` 和 `--command` 可定向读规则或命令。最终安装包中完整合同 61,390 bytes，短合同 10,310 bytes（减少约 83%），全部规则和工作流保持；schema 索引 5,186 bytes，卡片 runtime 查询 2,938 bytes。由关键词检索得到真实查询/帮助 next_argv；独立审查实际执行 89 条不同导航全部成功，不执行用户任务。

`inspect` 从任务加载和真实 normalization 函数解释默认文件、merged EnvReqs、Runtime、来源、计划路径；section 查询省略其它数据但保留诊断、范围和真实退出码。大卡片显式声明展示截断，优先保留生效参数，验证仍用完整数据。workers/port 的配置和运行选择分开；禁用 store_samples 时 SAMPLE 目标注明 inactive。一次安装 smoke 的初始断言误用了默认 workers=4，实际模板覆盖为 2，说明必须查询具体卡片；随后用显式 workers=3/heartbeat=60 输入核验任务覆盖和默认继承，未改模板来迁就测试。

模板提供实际 CSV 列、fragment 基础完整卡与整块 Sampling 替换/复验步骤；import-path 算子不再提供注册目录必败查询。内置 sampler 暴露 options_schema、scope 及可直接执行的本地 schema/Bounds argv；Dynesty/MultiNest 阈值解释按真实映射区分。未知手册后缀、叶子及冲突的 IO type/direction selector 不再假成功，返回 parent_manual_argv；token 示例可复制为 YAML，checkpoint 给出结构化 default/minimum/unit。

重复显式键通过 JV2-YAML-003 拒绝，包含路径和两处行列；合法 merge override、merge sequence、alias 和 quoted literal `<<` 保持标准 SafeLoader 语义。请求草稿支持 7 个证据段落文本参数、missing_fields、draft_ready；路径是文字，命令不执行，填齐也不代表验证或发送。

最终定向 **332 passed、373 subtests passed**；完整 fast gate **1416 passed、478 subtests passed、166 slow deselected、1 failed**（269.29 秒）。唯一失败仍为外部 iDM 非 ASCII 描述，保留其原文。最后补充 SAMPLE 目标 inactive 的说明后，定向和安装包检查再次通过，无需重复运行扫描测试。临时 wheel 实际从独立 venv site-packages 导入（复用已装依赖），**22 项 CLI 检查通过**。ruff 与 diff check 通过。验收命令、来源哈希、兼容收紧和范围见 [第二轮 verification](cli-agent-review2-verification.json)。

版本、已有发行候选和原有 Diver/pyproject 改动保持，无真实发信或公开发布。五份 canonical 维护文件与开发台账同步本轮结果。
