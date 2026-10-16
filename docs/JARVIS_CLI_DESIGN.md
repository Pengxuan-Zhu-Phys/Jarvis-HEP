# Jarvis 生态 CLI 设计规范

2026-10-17 定稿。用户决定：Jarvis 生态的所有命令行工具使用同一种 Jarvis CLI 风格，`Jarvis` 与各项目自己的命令（如 `jsampler`）在帮助、错误、退出码和 JSON 上保持一致。本文件是唯一来源；各项目的设计文档引用本文件，不复制规则。§13 保留 2026-10-17 修复前的实测差距，§14 记录其后落实状态。

## 1. 适用范围

| 命令 | 项目（PyPI 发行名） | 经 `Jarvis` 转发 | 当前渲染依赖 |
| --- | --- | --- | --- |
| `Jarvis` | Jarvis-HEP | — | typer、click、rich |
| `jsampler` | Jarvis-HEP-Sampler | `Jarvis sampler` | rich |
| `jopera` | Jarvis-Operas | `Jarvis operas` | typer、click、rich |
| `jplot` | JarvisPLOT | `Jarvis plot` | rich |
| `jportal` | Jarvis-HEP-Portal | `Jarvis portal` | rich |
| `jarvis-agent` | jarvis-agent | — | argparse 原样（MVP） |

发布脚本不是用户 CLI，不在本规范内。

## 2. 命令结构

- 形式为 `<prog> <command> [<subcommand>] [参数] [选项]`。命令名小写，多个词用连字符，如 `gen-plot-yaml`。
- `-h`/`--help` 与 `--json` 在任何位置都有效，子命令前后均可；其他选项属于各自的命令。
- 保留文件简写：`Jarvis task.yaml`、`jplot <file>`、`jportal <file>`。第一个参数既不是命令、也不是已存在的文件时：
  - 看起来像路径（含 `/` 或有扩展名）：报"文件不存在"；
  - 否则报"未知命令"，给出相近命令和 `<prog> -h`。

  两种情况都 exit 2。
- 旧选项保留，在帮助中单列 `Legacy options` 面板，说明对应的新命令。
- 删除类操作需要确认；非交互且未给 `--yes` 时拒绝执行并 exit 2，如 `jsampler remove`。

## 3. 帮助

版式采用 Typer 的 rich 帮助格式（`typer.rich_utils`），以 `Jarvis run -h` 为准：

```text
 Usage: Jarvis run [OPTIONS] <task_yaml>

╭─ Arguments ──────────────────────────────────────────────────────────────╮
│ task_yaml                 <str>   Path to scan task YAML [required]      │
╰──────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────╮
│ --resume                          Resume from checkpoint                 │
│ --help                    -h      Show this message and exit.            │
╰──────────────────────────────────────────────────────────────────────────╯
```

- 帮助输出到 stdout，exit 0。每一级命令都有帮助。
- 首行为 `Usage: <完整命令> [OPTIONS] ...`，不放进外框面板；有说明时紧跟其后，然后依次是 Arguments、Options 面板。
- 有子命令时，命令按分组面板列出；根命令在说明之后加一行 `Command help: <prog> COMMAND -h`。
- 参数占位写成 `<str>`、`<int>`、`<float>`；有固定取值时写 `<a|b>`；也可以用语义名，如 `<task_yaml>`。必填项标 `[required]`。
- 说明文字句首大写。`--help` 一律写 `Show this message and exit.`。
- 面板标题固定用英文：Arguments、Options、Commands（或项目自定义的命令分组名）、Legacy options。
- 宽度跟随终端；输出不是终端时不带颜色（§5）。
- 帮助和版本信息不写依赖的版本约束；版本要求只写在 `pyproject.toml`。
- 实现方式：已依赖 typer 的项目（HEP、Operas）用 `typer.rich_utils`；只依赖 rich 的项目（SDK、PLOT、Portal）不为帮助新增依赖，用 rich 按同一版式渲染。各项目自带一份渲染代码，与现在 HEP、Operas 各带一份相同。验收以输出为准（§12），不以实现方式为准。

## 4. 版本

- `-v`/`--version` 向 stdout 输出一行 `<PyPI 发行名> <版本>`，exit 0，如 `Jarvis-HEP-Sampler 1.0.3`。
- HEP 的 logo 横幅只在 stdout 是终端时显示；管道或重定向时只输出这一行。
- 与 `--json` 同时使用时（顺序不限）向 stdout 输出一个 JSON 对象，exit 0，不输出 logo：

  ```json
  {"ok": true, "command": "version", "name": "Jarvis-HEP-Sampler", "version": "1.0.3"}
  ```

  `name` 为 PyPI 发行名。转发命令（如 `Jarvis sampler --version --json`）输出子项目的发行名与版本。

## 5. 输出流与颜色

- 命令结果写 stdout；错误、警告、进度写 stderr。
- `--json` 时 stdout 只有一个 JSON 对象，出错时也是；日志不得写入 stdout。
- 只有对应的流是终端、且 `NO_COLOR` 未设置或为空时才输出颜色。`--json` 永不带颜色。不直接写 ANSI 转义码，由 rich 按终端决定。

## 6. 错误

- 文本模式在 stderr 输出一个红色圆角面板，标题为程序名（如 `jsampler`、`Jarvis sampler`）。内容依次为：
  1. 出了什么错，带上相关的路径、名字或值；
  2. 怎么修，给一条能直接运行的命令，如 `Jarvis sampler remove grid --yes`；
  3. 可选：相近的候选（did you mean）。
- 不打印 Python traceback；`--debug` 时附上 traceback。
- 未知命令、未知选项、参数值非法也用同一面板，不输出 argparse 原样的 `usage: … error: …`。
- 命令自身的结构化报告不受此条约束，如 `Jarvis validate` 的诊断列表（code、YAML 路径、建议、示例）。它们属于命令的结果。

## 7. 退出码

| 码 | 含义 | 例 |
| --- | --- | --- |
| 0 | 成功 | |
| 1 | 工作已经执行，但结果是失败 | 运行错误或样本失败；`jsampler check` 有检查项不通过；calculator 执行失败 |
| 2 | 用法或输入错误，工作没有开始 | 未知命令或选项；参数值非法（如 `--dimension 0`）；文件不存在；卡片配置错误（`Jarvis validate`）；名字未登记；登记文件损坏；与内置方法重名；非交互且缺 `--yes` |
| 130 | 被中断（Ctrl-C） | |

转发命令原样透传子命令的退出码。

## 8. JSON 模式

- 选项名统一为 `--json`。所有输出数据的命令都支持，子命令前后均可。
- stdout 只输出一个 JSON 对象，UTF-8，缩进 2，键名 snake_case。
- 每个 JSON 输出都必须带以下公共字段：
  - `ok`：布尔值；
  - `command`：完整命令路径，如 `"sampler info"`；
  - 失败时 `error`：`{"type": "...", "message": "...", "hint": "..." 或 null}`。
- 其余字段由各命令定义。现有字段只增不删，以免破坏已有消费方（`Jarvis agent`、`jplot man --json` 等）；已有的 `api_version`、`kind`、`diagnostics` 保留。
- JSON 字段或退出码的变化属于公共接口变化，须在各项目的变更说明中写明。

## 9. 日志与调试

- `--debug`（短选项 `-d`）在控制台显示 DEBUG 日志，出错时附 traceback。HEP `run` 与 `jopera` 已有此选项。
- 项目自己的日志选项保留，如 `jopera --log-mode`。

## 10. 转发命令

- `Jarvis sampler`、`Jarvis plot`、`Jarvis operas`、`Jarvis portal` 与 `jsampler`、`jplot`、`jopera`、`jportal` 对同样的参数保持子命令的解析、错误与退出码一致；帮助、错误和用法中显示 `Jarvis <command>`。`Jarvis sampler list`/`info` 通过 SDK 扩展点增加 HEP 内置方法说明与结果。
- 解析器只有一份，由子项目拥有。子项目的入口接受程序名参数（SDK 已有 `main(argv, prog=...)`）。HEP 不复制子项目的解析器，补充内容只在子项目提供的扩展点上添加，如 `Jarvis sampler list` 额外列出内置方法。
- 退出码和 JSON 原样透传。

## 11. 手册与 agent 发现

- 提供手册的 CLI 统一用 `man [topic] [--json]` 命名（`Jarvis man`、`jplot man`、`jportal man` 已是如此）。`jopera help` 改名为 `man`，`help` 保留为别名。
- 面向 coding agent 的生态入口是 `Jarvis agent --brief --json`。

## 12. 一致性测试

每个项目各加一份 CLI 一致性测试，至少覆盖：

1. 根命令和每个子命令的 `-h`：exit 0，首行为 `Usage: <prog> …`，`--help` 一行为 `Show this message and exit.`。
2. `--version` 在管道下输出一行 `<发行名> <版本>`；`--version --json` 输出 §4 的 JSON 对象。
3. 未知命令：exit 2，stderr 中没有 `Traceback`，也没有 argparse 原样的 `usage:`。
4. `--json` 下的一个错误：stdout 可解析为单个 JSON 对象，`ok` 为 false，并含 `error.type` 与 `error.message`。
5. `NO_COLOR=1` 且输出为管道时，输出中没有 `\x1b[`。
6. 仅 HEP：无宿主扩展的转发子命令 `-h` 与子项目命令的 `-h` 在替换程序名后相同；有扩展的 `Jarvis sampler list/info` 只增加 HEP 内置方法说明。

## 13. 修复前现状差距（2026-10-17 实测基线）

| 命令 | 与本规范的差距 |
| --- | --- |
| `Jarvis` | `--version --json` 被拒绝，输出 argparse 原样的 `usage:` 并 exit 2（未知选项同样输出原样 `usage:`）；未知命令被当作卡片路径，报 "task YAML not found" 并 exit 1；`validate` 文件不存在 exit 1（卡片配置错误为 2）；`-v` 在管道下、`NO_COLOR=1` 时仍输出彩色 logo；`Jarvis operas --help` 显示 `Usage: jopera`，`Jarvis plot --help` 显示 `jplot <file>`；`Jarvis portal --help` 比 `jportal --help` 多出格式行，其中 SLHA 重复一行；`Jarvis sampler` 帮助用的是 HEP 自带的一份 SDK 解析器拷贝（`cli/help.py` `build_sampler_help_parser`）；`cli/sampler.py` `_command_error` 输出纯文本 `[Jarvis] …` 并一律 exit 2 |
| `jsampler` | 帮助为 argparse 原样用法行（`usage: jsampler [-h] …`、`TARGET`、`show this help message and exit`）；登记损坏时 `list` exit 1、`remove` exit 2、`info` 抛 traceback；JSON 错误缺 `error.type`；`--version --json` 忽略 `--json`，输出纯文本；帮助中写了 "before 2" 的版本约束 |
| `jopera` | 用法行包在标题面板里；根帮助的 Arguments 面板重复列出 `--help`；未知命令输出 argparse 原样错误；`info` 找不到函数 exit 1，加 `--json` 仍输出彩色文本；`NO_COLOR` 无效；`--version` 输出 `jopera 1.5.0`，不是发行名，加 `--json` 仍为纯文本；手册命令名为 `help` |
| `jplot` | 用法为手写的多行列表；未知命令为一行纯文本（退出码 2 正确）；`data describe` 文件不存在 exit 1，且不加 `--json` 也输出 JSON；JSON 用 `kind`，没有 `command`；`--version --json` 输出纯文本 |
| `jportal` | 用法为手写的多行列表；`man -h` 中 `show this help message and exit` 为小写；卡片不存在输出纯文本 `[jportal] …` 并 exit 1；未知命令被当作文件；`--version` 输出 `Jarvis-Portal`，发行名为 Jarvis-HEP-Portal；`--version --json` 报 unrecognized arguments 并 exit 2 |
| `jarvis-agent` | 全部为 argparse 原样输出；没有 `--version` |

## 14. 落地顺序

2026-10-17 实现与发布更新：上表的 `jsampler` 差距及 `Jarvis sampler` 自带解析器、纯文本错误差距已按 N2/N4 修复。SDK 本地测试 234 passed；已公开的 SDK 1.0.3 经公共索引干净安装，已安装命令的帮助、JSON 版本与错误和 file-only check 均通过。HEP 当前源码的转发定向测试覆盖帮助和 JSON；已发布的 HEP 2.0.17 仍是 SDK-free 历史版本。`Jarvis` 根命令及 Operas/PLOT/Portal/agent 的差距仍待各自阶段处理。`Jarvis sampler list`/`info` 通过 SDK 解析器扩展点展示 HEP 内置方法，属于 §10 允许的补充内容。

1. **SDK**（HEP 实施路线 N2）：帮助版式、错误面板、退出码、JSON 错误字段。
2. **HEP**（N4）：删除 `cli/help.py` 中的 SDK 解析器拷贝，`Jarvis sampler` 的帮助与 agent 命令目录改用 SDK 的 `build_parser("Jarvis sampler")`；`_command_error` 改用同一错误面板。`Jarvis` 自身的其他差距另立步骤。
3. **Operas、PLOT、Portal**：各做一轮，按 §13 的差距清单和 §12 的测试。
4. **jarvis-agent**：MVP 阶段结束后再对齐。

各项目按自己的变更流程实施，例如 HEP 的变更门禁与 PR 开头约定。
