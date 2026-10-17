# Jarvis-HEP-Sampler — Technical Design

版本：**2.5**（2026-10-16，Australia/Adelaide）。用户合约是 `Jarvis sampler load <file>` 后按内置方法使用，外部方法统一由 `Jarvis sampler` 管理。

1. **Roadmap phase**：M1 可靠性修复与 M2 本机联合故障验收完成；执行台账位于 SDK `docs/development/DEVELOPMENT_LEDGER.md`。
2. **Affected modules**：SDK registry/metadata/session/output，HEP plugin identity/preflight/checkpoint、ResumeService、supervisor、core 与 FixedSet 提交去重，以及定向和独立 CLI 回归。
3. **Public interfaces**：CLI/YAML、插件 API v1、STREAM_PROTOCOL=1、registry format=1 不变；checkpoint/source 路径 helper 增加兼容的可选 `check=False`，SDK 私有待清理文件不进入插件合约。
4. **Main risks**：单入口不含第三方依赖/资源，用户函数无强制超时；远程 OS/Python、最低依赖组合和性能目标仍开放。check 的通用 runtime/log/summary 仍遵循原生共享布局。
5. **Verification performed**：SDK 全测试/完整构件，HEP 确定性回归、独立 CLI 管理/迁移、真实 SIGINT、联合故障、运行中登记和启动失败恢复及完整依赖安装；命令、结果、源码 hash 与例外见 [M1 验收](native-parity-m1-verification.json)和[M2 验收](native-parity-m2-verification.json)。

当前 M2 验收基线为 SDK `cf924cb` / 当时的包版本 0.1.0 与 HEP `569d3865` / 2.0.17 加 sampler 工作区修改；M1 的 `ef6c223a` / 2.0.16 结果保留为历史。SDK 1.0.1 起提供 `jsampler`。当前 HEP 源码要求 `Jarvis-HEP-Sampler>=1.0.3`，不设上限；API、STREAM_PROTOCOL 与 registry 仍是 1。SDK 1.0.3 已由维护者发布，当前依赖下限最初对应其 N1/N2 行为；本轮入口分层与共享渲染需要尚未发布的 HOST_CONTRACT_VERSION 2。源码版本号不变，HEP 明确校验契约，维护者后续发布新 SDK 后再同步新的发行下限。S0/S1、S2a 默认严格恢复与 S2b 源码打包已有实现；本轮修复 check 污染正式绑定、indexed 中断保存超前游标、恢复后重置游标，以及 prefix 外已落盘行重复计算和归档目标超计。S2c 内置推广、S3 generational 尚未交付。

**更正旧维护描述**：AGENTS 指定的五份 canonical 文档现已重建并维护 HEP 总体上下文、PRD、架构、代码地图和路线。本文维护 sampler 技术决定，SDK 开发台账维护本轮任务与验收；历史章节的旧版本/推送/发行状态以其记录日期理解，不作为当前公开发行事实。

---

## 1. 目标与对等契约

### 1.1 目标

用户不修改 Jarvis-HEP 源码，只通过外部载入定义自己的 sampler。载入后，这个 sampler 的使用方式和行为与内置 sampler 完全相同：task card 写法、校验、Redis 发布、`sampler.log`、monitor、checkpoint/resume、手册和打包都一样。

与内置 sampler 唯一的区别：用户可以重新 `load` 同名插件，新版本覆盖旧版本。

这是发布目标，不是当前 0.1.0 / S1 工作区已经全部具备的能力。首个对外版本支持 independent profile；恢复与打包必须分别通过 S2a、S2b。重新 load 影响后续新 scan，不能改变已经启动或待续跑 scan 的绑定源码。

责任边界：用户负责算法、依赖、初始化及科学正确性，可以不写 schema、`validate_options` 或观测 hook。HEP/SDK 仍负责选对源码、校验恢复身份、如实记录来源、正确发布和归档点、报告失败与清理运行资源。不能把这些接入层责任归入“用户插件有 bug”。

2026-10-16：首发只交付上述使用体验。内置 Sampler 没有的行为，外部 Sampler 也不要有。`Jarvis sampler load` 之后，validate、check、run、resume、man 与内置方法使用同一套命令和参数，不再新增任何 CLI。已有的 `Jarvis sampler init|check|load|list|info|remove` 只负责把文件登记成 Method，不再扩展。`--accept-plugin-change` 不进入首发，不实现；续跑使用该 scan 已绑定的源码。已绑定 scan 上更换 NumPy、Python 或 SDK 后继续运行，只警告并按用户当前环境接着跑，与内置 Sampler 相同。

### 1.2 对等契约

下表是 S1/S2 的验收标准，§9 的每一条测试都对应这里的一行。

| # | 用户可观察面 | 要求 |
| --- | --- | --- |
| P1 | task card | `Sampling.Method: <插件名>`，插件参数直接放在 `Sampling.Bounds`。没有任何额外字段 |
| P2 | `Jarvis validate` / `Jarvis man` | 未知方法的报错列出插件名；有 schema 时 Bounds 的报错和内置一样定位到 `$.Sampling.Bounds.<key>`；`Jarvis man sampler.<插件名>` 有自己的手册页。两条命令都不执行用户插件代码 |
| P3 | 采样执行 | 点经 Redis 发布给 Worker，backpressure、batch、selection、确定性 UUID 都走内置同一路径 |
| P4 | `sampler.log` | 标签 `Jarvis-HEP.Sampler.<插件名>`；S1–S11 的措辞、级别、表格格式与内置相同；没有 schema 或 hook 时也包含 card 显式提供的插件参数；插件自己应用的默认值由可选 hook 补充 |
| P5 | monitor | 显示进度、accepted、seed 和插件自报的指标，不出现 `CUSTOM METHOD` |
| P6 | run metadata | 记录实际执行的插件名、版本、源码 hash，以及 allowlist 内的 seed、预算；完整科学配置进入本地身份记录。续跑副本的来源不能被当前登记覆盖 |
| P7 | checkpoint / resume | 中断后续跑，坐标和 UUID 与不中断时逐位一致；改了 card 的科学配置时拒绝续跑 |
| P8 | 打包 | `Jarvis project pack` 带上 card 用到的插件 |
| P9 | 失败 | 插件出错时，S8 日志带插件自己的 traceback 和定位信息，scan 正常停止，不回退到其他 sampler |

### 1.3 非目标

- 不提供沙箱，也不对插件函数做强制超时。插件是用户自己的可信代码，由用户负责。
- 不承诺跨 Python/NumPy 版本逐位复现。只保证同一环境下的复现。
- 插件不能覆盖内置方法，也不能实现 CSV 这类不使用 `Variables` 的方法。
- 插件生成的点不会因为名字或标记就被当成后验样本（§8）。

---

## 2. 端到端示例

以用户方法 `MyLatinHypercube` 为例。官方 `LatinHypercube` 已是保留名，不能用于用户的 init/load；两者可以使用相同算法，但名称必须不同。以下日志与恢复流程描述目标体验，当前完成状态见 §10。

### 2.1 生成模板并修改

```text
$ Jarvis sampler init MyLatinHypercube
Created MyLatinHypercube.py from the independent sampler template.
Edit draw(), then run: Jarvis sampler load MyLatinHypercube.py
```

用户把模板改成：

```python
"""Latin hypercube design: one point per stratum in every variable."""

import numpy as np

PLUGIN = {
    "api_version": 1,
    "method": "MyLatinHypercube",
    "version": "0.1.0",
    "profile": "independent",
    "description": "Latin hypercube design with point_number strata per variable.",
}

OPTIONS_SCHEMA = {
    "type": "object",
    "properties": {
        "centered": {
            "type": "boolean",
            "description": "Put each point at the centre of its stratum.",
        },
    },
    "additionalProperties": False,
}


def prepare(setup):
    rng = setup.rng()
    perms = np.stack([rng.permutation(setup.point_number) for _ in range(setup.dimension)])
    return {"perms": perms, "centered": bool(setup.options.get("centered", False))}


def draw(context, rng, plan):
    strata = plan["perms"][:, context.candidate_index]
    offset = 0.5 if plan["centered"] else rng.random(context.dimension)
    return (strata + offset) / context.point_number


def settings_rows(setup, plan):
    return [("centered", plan["centered"])]
```

### 2.2 离线检查（可选）

```text
$ Jarvis sampler check MyLatinHypercube.py --dimension 2 --points 1000 --seed 7
MyLatinHypercube check ->
  prepare           ok   (0.4 ms)
  draw              1000 points, shape (2,), all in [0, 1]
  repeatability     same seed gives identical points
  options schema    valid
```

`check` 不需要 Redis，也不需要先 load。它只检查形状、范围和可重复性，不评价算法本身对不对。

### 2.3 载入

```text
$ Jarvis sampler load MyLatinHypercube.py
Loaded sampler MyLatinHypercube ->
  version          0.1.0
  profile          independent
  options schema   yes
  source           /home/alice/proj/bin/sampling/MyLatinHypercube.py
  snapshot         ~/.jarvis/samplers/user_sources/MyLatinHypercube-3f9a1c0d2e4b.py
```

载入之后，新开的终端也能直接用，不需要再 load。再次 load 同名插件会覆盖，输出里多一行 `replaced  0.1.0 (3f9a1c0d2e4b) -> 0.2.0 (8b7e…)`。

### 2.4 task card 与运行

```yaml
Sampling:
  Method: MyLatinHypercube
  Variables:
    - name: mass
      distribution: {type: Flat, parameters: {min: 100.0, max: 1000.0}}
    - name: coupling
      distribution: {type: Log, parameters: {min: 0.001, max: 1.0}}
  Bounds:
    seed: 7
    point_number: 10000
    centered: false
```

写法与 `Method: Random` 完全一样。`Jarvis validate` 和 `Jarvis run` 的行为与内置方法相同。`sampler.log` 摘录如下。第一条记录给出了完整的 V1 记录头，后面两条只列出正文：

```text

·•· Jarvis-HEP.Sampler.MyLatinHypercube
	-> 10-06 23:15:02.114 - [WARNING] >>>
Initializing the MyLatinHypercube Sampling

MyLatinHypercube Sampler Settings ->
Field        | Value
------------ | -----------------------------------------------------------------
method       | MyLatinHypercube
seed         | 7
variables    | 2 (mass, coupling)
workers      | 8
resume       | no
point number | 10000
plugin       | MyLatinHypercube 0.1.0 (3f9a1c0d2e4b)
centered     | False
...
MyLatinHypercube Sampler Summary ->
Field            | Value
---------------- | -------------------------------------------------------------
stop reason      | point_number reached
submitted        | 10000
completed        | 10000
failed           | 0
elapsed          | 2 h 13 min
candidates drawn | 10000
points submitted | 10000
```

monitor 的 sampler 面板显示 `PLUGIN 0.1.0 · ACCEPTED 4213 · SEED 7 · SAVED …` 和进度条。中断后用 `Jarvis run task.yaml --resume` 续跑，行为与 `Random` 相同（§7.4）。

---

## 3. 职责边界

```mermaid
flowchart LR
    T[Jarvis sampler init] -->|SDK 固定模板| F[用户的 .py]
    F -->|Jarvis sampler load| R[(~/.jarvis/samplers<br/>sources.json + 快照)]
    R --> C[HEP 分层 catalog]
    Y[task card] --> V[Jarvis validate]
    C --> V
    C --> P[run preflight: 导入快照 + prepare]
    P --> A[PluginSampler adapter]
    A -->|u| M[Variables + Mapper + selection]
    M --> Q[Sample → Redis → Worker]
    Q --> D[Archiver → DATABASE]
    A --> L[lifecycle log / telemetry / checkpoint]
```

**SDK（Jarvis-HEP-Sampler，import 名 `jarvis_hep_sampler`）负责插件格式本身：**

- 模板；
- `PLUGIN` 元数据规则和方法名规则；
- 从文件加载插件；
- `~/.jarvis/samplers/` 登记的读写；
- 随机流派生规则；
- 输出检查；
- 日志接口（只定义接口，后端由 HEP 注入）；
- independent profile 的唯一实现 `IndependentSession`：负责 options 检查、`validate_options`、`prepare`、按索引派生随机流、`draw_block` 块缓存和输出检查。HEP adapter 和离线驱动都调用它，所以 `check` 与真实 scan 产生同样的点；
- 离线参考驱动和一致性测试套件。

SDK 只依赖 Python 标准库、NumPy、jsonschema，不导入 HEP、Redis、Portal。

**HEP 负责运行时的一切：**

- `Jarvis sampler` 命令行；
- catalog 集成与 task card 校验；
- run preflight；
- u→物理量映射和 selection；
- Sample、UUID、Redis 调度和 backpressure；
- 归档、checkpoint/resume；
- lifecycle 日志、telemetry、monitor；
- 手册和打包。

**插件只负责生成单位立方体里的坐标**，不接触 Sample、UUID、Redis、Worker 或 DATABASE。

**依赖方向：** 当前 HEP 源码硬依赖已发布的 SDK（`Jarvis-HEP-Sampler>=1.0.3`，无上限），与 Portal 的处理方式一致；SDK 绝不导入 HEP。SDK 是单独的仓库，放在 `Jarvis-Workshop/Jarvis-HEP-Sampler/`，PyPI 名为 `Jarvis-HEP-Sampler`。

**命令行入口（更正）：** SDK 自带 `jsampler`，HEP 的 `Jarvis sampler ...` 转发其模板、检查、加载等命令，并扩展内置方法的 list/info 展示。旧文中“只有一个入口、SDK 不提供命令行”的说法不再符合已发布的 1.0.3。

---

## 4. 插件契约（SDK API v1）

### 4.1 文件约定

一个插件是一个 `.py` 文件，包含以下模块级对象：

| 名字 | 必需 | 说明 |
| --- | --- | --- |
| `PLUGIN` | 是 | 元数据 dict（§4.2），必须可以 JSON 序列化 |
| `OPTIONS_SCHEMA` | 否 | Bounds 中插件参数的 JSON Schema（draft 2020-12）。缺失或为 `None` 时不检查这些参数 |
| profile 规定的函数 | 是 | independent profile 见 §4.3，generational profile 见 §4.6 |
| `validate_options(options)` | 否 | 语义检查，抛出异常表示不合法。在 run preflight 和 `check` 中调用 |
| `settings_rows` / `summary_rows` / `stop_reason` / `status` | 否 | 可观测性 hook（§4.5） |

用户不需要构造任何 descriptor 对象。插件可以导入当前 Python 环境里装好的任何库，依赖缺失时会在 load、check 或运行时报错。

插件也可以有模块级状态，HEP 不禁止。但 §7 的复现和续跑保证只在一个前提下成立：`draw` 的结果只依赖它的参数 `(context, rng, plan)`。

v1 的载入单位是单个 Python 文件，不复制同目录模块、数据文件或虚拟环境，也不修改 `sys.path`。复杂算法可以封装在当前环境已安装的包中，由该文件导入；源码快照只能固定入口文件。发布文档必须说明依赖安装和 `__file__` 随快照位置改变的行为，并提供一个导入外部包的示例。资源打包扩展不作为 v1 的隐式承诺。

### 4.2 `PLUGIN` 元数据

| 键 | 必需 | 取值 |
| --- | --- | --- |
| `api_version` | 是 | `1`。SDK 遇到不认识的 major 版本时拒绝加载 |
| `method` | 是 | 方法名，必须匹配 `^[A-Z][A-Za-z0-9]{1,63}$`，与全部内置方法名相同 |
| `version` | 是 | 非空字符串。作者用它声明"算法语义是否改变"，属于硬身份（§7.3） |
| `profile` | 是 | `"independent"`（S1）或 `"generational"`（S3） |
| `description` | 否 | 一句话说明，用于 `Jarvis man` 和 `Jarvis sampler list` |
| `statistical` | 否 | 默认 `False`。为 `True` 时，mapper 对维度扩张和非线性重参数化给出与统计类内置方法相同的警告（JV2-MAP-050/051） |
| `block_size` | 否 | 正整数。声明了它就必须提供 `draw_block`，不再要求 `draw`（§4.3） |
| `state_version` | generational 必需 | 正整数，状态格式的版本 |
| `feedback_return` | 否，仅 generational | 与 `Sampling.FeedbackReturn` 同构（`mode` 为 `minimal` / `fields` / `all_flat`，另有 `fields`）。缺省为 `{uuid, logL}` |

方法名规则：

- 与内置方法同名（不分大小写）时拒绝加载。
- 与已登记的插件同名（拼写完全一致）时，覆盖旧版本。
- 只差大小写时报错，要求先 `remove` 旧插件。

### 4.3 independent profile

候选点或块可按索引重建，生成结果与 Worker 数量、batch 大小、完成顺序都无关。`independent` 指不依赖 Worker 的评估反馈，不要求各点在统计上相互独立；LHS 和低差异序列可以有相关性。

```python
def prepare(setup: SetupContext) -> Any: ...                          # 可选
def draw(context: DrawContext, rng: np.random.Generator, plan) -> ArrayLike: ...
def draw_block(context: BlockContext, rng: np.random.Generator, plan) -> ArrayLike: ...
```

`draw` 和 `draw_block` 二选一：声明了 `block_size` 就用 `draw_block`，否则用 `draw`。

**`prepare` 在每次 scan 开始时调用一次，resume 时再调用一次。** 返回的 `plan` 可以是任意 Python 对象，不写进 checkpoint，所以必须完全由 `setup` 决定。没有 `prepare` 时，`plan` 为 `None`。`plan` 的内存开销由用户负责：例如 10⁶ 个点、20 维的置换表约 160 MB，可以改用 `int32` 或按索引计算的置换来节省。

**`draw` 对每个候选索引 `i ∈ [0, point_number)` 调用一次。** 返回长度为 `dimension` 的一维数组。

**`draw_block` 对每个块 `b` 调用一次。** 返回 `(block_size, dimension)` 的数组，覆盖候选索引 `[b·block_size, (b+1)·block_size)`；超出 `point_number` 的行由 adapter 丢弃。块边界只取决于插件声明的 `block_size`，与运行时的 `batch_size` 无关，所以生成结果与 batch 无关。适合 Sobol 这类按段顺序生成更高效的序列。

三种 context 都是 frozen dataclass：

| 字段 | SetupContext | DrawContext | BlockContext |
| --- | --- | --- | --- |
| `seed`、`point_number`、`dimension` | ✓ | `point_number`、`dimension` | `point_number`、`dimension` |
| `variable_names: tuple[str, ...]` | ✓ | ✓ | ✓ |
| `options: Mapping`（只读） | ✓ | ✓ | ✓ |
| `log`（§4.5） | ✓ | ✓ | ✓ |
| `rng()` | ✓（prepare 专用的随机流） | — | — |
| `candidate_index` | — | ✓ | — |
| `block_index`、`block_size` | — | — | ✓ |

- `options` 是 `Bounds` 去掉 HEP 保留键（§5.3）后的只读视图：list 变成 tuple，dict 变成只读 Mapping。options 进入 scan 身份，所以只接受 JSON 数据（字符串键、None、bool、int、有限 float、str、list、mapping），日期等其他类型报 `PluginConfigError`。
- `variable_names` 的顺序就是 `Variables` 的书写顺序。
- schema 的 `default` 不会自动填入。插件在代码里用 `options.get(key, default)` 自己给默认值。

### 4.4 输出契约

adapter 在建立 Sample 之前检查每个点：

1. 必须是一维（`draw`）或二维（`draw_block`）的数组或序列，形状完全匹配，不自动 flatten；
2. 元素必须是实数，不接受 bool、complex、字符串或 object；
3. 每个值都是有限数，且 `0 ≤ u ≤ 1`。闭区间与内置 Grid 一致；mapper 本身已经把无界分布的 u 夹到 `[tiny, 1−eps]`（`sampling/variables.py:14-15, 103-133`）；
4. 经 Variables 和 Mapper 映射后，每个物理量都是有限数。

任何一条不满足，就抛出 `PluginOutputError`，信息包含方法名、candidate_index、违规的分量和值，scan 停止。HEP 不会 clip，也不会替用户重抽。

### 4.5 可观测性

**日志。** 三种 context 都带 `log`。在 HEP 内，它写到 `Jarvis-HEP.Sampler.<Method>`，格式与内置相同；在 `check` 和离线驱动里，它写到标准 `logging`。可用的方法：

| 方法 | 对应 logging-spec | 级别 |
| --- | --- | --- |
| `log.step(title, rows)` | S5，每代或每个阶段一张表 | INFO |
| `log.adjust(name, old, new, why)` | S7，参数被修改、忽略或夹紧 | WARNING |
| `log.saved(what, path)` | S11，`<What> saved to <path>` | WARNING |
| `log.info(msg)` / `log.warning(msg)` / `log.debug(msg)` | 其他 | 同名 |

S1/S2/S3/S6/S8/S9/S10 由 runtime 独占写出，插件只通过下面的 hook 提供内容。

**hook（全部可选）：**

```python
def settings_rows(setup, plan) -> list[tuple[str, object]]: ...   # 追加到 S2 设置表
def summary_rows(plan, stats) -> list[tuple[str, object]]: ...    # 追加到 S10 总结表
def stop_reason(plan, stats) -> str | None: ...                   # S10 的停止原因
def status(plan, stats) -> dict[str, float | int | str]: ...      # monitor 指标
```

- `stats` 是只读的计数：`candidates`、`accepted`、`point_number`。
- `status` 最多 8 个标量键，字符串不超过 32 字符，整体受 telemetry 的 4 KiB 上限约束。超出的部分会被丢弃，不会让 telemetry 失败。
- hook 抛出的异常只记日志，不影响采样，这与 `lifecycle_log._hook_rows` 现有的处理一致。
- **S1 已实现的宿主约束（R03/R05）**：基础设置表由 adapter 记录显式 options，hook 补充默认/派生值，匹配的参数只显示一行；插件行最多 32 行，文本值最多 256 字符。SDK 每个 hook 每个 session 只告警一次，错误列表最多保留 16 条、每条 512 字符。`status` 在采样线程初始化及候选/等待安全位置计算，更新间隔至少 1 秒；租约/心跳线程只读快照，不执行用户代码或等待 hook。状态失败返回空插件指标，宿主计数继续有效；不要求插件线程安全，不引入强制超时。
- 停止原因必须来自真实状态。没有 hook 时，independent profile 的停止原因是 `point_number reached`；被中断时是 `interrupted`，由 runtime 判定。

### 4.6 generational profile（名称与元数据预留，运行语义在 S3 验收前冻结）

用于 MCMC 变体、优化器、自适应网格等需要反馈的算法。状态以纯函数方式传递：

```python
def init_state(setup) -> dict: ...                             # JSON 状态
def ask(state, context, rng) -> list[dict] | None: ...         # None 表示结束；[] 表示只记账、不发点的一代
def tell(state, context, results) -> dict: ...                 # 返回新状态
```

- **proposal** 的格式是 `{"u": [...], "tag": <小型 JSON，可选>}`。HEP 负责分配 id 和 UUID，插件不能自己生成 UUID。
- **results** 的顺序与 `ask` 返回的顺序一致，不按 Worker 完成顺序排列。每条结果是 `{"u", "tag", "status": "ok" | "failed", "logL": float | None, "fields": {...}}`，失败用 `status` 明确标出，不用 `None` 冒充 0。
- **状态**：`state` 必须是纯 JSON。HEP 在每一代结束（barrier）时把 `state` 和 generation 编号写进 checkpoint。不允许 pickle 外部对象；大型数值状态的格式留到 S3 再定。
- **随机流**：第 `g` 代的 rng 由 §7.1 派生，与运行时的并行度无关。
- **context**：包含 `generation`、`dimension`、`variable_names`、`options`、`log`。
- **HEP 侧实现**：adapter 复用 `FeedbackSampler` 的 barrier 循环（`feedback_sampler.py:438-469`）。`propose_generation` 对应 `ask`，`absorb_generation` 对应 `tell`，`feedback_return` 来自元数据。
- **预留**：async_chains profile（每条链异步推进）的名字已保留，v1 不开放。
- **S0 现状**：SDK 0.1 能校验 generational 的元数据（`state_version`、`feedback_return`），但加载时拒绝运行，提示要等后续版本。

以上是接口草案，不能作为已稳定的 ask/tell API。S3 必须先定义：`ask` 使用随机数后的新状态如何交给 `tell`；selection 拒绝、计算失败、缺失 score 如何逐 proposal 回传；空代的预算与停止规则；pending proposal、generation RNG 和 barrier 的 checkpoint 原子边界；`feedback_return` 与 card 的优先级。再用有拒绝、失败和中断的 toy 驱动验证。现有元数据校验不能证明任意 MCMC/优化器可接入。

### 4.7 错误类型

| 类型 | 触发场合 |
| --- | --- |
| `PluginLoadError` | 导入失败、`PLUGIN` 不合法、缺少必需函数、方法名冲突（load、check、preflight） |
| `PluginConfigError` | options 不符合 schema、`validate_options` 抛异常、`prepare` 抛异常 |
| `PluginOutputError` | `draw` / `draw_block` 抛出异常，或输出违反 §4.4 |
| `PluginResumeMismatch` | 续跑时身份不一致（§7.4） |
| `RegistryError` | `sources.json` 无法读取或格式不对、快照缺失或 hash 不符（HEP 侧报 JV2-PLG-002） |

所有错误都带方法名、插件快照路径、YAML 路径（例如 `$.Sampling.Bounds.centered`）、候选索引（如果适用），以及插件自己的 traceback。插件出错时 scan 停止，不回退到其他 sampler。

---

## 5. 登记、命名与校验

### 5.1 `Jarvis sampler` 命令

| 命令 | 行为 |
| --- | --- |
| `init <Method> [--output PATH]` | 用 SDK 的 independent 固定模板生成 `<Method>.py`。检查方法名格式，检查是否与内置方法同名；目标文件已存在时拒绝。模板见 §5.5 |
| `check <file> [--dimension D \| --card TASK.yaml] [--options OPTIONS.yaml] [--points N] [--seed S]` | 只接受文件路径。第一项检查方法名不与内置方法重名（§5.2 保留名），重名 exit 2，后续项不执行；算法检查失败 exit 1。然后离线执行参数/输出检查，重新导入后比较同 seed 流，并在新 Python 进程从中间候选重放尾段；覆盖 draw_block 的块边界。不需要 load；不能证明任意算法的科学正确性 |
| `load <file>` | 导入一次，先比对保留名，再按 §4.1–§4.2 检查结构；有 `OPTIONS_SCHEMA` 时运行 `check_schema`；然后写快照和登记。**不调用** `prepare` / `draw`。与 `jsampler load` 完全相同 |
| `list` | 内置方法与用户登记合并列出（§5.4 注册表报告）：方法、来源、版本、profile、源文件路径、状态（`ok` / `source changed since load` / `source missing` / `name-conflict` / `snapshot-broken`）。载入时间可由 `info` 查询；重名条目照常列出，不隐藏 |
| `info <Method>` | 显示登记的全部字段和插件提供了哪些 hook；重名或快照有问题时给出原因和修复命令 |
| `remove <Method> [--yes]` | 删除外部登记并清理其不再引用的托管快照；用户的源文件与 HEP 内置方法不受影响；按登记中的确切名字匹配，不做保留名检查，与内置方法重名的外部条目可以正常删除 |

**`check` 的 S1 语义（R06/R07 已验收）**：`--card` 与 `--dimension` 互斥；提供 card 时沿用变量顺序、seed、options 和 `point_number`，没有 card 时默认预算为 100。显式 `--points` 改变本次检查的总预算，输出显示有效配置与覆盖项；`--options` 替换插件 options。可重复性检查重新导入同一源码字节，并由 SDK 私有 worker 在新 Python 进程重放中间尾段。此处重启的是离线 session，不是 HEP 的 DATABASE/checkpoint 恢复验收；将来若增加部分检查上限，须另设参数并保留 context 原预算。

**CLI 风格。** 两个入口的帮助、错误面板、退出码与 JSON 遵循 [Jarvis 生态 CLI 设计规范](../JARVIS_CLI_DESIGN.md)；`Jarvis sampler` 直接使用 SDK 的解析器与渲染，通过 SDK `CatalogExtension` 数据回调为 `list`/`info` 补充 HEP 内置方法及报告；文本、JSON 均由 SDK 渲染，HEP 不再捕获 JSON 后重绘。内置 info 的 `conflicts` 字段列出大小写无关的旧冲突。

**`check` 只接受文件（2026-10-17 已在本地源码实现）。** check 用于用户写完 sampler 之后、load 之前，确认它能被 Jarvis 原生调用。已登记方法的快照就是 load 时那份源文件的字节，按名字 check 不提供新信息，还要求先校验快照 hash，因此删除按名字 check：参数不是文件时 exit 2，提示用 `Jarvis sampler info <Method>` 找到源文件路径再 check。`jsampler check` 与 `Jarvis sampler check` 行为相同。报告的第一项是方法名，单独列出；重名时后续项标为未执行。

### 5.2 用户登记 `~/.jarvis/samplers/`

**位置。** 默认目录是 `~/.jarvis/samplers/`，可以用环境变量 `JARVIS_SAMPLERS_HOME` 覆盖（测试和 CI 用）。`~/.jarvis/` 已经是 HEP 的用户目录，`scan_refs.json` 和 `cache/` 都在这里。这套做法与 Jarvis-Operas 的 `~/.jarvis-operas/sources.json` + `user_sources/` 一致。

```text
~/.jarvis/samplers/
  sources.json
  user_sources/MyLatinHypercube-3f9a1c0d2e4b.py
  .lock
```

**`sources.json` 的格式：**

```json
{
  "version": 1,
  "entries": {
    "MyLatinHypercube": {
      "source_path": "/home/alice/proj/bin/sampling/MyLatinHypercube.py",
      "snapshot": "user_sources/MyLatinHypercube-3f9a1c0d2e4b.py",
      "source_sha256": "3f9a1c0d2e4b…",
      "loaded_at": "2026-10-06T23:10:00+10:30",
      "sdk_version": "0.1.0",
      "plugin": {"api_version": 1, "method": "MyLatinHypercube", "version": "0.1.0", "profile": "independent", "description": "…"},
      "options_schema": {"type": "object", "properties": {"centered": {"type": "boolean"}}},
      "hooks": ["prepare", "draw", "settings_rows"]
    }
  }
}
```

**写入规则：**

- 写之前用 `.lock` 加文件锁（`fcntl.flock`），把多个终端同时 load 串行化。
- 先写同目录的临时文件，再 `os.replace` 替换。
- 覆盖同名插件时，只清理明确退役且不再引用的旧快照，保护工作源文件。提交后清理失败不回报登记失败；私有 `.pending_cleanup.json` 记录快照相对路径和完整 hash，后续成功 load/remove 在同一锁内重试。内容变更者保留，不扫描目录删除未知文件。提交前失败遗留的无登记文件也保守保留。

**快照语义。** HEP 运行时只用快照，不用用户正在改的工作文件。改了文件但没有重新 load，修改不会生效；`Jarvis sampler list` 会把这种情况标成 `source changed since load`。

**一致性检查。** 读取登记时，如果快照文件不存在，或快照的 hash 与 `source_sha256` 不符，就报 JV2-PLG-002，提示重新 load。

**保留名（硬约束；入口分层为本轮未发布源码修复）。** 外部 sampler 不能与内置方法重名，大小写不敏感，没有绕过选项。保留名单只有一份，由 SDK 定义为 `RESERVED_METHOD_NAMES`（当前 16 个：AMMCMC、AdaptiveBridson、Bridson、CSV、DEMCMC、DRAM、Dynesty、EnsembleMCMC、Grid、LatinHypercube、MCMC、MultiNest、PTEnsemble、PTMCMC、Random、ToyMCMC），`check_method_name` 只检查语法，`check_external_method_name` 在 init/load/check 入口默认使用该名单。通用 metadata/loader 不执行外部保留名策略。因此 `init`、`load`、`check` 在 `jsampler` 和 `Jarvis sampler` 两个入口下拒绝同样的名字。HEP 的 `reserved_method_names()` 为 SDK 名单与 HEP catalog 的并集，并由 HEP 测试断言 HEP 内置方法都在 SDK 名单内。新增内置方法时先更新 SDK 名单并发布，HEP 再提高 SDK 版本下限（不加上限）。旧版 SDK 或手工编辑仍可能在登记中留下重名条目，由 §5.4 的注册表报告在运行时拦截。

**可移植性。** 登记是"每个用户、每台机器"一份。在另一台机器上运行同一张 card，需要先 load，或者从打包的项目里拿到插件（§6.7）。

SDK 的 `RESERVED_METHOD_NAMES` 对外部 `init`、`load`、`check` 生效。先前 metadata 层默认拒绝保留名再靠 `trusted_builtin_method` 豁免官方插件的设计已更正：豁免参数与路径尾缀判断全部删除。registry.load/check_plugin 通过 loader 可选 metadata_validator 在 options schema 与 hooks 之前执行入口检查；运行时通用加载不传策略。

### 5.3 Bounds：保留键与插件参数

| 键 | 归属 | 规则 |
| --- | --- | --- |
| `seed` | HEP，所有 profile | 非负整数，`0` 有效，不接受 bool。默认 `0`，与内置方法相同 |
| `point_number` | HEP，independent 专用 | 正整数，表示候选点的数量（被 selection 拒绝的点也计入），必需 |
| `max_worker` | HEP | 与 FixedSetSampler 相同，是运行时的 inflight 上限，不进入身份 |
| 其余所有键 | 插件 | 组成 `options` |

插件的 `OPTIONS_SCHEMA` 只描述自己的参数，不需要声明保留键。HEP 先检查保留键，再把剩下的键交给插件 schema。

### 5.4 校验与 preflight

**`Jarvis validate` 和 `Jarvis man`** 对用户插件只读登记元数据及必要的快照完整性信息，不导入用户插件：

- 方法不存在：JV2-MTH-003，可选方法列表里包含已登记的插件，并提示 `Jarvis sampler list` 和 `Jarvis sampler load <file>`。
- 保留键错误：JV2-PLG-003。
- 有 schema：用它检查插件参数，违反时报 JV2-PLG-004，路径为 `Sampling.Bounds.<key>`。
- 没有 schema：给一条 info 级的 JV2-PLG-001，说明插件参数原样传递、未检查。不是 warning，因为不写 schema 是用户的合法选择；`--strict` 也不会让它失败。
- 登记损坏或快照不符：JV2-PLG-002。

插件参数的检查直接用登记里的 `options_schema` 调 SDK 的 `schema_issues`，不注入 bundled schema 的 `lru_cache` registry（`cardload/task_schema.py:86-107`）。

schema 仍然完全可选。S1 已修正 man 的解析根，本地 `$defs` / `$ref` 的类型、说明和默认值与校验一致。SDK 禁止引用解析时自动访问网络：无法解析的本地引用或外部资源引用返回定位明确的离线诊断。需要复用 schema 时把定义放进同一文档；schema 默认值只供展示，不自动写入 options，插件自行应用默认值（R09 已验收）。

**注册表报告（2026-10-17 已在本地源码实现）。** `sampler_catalog.registry_report()` 在内存里合并 builtin、builtin-plugin 与 user 三层，不把内置方法写进 `sources.json`。每个 user 条目带一个状态：`ok`、与内置方法重名、快照缺失或被改（复用 SDK `SamplerRegistry.problems`）；`sources.json` 整体无法读取时报告记一个全局错误。报告只读文件，不导入插件代码。登记部分沿用 `_user_layer` 的路径与签名缓存；快照校验每次都做，因为快照被改不会改变 `sources.json` 的签名。

`Jarvis run` 与 `Jarvis check` 每次启动都生成这份报告，位置在 preflight 第 1 步之前，先于 Redis、Worker、Archiver 和 DATABASE 整理。按报告处理：

| 报告内容 | 处理 |
| --- | --- |
| 任一 user 条目与内置方法重名 | 退出（配置错误，exit 2），不论卡片用的是哪个方法。列出全部重名条目，提示 `Jarvis sampler remove <X> --yes` 后改名重新 load，`<X>` 为登记中的确切名字（如 `grid`） |
| 卡片所用插件的快照缺失或被改 | 退出，不导入快照 |
| 其他条目的快照问题 | 警告，写入日志，运行继续 |
| 登记整体无法读取 | 卡片用内置方法时警告并继续；用插件时退出 |

续跑使用 scan 绑定的源码（§7.4），重名检查照样执行。`Jarvis sampler list` / `info` 读取同一份报告。原先分散在 Core preflight、Distributor 和 `info` 的 `assert_no_reserved_name_collision` 拦截由这份报告取代。

**`Jarvis run` 的 preflight** 在恢复服务中先于 DATABASE 整理/恢复和旧游标导入；fresh 路径在 supervisor 中先于 scan 清理、HDF5 写入准备和运行资源启动。导入、参数及初始化阶段的问题会在 Redis、Worker、Archiver 启动之前暴露；draw 的错误在执行时报告。目标顺序为：

0. 生成注册表报告，按上表处理；
1. 确定这次用哪份源码（§7.4）；
2. 核对将执行的源码字节与绑定 hash，恢复时先完成身份与 checkpoint 检查；
3. 从同一份字节导入；
4. 调用 `validate_options`；
5. 调用 `prepare`。

内置方法的 preflight 什么都不做，行为不变。之后 `init_sampler_from_config()`（`:114`）直接使用 preflight 已经导入的模块和 `plan`。插件只在 control 进程里导入，Worker 永远不导入插件。

实现为 `core.preflight_sampler_plugin()`，调用 `sampling/plugin_sampler.preflight_plugin`。失败时在 `sampler.log` 写 `<Method> Sampler meets error when preparing the plug-in -> …`（带插件自己的 traceback），并标记为已记录，避免 core 重复写一遍。

### 5.5 固定模板（independent）

`Jarvis sampler init` 生成的文件如下，`<Method>` 会被替换成用户给的名字：

```python
"""<Method>: one-line description (shown by Jarvis man sampler.<Method>)."""

import numpy as np
from jarvis_hep_sampler import DrawContext, SetupContext

PLUGIN = {
    "api_version": 1,
    "method": "<Method>",
    "version": "0.1.0",
    "profile": "independent",
    "description": "One-line description.",
    # "statistical": False,
    # "block_size": 64,          # declare it to use draw_block() instead of draw()
}

# Optional. Delete it, or set it to None, to leave the Bounds options unchecked.
OPTIONS_SCHEMA = {
    "type": "object",
    "properties": {
        "upper": {
            "type": "number",
            "exclusiveMinimum": 0,
            "maximum": 1,
            "description": "Upper edge of the scanned unit cube.",
        },
    },
    "additionalProperties": False,
}


def prepare(setup: SetupContext):
    """Run once per scan, and again on resume. Return what draw() needs."""
    return {"upper": float(setup.options.get("upper", 1.0))}


def draw(context: DrawContext, rng: np.random.Generator, plan) -> np.ndarray:
    """Return one point: shape (context.dimension,), values in [0, 1]."""
    return rng.random(context.dimension) * plan["upper"]


# Optional hooks. Delete the ones you do not need.
# def validate_options(options): ...
# def settings_rows(setup, plan): return [("upper", plan["upper"])]
# def summary_rows(plan, stats): return []
# def stop_reason(plan, stats): return None
# def status(plan, stats): return {}
```

---

## 6. HEP 集成

### 6.1 分层 catalog

`sampler_catalog.py` 成为唯一的查询入口，由三层组成：

| 层 | 来源 | 能否覆盖 |
| --- | --- | --- |
| builtin | 现有的 `_SPECS` | 不能 |
| builtin-plugin | 随 HEP 发布、用插件格式写的方法（§6.6） | 不能 |
| user | `~/.jarvis/samplers/sources.json` | 重新 load 时覆盖 |

**恢复上下文（S2a）。** `scan_plugin_context` 用 ContextVar 暂时提供当前 scan 的元数据，退出后恢复原 catalog；不覆盖用户登记，也不改变其他方法。`SOURCE_SCAN` 可用于该次 card 校验、方法解析、日志和 run metadata；`registered_plugin_record` 专门查询实时登记。扫描源码只在 preflight 核对 hash 后导入。

**`SamplerSpec` 的变化。** 新增字段 `source`（取值 `builtin` / `builtin-plugin` / `user` / 恢复时的 `scan`）和属性 `is_plugin`。插件的登记内容不放进 spec，而是由 `plugin_record(method)` 单独返回一个 `PluginRecord`（版本、profile、description、`options_schema`、hooks、源文件、快照路径、hash）。user 层从登记读出，不执行插件代码；builtin-plugin 层在第一次需要时加载包内文件一次。插件的 spec 由元数据生成：

- `stateless = (profile == "independent")`；
- `statistical` 取自元数据；
- `needs_variables = True`；
- `nested = mcmc = False`。插件不会触发 nested 结果 CSV 导出，也不会触发 MCMC 诊断；
- `family = "plugin"`。

**user 层的读取。** 按 `sources.json` 的 mtime 和 size 缓存，只读 JSON，因此 catalog 仍然可以在 validate/man 阶段安全导入。与内置方法重名的条目不进入 catalog 的方法解析，记入注册表报告（§5.4），由 run/check 启动时退出、由 `list` 显示。

**查询函数。** `spec_for`、`names`、`available_methods`、`is_stateless`、`statistical_names`、`family_of` 等都改成在调用时合并三层；新增 `canonical_name`（大小写不敏感）、`is_plugin`、`is_statistical`、`plugin_record`、`reserved_method_names`、`user_registry_error`、`registry_report`（本地已实现）。`builtin_specs()` 和 `schema_ids()` 仍只覆盖 builtin 层：只有内置方法有 bundled schema 和 Distributor factory。

**去掉 import 期冻结的名字集合（实际改法）：**

| 位置 | 现状 | 改为 |
| --- | --- | --- |
| `runtime/core.py:879` | `method in STATELESS_METHODS` | 调用处不变；`distributor._StatelessMethodsView` 对未注册 factory 的名字改查 catalog |
| `cardload/contracts/mapper.py:118` | `STATISTICAL_METHODS` | 调用时 `is_statistical(method)`；`mapper.STATISTICAL_METHODS` 保留作兼容导出 |
| `cardload/contracts/methods.py` | `_METHODS_NEED_VARIABLES` 冻结集合 | `method_needs_variables` 调用时查 spec |
| `logging/toplevel.py:132-152` | 日志标签的规范拼写在 import 时生成 | 表内查不到的全小写/全大写段，调用 `canonical_name`（带 LRU 缓存） |
| `monitor/queues.py:23-33` | 硬编码的反馈族集合 | 内置名单之外查 spec：independent 插件为 `inactive`，generational 为 `shared` |

`runtime/worker_config.py:126-130` 已经调用 `is_stateless(method)`，catalog 分层之后自动适用于插件。

**测试隔离。** `tests/conftest.py` 在测试会话开始时把 `JARVIS_SAMPLERS_HOME` 指向临时目录，测试永远不读开发者自己的 `~/.jarvis/samplers`。

### 6.2 Distributor 与 adapter

**`Distributor.set_method`。** 内置方法的 factory 不变。没有 factory 的名字查 catalog：插件方法由 `create_plugin_sampler(spec)` 构造 `PluginSampler`；未知方法的报错提示 `Jarvis sampler load <file>`。`Distributor.available_methods()` 是注册的 factory 加上 catalog 里的插件方法。`register_builtin_samplers` 的漂移检查只覆盖 builtin 层。

**preflight 的交接。** `core.init_sampler_from_config` 在 `set_config` 之前调用 `sampler.attach_preflight(core._plugin_preflight)`，adapter 直接使用 preflight 已导入的模块、已执行的 `prepare` 和身份。没有 preflight 时（直接在 Python 里用、测试），`set_config` 自己调用 `preflight_plugin(config)`，只是不写 scan 级副本。

**`jarvishep2/sampling/plugin_sampler.py`（新文件）** 里的 `PluginSampler(FixedSetSampler)` 负责 independent profile：

- **基本属性**：`method` 取元数据里的方法名。`uuid_prefix = f"plugin:{method}:{hard_identity[:16]}"`，再交给已有的 `deterministic_sampler_uuid`。
- **`set_config`**：沿用 FixedSetSampler 的逻辑（变量、mapper、selection、seed、batch、inflight），再读 `point_number`，并组装只读的 `options`。
- **`propose_next`**：与 `RandomS.propose_next` 结构相同。差别只在 u 的来源：调用 SDK 的 `IndependentSession.point(i)`，由它完成 `draw` 或块缓存，以及 §4.4 的输出检查。adapter 之后依次执行映射、物理量有限检查、selection、accepted 编号和确定性 UUID。进度条用 `PermilleProgress`，标签是 `<method> candidates sampled`。
- **checkpoint**：只保存两个游标 `_index` 和 `_accepted_index`，外加身份。不需要保存随机数状态，因为每个点的 rng 都由索引派生（§7.1）。续跑时，`prepare` 由 preflight 重新计算；`advance_to_persisted_prefix` 沿用 FixedSetSampler 的重放逻辑。
- **日志 hook**：
  - `log_settings_rows` 返回预算、`plugin <method> <version> (<hash12>)`、selection 和 card 显式 options，再用可选 `settings_rows` 补充默认/派生值；重复参数合并，展示受 §4.5 的大小限制。
  - `log_summary_rows` 返回 FixedSetSampler 的行，再加上插件的 `summary_rows`。
  - `log_stop_reason` 优先用插件的 `stop_reason`，否则是 `point_number reached`。
  - 三个 hook 经 session 的同名方法转发；异常每个 hook/session 只记一次 WARNING。`status` 由采样线程更新缓存，telemetry 只读宿主快照。
- **日志后端**：`JarvisLogBackend` 把 SDK `log` 绑定到 `Jarvis-HEP.Sampler.<Method>` 和 `format_two_column_log`。
- **checkpoint 身份核对（S2a）**：导入前核对 identity、源码 hash、预算、seed、selection 与合法游标；旧 state 不能覆盖当前 session 的科学配置。身份不符抛 PluginResumeMismatch，不回退 fresh、不删除 checkpoint。
- **运行总结**：`[Scan Performance]` 的 sampler 一栏对插件写方法名（内置方法仍是类名）。

generational profile 的 adapter 放在同一个文件里，在 S3 实现（§4.6）。

### 6.3 校验与 issue code

- **`cardload/contracts/plugin.py`（新）**：`validate_method_sampling` 遇到插件方法时，在 Variables 检查之后改走这里，不再进入内置方法的分支。插件方法没有 bundled schema，`task_schema` 的 method schema 分派自然跳过它们。
- **`cardload/task_validation.py`**：未知方法报 **JV2-MTH-003**（JV2-MTH-002 是"没写 Method"），suggestion 提示 `Jarvis sampler load`，大小写写错时给出正确拼写，登记无法读取时附上原因。
- **新的 `info` 级别**：`Level` 增加 `"info"`，`ValidationReport.notes()` 返回它们。info 不计入 errors/warnings，不影响 `ok`，`--strict` 也不会把它升级为 error。`Jarvis validate` 成功时打印 notes；run 路径把 notes 写进 `core.log`（INFO）。
- **`man_codes.py`**：登记四个新 code，都指向 sampler 手册：
  - JV2-PLG-001（info）：插件没有 `OPTIONS_SCHEMA`，它的参数原样传递、不做检查；
  - JV2-PLG-002（error）：快照缺失或 hash 不符，或登记无法读取；
  - JV2-PLG-003（error）：HEP 保留键（`seed`、`point_number`、`max_worker`）不合法；
  - JV2-PLG-004（error）：插件参数违反 `OPTIONS_SCHEMA`，或不是 JSON 数据。未知参数报在它自己的路径上，并给出最接近的已知参数名（SDK `schema_issues`）。

### 6.4 日志、telemetry、monitor、run metadata

- **`sampling/lifecycle_log.py`**：不需要改。插件经 adapter 提供三个 hook，标签由 `sampler.method` 得到，即 `Jarvis-HEP.Sampler.<Method>`。
- **`sampling/telemetry.py`**：有 `plugin_status()` 的 sampler 走插件分支。
  - progress 为 `finite`（`_index` / `_maxp`）；
  - metrics 包含 `accepted`、`seed`；插件 `status()` 的结果放在 `metrics["plugin_status"]` 下，与 adapter 自己的字段分开；
  - payload schema 不变（仍为 1）。
- **`monitor/samples.py`**：按 run metadata 的 `family == "plugin"` 选择通用插件布局，不按方法名。
  - 前两行与 Random 相同：`CANDIDATES n / point_number` 和进度条；
  - 副标题为 `PLUGIN <version> · ACCEPTED n · SEED s · <插件指标，最多两个> · SAVED …`；
  - `CUSTOM METHOD` 只留给既不是内置、也不是插件的方法。
- **`runtime_metadata.py`**：插件 knobs 的 allowlist 为预算、seed、插件版本、源码 hash12 和完整 `source_sha256`。运行 supervisor 从当次 preflight 取得实际执行来源，重新 load 不会改写这次运行的版本/hash（R04 已验收）。插件 options 不进入 run metadata；完整科学身份持久化在 DATABASE 的独立绑定记录。

### 6.5 CLI 与手册

- **新文件 `jarvishep2/cli/sampler.py`**：实现 `dispatch_sampler`。在 `cli/client.py` 的意图路由中注册 `sampler`，在 `cli/help.py` 中加帮助，归到 "Extensions & projects" 分组，与 `project` 放在一起。
- **`cli/man.py`**：`Jarvis man sampler.<插件名>`（大小写不敏感）由登记内容生成手册页，不执行插件代码。内容包括：
  - 说明，取自 `description`；没有 `OPTIONS_SCHEMA` 时注明参数不做检查；
  - Capabilities：type plugin、profile、版本、来源（用户文件路径，或 "shipped with Jarvis-HEP"）；
  - 保留键 `seed`、`point_number`，以及由 `OPTIONS_SCHEMA` 生成的参数表（类型、description、default）；
  - 自动生成的 `Sampling` 示例（带上 schema 里有 default 的参数）；
  - 相关诊断 JV2-PLG-002/003/004。

  `Jarvis man sampler` 的索引按 (type, method) 排序，插件标为 plugin 族；插件 description 超过 50 字符时截断，保持索引一行一项。
- **`--accept-plugin-change`**：2026-10-16 决定不进入首发，不实现。续跑使用该 scan 已绑定的源码。

### 6.6 Dogfooding：官方的 `LatinHypercube`

与 §2.1 用户示例算法相同的官方 `LatinHypercube` 随 HEP 发布，文件位于 `jarvishep2/sampling/plugins/latin_hypercube.py`，在 catalog 的 builtin-plugin 层登记，名字受保护。它只使用公开的 SDK API，经由 `PluginSampler` 执行。它用于持续覆盖插件路径，但不能代替 user 层的加载、覆盖、无 schema、无 hook 和恢复故障测试。

`Random` 不迁移。它依赖全局 NumPy 随机状态，迁移会改变已有 scan 的轨迹。

### 6.7 打包

**S2b 已实现；验收见 §10.7。** `project_packager` 在打包 share、repro、full 三种 profile 时：

- 找出 card 用到的 user 层插件，把它们的快照放进包内的 `.jarvis-pack/samplers/<Method>-<完整 SHA-256>.py`，清单为 `.jarvis-pack/samplers.json`；
- repro/full 导出已有 scan 时，以 scan 绑定的源码为准，不取用户后来 load 的版本；包内清单记录完整 hash、插件/API 版本和已知环境信息，fetch 时核对完整性；
- `Jarvis project fetch` 或解包之后，列出所需插件。未登记或同名但 hash 不同时，打印包内来源与 `Jarvis sampler load <path>` 提示；不能因为本机已有同名就静默跳过。

同名不同 hash 在清单中分别保留；重新打包替换旧生成清单/快照，避免重复 tar 成员。fetch 在复制到目标目录前核对来源路径与完整 hash；损坏时失败且不留下项目目录。不会自动 load。登记属于用户本机的状态。share 仍可带已有结果；新采样使用新 scan 名，恢复既有 scan 则使用 `--resume`。

验收必须在空登记目录、无 sibling checkout 的干净环境中完成，并覆盖同名不同版本；提示安装插件导入的依赖。源码入口与依赖环境是两项不同的可移植性要求，单个 `.py` 快照不代表环境已打包。

---

## 7. 随机性、预算与恢复

### 7.1 随机流

所有随机流都从 `Bounds.seed` 派生，使用显式的 PCG64，不依赖 `default_rng` 将来的默认算法：

```python
STREAM_PROTOCOL = 1

def stream(seed: int, *key: int) -> np.random.Generator:
    seq = np.random.SeedSequence(seed, spawn_key=(STREAM_PROTOCOL, *key))
    return np.random.Generator(np.random.PCG64(seq))

# prepare:            stream(seed, 0)
# 第 i 个候选点（draw）: stream(seed, 1, i)
# 第 b 块（draw_block）: stream(seed, 2, b)
# 第 g 代（generational）: stream(seed, 3, g)
```

每个点的坐标只由 `(seed, 索引, plan)` 决定，与 Worker 数、batch、完成顺序都无关。代价是每个候选点都要构造一次 SeedSequence 和 PCG64，§9 的性能 gate 会测量这部分开销。

保证范围：在同一个 Python/NumPy/插件依赖环境下逐位复现；跨版本只记录差异，不做保证。

### 7.2 预算与两种索引

- `candidate_index` 统计所有候选点，包括被 selection 拒绝的；`sample_index`（accepted）只统计进入 Worker 的点。
- 候选预算用完就停止，不会为了凑够 accepted 数量而继续抽点。
- selection 全部拒绝时正常结束：candidates 等于 `point_number`，submitted 等于 0，日志里照常有 summary。
- 被 selection 拒绝的点没有 UUID，也不进入 Worker。selection 表达式本身出错属于配置错误，不会被当作拒绝。
- 物理计算失败的点照样占用 accepted 编号，归档为 Failed。

### 7.3 身份

**硬身份（`hard_identity`，identity_version=2）** 是以下内容 canonical JSON 的 SHA-256，已包含生成模式与块长：

- 方法名、插件 `version`、profile、`api_version`、`state_version`、`STREAM_PROTOCOL`；
- profile 的生成模式（`draw` / `draw_block`）及 `block_size`；后者决定块随机流，不能被当作仅影响性能的 batch 参数；
- `Variables`（按书写顺序）、`Mapper`、`selection`；
- `seed`、`point_number`、`options`。

身份记在两个地方：

- **DATABASE 目录**：`<database_dir>/sampler_identity.json`，内容包括硬身份 hash、它的全部输入字段、源码 hash、SDK/Python/NumPy 版本。全新 scan 在 Archiver 启动之前，由 control 进程原子写入。之所以用独立文件而不是 HDF5 属性，是因为 control 进程不能以追加模式打开 `samples.hdf5`，否则会和 Archiver 的写锁冲突（`_runtime_supervisor.py:76-90`）。
- **checkpoint payload**：`sampler_state.plugin` 保存硬身份 hash 与来源 hash；完整硬身份输入来自独立绑定，checkpoint 另有预算/seed/selection/游标。

**软身份** 是源码 hash（插件 `version` 不变的情况下）、SDK 版本、Python 版本、NumPy 版本。

UUID 前缀只取自硬身份。相同 UUID 只能证明编号相同，不能证明换过源码字节后的坐标相同。2026-10-16 起首发不提供 `--accept-plugin-change`。若以后另做源码切换，仍不得绕过硬身份，并须先核对全部已持久化点（包括 prefix 之外的点）的坐标与 accepted 对应；对不上就要求新 scan。源文件字节改变后不能继续宣称原源码的逐位复现保证。

### 7.4 续跑规则

每个新 formal scan 在插件参数检查和 prepare 成功后，preflight 将源码原子复制到 `<task_root>/checkpoints/<scan>/<Method>/plugin/<Method>-<完整 SHA-256>.py`，随后原子发布 `<task_result_dir>/DATABASE/sampler_identity.json`（binding_version=1）。记录绑定确切路径、完整 hash、完整硬身份字段、元数据、schema/hooks 和 SDK/Python/NumPy 版本。插件 prepare 失败不发布绑定；后续 Redis/Worker 等资源启动失败可能保留已经发布的绑定，不能笼统声称所有初始化失败均不留绑定。此时仍按绑定和严格 resume 规则处理，M2 已通过真实 Redis/Worker 启动失败后的清理与 `--resume` 重试，fresh 重跑明确提示使用 `--resume`；多份副本或 mtime 变化不影响来源选择。

**check 隔离（SAM-HEP-01）**：`Jarvis check` 的身份写入 `<result>/DATABASE/test/sampler_identity.json`，源码和 checkpoint 写入 `checkpoints/<scan>/<Method>/check/{plugin/,state.pkl}`。它是每次独立的 calculator smoke，允许重跑，不读取 formal binding/checkpoint，也不改变正式科学 DATABASE、源码和恢复 checkpoint。通用 runtime.json、日志、run_summary 沿用现有原生共享布局。

`--resume` 先只读绑定，使用该 scan 元数据校验 card，不依赖实时登记。随后执行：

1. **核对硬身份和 checkpoint。** 在导入插件源码、旧游标或启动写入资源前比较科学字段；checkpoint 还须匹配来源、预算、seed、selection 和游标范围。不一致列出字段并拒绝，不能用旧 state 覆盖当前配置。损坏 checkpoint 不删除、不降级为 fresh。
2. **身份缺失或损坏。** 已有插件 DATABASE/checkpoint 不得从当前登记猜测来源。旧 S1 格式不能自动“升级”成已验证身份；使用新 scan，或另行设计可验证的迁移。
3. **选择源码。** 依据绑定路径核对完整 hash 后，从同一字节导入。缺失时只可用 byte-identical 登记/内置快照修复，并在初始化成功后补回；不同 hash 或损坏副本拒绝。删除登记、覆盖版本和较新的孤立副本均不影响原 scan。
4. **环境变化。** SDK/Python/NumPy 差异写 WARNING，绑定不变；不承诺跨环境逐位相同。
5. **prefix 重放。** 重放不足或旧游标超过 durable prefix 不得静默继续；adapter 抛 PluginResumeMismatch。M2 四组真实联合故障矩阵和两组全 durable 拒绝尾段已验收，见 §10.9；科学配置、源码或 checkpoint 损坏仍严格拒绝。

**首发不做：** `--accept-plugin-change`（2026-10-16）。续跑始终使用绑定源码。SDK、Python、NumPy 变化只警告并继续，与内置 Sampler 相同。

官方 `LatinHypercube` 的续跑只依赖 scan 自己保存的源码副本、元数据、科学设置和完整 SHA-256；原路径不再被当成来源认证或豁免依据。当前安装包的绝对路径或字节变化不会令完好的副本失效。`sampler_identity.json` 与副本都是本地可编辑的项目状态，哈希校验不认证蓄意同时篡改它们的行为。

续跑的机制本身沿用现有实现：

- 只从 `accepted ≤ DATABASE 连续 prefix` 的安全快照恢复；
- 从快照出发，重放候选点和 selection，直到 accepted 等于 prefix，这样就能对齐候选游标；
- 外部 PluginSampler 及 FixedSetSampler 在提交边界跳过 prefix 与 pending sample_index 中已 durable 的终态行，Failed 与 Completed 相同；
- pending 索引由 DATABASE 重建为内存副本，不写 checkpoint，也不累计全历史 UUID；legacy UUID 与其它 sampler 路径保持各自现有实现；
- candidate/accepted 和 selection 仍照常推进，flush 返回实际新提交数，backpressure 与归档目标据此计算；
- 重放过程中不输出进度刷屏。

续跑是 at-least-once 语义：物理计算可能重复执行，但归档按 UUID 去重。

### 7.5 推广到内置方法（S2c，外部插件恢复闭环之后）

`sampler_identity.json` 和硬身份检查再推广到内置方法；逐方法列出科学键、运行时键及旧数据库兼容策略，不能用一个泛化规则视为全部验证完成。内置方法的硬身份包括方法名、Variables、Mapper、selection、seed，以及 Bounds 去掉运行时键后的部分。运行时键由 catalog 的 `identity_exclude` 声明，例如 `max_worker`。旧 DATABASE 没有身份文件时只能标为历史身份未知；若允许告警后建立新基线，不能宣称已验证先前配置。内置方法的 `uuid_prefix` 保持不变。此推广不阻塞外部 SDK 的首个可用版本。

### 7.6 中断与开销

- **中断（S1 已实现，R10）**：FixedSetSampler 在批次、候选点、prefix 重放和 backpressure 等待中检查中断；runtime 绑定自己的中断标记。Redis proc board 的 paused/stopping/draining 状态最多每 0.1 秒读取一次。PluginSampler 也在每次 draw/status 返回后检查控制状态，selection 连续拒绝不会绕过检查。单次用户函数仍是同步调用，不承诺强制超时；Ctrl+C 清理与协作停止的证据分别记录在 §10.6，严格恢复另属 S2a。
- **内存**：上界为 `O(max_inflight · dimension + block_size · dimension + 快照窗口)`，再加上插件的 `plan`。adapter 不保存候选点到 accepted 的终身映射。
- **重放时间**：续跑的重放时间是 `O(候选 prefix)`。selection 很稀疏时会比较长，但不能为了加速跳过 selection。

---

## 8. 科学语义

- independent profile 是扫描 API：生成待评估的点。用户生成的偏置分布、相关点、低差异序列，都不会自动成为后验样本。`statistical=True` 只影响 mapper 的警告，不让 HEP 把结果当作后验处理。
- generational profile 的算法正确性由插件作者负责，包括 MCMC 的接受率与 Hastings 修正、越界拒绝、优化目标的方向等。HEP 只保证以下几点：
  - 结果按 `ask` 的顺序交给 `tell`；
  - 失败有明确的状态；
  - 状态在 barrier 上保存；
  - 恢复后不跳点。
- 停止原因必须来自实际状态，不能根据 finished 标志推测"已收敛"。`recent-code-review.md` 的 F3 是现有 nested sampler 的反例。

---

## 9. 测试与验收

### 9.1 SDK（S0）

- `PLUGIN` 元数据校验和方法名规则；
- 模板生成，并能通过一致性测试；
- 输出契约：`(1,d)`、bool、complex、字符串、NaN、Inf、负数、大于 1 都要拒绝；0 和 1 要接受；
- 随机流：同一 seed 逐位相同；`draw` 与 `draw_block` 的块边界不随运行时 batch 变化；
- 登记：原子写、并发 load 时的锁、同名覆盖时删除旧快照、只差大小写时报错、`remove`、快照损坏时报错；
- `testing.check_plugin`（`Jarvis sampler check` 使用）与 `assert_conforms`（pytest 使用）；`run_offline`；`reference_stream`（HEP adapter 测试的参照序列）；
- SDK 内不出现任何 `jarvishep2` 导入。

### 9.2 HEP 对等验收（S1/S2）

| 对等项 | 测试 | 判定标准 |
| --- | --- | --- |
| P1/P2 | 用临时 `JARVIS_SAMPLERS_HOME` load 测试插件，再 validate 一张只写 `Method` 的 card | 通过；拼错的参数在有 schema 时报到 `$.Sampling.Bounds.<key>`，没有 schema 时只有 JV2-PLG-001 |
| P2 | 插件顶层代码写一个 marker 文件 | `validate`、`man`、`list` 都不产生 marker |
| P3 | seed 固定，batch 取 1/16/64，workers 取 1/4 | 按 sample_index 比较，u 和 UUID 逐位一致 |
| P3 | selection 拒绝前两个点和中间的点；selection 全部拒绝 | 候选和 accepted 序列与参考枚举一致；全部拒绝时在有限时间内结束，并写出 summary |
| P4 | 现有 `tests/test_sampler_logging_contract.py` 遍历 catalog，覆盖 builtin-plugin 的 LatinHypercube 和一个 user 插件 | S1/S2/S3/S6/S9/S10 全部通过，设置表包含插件参数 |
| P5/P6 | telemetry payload 快照；monitor 渲染快照 | 出现 plugin 布局和指标，没有 `CUSTOM METHOD` |
| P7 | 中断后续跑 vs 不中断；DATABASE 有空洞；部分 batch 失败；checkpoint 之后只剩 selection 拒绝的尾巴 | 坐标和 UUID 逐位一致；prefix 不越过空洞；每个点只归档一行；尾部不重跑物理计算 |
| P7 | 两次运行之间重新 load/remove 插件；副本缺失、损坏或有多个副本 | 根据身份绑定恢复，完整 hash 不符则拒绝；登记更新不影响旧 scan；改 options、Variables、生成模式或 block_size 时拒绝。首发不提供源码显式切换 |
| P8 | 打包一个使用 user 插件的项目 | 包里有快照；fetch 之后打印 load 提示 |
| P9 | `prepare` 抛异常；`draw` 返回 NaN 或越界 | 在 Redis 启动前，或在第一个违规点处停止；日志有插件 traceback 和定位信息 |
| Worker | 真实 spawn 的 Worker | Worker 进程不导入插件模块，物理结果正常 |

### 9.3 运行环境与性能

- **CI 矩阵**：Python 3.10–3.13，Linux 和 macOS；NumPy 1.24 与 2.x 的合法组合。逐位比较只在同一环境内进行。真实 Redis 和 spawn 的测试标记为 `slow`，不能用 fakeredis 单测代替。
- **性能 gate（本机证据已记录，20 μs 不是通过门槛）**：分别测 2、20、200 维，各 10⁵ 次 draw。adapter 的额外开销（随机流、检查、映射，不含插件自己的 `draw`）中位数 20 μs/点是待校准目标。2026-10-15 的 pyenv 3.13.13 / NumPy 2.5.4 / Apple M4 Pro 样本、1 ms / 100 ms / 1 s 物理吞吐，以及 4000 次 hook 报错后日志仍为 4 行，写在 SDK 台账「性能证据」。2 维和 20 维低于 20 μs，200 维高于它。远程矩阵不在这组证据里。

---

## 10. 路线与开发台账

### 10.1 阶段状态与退出条件

“已实现”指代码存在；“已验收”必须有对应测试结果；“已发布”必须有用户可安装的构件与版本。三者分别记录。当前不能把 S1 的路径覆盖概括为 P1–P9 全部通过。

| 阶段 | 当前状态 | 退出条件 |
| --- | --- | --- |
| **S0 SDK 初版** | 初版 `0e841bd` / `78a2cd2`，S1 补齐 `79972c7` 已提交；发布准备 `cf924cb`；未发布 | 146 项 SDK 测试通过；最终构件安装与 smoke 再次通过。远程矩阵仍未执行 |
| **S1 HEP 接入** | 收尾验收完成并提交 `0d85950`；未发布 | R03–R10 已关闭；TeX 第 38 章已接入。本次不修改/编译 TeX；全仓已有失败单列记录 |
| **S2a 外部插件恢复安全** | 默认严格恢复及本机联合故障验收完成，R01/R02及R13本机P7部分关闭；2026-10-16 起 `--accept-plugin-change` 不进入首发 | 真实SIGINT、配置/来源/checkpoint拒绝、空洞/混合失败/拒绝尾段、运行中登记管理与资源启动重试均有本机证据。续跑使用绑定源码 |
| **S2b 可移植性与发布** | 项目来源打包/fetch 已验收；完整构件验证见 §10.7；pyenv 3.13 性能证据已写入 SDK 台账；未公开发布 | P8、无 sibling/PYTHONPATH 的完整 wheel 安装与空/冲突登记运行；远程兼容矩阵留到开发结束后。20 μs 仍不是通过门槛。SDK 先发布，再发布依赖它的 HEP |
| **S2c 内置方法推广** | 未开始，后续独立范围 | 逐方法验收身份与历史 DATABASE 兼容；不把所有内置方法改造混入 SDK 首发关键路径 |
| **S3 generational** | 仅有名称/元数据预留，SDK 明确拒绝运行 | 先关闭 R11 的状态与反馈语义缺口，再交付 adapter、模板、离线驱动和 barrier 故障恢复测试 |

本机 P7 故障矩阵已补齐。pyenv 3.13.13 的性能证据已写入 SDK 台账，20 μs 仍是待校准目标；远程兼容矩阵留到开发结束后。2026-10-16 决定 `--accept-plugin-change` 不进入首发。默认严格绑定、项目来源打包和本机完整构件验收是本次落实范围。S2c、S3 后置。当前不宜优先增加更多算法模板。

### 10.2 审阅发现与待办（2026-10-08）

严重性用“高/中”表示，与对等契约 P1–P9 的编号无关。第二列保留 2026-10-08 审阅时的事实，第三列记录当前状态：R01–R10 的默认 independent 路径已关闭；R11、R12 的发布/跨平台部分仍开放；R13 本机 P7 联合故障部分已在 M2 补齐，不能据此推断跨环境恢复一致性。

| ID / 严重性 / 阶段 | 审阅时的事实与影响 | 当前状态、责任模块与验收 |
| --- | --- | --- |
| **R01 高 · S2a** 恢复会混合旧游标与新科学配置 | 原实现把 checkpoint `_maxp` 等导入新 session，block_size 不进入身份 | **已关闭（默认严格恢复）**。identity_version=2 包含生成模式/块长；metadata-only 比较先于源码导入；checkpoint 不得覆盖预算/seed/selection。回归覆盖科学字段、游标与真实 CLI 拒绝，见 §10.7 |
| **R02 高 · S2a** scan 源码没有稳定绑定 | 原实现按 mtime 取副本，失败 prepare 留副本；恢复入口依赖实时登记 | **已关闭（默认绑定）**。成功 prepare 后原子记录确切路径/full hash；read-only scan 元数据贯通解析与运行。真实 reload/remove、腐坏源码、不同科学配置和坏 checkpoint 覆盖；单测覆盖失败初始化、孤立较新副本、缺失仅同字节修复 |
| **R03 高 · S1 + SDK 补齐** 观测 hook 能阻塞租约 | 静态调用链：`_control_lease_loop → _ensure_managed_redis_alive → _publish_lease_proc_boards → sampler_status_json → plugin_status → session.status`。用户 hook 同步运行在租约线程，可能阻止下一次续租，并与 draw 并发访问 plan；SDK `_hook_failed` 每次追加错误和写日志，无上限 | **已关闭**。PluginSampler 在采样线程更新 status，telemetry/lease 只读缓存；SDK 错误限频、限量。HEP `test_slow_status_cannot_block_control_lease` 用事件阻塞 hook，实际 control loop 仍调用续租且不并发执行 hook；SDK 连续 1000 次变更错误消息只保留/告警一次 |
| **R04 中 · S1** provenance 记录了错误版本 | `runtime_metadata._plugin_knobs` 查当前登记。加载 v2 后恢复旧 scan：实际执行 v1，metadata 却写 v2/hash(v2) | **已关闭（运行来源）**。supervisor 传入 preflight 的实际 version/full hash/N/seed；HEP `test_metadata_uses_the_prepared_source_after_reload` 验证已登记 v2、实际执行 v1 时 metadata 仍为 v1/full hash。严格 scan 身份绑定与恢复拒绝规则继续归 R01/R02 |
| **R05 中 · S1 + SDK 补齐** 可选 hook 实际成为参数日志的前提 | 没有 `settings_rows` 时 `log_settings_rows` 只有预算和 plugin 行；显式 `upper=0.5` 消失。当前测试插件实现了 hook，掩盖这一缺口 | **已关闭（参数日志）**。adapter 总是记录显式 options，hook 只补充/合并展示；32 行与 256 字符文本限制。HEP 三种参数化用例覆盖无 schema/无 hook、有 schema/无 hook、hook 抛异常；完整科学配置持久绑定继续属 S2a |
| **R06 中 · S1** `check --card` 检查另一份预算 | `cli.sampler._run_check` 读取 card 后仍将 `args.points` 传给 checker；探针 card N=128，传入 N=100。LHS 等依赖总点数的算法检查结果与 run 不同 | **已关闭**。CLI 沿用 card 配置并打印覆盖项；HEP 测试验证原预算 128、显式覆盖 16 与 seed 3、变量顺序和 options。TeX LHS 示例按 card N=1000 执行真实 check/run |
| **R07 中 · SDK 补齐 / S1 gate** 可重复性检查存在假通过 | `testing.check_plugin` 的 repeat/resume 共用一次导入的 `LoadedPlugin`。模块导入时生成随机常数的插件通过全部 check，但两次独立导入的流不同。`reference_stream` 与 adapter 共用 session，只能验证接线，不能充当独立 RNG 正确性证据 | **已关闭（离线检查范围）**。repeat 重新导入源码，新 Python 进程重放中间尾段；SDK 测试能拒绝 import-time 随机量与进程 PID 依赖。SDK/HEP 均以解析 indexed 序列验证 draw_block 非整块尾段和中间开始，不以共享 session 自证科学正确性 |
| **R08 中 · S1 文档** 新用户教程第一步失败 | SDK README 与原设计示例使用受保护的 `LatinHypercube`；实际 `Jarvis sampler init LatinHypercube` 返回 exit 2 | **已关闭**。用户示例、SDK README/authoring/example 和 CLI help 统一为 MyLatinHypercube；主教程进入 Jarvis-Books TeX 第 38 章。逐条执行该章的源码与 card，在独立 CLI 进程完成加载、1000 点运行、覆盖和移除，源文件保留 |
| **R09 中 · S1 + SDK 补齐** schema 的展示与校验不一致 | 合法的局部 `$defs/$ref` 在 check 中通过；`cli.man._plugin_method_page` 调 `_key_entries` 未传 schema 根，手册中的 number/说明变成 any/空白。外部引用的离线策略未定 | **已关闭**。HEP man/example 传入 schema 根；SDK 引用解析不自动联网。回归覆盖本地引用的类型/说明/default、无 schema、无效本地与外部引用的离线诊断；schema 仍可省略 |
| **R10 中 · S1** 协作停止承诺没有实现证据 | `PluginSampler.propose_next` 在 selection 持续拒绝时一直循环；`FixedSetSampler._wait_for_inflight_room` 只 sleep/poll；都没有 §7.6 宣称的中断/scan-mode 检查 | **已关闭（协作停止）**。HEP 回归覆盖全拒绝候选、prefix 重放、inflight 无空位时的 interrupt/paused 状态；真实 Redis/spawn 的运行收到 SIGINT 后 exit 130、保存 checkpoint、清理子进程和 Redis 端口，见 §10.6。无强制用户函数超时；中断后的严格恢复仍属 S2a |
| **R11 中 · S3 前置** generational 过早声称接口冻结 | `ask(state, context, rng)` 只返回 proposals，随机 proposal 之后的状态归属不清；result 没有 selection 拒绝状态；空代上限、score 缺失、pending 状态恢复未定义。SDK 当前只校验元数据，拒绝载入执行 | **SDK 契约 + HEP FeedbackSampler adapter**：先完成 §4.6 列出的状态机与逐 proposal 终态；toy 优化/反馈算法跑乱序、失败、全拒绝、空代和 barrier 前后中断。通过后再冻结可运行接口，不在首发说明中声称已支持任意 MCMC |
| **R12 中 · S2b** 可分发性需要构件证据 | 历史 editable/sibling 安装不能证明用户安装；单入口不含其依赖环境 | **源码打包/fetch 已验收，构件结果见 §10.7；公开发布/远程矩阵仍开放**。清单保存 full hash 和环境，同名多版本保留；fetch 拒绝损坏来源而不导入源码。两仓统一完整 wheel/sdist 安装和可选显式上传；HEP API 依赖范围 `>1`。用户依赖仍由用户安装 |
| **R13 中 · S1/S2a gate** 测试措辞失真 | 历史完成后 resume 冒充中断验收，slow 参数解析会改变 gate 范围 | **默认严格恢复与本机 P7 联合故障证据已完成，远程矩阵仍开放**。真实 SIGINT→覆盖/移除→恢复逐坐标/UUID 核对；full gate 初跑 CLI 两个新增失败已修复并定向回归，最终结果如实见 §10.7 |

### 10.3 历史验证记录与边界

以下是前次开发记录，不是本次重新执行的结果：

- S0 初始记录：140 tests、ruff、wheel 干净 venv 安装；SDK 后续记录为 141 tests。仓库有 Linux/macOS、Python 3.10–3.13 及最低依赖 CI 配置，尚未取得本次审阅基线的矩阵执行结果。
- S1 记录：`tests/test_sampler_plugins.py` 25 tests、官方 LatinHypercube 的 logging 契约、真实 Redis 100 点 scan 及完成后 resume。它们证明主路径可用，不足以关闭 R01–R10 或完整 P7/P9 故障矩阵。
- 历史总 gate：1368 passed、8 failed、1 deselected。**非全绿；不再写“全量 fast gate 已通过”。** 历史失败包括 corpus ASCII、variable distributions 的 3 个卡配置用例、AdaptiveBridson 的 3 个用例、worker scaling 性能用例。部分原因曾归于 fixture、遗留 checkpoint 或机器负载，但缺少完整的逐项同环境对照，不能笼统认定八项全与 S1 无关。
- `tests/conftest.py` 的选择参数问题意味着带 `--deselect <nodeid>` 的命令可能包含 slow；应按实际收集范围归档，不能只看任务名称。后台通知的 exit 0 也不能替代 pytest 的测试摘要与退出码。
- 现有 `test-summary.json` 日期是 2026-10-05，记录另一次基线检查（1205 passed、2 failed、163 deselected）；不要当作 2026-10-08 的 S1 验收附件。

### 10.4 历史审阅探针（2026-10-08，修复前）

环境：本机 Python 3.13.13、当时 HEP 工作区、SDK `78a2cd2`；`JARVIS_SAMPLERS_HOME` 指向临时目录，未触碰用户真实登记。以下保留当次 10 项轻量探针的诊断结果；修复后的长期回归与实际结果见 §10.6，本文不另建 review 报告。

| 探针 | 实际结果 | 关联 |
| --- | --- | --- |
| 新配置导入旧 checkpoint state | 未抛错；adapter N=8 / session N=2，options 使用新值 | R01 |
| block_size 4 与 8 | 硬身份相同，`candidate_index=4` 的坐标不同 | R01 |
| 无 settings hook，传入 upper=0.5 | 设置表没有 upper | R05 |
| init LatinHypercube | exit 2：内置保留名 | R08 |
| import-time 随机常数插件 | check 返回 ok；两次重新导入的首点不同（本次为同进程重新导入探针，未冒充 spawn 测试） | R07 |
| check --card 参数路由 | card N=128；传给 checker N=100（mock 截获调用参数） | R06 |
| 本地 $defs/$ref | check 通过；man 的 upper 类型 any，说明为空 | R09 |
| 登记 v2，resume v1 副本 | 执行 v1；metadata 报 v2 | R04 |
| 同 scan 初始化失败后 resume | 选中了失败 prepare 留下的较新副本，再次失败 | R02 |
| 默认 fast marker 参数解析 | 普通 tests/ 保留 not slow；加 --deselect nodeid 后变为空；显式 -m 可保留 | R13 |

当次 R03/R10 的线程和停止风险来自静态调用链审阅，未模拟租约过期或跑真实中断；当次只更新本文，没有执行全量 pytest、Redis 集成、wheel 安装或性能基准。此段是修复前的历史记录。

### 10.5 修改入口

| 范围 | 主要文件 |
| --- | --- |
| SDK 修正 | `Jarvis-HEP-Sampler/src/jarvis_hep_sampler/{session,testing,options,metadata}.py`；`README.md`、`docs/development/PLUGIN_AUTHORING.md`、模板/示例与对应测试 |
| 加载与恢复 | `jarvishep2/sampling/plugin_sampler.py`、`sampler_catalog.py`、`runtime/{core,_runtime_supervisor,_resume_service}.py`、`cardload/contracts/plugin.py` |
| 用户入口与观测 | `jarvishep2/cli/{sampler,help,man}.py`、`runtime_metadata.py`、`sampling/{telemetry,fixed_set_sampler}.py`、`docs/logging-spec.md` |
| 打包与发布 | `jarvishep2/project_packager.py`、两仓 `pyproject.toml`/CI/release 文档 |
| 验收 | SDK tests、HEP `tests/test_sampler_plugins.py`、`tests/test_sampler_logging_contract.py`、`tests/conftest.py`；新增真实恢复故障用例 |

### 10.6 S1 收尾验收（2026-10-08–09，Australia/Adelaide）

**当时范围与环境（提交状态保留为历史）。** macOS arm64、Python 3.13.13、NumPy 2.5.0；HEP `0bbcf675` 加插件接入工作区，SDK `78a2cd2` 加本次补齐。改动尚未提交、未发布。pytest 使用隔离登记；初次受 sandbox 限制的 loopback 测试不计为通过，随后以允许本地 TCP/spawn 的执行环境重跑。SDK 在 `/private/tmp/jhs-s1-work` 验证后写回实际仓库，逐文件核对 src/tests 字节一致；真实教程去掉 PYTHONPATH，使用实际 editable SDK。

| 验证 | 结果与退出码 | 可检查证据 |
| --- | --- | --- |
| SDK 全测试与 ruff | **146 passed**，pytest exit 0；实际仓库 `ruff check --no-cache .` exit 0 | SDK `tests/test_session.py`、`test_options.py`、`test_template_and_testing.py` 等；包括重复 hook 错误、本地引用、重新导入/新进程与 indexed block 独立预期 |
| HEP 定向 gate | **56 passed / 24 subtests**，exit 0，15.63 秒；1 条既有除零警告 | `/private/tmp/jarvis-s1-focused.xml`；下面列出的 5 个测试文件 |
| HEP 全量 fast | **1241 passed、2 failed、163 deselected / 141 subtests**，exit 1，279.98 秒 | `/private/tmp/jarvis-s1-fast.log`、`/private/tmp/jarvis-s1-fast.xml`；两项失败逐项见下文 |
| TeX 教程完整用户流程 | init→edit→check→load→独立进程 validate/man/run→reload→info→remove 均 exit 0；**1000 submitted/completed/archived，0 failed**；重新登记版本为 0.2.0，remove 后工作源码仍存在 | `/private/tmp/test_s1_tutorial.py`、`/private/tmp/jarvis-s1-tutorial/` 下各命令日志、run_summary.json、DATABASE 与 sampler.log；真实 managed Redis、2 个 spawn Worker 和 Archiver |
| 实际 SIGINT 清理 | 预算 200000 的插件运行就绪后发送 SIGINT；**exit 130**，保存 state.pkl，剩余测试子进程 0，Redis 端口关闭，sampler.log 记录 interrupted | `/private/tmp/test_s1_interrupt.py`、`/private/tmp/jarvis-s1-interrupt/run.log` 与 scan 日志/checkpoint；此项没有继续 resume |
| TeX 接入与内容 | 第 38 章共 8 页，沿用手册的 terminal/YAML/代码样式，源码含完整 LHS 与 card；本次整本编译 exit 0、新章节无 overfull、无未解析引用；章节源码与编译快照一致 | `Jarvis-Books/TeX/Jarvis-HEP-V2 manual/chapters/38-jarvis-hep-sampler.tex`，已在 main.tex 的第 37 章之后 include |

HEP 定向命令（从 HEP 仓库执行）：

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=/private/tmp/jhs-s1-work/src \
python -m pytest -q -m 'not slow' \
  tests/test_sampler_plugins.py tests/test_sampler_telemetry.py \
  tests/test_sampler_logging_contract.py tests/test_runtime_metadata.py \
  tests/test_pytest_slow_filter.py --tb=short \
  --junitxml=/private/tmp/jarvis-s1-focused.xml \
  -o cache_dir=/private/tmp/jarvis-s1-pytest-cache
```

全量命令：同一环境执行 `python -m pytest -q -m 'not slow' --tb=short --junitxml=/private/tmp/jarvis-s1-fast.xml -o cache_dir=/private/tmp/jarvis-s1-pytest-cache`，stdout/stderr 写入上述 fast.log。

全量失败与既有 `test-summary.json`（2026-10-05）的逐项对照：

- `tests/test_monitor_resources.py::ScanTreeCpuTests::test_short_lived_children_are_counted`：`min(cores)` 收到 None，TypeError；测试名与失败位置/原因相同。
- `tests/test_task_schema_corpus.py::TaskSchemaCorpusTests::test_shipped_yaml_corpus_is_ascii_only`：`Jarvis-Examples/iDM/bin/iDM_Vector_LZ2026_HighER.yaml` 的变量 description 非 ASCII；文件/字段与失败原因相同。

这支持“与已记录基线相同的失败”，不冒充重新 checkout 基线后做过同日 A/B。未改这两项，也不能把 full fast 标为全绿。core/supervisor 的既有未使用 import lint 未作无关清理；本次 adapter、CLI、metadata、conftest 与新增测试的 `ruff --select F` 通过。

**文档归属。** 用户主教程为 Jarvis-Books TeX 第 38 章 **Jarvis-HEP-Sampler: Implementing Your Own Sampler**；SDK README 指向该章，authoring 文档保留 API 参考和能力边界。本设计台账不再添加另一份用户教程。按维护者要求，后续 TeX 修改交给其自动化工具编译，不默认运行手工编译。本次编译曾绕过本机 biber universal 二进制的 lipo 兼容问题，未修改系统安装。

**S1 当时结论（后续进展见 §10.7）。** 外部 independent sampler 的生成、可选校验、持久载入、普通 card 执行、Redis 发布/归档、日志/monitor 与协作停止收尾完成。尚未执行 P7 的严格恢复故障矩阵、P8 插件打包、干净构件安装/跨平台 CI 或性能基准；R01/R02/R12 与 R13 的 S2a 部分继续开放，generational 的 R11 继续留在 S3。

### 10.7 S2 落地与发布准备（2026-10-09）

**提交与范围。** S1 SDK `79972c7`、HEP `0d85950`；SDK 发布准备 `cf924cb`。本次新增严格绑定、checkpoint 配置保护、扫描上下文、源码项目包和完整构件 shell 流程。未写 TeX、未编译、未 push/tag/upload。

| 验证 | 实际结果 | 证据 |
| --- | --- | --- |
| 最终 HEP 定向（含 CLI 兼容） | **158 passed / 24 subtests**，exit 0；1 条除零警告，36.31 秒 | `/private/tmp/jarvis-s2-focused-final.log`、对应 XML；CLI、runtime、插件身份、pack/fetch、logging/metadata、slow 参数选择 |
| 真实 Redis/spawn 中断与恢复 | **2 passed**，exit 0，122.40 秒 | `/private/tmp/jarvis-s2-runtime-final.log`、XML；211 候选、block=8、selection 后 169 个 accepted，最终坐标与解析 indexed 序列/UUID 核对，单行归档 |
| 拒绝时副作用 | seed 改变、损坏副本、坏 checkpoint 均 exit 1；插件导入计数不变、Redis 端口关闭、DATABASE/checkpoint 字节不变 | 同一真实测试；checkpoint 与 DATABASE-only 两条路径都在 reload/remove 登记之后恢复，Worker 不导入入口源码 |
| 项目打包/fetch | share/repro/full/manifest 的来源与完整 hash 正确；两份同名版本保留；坏来源 fetch 不保留目标、不导入、不登记；重复打包无重复 tar 成员 | `tests/test_project_samplers.py`；与 project tools / plugin identity 的组合 **60 passed**，4.09 秒；已纳入最终 158 项 |
| 最终 SDK 构件安装 | wheel/sdist + Twine check；全新 venv 安装全部依赖，pip check；模板/跨进程 check/登记；**146 项已安装 wheel 测试通过**，2.34 秒 | `/private/tmp/jarvis-sdk-release-final.log`；构件、SHA256SUMS、release-verification.json 在 `/private/tmp/jarvis-sdk-release-final/` |
| HEP 构件安装 | **exit 0**：wheel/sdist + Twine、完整依赖安装/pip check、CLI init/check/load/validate/man、真实 Redis/spawn 16 点、项目包冲突/空登记提示→显式 load→新 scan 16 点均通过 | `/private/tmp/jarvis-hep-release.log`；构件、SHA256SUMS、release-verification.json 在 `/private/tmp/jarvis-hep-release/` |
| 全量 fast | 最终 **1267 passed / 4 failed / 165 deselected / 141 subtests**，exit 1，369.46 秒；初跑的两项 CLI 新增失败已修复。最终失败逐项见下文 | `/private/tmp/jarvis-s2-fast-final.log`、XML；初跑为 `/private/tmp/jarvis-s2-fast.log`、XML |
| 全仓额外两项单独复核 | layer force-serial 相对耗时、Worker monitor 就绪两项 **2 passed**，13.78 秒；不能据此把原全量结果改成通过 | `/private/tmp/jarvis-s2-timing-recheck.log`、XML |
| lint 与 shell | SDK `ruff check --no-cache .` 通过；HEP 新增模块/测试/scripts 的 `ruff --select F` 通过；两仓 bash -n、Bash/Zsh source/direct help 通过 | 不作无关的全仓 lint 清理；远程 CI 尚未执行 |

最终 full fast 的失败与复核：

- `test_monitor_resources.py::ScanTreeCpuTests::test_short_lived_children_are_counted`：`min(cores)` 收到 None，TypeError；与 §10.6 / 既有台账记录相同。
- `test_task_schema_corpus.py::TaskSchemaCorpusTests::test_shipped_yaml_corpus_is_ascii_only`：同一外部 Jarvis-Examples iDM card 的 description 非 ASCII，与 §10.6 相同。
- `test_layer_concurrency.py::LayerConcurrencyIntegrationTests::test_force_serial_layers_matches_explicit_serial_layers`：两次 wall time 差 1.169 秒，超过 0.56 秒阈值；独立复核通过。
- `test_worker_mvp.py::WorkerMVPTests::test_task_factory_starts_worker_and_monitor_snapshot`：固定等待 0.3 秒后 snapshot 仍为空；独立复核通过。

后两项的单独复核支持时序敏感的判断，但没有同日旧 checkout A/B，不能笼统认定“全部失败都与实现无关”。没有修改这四项用例或其实现来掩盖失败；GitHub 的全仓发布 gate 继续保留，不绕过红灯。

本机 runtime 测试仍为 macOS arm64 / Python 3.13.13 / NumPy 2.5.0。全新完整安装使用 SDK 0.1.0 与 HEP 2.0.14，NumPy 2.5.3，jsonschema 4.26.0；完整依赖版本由构件目录的验证 JSON 记录。运行流程从安装 venv 取 HEP/SDK，清除 PYTHONPATH，隔离 HOME/登记，运行目录没有 sibling 源码。

**发布入口。** 每个项目各有一个独立发布脚本：HEP 是 `./jhrel <version> [--dry]`，SDK 是 `jhsrel <version> [--dry]`（先 `source ./jhsrel`）。脚本更新 pyproject 和 runtime 版本，清理旧构件，构建 wheel/sdist 并 `twine check`；不带 `--dry` 直接上传 PyPI，`--dry` 停在上传之前。`--dry` 是唯一选项。

正式发布前选尚未使用的版本，确认远程 gate，先发布 SDK，再发布 HEP。

**本节 S2 基线当时未验收的范围（本机 P7 后续见 §10.9）。** 远程 Linux/macOS×Python 3.10–3.13 矩阵、最低依赖组合、§9.3 性能测量，以及 P7 的 DATABASE 空洞/部分 batch 失败/拒绝尾段联合故障。显式源码切换、S2c、S3 尚未实现。这些不以本机 smoke 或 146/158 条测试概括为已完成。

---

## 11. 决策记录

2026-10-06 由维护者确认。以下保留当时决策，取代 v1.0 的对应设计；2026-10-08 的验收澄清列在表后，不把已有决定当作实现完成证明。

| # | 决策 | 取代的 v1.0 设计 |
| --- | --- | --- |
| D1 | 插件是一等 Method，与内置方法共用命名空间，参数放在 Bounds | `Method: External` + `Bounds.plugin` |
| D2 | `Jarvis sampler init/check/load/list/info/remove`；`~/.jarvis/samplers/` 登记 + 快照；重新 load 覆盖；只有一个 CLI 入口 | 文件路径写在 card 里；sidecar manifest；entry-point 发现；SDK 自带 CLI |
| D3 | `OPTIONS_SCHEMA` 可选；插件的正确性由用户负责；HEP 只守自己的不变量（§4.4、§5.3、§7） | sidecar 强校验、strict validate、资源依赖 SHA 清单 |
| D4 | profile + 完整元数据在 v1 冻结；independent 在 S1 实现，generational 在 S3 实现 | 只开放 `IndependentSampler`，S3 另起一套 `FeedbackEngine` |
| D5 | `prepare` / `draw` / `draw_block`；输出区间为 `[0,1]`，并检查映射后的物理量是否有限 | 单点纯函数 `draw`；`[0,1)` |
| D6 | SDK 提供 `log` 接口和 4 个 hook；telemetry 和 monitor 有通用的插件布局 | 插件没有 logger；monitor 回退成 `CUSTOM METHOD` |
| D7 | 硬身份/软身份两层；续跑默认用 scan 自己的源码副本；`--accept-plugin-change`；S2 推广到内置方法 | 源码或环境版本不同就拒绝，且没有覆盖手段；只给 External 记录身份 |
| D8 | SDK 是 HEP 的硬依赖 | lazy 可选依赖 |
| D9 | 官方 `LatinHypercube` 走插件路径（dogfooding） | 不迁移任何方法 |

v1.0 中仍然保留的设计：

- SDK 与 HEP 单向依赖；插件不接触 Sample、UUID、Redis、Worker；
- 每个候选点派生独立随机流；
- 候选/accepted 双游标，从 DATABASE 的 durable prefix 重放；
- preflight 在 Redis、Worker 启动之前；
- 只导入 card 选中的插件，Worker 不导入插件；
- 插件出错时不回退到其他 sampler；
- 不承诺沙箱和强制超时。

**2026-10-08 审阅澄清：**

- D1/D2/D3/D8/D9 保持：一等 Method、SDK 模板、单一 HEP CLI、持久登记与覆盖、schema 可选、SDK 单向依赖、官方插件 dogfood。此次不引入 manifest 编辑、强制 schema 或修改 HEP 源码的用户步骤。
- D4 的“冻结”限于预留 profile 名称和已有元数据；运行中的 generational 状态/反馈语义尚需 R11 验收。独立于 S1/S2 首发。
- D6 的 optional hooks 必须真的可选：基础参数日志、宿主计数和租约维护不能依赖它们。实现修正见 R03/R05。
- D7 保留硬/软身份与默认 scan 副本的方向，补上确切 hash 绑定、切换源码前的已持久化点一致性检查；相同 UUID 不构成科学一致性证明。内置方法推广调整到 S2c，不抢占外部插件闭环。
- 2026-10-16：`--accept-plugin-change` 不进入首发，不实现。首发只提供与内置 Sampler 相同的使用体验。已绑定 scan 上 NumPy、Python 或 SDK 变化维持警告后继续。
- 台账不再把测试条数、源码文件存在、发布 workflow 存在等同于全契约验收或公开发布；今后按 §10.1 的退出条件更新同一份文档。

`recent-code-review.md` 的 F1–F5 属于更早的 HEP bug 审阅，与本台账 R01–R13 分开跟踪；本次未重新核验其修复状态。


### 10.8 原生使用 M1 修复与验收

用户再次确认 load 后原生使用与 sampler 统一管理。SDK 在共享锁内取得 entry 和源码并核对 hash，锁外执行同一 bytes；受管快照拒绝绝对路径、越界和 symlink；sources.json 替换是提交点，退役清理失败只告警并保守重试。报告 hook 调用/迭代/转换统一隔离，控制流异常仍传播；NumPy 常见标量正常归一，float64 溢出转为带 cause 和候选位置的 PluginOutputError。五组固定 PCG64 raw 向量补充协议回归，不修改随机流。

SAM-HEP-02 的真实中断先暴露 live cursor 超出 durable prefix，修复为 indexed 只保存已认可的批次安全快照；随后发现 supervisor 在 prefix 重放后再次导入旧 payload，造成重复提交和归档等待，现只给非 indexed 保留晚导入。三组确定性失败回归（check 隔离、安全保存、恢复顺序）已转为通过，正式身份和 prefix 严格拒绝保持，feedback/nested 兼容行为有回归。

本机独立 CLI 跑通 load→query→check×2→run→reload→新 scan→remove→旧绑定 resume，以及空 registry 的 pack/fetch→显式 load→原生运行。真实 SIGINT 的 checkpoint 和 DATABASE-only 两条恢复核对169个 accepted 点的坐标、UUID和来源。完整安装包还复核主/子进程同一 venv、包从 site-packages 导入、五个 HEP 实现文件与工作源码 hash 一致。SDK 新增完整 `examples/MyLatinHypercube.task.yaml`，使用原生命令验证32点LHS设计。

有效测试命令、实际数量、环境、构件与文件 hash 见 [验收记录](native-parity-m1-verification.json)。M1 结束时联合矩阵、并发移除和资源启动重试仍开放，后续本机结果见 §10.9；远程 Linux/macOS/Python、Books 最终教程、显式改源恢复和 generational 保持各自后续范围。


### 10.9 本机 M2 联合故障验收

新增 SAM-HEP-03 来自真实稀疏 DATABASE：prefix=2，保留索引0/1/3/6/8（其中1/6为 Failed）。原恢复重复计算3/6/8，最终科学行22却等待25，报归档停滞。现保留 pending={3,6,8} 并在 FixedSet 提交时去重，只新计算17个缺失点；5个 Failed 保持终态。候选37/accepted22、坐标/UUID与不中断参考逐行一致。checkpoint/DB-only × workers/batch=(1,3)/(2,5) 四组真实 Redis/spawn 均通过；稀疏 crash-state 是明确构造的故障夹具，真实 SIGINT 证据独立保留。

补齐两条“22 accepted 已全部 durable，只剩末11候选全拒绝”的恢复：新增提交与 Worker 计算均0，DATABASE不增行，无超时；checkpoint只draw块3/4，DB-only在控制端重放0..4。既有 RunOutcome 在零新任务时 submitted 回退22个历史终态，未改其公共语义；以实际提交返回值和 Worker 调用记录证明零物理计算。

运行中 v1 经 reload v2/remove 后继续完成或 SIGINT→无登记resume，两组各37点全量坐标/UUID一致。资源故障测试记录实际启动的 managed Redis/两个 Worker，失败后全部清理、0科学行；fresh提示 `--resume`，reload/remove与配置/源码拒绝后仍按v1成功重试8点。

最终 fast gate：1519 passed、478 subtests passed、3 skipped、179 slow deselected。安装包23个不同定向用例通过（19项M2新增及4项既有真实check/SIGINT回归，重复复验单列），SDK193项及Ruff通过，wheel/sdist与依赖检查通过。7个 HEP 实现文件与工作源码一致，SDK仍0.1.0/API1/STREAM1/registry1。原始产品失败、夹具隔离修正、实际命令、构件/hash及范围见 [M2验收](native-parity-m2-verification.json)。M2只关闭本机 independent 范围，不宣称CSV/所有内置方法、远程兼容、性能、Books、显式改源或generational已交付。
