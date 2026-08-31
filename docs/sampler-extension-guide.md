# Adding `RLTPMCMCSampler`

这份说明只覆盖一个扩展点：如何在 Jarvis-HEP 中添加
`RLTPMCMCSampler`。它的职责是运行 TPMCMC，并使用 RL 控制 TPMCMC 的
超参数；它不是一个与 TPMCMC 并列的普通 sampler。

## 1. 继承关系

推荐的继承关系是：

```
MCMCBaseSampler
└── PTMCMCBase
    └── PTMCMCSampler
        └── RLTPMCMCSampler
```

因此，RL-TPMCMC 应该继承 `PTMCMCSampler`，而不是直接继承
`MCMCBaseSampler`：

```
from jarvishep2.sampling.ptmcmc import PTMCMCSampler


class RLTPMCMCSampler(PTMCMCSampler):
    """TPMCMC whose ladder/proposal controls are updated by an RL policy."""
```

原因是 `PTMCMCSampler` 已经提供了 TPMCMC 所需的：

- temperature ladder；
- parallel-chain exchange；
- exchange acceptance statistics；
- PT 的 barrier/generation 调度；
- method-state 的 checkpoint/resume 支持。

RL 层只负责观察运行状态、产生 action，并在安全的采样屏障上应用
TPMCMC 超参数更新。`RLController` 应该作为组合对象保存，不要把它作为
第二个父类。

## 2. RL 应该接入哪里

不要一开始重写 `run_adaptive()`。直接复用父类的 barrier loop。其关键顺序是：

```
run generation
    ↓
exchange chains
    ↓
_on_generation_completed()
    ↓
checkpoint at sampling barrier
```

`RLTPMCMCSampler` 应该在 `_on_generation_completed()` 中做控制。这样每次
RL action 执行时，上一代的 chain state 和 exchange 结果已经完整，且没有
未完成的 generation。

```
from collections.abc import Mapping
from typing import Any


class RLTPMCMCSampler(PTMCMCSampler):
    """TPMCMC sampler with RL-controlled runtime hyperparameters."""

    def __init__(self) -> None:
        super().__init__()
        self._rl_controller = None
        self._rl_state: dict[str, Any] = {}
        self._rl_effective_config: dict[str, Any] = {}

    def _configure_method(self, bounds: Mapping[str, Any]) -> None:
        # Keeps PTMCMC's ladder, exchange, scale, and bounds validation.
        super()._configure_method(bounds)

        rl_config = bounds.get("rl", {})
        self._rl_controller = RLController()
        self._rl_controller.configure(rl_config)

    def _on_generation_completed(self) -> None:
        # PTMCMCSampler currently uses this hook for progress reporting.
        super()._on_generation_completed()

        if self._rl_controller is None:
            return
        if not self._rl_controller.due(self._generation):
            return

        observation = self._rl_controller.observe(self._ensure_registry())
        action = self._rl_controller.act(observation)
        self._apply_rl_action(action)
        self._rl_state = dict(self._rl_controller.export_state())

    def _apply_rl_action(self, action: Mapping[str, Any]) -> None:
        """Validate and apply only supported, runtime-safe controls."""
        # See the next section for the required synchronization rules.
        ...
```

`_ensure_registry()` 只是示意：实现时应使用该类现有的 chain registry
访问方式，不要复制一份独立的 chain state。`RLController` 的最小职责可以是：

```
class RLController:
    def configure(self, config: Mapping[str, Any]) -> None: ...
    def due(self, generation: int) -> bool: ...
    def observe(self, registry: Any) -> Mapping[str, Any]: ...
    def act(self, observation: Mapping[str, Any]) -> Mapping[str, Any]: ...
    def export_state(self) -> Mapping[str, Any]: ...
    def import_state(self, state: Mapping[str, Any]) -> None: ...
```

## 3. action 可以控制什么

第一版只开放 TPMCMC 中确实可安全更新的参数，例如：

```
rl:
  enabled: true
  algorithm: ppo
  control_interval: 50
  warmup_iters: 2000
  action_space:
    - temperature_ladder
    - proposal_scale
    - exchange_interval
```

### Temperature ladder

如果 action 更新 ladder，必须同时更新：

1. sampler 的 `_temperature_ladder`；
2. registry 中每一条 `ChainRuntime.temperature`；
3. RL checkpoint 中的 effective ladder。

只改 `_temperature_ladder` 不够，因为运行中的 swap acceptance 和
generation absorb 使用的是每条 chain 的 runtime temperature。

必须验证：

- ladder 长度等于 `num_chains`；
- 第一项为 `1.0`；
- 所有温度为正；
- 温度严格递增；
- action 不改变 `num_chains`。

### Proposal scale

如果 action 更新 proposal scale，必须更新实际运行 engine 使用的
`proposal_scale`；只改 YAML/config 副本不会影响当前 chain。

所有 scale 必须为有限正数，并且不能改变参数维度或 sampler engine 类型。

### Exchange interval

`exchange_interval` 必须是正整数。它可以在 generation barrier 更新，但不应
让 RL action 改变当前 generation 的执行方式。

## 4. 不允许动态改变的内容

为了保持 chain registry、checkpoint 和统计量的一致性，第一版不要让 RL
动态改变：

- `num_chains`；
- 参数维度或 parameter bounds；
- chain engine 类型；
- sampler 的 method name；
- 已经完成的 sample/generation 编号；
- checkpoint 文件格式。

RL action 必须是一个经过边界投影或拒绝的有限值；遇到非法 action 时应记录
原因并使用上一组有效参数，而不是让 sampler 进入半更新状态。

## 5. checkpoint / resume

RL controller 的可恢复状态必须是纯数据，不能直接 pickle controller 对象。
至少保存：

```
{
    "policy_state": ...,
    "rng_state": ...,
    "last_action": ...,
    "effective_temperature_ladder": ...,
    "effective_proposal_scales": ...,
    "effective_exchange_interval": ...,
}
```

实现 method-state hook 时必须保留父类状态：

```
def _export_method_state(self) -> dict[str, Any]:
    state = dict(super()._export_method_state())
    state["rl_state"] = dict(self._rl_state)
    state["rl_effective_config"] = dict(self._rl_effective_config)
    return state


def _import_method_state(self, state: Mapping[str, Any]) -> None:
    super()._import_method_state(state)
    self._rl_state = dict(state.get("rl_state", {}))
    self._rl_effective_config = dict(
        state.get("rl_effective_config", {})
    )
    if self._rl_controller is not None:
        self._rl_controller.import_state(self._rl_state)
```

新增的 instance attributes 也要按 checkpoint 语义分类：

```
_checkpoint_saved_attributes = frozenset({"_rl_state", "_rl_effective_config"})
_checkpoint_excluded_attributes = frozenset(
    set(PTMCMCSampler._checkpoint_excluded_attributes)
    | {"_rl_controller"}
)
```

controller 对象本身排除在 checkpoint 外；它的 policy、RNG 和当前 action
必须通过 `_rl_state` 保存。若新增 observation cache 或临时 buffer，也要
明确标记为 excluded，或把它作为纯数据保存。

动态 proposal scale 需要特别注意：当前 checkpoint validator 会把
`Bounds.proposal_scale` 当作配置基线进行检查。建议保留 YAML 中的初始
proposal scale，把运行中的动态值放在 engine/native method state 中；如果
要覆盖配置字段，就必须同步修改 validator 的语义和对应测试。

## 6. task-card 配置接口

新 sampler 的目标配置可以设计为：

```
Sampling:
  Method: RLTPMCMC
  Bounds:
    num_chains: 4
    num_iters: 10000
    proposal_scale: 0.2
    temperature_ladder: [1.0, 2.0, 4.0, 8.0]
    exchange_interval: 2
    seed: 7
    rl:
      enabled: true
      algorithm: ppo
      control_interval: 50
      warmup_iters: 2000
```

`rl` 应放在 `Sampling.Bounds` 下，由 sampler 自己消费。不要重新启用旧的
顶层 `Sampling.Control`、`Sampling.Diagnostics`、`Sampling.PPO` 或
`Sampling.Reward` block；V2 schema 当前会拒绝这些旧接口。

## 7. 注册到 Jarvis-HEP

除了实现类，还必须完成以下注册，否则 YAML 中的 `Method: RLTPMCMC`
不能通过 catalog/distributor。

### Distributor factory

在
`jarvishep2/distributor.py` 中添加 lazy factory，并放入
`_builtin_factories`：

```
def _factory_rltpmcmc():
    from jarvishep2.sampling.rl_tpmcmc import RLTPMCMCSampler

    return RLTPMCMCSampler


_builtin_factories["RLTPMCMC"] = _factory_rltpmcmc
```

factory 必须是零参数、返回 class 的 callable。不要在 module import 时
创建 controller、读取配置或启动 worker。

### Sampler catalog

在 `jarvishep2/sampler_catalog.py` 的 `_SPECS` 中加入：

```
_method(
    "RLTPMCMC",
    schema="rltpmcmc.json",
    stateless=False,
    statistical=True,
    mcmc=True,
    mcmc_pipeline="barrier_coupled",
    family="MCMC",
)
```

`stateless=False` 是必须的，因为 sampler 持有 chain、exchange、RL policy
和 checkpoint state。`mcmc_pipeline="barrier_coupled"` 表示 RL action 与
MCMC generation/exchange 共用 barrier 语义。

### Schema 和 manifest

新增：

```
jarvishep2/cardload/schema/sampling/methods/rltpmcmc.json
```

并在
`jarvishep2/cardload/schema/sampling/methods/manifest.json` 中登记 schema
文件和 `RLTPMCMC` method。schema 至少要声明：

- `Sampling.Method == "RLTPMCMC"`；
- 必需的 TPMCMC bounds；
- `rl` 的类型和 action 配置；
- `num_chains`、ladder 长度和正值约束；
- 不允许旧顶层 RL blocks。

只有当存在 schema 静态类型无法表达的跨字段规则时，才在 semantic contract
层补充检查，例如 `len(temperature_ladder) == num_chains`。

## 8. 最小验收清单

提交前至少检查：

- `RLTPMCMCSampler` 继承 `PTMCMCSampler`；
- 没有重复实现 PT ladder、exchange 或 barrier loop；
- RL action 只在 generation barrier 生效；
- ladder 更新同步到 sampler 和所有 chain runtime；
- proposal scale 更新同步到实际 engine；
- 非法 action 会被拒绝或投影；
- `num_chains`、维度和 engine 类型保持静态；
- `_export_method_state()` / `_import_method_state()` 都调用 `super()`；
- controller 对象不直接进入 checkpoint；
- distributor factory 为零参数 lazy factory；
- catalog、schema、manifest 三处都已注册；
- fresh run、checkpoint/resume、非法 action、ladder 更新和旧 RL block
  rejection 都有测试。

一个最小的 smoke test 应验证：使用同一 seed 运行到 checkpoint，resume 后
能够恢复 RL policy state、effective ladder、proposal scales 和 exchange
interval，而不出现 schema mismatch 或 chain-count mismatch。
