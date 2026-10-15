# 近期代码 Review — 2026-10-05

1. **Roadmap phase**：现有采样器/monitor 正确性回归检查；新 S0–S3 外部 sampler 设计准备。
2. **Affected modules**：sampling lifecycle、Dynesty/MultiNest、monitor collector、未提交 Diver engine，以及对应测试；静态检查覆盖近期 process cleanup、crypto、CLI help、runtime logging 等变更。
3. **Public interfaces**：本次不修改 Python/YAML 运行接口；新增技术设计与审查文档，`.gitignore` 放行本目录。
4. **Main risks**：错误收敛标记、跨平台 CPU 测量失效、Diver 除零与 checkpoint 恢复/保存缺陷。
5. **Verification performed**：全仓 fast gate、失败项沙箱外重跑、确定性最小复现、实际 nested budget-stop 与 checkpoint 写盘故障注入。

## 结论与审查范围

发现 **5 个可采取行动的问题**：已合并近期代码中 2 个，工作区未提交 Diver port 中 3 个。没有实施 bug 修复；下列建议是后续 patch 的要求。外部 sampler 的设计见 [技术设计](technical-design.md)（v2.0，已定型）。

比较基线为 `acb2d9fc..0bbcf675`：2026-09-22 后的 48 个 tracked 文件变更，重点检查最近两个合并 PR 的监控、日志和采样器增量；同时审查开始任务时已有的 `pyproject.toml` package-data 增量、`sampling/Source/Diver/` 和 `tests/test_diver_engine.py`。没有把工作区用户改动当成已发布版本，也没有改写它们。

这不是所有物理算法或 vendored 数学库的完整科学认证，也不是所有文件逐行审查。既有 MCMC/selection/checkpoint 单测参与 fast gate，但没有重做高维多峰分布验证、所有慢测试或所有 Python/OS 组合。

## 已确认问题（按严重程度）

### F1 — [P1] Diver 的零目标均值会触发除零崩溃（未提交）

位置：[engine.py](../../jarvishep2/sampling/Source/Diver/engine.py)，`_update_convergence`，行 713–719，尤其行 718。

合法目标函数恒返回 0 时，第一代把 meanlike 设为 0；下一代执行 Python float 的 `current / self.meanlike`，抛出 `ZeroDivisionError`。全零 χ²、平坦/离散 plateau 都可能出现。不是仅有 NaN 数值警告，而是直接中止运行。

最小复现使用 seed=7、1D、maxgen=5、convsteps=2、disableIO=True；结果：`ZeroDivisionError: float division by zero`。上游 Fortran 也使用同样比值，但 IEEE 浮点除零与 Python float 抛异常不具有相同行为；移植必须明确处理零均值和非有限值，不能照搬表达式后假定等价。

建议：定义 meanlike/current 均为零时的稳定停止判据，以及符号变化/负目标时的适用范围；统一处理 NaN/Inf，不让错误数值进入 best/convergence。回归覆盖恒 0、接近 0、平台和负目标的预期行为。

### F2 — [P1] Diver 覆盖唯一 checkpoint 时不具备原子性（未提交）

位置：[engine.py](../../jarvishep2/sampling/Source/Diver/engine.py)，`_write_checkpoint`，行 784。

`np.savez_compressed(self._checkpoint_path(), ...)` 直接覆盖上次的唯一 `.diver.npz`。进程被杀、磁盘写失败或序列化异常发生在覆盖期间，会同时损失新 checkpoint 和上一份可恢复状态。与 HEP runtime 的临时文件+replace 保存方式不同。

故障注入：先生成可读 checkpoint，再在第二代保存时让 `numpy.lib.format.write_array` 抛 OSError。上一份 **2813 bytes** 文件变为 **136 bytes**，内容已变，重新 resume 报错。复现只使用临时目录，没有损坏真实 scan。

建议：同目录临时文件保存完整 archive，关闭并校验后 `os.replace`；异常保留旧 checkpoint。`.raw/.sam` append 与 checkpoint generation 也应明确 crash 后重复/截尾策略。回归要求注入写盘失败后旧文件仍可恢复；不是只测试正常保存。

### F3 — [P2] nested 正常返回被无条件报告为已收敛（已合并）

位置：[dynesty_sampler.py](../../jarvishep2/sampling/dynesty_sampler.py)，行 1541–1543；由 `b5ab0729` 新增。`run_adaptive` 行 1380–1381 在任何正常 `run_nested` 返回后设 `_finished=True`。

`run_nested` 可能因为 maxiter/maxcall/maxbatch 耗尽或用户 stopping function 返回而正常结束。`log_stop_reason` 只检查 `_finished`，因此统一打印 `converged (dlogz < ...)`。动态 Dynesty 的实际停止规则还可能是 ESS/custom stop，不能用静态 dlogz 描述。

实际复现：MultiNest profile、nlive=10、seed=7、`run_nested.maxiter=1`、`dlogz=1e-8`，使用 in-memory Redis 和反馈线程。引擎明确警告 maxiter/maxcall 提前停止、证据阈值未达到；但 hook 输出 `converged (dlogz < 1e-08)`。这会误导用户把预算耗尽的结果当成科学收敛结果。

建议：保存明确 stop reason，区分 convergence / budget_exhausted / custom_stop / interrupted；未知时写“engine returned”，不要猜测已满足阈值。static/dynamic 分别测试，加入 budget-stop summary contract test。

### F4 — [P2] macOS 短命 calculator 的 scan CPU 会持续未知/低估/冻结（已合并）

位置：[collector.py](../../jarvishep2/monitor/collector.py)，行 763–765、scan aggregate 行 647–649；由 `cdd9859a` 引入的累计 CPU 方案。

算法把所有活进程的 own CPU 与 reaped children CPU 相加后求差，假设父进程的 children_user/system 能补上已经退出的 calculator。macOS 的这两个 psutil 字段总为 0。活子进程退出后 aggregate 下降，`_CpuMeter.rate` 返回旧值；尚无旧值时返回 None。CPU 忙碌也可以多次显示未知或沿用旧的低读数。

沙箱外 `test_short_lived_children_are_counted` 在全仓和单独重跑都失败：三次测量含 None，`min(cores)` 抛 TypeError。独立子进程诊断连续六次读父进程 CPU：children_user/system 一直 0，前三次 scan_cores=None，后三次同一个约 0.417 的旧读数。此问题不归因于 TCP/sysctl 沙箱限制。[psutil 官方 Process.cpu_times 文档](https://psutil.io/api/#psutil.Process.cpu_times) 明确标注 Windows/macOS children counters 为 0。

建议：能力检测并区分 available/live-tree-estimate/full-accounted；在不支持累计子进程 CPU 的平台按 PID+create_time 计算活进程增量，标明短命进程遗漏，或由启动/回收 calculator 的 Worker 发布 wait/rusage 累计 CPU。没有数据的部分不能通过保留旧值伪装成最新测量；CI 分平台验证精确/估算语义。

### F5 — [P2] 恢复已收敛的 Diver checkpoint 会继续演化（未提交）

位置：[engine.py](../../jarvishep2/sampling/Source/Diver/engine.py)，resume 分支行 381–383，checkpoint metadata 行 762–768 与 `_load` 行 806–811。

checkpoint 没保存/恢复 converged、quit 或 terminal reason；加载后 `_Engine` 默认 converged=False/quit=False，直接进入下一代。与 fresh 分支的停止检查不同，resume 已收敛的结果会再次评估目标、修改 population 和输出。

复现：常值目标 1、seed=7、convsteps=1，第一轮 generation=2 / calls=17 / converged=True；相同配置 resume 后 generation=3 / calls=27，再次 converged=True。每次 resume 都会新增评估直到 maxgen，不能维持完成状态的语义。

建议：保存 terminal flags/reason，恢复后保持结果或明确提供新一轮 continue 的显式操作；用户 quit 与数学 convergence 区分，旧 checkpoint 缺标记需明确兼容规则。回归包括已收敛、已达预算、quit 和中途恢复，并比较 calls/population/output 未增加。

## 测试与证据

运行环境：macOS，Python **3.13.13**、NumPy **2.5.0**、psutil **7.2.2**。日期使用 Adelaide client date。

| 验证 | 结果 | 解释 |
| --- | --- | --- |
| 沙箱内 `python -m pytest -q -m 'not slow'` | 1166 passed / 41 failed / 163 deselected；140 subtests passed | 多个 loopback TCP、psutil/sysctl 权限失败；不能直接判为 41 个回归 |
| 沙箱外同一 fast gate，`--tb=short --junitxml=...` | **1205 passed / 2 failed / 163 deselected；140 subtests passed**，240.24 s | 完整回归仍未全绿 |
| CPU 失败项沙箱外单独重跑 | 1 failed，3.60 s | 确认 F4 |
| 第二个 fast gate failure | sibling `Jarvis-Examples/iDM/bin/iDM_Vector_LZ2026_HighER.yaml` 的 Variables[1].description 非 ASCII | 外部 fixture 状态；不是本仓近期变更新引入，不改兄弟项目 |
| [最小复现脚本](reproduce-findings.py) | F1/F2/F3/F5 均观察到描述的错误行为 | exit=0 只表示 evidence script 完成，不代表算法通过 |

完整本地 pytest 日志和 JUnit 位于 `/private/tmp/jarvis-hep-review-2026-10-05.log` 与 `.xml`；简明可持久证据存于本目录 [test-summary.json](test-summary.json) 和 [reproduction-results.json](reproduction-results.json)。沙箱外执行已获本次工具授权，没有外部消息或服务发布。

## 对新 Sampler 设计的直接影响

- stop reason 必须来自实际算法/预算，不能从 finished flag 推断。
- 用户 Python 的数值与返回 shape 必须在 Sample 提交前检查；checkpoint 数值同样需要验证。
- 恢复身份和完成状态是协议组成部分，不能只保存一份 RNG seed。
- HEP 继续拥有已有原子 checkpoint 基础设施；不要让外部 draw 函数自行保存运行状态。
- worker/spawn 与 Linux/macOS 都在 CI 矩阵中；单平台通过不足以确认运行时兼容。

## 文档与未完成项

上级 AGENTS 所指的五份 canonical planning 文件缺失；现有本地 sampler extension guide 还有 factory 返回 class 的过时描述。本设计明确记录当前实例契约与新的 S0–S3 路线，未伪造历史规划文件。下一实现阶段应将这些接口决定纳入实际维护层并同步用户 guide。

外部 sampler 设计已定型为 technical-design.md v2.0。以上 5 个代码问题仍为 **open**，未提交修复 PR，也未对外创建 Jarvis-HEP-Sampler repository。此报告足以作为下一批 bugfix 与 SDK 实现的输入。
