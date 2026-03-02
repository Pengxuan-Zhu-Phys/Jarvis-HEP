# Jarvis-HEP 历史瘦身总结

Last updated: 2026-03-02
对应计划: [HISTORY_SLIMMING_TODO.md](./HISTORY_SLIMMING_TODO.md)

## 多 Round 总览（倒叙）

### Round 5（Phase 3 元数据合规修复）

- 范围: `pyproject.toml`
- 完成:
  - 将 `project.license` 从旧式 table 写法 `{"file" = "LICENSE"}` 切换为 SPDX 字符串 `"MIT"`
- 验证:
  - `python -m build` 通过
  - 构建日志中已不再出现 `project.license` deprecation warning

### Round 4（Phase 2 第三方快照瘦身收尾）

- 范围: `src/Sampling/Source/Diver`, `src/Sampling/Source/Dynesty`
- 完成:
  - 删除失效第三方 `gitlink`：`src/Sampling/Source/Diver/Diver`
  - 删除 Dynesty 冗余文档/打包元文件（`MANIFEST.in`, `README.md`, `RELEASE.md`, `TESTING.md`, `pytest.ini`, `requirements*.txt`, `setup.cfg`）
  - 清理 `Source` 目录下缓存目录（`__pycache__`，工作区清理）
- 体积结果（源码目录）:
  - `src/Sampling/Source`：`~27.9 MB -> 664 KB`
  - 主要来源：移除 `Diver/Diver` Fortran 快照与历史文档示例
- 验证:
  - `python -m unittest discover -s tests -p 'test_*.py'` 通过
  - `python -m py_compile`（核心相关文件）通过
  - `python -m build` 通过

### Round 3（Phase 3 发布体积优化）

- 范围: `pyproject.toml`, `dist/*`（重建验证）
- 完成:
  - 将 `tool.setuptools.package-data.src` 从 `Sampling/Source/**/*` 收敛为 `Sampling/Source/**/*.py`
  - 移除 wheel/sdist 中与运行无关的 `pdf`、示例文件与 `__pycache__/*.pyc`
  - 保留 `card/**/*`、`icons/**/*` 以及 `Sampling/Source` 运行所需 Python 源码
- 体积结果（同版本号重新构建对比）:
  - wheel: `8,131,499 -> 3,453,838 bytes`（减少 `4,677,661 bytes`, `57.53%`）
  - sdist: `8,055,969 -> 3,389,528 bytes`（减少 `4,666,441 bytes`, `57.93%`）
- 验证:
  - 包内容检查中 `/.git/`、`__pycache__/`、`*.pdf`、`example_*` 均为 0
  - `python -m unittest discover -s tests -p 'test_*.py'` 通过
  - `python -m py_compile`（核心改动相关文件）通过

### Round 2（Phase 2 局部代码清理）

- 范围: `src/dataconvert.py`, `src/Sampling/dnn.py`
- 完成:
  - 清理 `dataconvert.py` 未完成旧 API 分支，补齐 `_ensure_dataframe` 统一转换入口
  - 清理 `dnn.py` 无用 import（历史残留）
  - 修复 `FeedforwardNN.forward` 历史残留实现，改为直接走 `self.model(x)`
- 变更规模:
  - `2 files changed, 24 insertions(+), 36 deletions(-)`
- 验证:
  - `python -m unittest discover -s tests -p 'test_*.py'` 通过
  - `python -m py_compile`（相关核心文件）通过

### Round 1（Phase 1 安全瘦身）

- 范围: 历史目录、失活配置、打包负担
- 完成:
  - 删除 `src/Sampling/Old/`
  - 删除 `src/Sampling/Source/rjdnn.zip`
  - 删除 `src/logger.py`, `src/logger_old.py`
  - 删除 `src/card/possion.json`
  - `src/card/preference.json` 移除失活键 `PoissonSchema`
  - `pyproject.toml` 移除 `Sampling/Old/**/*` package-data 打包项
- 变更规模:
  - 删除文件数: 16（索引统计）
  - 文本删除: 2559 行
  - 删除大文件: 29,104,709 bytes + 13,652,498 bytes
- 验证:
  - `python -m unittest discover -s tests -p 'test_*.py'` 通过
  - `py_compile` 与 import smoke 通过

## Round 5（Phase 3 元数据合规修复）

### 目标

- 修复 `setuptools` 对 `project.license` 旧写法的弃用警告。
- 保持构建行为不变，仅提升发布元数据合规性。

### 已完成变更

1. `pyproject.toml`
- 修改:
  - `license = { file = "LICENSE" }`
  - -> `license = "MIT"`

### 验证结果

1. 构建检查
- 命令: `python -m build`
- 结果: 构建通过，且不再出现 `project.license` 弃用警告。

### 风险评估

- 该变更仅作用于打包元数据，不影响运行时代码路径与计算语义。

## Round 4（Phase 2 第三方快照瘦身收尾）

### 目标

- 完成 Phase 2 最后一项：移除第三方源码快照中的非运行冗余内容。
- 在不改变采样主链路的前提下，减少仓库历史包袱。

### 已完成变更

1. 删除失效 gitlink（第三方快照）
- 删除: `src/Sampling/Source/Diver/Diver`
- 说明: 该路径为历史遗留 submodule gitlink，当前运行链路不依赖此 Fortran 快照。

2. 删除 Dynesty 冗余文件
- 删除:
  - `src/Sampling/Source/Dynesty/MANIFEST.in`
  - `src/Sampling/Source/Dynesty/README.md`
  - `src/Sampling/Source/Dynesty/RELEASE.md`
  - `src/Sampling/Source/Dynesty/TESTING.md`
  - `src/Sampling/Source/Dynesty/pytest.ini`
  - `src/Sampling/Source/Dynesty/requirements-dev.txt`
  - `src/Sampling/Source/Dynesty/requirements.txt`
  - `src/Sampling/Source/Dynesty/setup.cfg`
- 保留:
  - Dynesty 运行代码（`py/dynesty/*.py`）
  - `LICENSE` 与本地运行相关 `.py` 文件（如 `setup.py`, `priors.py`）

3. 工作区清理
- 删除缓存目录（未追踪）:
  - `src/Sampling/Source/**/__pycache__/`

### 量化结果

1. 源码目录体积
- `src/Sampling/Source`: `~27.9 MB -> 664 KB`
- `src/Sampling/Source/Diver`: `~27 MB -> 96 KB`
- `src/Sampling/Source/Dynesty`: `920 KB -> 536 KB`

### 验证结果

1. 回归测试
- 命令: `python -m unittest discover -s tests -p 'test_*.py'`
- 结果: `Ran 4 tests ... OK`

2. 语法检查
- `python -m py_compile src/Sampling/diver.py src/Sampling/dynesty.py src/distributor.py src/client.py src/hdf5writer.py src/observable_io.py` 通过

3. 构建检查
- `python -m build` 通过

### 风险评估

- 删除目标均为第三方快照文档/示例或失效 gitlink，不在主运行导入路径上。
- 当前采样主链路仍使用 `Source/Diver/*.py` 与 `Source/Dynesty/py/dynesty/*.py`，功能不变。

## Round 3（Phase 3 发布体积优化）

### 目标

- 对 PyPI 构建产物做“可量化”的体积基线与瘦身闭环。
- 在不影响运行链路的前提下，仅保留发布包所需资源。

### 已完成变更

1. `pyproject.toml`
- `tool.setuptools.package-data.src`:
  - `Sampling/Source/**/*` -> `Sampling/Source/**/*.py`
- 保留:
  - `card/**/*`
  - `icons/**/*`

2. 构建产物内容收敛（`dist/*`）
- 删除了发布包内的非运行项:
  - `src/Sampling/Source/Diver/Diver/*.pdf`
  - `src/Sampling/Source/Diver/Diver/example_*/*`
  - 全部 `__pycache__/*.pyc`
- `Source` 目录仅保留运行所需 Python 代码。

### 量化结果

1. 构建产物体积
- wheel: `8,131,499 -> 3,453,838 bytes`（-`57.53%`）
- sdist: `8,055,969 -> 3,389,528 bytes`（-`57.93%`）

2. 风险项扫描
- wheel/sdist 中:
  - `/.git/` 命中: 0
  - `__pycache__/` 命中: 0
  - `*.pdf` 命中: 0
  - `example_*` 命中: 0

### 验证结果

1. 回归测试
- 命令: `python -m unittest discover -s tests -p 'test_*.py'`
- 结果: `Ran 4 tests ... OK`

2. 语法检查
- `python -m py_compile src/hdf5writer.py src/observable_io.py src/client.py src/distributor.py src/Sampling/dynesty.py src/Sampling/diver.py` 通过

### 风险评估

- 本轮只改发布白名单，不改计算逻辑与流程控制。
- 已确认 `Source` 目录运行所需 `.py` 文件仍在包内。
- 仍有一个发布规范项待处理: `project.license` 旧写法会触发 setuptools 弃用警告（不影响当前运行）。

## Round 2（Phase 2 局部代码清理）

### 目标

- 清理陈年残留实现，减少“看起来可用但实际不完整”的代码分支。
- 在不改变主流程行为的前提下，提高后续维护可读性。

### 已完成变更

1. `src/dataconvert.py`
- 删除重复/无意义分支（`o(idx=None)` 中重复判断）。
- 移除未使用 import。
- 为历史转换接口补齐统一入口 `_ensure_dataframe`。
- `to_numpy` / `to_listdict` 改为基于统一入口的稳定实现，移除未定义成员依赖。

2. `src/Sampling/dnn.py`
- 删除无用 import: `ABCMeta`, `abstractmethod`, `re`, `time`, `BoolConversionError`, `sympy`, `DataLoader`, `TensorDataset`, `is_solenoidal`, `futils`。
- 修复历史残留 `FeedforwardNN.forward`，避免访问不存在的 `self.layers/self.dropout/self.output_layer`。

### 验证结果

1. 回归测试
- 命令: `python -m unittest discover -s tests -p 'test_*.py'`
- 结果: `Ran 4 tests ... OK`

2. 语法检查
- `python -m py_compile src/dataconvert.py src/Sampling/dnn.py ...` 通过

### 风险评估

- 变更集中在历史残留和无用代码，不涉及核心 workflow 调度逻辑。
- 仍建议在后续真实 DNN 任务卡上做一次端到端烟测。

## Round 1（Phase 1 安全瘦身）

### 目标

- 删除已确认无主运行链路引用的历史代码与资产。
- 减少打包负担与维护噪音。
- 不改变当前主流程行为。

### 已完成变更

1. 删除历史目录与无引用模块
- `src/Sampling/Old/`（整目录）
- `src/Sampling/Source/rjdnn.zip`
- `src/logger.py`
- `src/logger_old.py`
- `src/card/possion.json`

2. 同步配置
- `src/card/preference.json`
- 删除失活配置键 `PoissonSchema`

3. 打包策略收敛
- `pyproject.toml`
- 移除 package-data 中的 `Sampling/Old/**/*`

### 变更规模（基于索引统计）

- 删除文件数: 16
- 文本删除行数: 2559 行
- 删除二进制历史包:
- `src/Sampling/Old/dynesty-a24...zip`（29,104,709 bytes）
- `src/Sampling/Source/rjdnn.zip`（13,652,498 bytes）
- 仅两份 zip 合计减少约 40.8 MB（不含其余源码删除体积）

### 验证结果

1. 引用与存在性检查
- `src/Sampling/Old`: 已移除
- `src/Sampling/Source/rjdnn.zip`: 已移除
- `src/logger.py`: 已移除
- `src/logger_old.py`: 已移除
- `src/card/possion.json`: 已移除

2. 回归测试
- 命令: `python -m unittest discover -s tests -p 'test_*.py'`
- 结果: `Ran 4 tests ... OK`

3. 语法/导入烟测
- 核心模块 `py_compile` 通过
- 关键模块 import smoke 通过

### 风险评估

- 本轮删除项均不在 `client -> core -> distributor -> Sampling/*` 主链路中。
- 属于“安全瘦身”，未触及核心采样、工作流调度、likelihood 计算语义。
