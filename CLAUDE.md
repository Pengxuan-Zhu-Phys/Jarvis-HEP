# CLAUDE.md — Jarvis-HEP Developer Context

## Project Overview

Jarvis-HEP is a YAML-driven orchestration framework for likelihood-based High Energy Physics (HEP) parameter space scans. It coordinates expensive external calculators, explores parameter spaces with various sampling algorithms, and persists structured outputs in HDF5/CSV format.

- **Version**: 1.6.11
- **Python**: 3.10+
- **License**: MIT
- **Entry point**: `Jarvis` CLI → `jarvishep.client:main`

## Repository Layout

```
jarvishep/
├── client.py              # CLI dispatch (Jarvis command)
├── core.py                # Main runtime orchestrator
├── config.py              # YAML config loader & validator
├── workflow.py            # Task DAG management, flowchart generation
├── factory.py             # WorkerFactory, thread pool management
├── distributor.py         # Sampler dispatch (match/case)
├── sample.py              # Sample data structure
├── base.py                # Base class with path management (&J/ resolution)
├── moduleManager.py       # Module lifecycle (singleton)
├── modulePool.py          # Module pooling & instance management
├── hdf5writer.py          # HDF5 storage backend (threaded writer)
├── observable_io.py       # Observable schema I/O, CSV conversion
├── io_manager.py          # I/O orchestration
├── async_subprocess.py    # Async subprocess scheduling
├── monitor.py             # Resource monitoring (CPU/memory)
├── plot.py                # Plotting engine
├── utils.py               # General utilities
├── project_scaffold.py    # `Jarvis project create` scaffolding
├── project_packager.py    # `Jarvis project pack` packaging
├── official_project_library.py  # Official project library manager
├── versioning.py          # Version branding
├── Module/
│   ├── module.py          # Base Module class
│   ├── calculator.py      # CalculatorModule (external program execution)
│   ├── operas.py          # OperasModule (Jarvis-Operas integration)
│   ├── likelihood.py      # Likelihood module
│   └── parameters.py      # Parameter management
├── Sampling/
│   ├── sampler.py         # SamplingVirtial base class (typo in name)
│   ├── variables.py       # Variable/distribution definitions
│   ├── randoms.py         # Random sampler
│   ├── grid.py            # Grid sampler
│   ├── csv_sampler.py     # CSV replay sampler
│   ├── bridson.py         # Bridson sampling
│   ├── dnn.py             # DNN-assisted iterative sampling
│   ├── dynesty.py         # Dynesty nested sampling
│   ├── multinest.py       # MultiNest nested sampling
│   ├── diver.py           # Diver optimizer
│   ├── mcmc_standard.py   # Standard MCMC
│   ├── tpmcmc.py          # Parallel Tempering MCMC
│   ├── ammcmc.py          # Adaptive Metropolis MCMC
│   ├── robustam.py        # Robust Adaptive Metropolis
│   ├── dram.py            # Delayed Rejection Adaptive Metropolis
│   ├── dream.py           # DREAM (DE-MCMC variant)
│   ├── dream_lite.py      # DREAM Lite
│   ├── ensemblemcmc.py    # Ensemble MCMC
│   ├── pt_ensemble.py     # Parallel Tempering Ensemble
│   ├── slicemcmc.py       # Slice MCMC
│   ├── ess.py             # Elliptical Slice Sampler
│   ├── mala.py            # MALA (gradient)
│   ├── hmc.py             # HMC (gradient)
│   ├── nuts.py            # NUTS (gradient)
│   ├── rltpmcmc.py        # RL-TPMCMC (experimental)
│   └── Source/MCMC/       # MCMC engine implementations
│       ├── state_machine_base.py           # Base state machine (~945 lines)
│       ├── state_machine_multistage_base.py
│       ├── runtime_checkpoint.py           # Checkpoint/resume infrastructure
│       ├── mcmc_chain.py                   # DEAD CODE - not used anywhere
│       ├── engine_standard.py              # Standard MH engine
│       ├── engine_hmc.py                   # PLACEHOLDER - not implemented
│       ├── engine_mala.py                  # PLACEHOLDER - not implemented
│       ├── engine_nuts.py                  # PLACEHOLDER - not implemented
│       └── engine_contract.py              # MCMCEngineProtocol
├── IOs/                   # Input/Output parameter system
├── monitoring/            # Run monitoring & run_summary generation
├── card/                  # Config files, schemas, templates, logos
└── project_template/      # Default project scaffold (bin/, data/, deps/)
tests/                     # 35 test files
```

## Execution Flow

```
client.py (CLI) → core.py (Core)
  → config.py (load YAML, validate against schema)
  → distributor.py (instantiate sampler by name)
  → workflow.py (build calculator DAG)
  → factory.py (WorkerFactory thread pool)
    → moduleManager.py → modulePool.py → calculator.py / operas.py
    → hdf5writer.py (threaded HDF5 writes)
  → monitoring/run_summary.py (end-of-run diagnostics)
```

## Key Conventions

- **Path marker** `&J/` resolves to standalone project root
- **Runtime tokens**: `@SampleID` (sample UUID), `@Sdir` (sample save dir), `@PackID` (calculator instance ID)
- **Checkpoint path**: `<task_root>/checkpoints/<scan_name>/<sampler>/state.pkl`
- **Output structure**: `outputs/<scan>/DATABASE/` (HDF5, CSV), `outputs/<scan>/SAMPLE/` (per-sample artifacts)

## Build & Test

```bash
python3 -m pip install -e ".[dev]"
pytest tests/
```

## Known Bugs (as of 2026-03-30 review)

### Critical (P0) — will crash or produce wrong results

1. **`sample.py:140,142`** — `self.info['status'] == "Running"` uses `==` (comparison) instead of `=` (assignment). Sample status is never set to "Running".

2. **`config.py:812`** — `self.schemablock` is undefined; should be `schema_block` (local variable). Config validation will `AttributeError`.

3. **`config.py:327`** — Catches `FileExistsError` but message says "not found"; should catch `FileNotFoundError`.

4. **`distributor.py:46`** — `def set_method(method)` missing `self` parameter. Every call will fail with wrong argument count.

5. **`calculator.py:136`** — `def custom_format(record)` missing `self` parameter; will `TypeError` when called as instance method.

6. **`sampler.py:415`** — `ax.axes("off")` should be `ax.axis("off")`. Matplotlib has no `axes()` method on Axes objects.

7. **`state_machine_base.py:379`** — `self._runtime_checkpoint_save_lock` is never initialized anywhere. All checkpoint saves will `AttributeError`. This breaks the entire resume feature.

8. **`state_machine_base.py:189,236`** — `_export_runtime_extras()` defined twice; second definition silently shadows the first.

### High (P1) — silent logic errors

9. **`config.py:125`** — Bare `except:` catches `KeyboardInterrupt`/`SystemExit`, prevents clean shutdown.

10. **`config.py:233-241`** — `subprocess.run()` without `check=True` never raises `CalledProcessError`. ROOT detection always reports True.

11. **`variables.py:91`** — Logit distribution mapping uses `scipy.special.logit()` but should use the inverse (logistic PPF). Math is wrong.

12. **`variables.py:48`** — `np.log(p / (1 - p))` will crash when p=0 or p=1 (division by zero / log of zero).

13. **`dnn.py:244-246`** — `Classifier.save()` uses `self.path` but it's never initialized in `__init__`.

14. **`bridson.py:397`** — `np.sum(..., axis=ndim)` should be `axis=-1`. Will `IndexError` when ndim >= array dimensions.

15. **`sample_archive.py:179-180`** — Non-atomic check-then-add pattern; race condition under multithreading.

16. **`dnn.py:556`** — `self.classifier.train_model(self.dataset.v, self.dataset.valid)` — `valid` is likely not the correct target labels.

### Spelling / Naming

- `SamplingVirtial` → `SamplingVirtual` (sampler.py:27)
- `"initializaing"` → `"initializing"` (appears in mcmc_standard.py, tpmcmc.py, ammcmc.py, dram.py, ensemblemcmc.py, slicemcmc.py, ess.py, robustam.py, pt_ensemble.py)
- `"Iegal"` → `"Illegal"`, `"founded"` → `"found"` (config.py:126)
- `cunsom_constants` → `custom_constants` (sample.py:130)

## Structural Design Issues

### God Classes (split recommended)

- **`CalculatorModule`** (~420 lines): mixes config parsing, logging, command execution, I/O, path resolution, index tracking
- **`ModuleManager`** (~270 lines): mixes module registry, workflow orchestration, likelihood calculation, nuisance handling, config propagation
- **`MCMCStateMachineBase`** (~945 lines): mixes state machine logic, checkpoint/resume, metrics, controller interface, chain validation

### Inheritance Problems

- `Module` base class does not use ABC/`@abstractmethod`; subclasses have inconsistent `execute()` signatures (`CalculatorModule.execute(input_data, sample_info)` vs `OperasModule.execute(observables, sample_info)`)
- `Module` base class has 7 unused imports (uuid, pprint, yaml, logging, subprocess, sleep, Parameter)
- `MCMCMultiStageStateMachineBase` uses explicit `MCMCStateMachineBase.__init__(self)` instead of `super().__init__()`, breaking MRO

### Configuration Propagation

Settings cascade manually through 3 layers: `Core → ModuleManager → ModulePool → CalculatorModule`, each with explicit `for` loops. Adding a new config field requires changes in all layers.

### Dead Code

- `Sampling/Source/MCMC/mcmc_chain.py` — entire file unused, contains infinite loop bug
- `engine_hmc.py`, `engine_mala.py`, `engine_nuts.py` — placeholder stubs with `gradient_contract_level = "placeholder"`
- `calculator.py:134` — `analyze_config_multi()` is `pass`, never called
- `bridson.py:5`, `grid.py:5` — `from re import S` unused
- `bridson.py:7`, `grid.py` — `from mpmath.functions.functions import re` shadows stdlib `re`

## Fix Priority

| Priority | Scope | Effort |
|----------|-------|--------|
| P0 | Critical bugs #1-8 (crashes, data corruption) | ~2 hours |
| P1 | High bugs #9-16 (silent logic errors) | ~3 hours |
| P2 | Remove dead code (mcmc_chain.py, unused imports, placeholder engines) | ~1 hour |
| P3 | Refactor god classes, extract responsibilities | ~2-3 days |
| P4 | Unify naming conventions, fix spelling | ~1 hour |
