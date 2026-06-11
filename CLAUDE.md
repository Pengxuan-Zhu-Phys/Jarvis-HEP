# CLAUDE.md — Jarvis-HEP Developer Context

## Source of Truth

| Layer | Authority | Role |
|-------|-----------|------|
| **Code** | Ground truth | Determines actual behavior |
| **CLAUDE.md** | Operational summary | Accurate snapshot for AI-assisted development |
| **docs/** | Explanatory layer | Design intent; may lag implementation |

When CLAUDE.md conflicts with code, trust the code. When docs/ conflicts with code, trust the code and update docs/.

---

## Project Overview

Jarvis-HEP is a YAML-driven orchestration framework for likelihood-based High Energy Physics (HEP) parameter space scans. It coordinates expensive external calculators, explores parameter spaces with various sampling algorithms, and persists structured outputs in HDF5/CSV format.

- **Version**: 1.6.11
- **Python**: 3.10+
- **License**: MIT
- **Entry point**: `Jarvis` CLI → `jarvishep.client:main`

---

## Repository Layout

```
jarvishep/
├── client.py              # CLI dispatch (Jarvis command)
├── core.py                # Main runtime orchestrator
├── config.py              # YAML config loader & validator
├── runtime_config.py      # Optional Runtime: block normalization (WP-1.1)
├── sample_logger.py       # SampleLogger (file sink) + BufferedSampleLogger (WP-1.2)
├── workflow.py            # Task DAG management, semantic flowchart JSON export
├── factory.py             # WorkerFactory, single scheduler thread pool (WP-1.3)
├── benchmark.py           # Additive throughput benchmark mode helpers
├── distributor.py         # Sampler dispatch (match/case)
├── sample.py              # Sample data structure
├── base.py                # Base class with path management (&J/ resolution)
├── moduleManager.py       # Module lifecycle (singleton)
├── modulePool.py          # Module pooling & instance management
├── hdf5writer.py          # HDF5 storage backend (threaded writer)
├── observable_io.py       # Observable schema I/O, CSV conversion
├── io_manager.py          # Async blocking-work executor (ThreadPoolExecutor wrapper)
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
│   ├── sampler.py         # SamplingVirtial base class (class name has typo, do not rename)
│   ├── variables.py       # Variable/distribution definitions
│   ├── randoms.py         # Random sampler
│   ├── grid.py            # Grid sampler
│   ├── csv_sampler.py     # CSV replay sampler
│   ├── bridson.py         # Bridson sampling
│   ├── dnn.py             # DNN-assisted iterative sampling
│   ├── dynesty.py         # Dynesty nested sampling
│   ├── multinest.py       # MultiNest nested sampling
│   ├── diver.py           # Diver optimizer
│   ├── mcmc_standard.py   # Standard MCMC (also registered as "ToyMCMC")
│   ├── tpmcmc.py          # Parallel Tempering MCMC
│   ├── ammcmc.py          # Adaptive Metropolis MCMC
│   ├── robustam.py        # Robust Adaptive Metropolis
│   ├── dram.py            # Delayed Rejection Adaptive Metropolis
│   ├── dream.py           # DREAM (DE-MCMC variant)
│   ├── dream_lite.py      # DREAM Lite (also registered as "DREAM-lite")
│   ├── demcmc.py          # Differential Evolution MCMC
│   ├── ensemblemcmc.py    # Ensemble MCMC
│   ├── pt_ensemble.py     # Parallel Tempering Ensemble
│   ├── slicemcmc.py       # Slice MCMC
│   ├── ess.py             # Elliptical Slice Sampler
│   ├── mala.py            # MALA (gradient)
│   ├── hmc.py             # HMC (gradient)
│   ├── nuts.py            # NUTS (gradient)
│   ├── rltpmcmc.py        # RL-TPMCMC (experimental)
│   └── Source/MCMC/       # MCMC engine implementations
│       ├── state_machine_base.py           # Base state machine (945 lines)
│       ├── state_machine_multistage_base.py
│       ├── runtime_checkpoint.py           # Checkpoint/resume infrastructure
│       ├── mcmc_chain.py                   # DEAD CODE — not imported anywhere
│       ├── engine_standard.py              # Standard MH engine
│       ├── engine_dream.py                 # DREAM engine
│       ├── engine_dream_lite.py            # DREAM Lite engine
│       ├── engine_demcmc.py                # DEMCMC engine
│       ├── engine_hmc.py                   # PLACEHOLDER — gradient_contract_level = "placeholder"
│       ├── engine_mala.py                  # PLACEHOLDER — gradient_contract_level = "placeholder"
│       ├── engine_nuts.py                  # PLACEHOLDER — gradient_contract_level = "placeholder"
│       └── engine_contract.py              # MCMCEngineProtocol
├── IOs/                   # Input/Output parameter system
├── monitoring/            # Run monitoring & run_summary generation
├── card/                  # Config files, schemas, templates, logos
└── project_template/      # Default project scaffold (bin/, data/, deps/)
tests/                     # Test suite
```

---

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

### Initialization sequence (Scan / 1PC modes)

1. `init_argparser()` — parse CLI args
2. `init_project()` — set up directory structure
3. `init_logger()` — configure logging
4. `init_configparser()` — load YAML, validate
5. `init_utils()` — load custom functions
6. `_preload_resume_checkpoint()` — check for resumable state
7. `init_sampler()` — instantiate sampling method
8. `init_StateSaver()` — set up checkpointing
9. `init_workflow()` — build DAG, export flowchart.json, optionally render PNG via JarvisPLOT
10. `init_librarys()` — install required libraries
11. `init_WorkerFactory()` — create thread pool
12. `init_likelihood()` — set up likelihood (if configured)
13. `init_database()` — open HDF5 writer
14. `init_run_summary()` — start metrics collection

---

## Key Conventions

- **Path marker** `&J/` resolves to standalone project root (identified by `.jarvis-project.json` or `jarvis.project.yaml`)
- **Runtime tokens**: `@SampleID` (sample UUID), `@Sdir` (sample save dir), `@PackID` (calculator instance ID)
- **Checkpoint path**: `<task_root>/checkpoints/<scan_name>/<sampler>/state.pkl`
- **Output structure**: `outputs/<scan>/DATABASE/` (HDF5, CSV), `outputs/<scan>/SAMPLE/` (per-sample artifacts; `outputs/<scan>/SAMPLE/tests/` in `--check-modules`; `DATABASE` stays at `outputs/<scan>/DATABASE/`)
- **Calculator YAML key**: `make_paraller` (intentional legacy spelling — do not rename)
- **Calculator timeout**: `Calculators.Modules[].timeout` is an optional total-second limit for one sample `execution` section; it starts after `initialization`.
- **Runtime block** (optional top-level YAML; WP-1.1): `Runtime.mode` (`auto|thread|process`), `Runtime.workers`, `Runtime.batch_size`, `Runtime.sample_artifacts` (`auto|always|never`). Only `sample_artifacts` is consumed at M1: `auto` (default) skips per-sample `SAMPLE/<bucket>/<uuid>/` creation for opera-only scans; `always` restores v1 behavior; `never` suppresses artifact creation (including failure replay). `@Sdir` resolution or any calculator workflow triggers `Sample.materialize()`. Lazy samples still expose `sample_info["logger_name"]` (metadata only) so Likelihood/Calculator child loggers bind correctly without filesystem materialization.
- **Failure-replay logging** (WP-1.2): non-materialized samples route module/likelihood logs through `BufferedSampleLogger` (bounded in-memory buffer). On success the buffer is discarded; on failure `materialize_failure_artifacts()` replays buffered events to `Sample_running.log` using the frozen `LOGGING_CONTRACT.md` format, then appends a failure footer.
- **Thread-pool flattening** (WP-1.3): `WorkerFactory` uses one named executor (`jarvis-hep-factory`); the legacy `log_executor` is removed (status/error logs are synchronous). Calculator installation runs on the calling factory worker thread (no per-pool installer executor). `IOManager` is created only for calculator workflows; opera-only scans skip the `jarvis-hep-io` pool entirely.

---

## CLI Reference

```
Jarvis <file.yaml>                  # Run scan
Jarvis <file.yaml> --resume         # Resume from checkpoint without prompt
Jarvis <file.yaml> --plot           # Generate plot config
Jarvis <file.yaml> --convert        # Convert HDF5 to CSV
Jarvis <file.yaml> --monitor        # Real-time resource monitor
Jarvis <file.yaml> --check-modules  # Test-run 10 samples (1PC mode)
Jarvis <file.yaml> --benchmark [s]  # Run throughput benchmark mode (default 30 s)
Jarvis project create <name>        # Create scaffold
Jarvis project pack [path] [--share|--repro|--full]
Jarvis project pack [path] [--share|--repro|--full] --man
Jarvis project pack <pack-manifest.yaml>
Jarvis project browse               # List official projects
Jarvis project fetch <name>         # Download official project
Jarvis project info <name>          # Show project metadata
Jarvis -v                           # Print version
Jarvis --refs                       # Print full reference information
```

---

## Checkpoint / Resume System

- **Lock**: `_runtime_checkpoint_save_lock` is initialized in `SamplingVirtial.__init__()` at `sampler.py:59` and inherited by all MCMC subclasses via `MCMCStateMachineBase`. The lock is correct and functional.
- **Interval**: 30-second heartbeat; saves only at safe barriers
- **Atomic write**: `os.replace` on a temporary file; latest checkpoint is always at same `state.pkl`
- **Resume prompt** (when `--resume` not passed and checkpoint exists):
  ```
  Detected checkpoint file. Re-run from scratch? [y/N] (default: resume in 30s):
  ```
  `y`/`yes` = fresh run; Enter / timeout (30s) = resume

### Actual checkpoint payload structure (`runtime_checkpoint.py`)

```python
{
    "format": "jarvis-hep.statesaver",
    "version": 1,
    "created_at_utc": "...",
    "jarvishep_version": "...",
    "run_spec": {
        "raw_yaml_text": "...",
        "normalized_config": {...},
        "scan_name": "...",
        "task_root": "...",
        "task_result_dir": "...",
        "workflow_modules": [...],
        "module_configs": {...},
        "sampler_method": "...",
    },
    "factory_blueprint": {
        "max_workers": <int>,
        "workflow": {...},
        "module_pools": {...},
    },
    "sampler_state": {
        "format": "jarvis-hep.mcmc-runtime",
        "version": 1,
        "timestamp_utc": "...",
        "sampler_signature": {...},
        "numpy_random_state": <tuple>,
        "bucket_allocator": {...} | null,
        "state_machine": {
            "state": "<MCMCState value>",
            "ready_queue": [...],
            "control_state": {...},
            "chains": [...],
            "extras": {"pending_samples": [...]},
        },
    },
    "integrity": {
        "config_hash": "...",
        "sampler_signature": {...},
        "variable_signature": {...},
        "checkpoint_reason": "...",
        "safe_barrier_confirmed": True,
    },
}
```

Note: There is no `resume_policy` key in the payload. Resume behavior (always rebuild factory, never restore inflight, always use checkpoint config) is hardcoded in `core.py`, not stored in the file.

---

## Build & Test

```bash
python3 -m pip install -e ".[dev]"
pytest tests/
```

---

## Known Bugs (validated against v1.6.11 source)

### Fixed during V2 M0

- **`sample.py:140,142`** — Fixed `==` vs `=` in `Sample.start()` status transitions; covered by sample status tests.
- **`config.py:812`** — Fixed schema block access so config validation no longer raises `AttributeError`.
- **`distributor.py:46`** — Made `Distributor.set_method(...)` a static method so CLI sampler initialization works.

### Critical (P0) — will crash or produce wrong results

1. **`config.py:327`** — Catches `FileExistsError` but the error message says "not found"; should catch `FileNotFoundError`. A missing default env YAML silently passes.

2. **`calculator.py:136`** — `def custom_format(record)` is missing `self`. Defined inside the class body but not as an instance method signature; will `TypeError` when called via `self.custom_format(...)`.

3. **`sampler.py:447`** — `ax.axes("off")` should be `ax.axis("off")`. Matplotlib `Axes` objects have no `axes()` method; this will raise `AttributeError` when drawing the summary plot.

4. **`state_machine_base.py:189`** — `_export_runtime_extras()` is defined twice in the same class body. The first definition at line 189 returns `{}` and is dead code — it is unconditionally shadowed by the second definition at line 236 (which returns `{"pending_samples": ...}`). The second definition is the correct one. Remove the dead definition at line 189.

### High (P1) — silent logic errors

8. **`config.py:125`** — Bare `except:` catches `KeyboardInterrupt` and `SystemExit`, preventing clean shutdown.

9. **`config.py:233-241`** — `subprocess.run()` is called without `check=True`, so it never raises `CalledProcessError`. The `except subprocess.CalledProcessError` block is unreachable. ROOT detection always reports success (`ROOT = True`) regardless of whether ROOT is actually installed.

10. **`variables.py:48`** — `np.log(p / (1 - p))` where `p = np.random.uniform(0, 1)`. When `p=0` (possible, inclusive lower bound), this evaluates to `log(0) = -inf`. Returns `-inf` silently rather than raising.

11. **`dnn.py:244-246`** — `Classifier.save()` writes to `self.path` but `self.path` is not initialized in `Classifier.__init__`. Will raise `AttributeError` on save.

12. **`bridson.py:397`** — `np.sum(..., axis=ndim)` should be `axis=-1`. When `ndim >= array.ndim`, numpy raises `AxisError`.

13. **`sample_archive.py:179-180`** — Non-atomic check-then-add pattern; race condition under multithreading.

14. **`dnn.py:556`** — `self.classifier.train_model(self.dataset.v, self.dataset.valid)` — `valid` is likely not the intended target labels for this training call.

### Previously Reported — NOT BUGS (removed from bug list)

- ~~`state_machine_base.py:379` — `_runtime_checkpoint_save_lock` never initialized~~ — **INCORRECT**. The lock is correctly initialized in the parent class `SamplingVirtial.__init__()` at `sampler.py:59` and inherited by all MCMC subclasses.

- ~~`variables.py:91` — `scipy.special.logit()` wrong; should use logistic PPF~~ — **INCORRECT**. `scipy.special.logit` IS the logistic PPF (quantile function). `logit(p) = log(p/(1-p))` is the correct inverse of the logistic CDF. The implementation is mathematically correct.

### Spelling / Naming

- `SamplingVirtial` → `SamplingVirtual` (class name in `sampler.py:27`) — do not rename; it would break all subclass imports
- `"initializaing"` → `"initializing"` (log strings in `mcmc_standard.py`, `tpmcmc.py`, `ammcmc.py`, `dram.py`, `ensemblemcmc.py`, `slicemcmc.py`, `ess.py`, `robustam.py`, `pt_ensemble.py`)
- `"Iegal"` → `"Illegal"`, `"founded"` → `"found"` (`config.py:126`)
- `cunsom_constants` → `custom_constants` (`sample.py:130`)

---

## Structural Design Issues

### God Classes (split recommended)

- **`CalculatorModule`** (~420 lines): mixes config parsing, logging, command execution, I/O, path resolution, instance tracking
- **`ModuleManager`** (~270 lines): mixes module registry, workflow orchestration, likelihood calculation, nuisance handling, config propagation
- **`MCMCStateMachineBase`** (945 lines): mixes state machine logic, checkpoint/resume, metrics, controller interface, chain validation

### Inheritance Problems

- `Module` base class does not use ABC/`@abstractmethod`; subclasses have inconsistent `execute()` signatures: `CalculatorModule.execute(input_data, sample_info)` vs `OperasModule.execute(observables, sample_info)`
- `Module` base class has 7 unused imports: `uuid`, `pprint`, `yaml`, `logging`, `subprocess`, `sleep`, `Parameter`
- `MCMCMultiStageStateMachineBase` uses explicit `MCMCStateMachineBase.__init__(self)` instead of `super().__init__()`, breaking MRO for future multiple inheritance

### Configuration Propagation

Settings cascade manually through 4 layers: `Core → ModuleManager → ModulePool → CalculatorModule`, each with explicit `for` loops. Adding a new config field requires changes in all layers.

### Dead Code

- `Sampling/Source/MCMC/mcmc_chain.py` — entire file unused; not imported anywhere; contains an infinite loop bug
- `engine_hmc.py`, `engine_mala.py`, `engine_nuts.py` — placeholder stubs with `gradient_contract_level = "placeholder"`; no real implementation
- `calculator.py:134` — `analyze_config_multi()` body is `pass`; never called
- `bridson.py:5`, `grid.py:5` — `from re import S` unused
- `bridson.py:7`, `grid.py` — `from mpmath.functions.functions import re` shadows stdlib `re`

---

## Fix Priority

| Priority | Scope | Effort |
|----------|-------|--------|
| P0 | Critical bugs #1–7 (crashes, wrong results) | ~2 hours |
| P1 | High bugs #8–14 (silent logic errors) | ~3 hours |
| P2 | Remove dead code (`mcmc_chain.py`, unused imports, placeholder engines) | ~1 hour |
| P3 | Refactor god classes, extract responsibilities | ~2–3 days |
| P4 | Fix spelling in log strings, variable names | ~1 hour |

---

## Undocumented Features (not in user-facing docs)

- `Sampling.ModuleFailurePolicy`: config key accepted by `ModuleManager`; values `"continue"` (default) or `"fail-fast"`.
- `Operas` module `call_mode: "acall"`: enables async operator dispatch. Not in docs.
- `ToyMCMC`: alias for `MCMC` in `distributor.py`. Not mentioned in README sampler list.
- `--resume` CLI flag triggers `checkpoint_policy = "resume"` in `core.py`. The complementary `"fresh"` policy (answering `y` to the resume prompt) is not exposed as a CLI flag.
