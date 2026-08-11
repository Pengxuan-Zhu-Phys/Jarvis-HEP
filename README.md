# Jarvis-HEP V2

Jarvis-HEP V2 is a distributed runtime for high-energy-physics parameter scans.
You describe a scan as a validated YAML task card; Jarvis generates parameter
points, evaluates external calculators and/or Python operators in parallel, and
archives observables, sample artifacts, logs, and run summaries.

The PyPI distribution is named `Jarvis-HEP` (current version `2.0.0`). Its
Python import package remains `jarvishep2`, and it exposes a single user-facing
command: `Jarvis`.

## What V2 provides

| Capability | V2 implementation |
| --- | --- |
| Distributed execution | Redis-backed queue with Python `spawn` Workers |
| Parameter sampling | Stable catalog: `Random`, `Grid`, `Bridson`, `AdaptiveBridson`, `CSV`, `Dynesty`, and `MultiNest` |
| HEP calculator integration | [Jarvis-HEP-Portal](https://github.com/Pengxuan-Zhu-Phys/Jarvis-Portal) for calculator file I/O |
| Python operators | [Jarvis-Operas](https://github.com/Pengxuan-Zhu-Phys/Jarvis-Operas) for registered operators and expressions |
| Reliability | Checkpoints, resume, graceful shutdown, and per-sample failure artifacts |
| Observability | Validation diagnostics, monitor snapshots, component logs, and performance summaries |
| Project workflow | Scaffold, validate, run, package, fetch, and inspect standalone projects |

The execution model is deliberately simple:

```text
Task YAML → sampler → Redis → Workers → Archiver
                                      ├─ DATABASE/samples.hdf5
                                      ├─ SAMPLE/<bucket>/<uuid>/
                                      └─ logs/<scan>/ + run_summary.*
```

## Installation

Requirements:

- Python 3.10 or newer
- Redis for distributed runs; V2 uses the local `127.0.0.1:6379` service

From PyPI:

```bash
python3 -m pip install 'Jarvis-HEP[distributed]'
```

If the short-lived `jarvishep2` distribution was installed previously, remove
it before installing V2 from `Jarvis-HEP`; both distributions provide the same
`jarvishep2` import package and should not be installed together:

```bash
python3 -m pip uninstall jarvishep2
python3 -m pip install --upgrade 'Jarvis-HEP[distributed]'
```

From a source checkout:

```bash
python3 -m pip install -e '.[distributed]'
```

For development and tests:

```bash
python3 -m pip install -e '.[distributed,dev]'
```

Start Redis before a real scan. For example, on macOS:

```bash
brew install redis
brew services start redis
redis-cli ping       # PONG
```

On the first `Jarvis` command after installation, Jarvis checks whether a
Redis-compatible server executable (`redis-server`, `redis6-server`, or
`valkey-server`) is available. If it is missing, the command's normal output
is followed by a one-time, OS-specific installation hint (written to stderr so
JSON/stdout output remains usable). The hint is advisory only; Jarvis does not
install operating-system packages automatically. Linux hints cover
Debian/Ubuntu, Fedora/RHEL, Amazon Linux, Arch, openSUSE, Alpine, Gentoo, Void,
NixOS and other package-manager families. The check marker is stored at
`~/.jarvis/redis-install-check-v1`.

The default install includes the current compatible releases of
`Jarvis-HEP-Portal`, `Jarvis-Operas`, and `Jarvis-PLOT`. The optional extras
are:

- `distributed`: Redis, `msgpack`, and `aiofiles`
- `plot`: compatibility alias for the default [Jarvis-PLOT](https://github.com/Pengxuan-Zhu-Phys/JarvisPLOT) dependency
- `dev`: pytest, fakeredis, and colorlog

For the complete installation guide, Redis options, and project-packaging
workflow, see [INSTALL.md](INSTALL.md).

## First run

The project scaffold includes a runnable Bridson + Operas example that does not
need an external calculator:

```bash
Jarvis project create MyScan
cd MyScan

Jarvis validate bin/quickstart_bridson_operas.yaml
Jarvis run bin/quickstart_bridson_operas.yaml
```

The scaffold also contains calculator, CSV, Dynesty, and MultiNest examples.
List the built-in examples with:

```bash
Jarvis man example
```

## Task-card contract

V2 task cards use one closed, consistent vocabulary across `validate`, `check`,
`run`, and `man`:

- Required top-level blocks: `Scan`, `Sampling`, and `EnvReqs`.
- Declare at least one execution backend: `Calculators` or `Operas`.
- Put method-specific sampler settings under `Sampling.Bounds` and use lower
  `snake_case` keys.
- Put V2 runtime settings such as worker count and checkpoint heartbeat under
  `EnvReqs.V2`.
- Use `Jarvis check TASK.yaml` for the fixed-point calculator smoke test.
- Outputs are written to the project output tree; `cleanup.strategy` and
  `archiver.handoff` are not V2 task-card interfaces.

A card starts with this shape:

```yaml
Scan:
  name: my_scan

Sampling:
  Method: Random
  Bounds:
    point_number: 100
    seed: 7
  Variables:
    - name: x
      distribution:
        type: Flat
        parameters: {min: 0.0, max: 1.0}

EnvReqs:
  V2:
    workers: 2

Operas:
  Modules:
    - name: my_operator
      operator: my_package.my_function
```

Use the generated manuals instead of guessing field names:

```bash
Jarvis man                         # interactive YAML authoring guide
Jarvis man --json                  # structured output for tooling and agents
Jarvis man sampler                 # sampler catalog
Jarvis man yaml.EnvReqs.V2         # one YAML section
Jarvis man calculator.execution.output --type JSON
```

`Jarvis validate` performs schema and semantic checks without starting Redis or
Workers. It is the recommended gate before `run`.

## CLI workflow

```bash
Jarvis validate TASK.yaml          # validate only
Jarvis check TASK.yaml             # fixed-point calculator smoke test
Jarvis run TASK.yaml               # execute a distributed scan
Jarvis run TASK.yaml --resume      # resume from a checkpoint
Jarvis convert TASK.yaml           # refresh DATABASE/samples.csv from HDF5
Jarvis monitor                     # read one live monitor snapshot
Jarvis ps                          # list Jarvis process groups
Jarvis kill --yes                  # terminate selected runtime processes
Jarvis --refs                      # print framework and sampler references
```

Additional integrations are available through native subcommands:

```bash
Jarvis project create MyScan
Jarvis project list
Jarvis project fetch Eggbox
Jarvis portal man
Jarvis operas list
Jarvis gen-plot-yaml TASK.yaml
Jarvis plot path/to/plot.yaml      # Jarvis-PLOT is installed by default
```

Use `Jarvis COMMAND -h` for command-specific options. Screen logging defaults
to `WARNING`; use `--console-level INFO`, `--debug`, or `--silence` as needed.
File logs remain available under `logs/<scan>/`.

## Outputs

For a scan named `my_scan`, the project root typically contains:

```text
outputs/my_scan/
├── DATABASE/
│   ├── samples.hdf5              # archived observables
│   └── samples.csv               # CSV snapshot after `Jarvis convert`
├── SAMPLE/
│   ├── 000001/<uuid>/             # per-sample files when retained
│   └── 000001.tar.gz              # packed sample bucket when enabled
├── run_summary.json
├── run_summary.csv
└── run_summary.txt

logs/my_scan/
├── core.log
├── sampler.log
├── archiver.log
└── worker-00.log ...
```

Checkpoints are stored under
`checkpoints/<scan>/<sampler>/state.pkl`. Press `Ctrl+C` to stop a scan cleanly;
Jarvis shuts down its Workers, Archiver, and any Redis process it manages.

## Development

```bash
python3 -m pip install -e '.[distributed,dev]'
python3 -m pytest -q
```

The long-running AdaptiveBridson integration tests are skipped by default;
run them explicitly with `python3 -m pytest -q tests/test_adaptive_bridson.py`
when changing that sampler.

The long-running feedback-loop coverage in `tests/test_ensemble_samplers.py`
is also skipped by default; run it explicitly for Ensemble/DEMCMC/PT changes:

```bash
python3 -m pytest -q tests/test_ensemble_samplers.py
```

The slow distributed acceptance gates in `tests/test_distributed_acceptance.py`
are also skipped by default; run them explicitly when validating Worker/
Archiver performance:

```bash
python3 -m pytest -q tests/test_distributed_acceptance.py
```

`tests/test_distributed_resume.py`, `tests/test_mcmc_sampler.py`, and
`tests/test_worker_pool.py` are also excluded from the default suite. They
exercise interruption/checkpoint resume, multi-process sampler execution, and
Worker calculator-pool concurrency; run the relevant file explicitly after
changing those paths:

```bash
python3 -m pytest -q tests/test_distributed_resume.py
python3 -m pytest -q tests/test_mcmc_sampler.py
python3 -m pytest -q tests/test_worker_pool.py
```

`tests/test_variable_distributions.py` and `tests/test_worker_failure.py` are
also excluded from the default suite while their V1-card/schema fixture and
calculator-pool SIGKILL fixture are reconciled. Run either explicitly when
working on those paths:

```bash
python3 -m pytest -q tests/test_variable_distributions.py
python3 -m pytest -q tests/test_worker_failure.py
```

Useful repository documentation:

- [INSTALL.md](INSTALL.md) — installation, CLI details, Redis, projects, and packaging
- [Task-card schema](docs/task-card-schema.md) — validation layers and V2 configuration rules
- [Validation diagnostics](docs/validation-diagnostics.md) — actionable `JV2-*` error codes
- [Project template](jarvishep2/project_template/README.md) — standalone project layout

V1 (`jarvishep`) is frozen and CLI-retired. V2 uses `jarvishep2` and the
`Jarvis` command; `Jarvis2` is not installed.
