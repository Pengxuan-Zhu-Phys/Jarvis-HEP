# Jarvis-HEP V2

Jarvis-HEP V2 is a distributed runtime for high-energy-physics parameter scans.
You describe a scan as a validated YAML task card; Jarvis generates parameter
points, evaluates external calculators and/or Python operators in parallel, and
archives observables, sample artifacts, logs, and run summaries.

The PyPI distribution is named `Jarvis-HEP` (current version `2.0.14`). Its
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
python3 -m pip install Jarvis-HEP
```

If the short-lived `jarvishep2` distribution was installed previously, remove
it before installing V2 from `Jarvis-HEP`; both distributions provide the same
`jarvishep2` import package and should not be installed together:

```bash
python3 -m pip uninstall jarvishep2
python3 -m pip install --upgrade Jarvis-HEP
```

From a source checkout:

```bash
python3 -m pip install -e .
```

For development and tests, install the test tools separately:

```bash
python3 -m pip install -e .
python3 -m pip install "pytest>=7.0" "fakeredis>=2.0" "colorlog>=6.0"
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

The default install includes every Jarvis runtime dependency: the Redis Python
client, `msgpack`, `aiofiles`, Textual Monitor, `Jarvis-HEP-Portal`,
`Jarvis-Operas`, and [JarvisPLOT](https://github.com/Pengxuan-Zhu-Phys/JarvisPLOT).
Jarvis does not publish optional dependency groups.

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

## Writing a task card

A task card is the YAML file that describes one scan. Every Jarvis command
(`validate`, `check`, `run`, `man`) reads the same set of fields:

- Every card needs three top-level sections: `Scan` (the scan name), `Sampling`
  (how parameter points are chosen), and `EnvReqs` (runtime settings).
- A card also needs at least one of `Calculators` (external programs) or
  `Operas` (Python functions) to compute results for each point.
- Settings for the chosen sampling method go under `Sampling.Bounds`, written
  in lowercase with underscores, e.g. `point_number`.
- Runtime settings go under `EnvReqs.V2`, for example the number of parallel
  `workers` and `checkpoint.heartbeat` (how often, in seconds, progress is
  saved so an interrupted scan can be resumed; minimum 30).
- To test your calculators on a few points before a full scan, use
  `Jarvis check TASK.yaml`.

A minimal card looks like this:

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

Jarvis has a built-in manual for every field, so you don't need to guess names:

```bash
Jarvis man                         # interactive YAML authoring guide
Jarvis man --json                  # structured output for tooling and agents
Jarvis man sampler                 # sampler catalog
Jarvis man sampler.ToyMCMC         # reference independent-chain MCMC Bounds
Jarvis man sampler.mcmc-runtime    # MCMC multi-chain Redis pipeline design
Jarvis man yaml.EnvReqs.V2         # one YAML section
Jarvis man calculator.execution.output --type JSON
```

`Jarvis validate` checks a card for mistakes without starting a scan (and
without needing Redis). Run it before every `Jarvis run`.

## CLI workflow

```bash
Jarvis validate TASK.yaml          # validate only
Jarvis check TASK.yaml             # quick test of calculators on a few points
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
Jarvis project browse
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

For a scan named `my_scan`, the project folder typically contains:

```text
outputs/my_scan/
├── DATABASE/
│   ├── samples.hdf5              # results for every evaluated point
│   └── samples.csv               # same results as CSV, written when the scan finishes
├── SAMPLE/
│   ├── 000001/<uuid>/             # files produced for each point (if kept)
│   └── 000001.tar.gz              # the same files, compressed in groups (if enabled)
├── run_summary.json              # counts, timing and speed of the run
├── run_summary.csv
└── run_summary.txt

logs/my_scan/
├── core.log                      # main process
├── sampler.log                   # choosing parameter points
├── archiver.log                  # writing results to disk
└── worker-00.log ...             # one log per parallel worker
```

If a scan was interrupted, or you want to regenerate `samples.csv`, run
`Jarvis convert TASK.yaml`.

Progress for resuming is saved under
`checkpoints/<scan>/<sampler>/state.pkl`. Press `Ctrl+C` to stop a scan cleanly;
Jarvis stops all of its worker processes and any Redis server it started
itself. Continue later with `Jarvis run TASK.yaml --resume`.

**Saving disk space.** By default Jarvis keeps a folder of files for every
point under `SAMPLE/`. For large scans where you only need the numbers in
`DATABASE/`, set:

```yaml
EnvReqs:
  V2:
    store_samples: false
```

Calculators still run normally, but their files are written to a temporary
folder and deleted after each point. Results, error status for failed points,
resume data and the main logs are kept; per-point files and per-point logs
are not. See [SAMPLE storage](docs/sample-storage.md) for details.

## Documentation

- [Online documentation](https://pengxuan-zhu-phys.github.io/Jarvis-Docs/)
- `Jarvis man` — built-in manual for every task-card field, with examples
- [INSTALL.md](INSTALL.md) — installation, Redis setup, full command list, and
  sharing projects
- [SAMPLE storage](docs/sample-storage.md) — keeping or discarding per-point files
- [Project template](jarvishep2/project_template/README.md) — what
  `Jarvis project create` puts in a new project

Coming from Jarvis-HEP 1.x? Version 2 replaces it and 1.x is no longer
developed. The command is still called `Jarvis`; check older task cards with
`Jarvis validate` before running them.

## Citation

If you use Jarvis-HEP in your research, please cite:

```bibtex
@article{Guo:2026jarvishep,
  title         = {Jarvis-HEP: A lightweight Python framework for workflow
                   composition and parameter scans in high-energy physics},
  author        = {Guo, Erdong and Jackson, Paul and Yang, Jin-Min and Zhu, Pengxuan},
  year          = {2026},
  eprint        = {2604.25557},
  archivePrefix = {arXiv}
}
```

Please also cite the sampling algorithms you use. `Jarvis --refs` prints the
references for Jarvis-HEP and for each built-in sampler.

## Development

```bash
python3 -m pip install -e .
python3 -m pip install "pytest>=7.0" "fakeredis>=2.0" "colorlog>=6.0"
python3 -m pytest -q
```

The long-running AdaptiveBridson integration tests are skipped by default;
run them explicitly with `python3 -m pytest -q tests/test_adaptive_bridson.py`
when changing that sampler.

The long-running feedback-loop coverage in `tests/test_ensemble_samplers.py`
is also skipped by default; run it explicitly for EnsembleMCMC/DEMCMC/PT changes:

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

## License

Jarvis-HEP is released under the [MIT License](LICENSE). Bundled third-party
code keeps its own license (for example
[Dynesty](jarvishep2/sampling/Source/Dynesty/LICENSE)).
