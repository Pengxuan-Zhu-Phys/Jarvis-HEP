# Jarvis-HEP V2

Jarvis-HEP V2 is a distributed runtime for high-energy-physics parameter scans.
You describe a scan as a validated YAML task card; Jarvis generates parameter
points, evaluates external calculators and/or Python operators in parallel, and
archives observables, sample artifacts, logs, and run summaries.

The distribution is named `Jarvis-HEP` (current source version `2.1.1`). Its
Python import package remains `jarvishep2`, and it exposes a single user-facing
command: `Jarvis`.

## What V2 provides

| Capability | V2 implementation |
| --- | --- |
| Distributed execution | Redis-backed queue with Python `spawn` Workers |
| Parameter sampling | Stable catalog: `Random`, `Grid`, `Bridson`, `AdaptiveBridson`, `CSV`, `Dynesty`, and `MultiNest` |
| Your own samplers | [Jarvis-HEP-Sampler](https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP-Sampler) plug-ins: write one Python file, `Jarvis sampler load` it, use it as `Sampling.Method` |
| HEP calculator integration | [Jarvis-HEP-Portal](https://github.com/Pengxuan-Zhu-Phys/Jarvis-Portal) for calculator file I/O |
| Python operators | [Jarvis-Operas](https://github.com/Pengxuan-Zhu-Phys/Jarvis-Operas) for registered operators and expressions |
| Reliability | Checkpoints, resume, graceful shutdown, and per-sample failure artifacts |
| Observability | Validation diagnostics, monitor snapshots, component logs, and performance summaries |
| Project workflow | Scaffold, validate, run, package, fetch, and inspect standalone projects |

`LatinHypercube` is an external sampler example in
[Jarvis-HEP-Sampler/examples](https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP-Sampler/tree/main/examples).
With the updated SDK, load `LatinHypercube.py` before using
`Sampling.Method: LatinHypercube` in a new scan. Existing scans resume from
their saved sampler source.

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
- A Redis-compatible server program (`redis-server`, `redis6-server`, or
  `valkey-server`) installed on the machine; Redis must be 6.2 or newer.
  You do not need to start a service: see below. On Windows, install and run
  both Jarvis-HEP and Redis inside WSL.

From PyPI:

```bash
python3 -m pip install Jarvis-HEP
```

From a source checkout:

```bash
python3 -m pip install -e .
```

Jarvis-HEP requires `Jarvis-HEP-Sampler>=1.1.2` (no upper bound); pip installs it with Jarvis-HEP. For development, install the sibling SDK checkout if you are changing it, then this checkout and the test tools:

```bash
python3 -m pip install -e ../Jarvis-HEP-Sampler
python3 -m pip install -e .
python3 -m pip install "pytest>=7.0" "fakeredis>=2.0" "colorlog>=6.0"
```

Install Redis, but do not start it as a service. Jarvis starts its own Redis
server for each scan and stops it when the scan ends. For example, on macOS:

```bash
brew install redis   # install only; no `brew services start redis`
```

If another Redis is already running on port 6379 (for example one started with
`brew services start redis`), Jarvis cannot use that port. For task cards that
load the project's default environment (the cards made by
`Jarvis project create` do), it picks the next free port and writes it into
`deps/environment_default.yaml`. Other task cards stop with a "port is already
used by another program" error instead.

After the first `Jarvis` command following installation, Jarvis checks whether a
Redis-compatible server executable (`redis-server`, `redis6-server`, or
`valkey-server`) is available. If it is missing, the command's normal output
is followed by a one-time, OS-specific installation hint (written to stderr so
JSON/stdout output remains usable). The hint is advisory only; Jarvis does not
install operating-system packages automatically. This is a PATH check on first
CLI use, not a hook run by `pip install`; it does not verify the server version.
The check marker is stored at `~/.jarvis/redis-install-check-v1`.

You can inspect the environment and read setup instructions at any time:

```bash
Jarvis redis                       # system, executable/version, install/start guidance
Jarvis redis check                 # repeat the PATH and bounded --version checks
Jarvis redis install               # print the installation commands for this OS
Jarvis redis start                 # explain automatic startup and optional manual debugging
Jarvis redis info --json           # machine-readable report; --json works on every topic
Jarvis redis --help                # list the info/check/install/start subcommands
Jarvis redis check --help          # help and options for this subcommand
```

Redis follows the CLI's group/subcommand help layout: the parent lists the
four commands, and each `Jarvis redis COMMAND --help` describes that command's
purpose and options. `--json` works before or after the subcommand, for example
`Jarvis redis --json check` and `Jarvis redis check --json`.

`Jarvis redis` and `Jarvis redis --json` retain the `info` default. `check` exits 1 if the
executable is absent, unusable, too old, or the runtime is native Windows.
These commands do not contact a server, test ports, start a scan, install
packages, or start services. `install` and `start` only print guidance; they do
not invoke the server even to read its version. A successful environment check
confirms the local executable, not a running service or scan readiness.

Installation guidance accounts for repository versions: Ubuntu 20.04/22.04 and
Debian 11 receive the official signed Redis APT repository steps; RHEL-family
8 systems receive the `redis:6` module instructions with the updated AppStream
requirement. Older releases without a verified recipe receive documentation
and upgrade guidance. Use `Jarvis redis install` for the full commands for
your detected system.

The default install includes every Jarvis runtime dependency: the Redis Python
client, `msgpack`, `aiofiles`, Textual Monitor, `Jarvis-HEP-Portal`,
`Jarvis-HEP-Sampler`, `Jarvis-Operas`, and [JarvisPLOT](https://github.com/Pengxuan-Zhu-Phys/JarvisPLOT).
Jarvis does not publish optional dependency groups.

For the complete installation guide, Redis options, and project-packaging
workflow, see the [installation guide](https://github.com/Pengxuan-Zhu-Phys/Jarvis-Books/blob/main/Docs/Jarvis-HEP/INSTALL.md).

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

`Jarvis -v` shows the framework banner followed by a table of locally installed
Jarvis ecosystem versions. `Jarvis --version --json` returns the same inventory
as JSON. Each `Jarvis sampler|operas|plot|portal -v` command prints the child
distribution's official name and version, matching its standalone CLI.

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
Jarvis sampler init MySampler      # write your own sampler from a template
Jarvis sampler load MySampler.py   # then use Sampling.Method: MySampler
Jarvis portal man
Jarvis operas list
Jarvis gen-plot-yaml TASK.yaml
Jarvis plot path/to/plot.yaml      # Jarvis-PLOT is installed by default
```

Use `Jarvis COMMAND -h` for command-specific options. Screen logging defaults
to `WARNING`; use `--console-level INFO`, `--debug`, or `--silence` as needed.
File logs remain available under `logs/<scan>/`.

## AI agents and YAML assistance

Codex, Claude and other agents should start with the installed CLI contract:

```bash
Jarvis agent --brief --json        # short contract retains all rules/workflows
Jarvis man agent-readme --brief --json  # same contract, through the manual index
Jarvis agent --section commands --command run --json
Jarvis agent search "checkpoint heartbeat" --json
Jarvis agent inspect bin/task.yaml --section runtime --json
Jarvis agent schema --list --json  # small schema index
Jarvis agent schema --json         # bundled schemas; all references resolve locally
Jarvis agent template --json       # list task cards and Sampling fragments
Jarvis agent template random-operas --json
Jarvis validate bin/task.yaml --json
```

The contract includes the installed version, question routes, argument arrays,
command effects, authorization boundaries and authoring workflows. Follow the
installed `man` pages for field names and registered sampler options, and the
owning Portal/Operas/PLOT CLI for runtime behavior. The schema bundle describes
the bundled structural vocabulary; semantic checks still require `validate`.
`Jarvis agent schema --id URI` prints one original schema document or local
JSONPointer fragment from that bundle without fetching the URI. Sampler manual
pages include `options_schema`, `schema_argv` and `bounds_schema_argv`. Omit
`--brief` to read the full original contract; `--section` reads one section.
Search performs keyword discovery and returns real manual/help `next_argv`;
absence of a match is not proof that a feature is unsupported.

`inspect` explains the selected card's loaded defaults, merged environment,
normalized runtime, planned paths and validation, with field sources. Use
`--section paths|runtime|sources|validation|card|policies` to limit output.
Configured worker counts and Redis ports are distinct from runtime selection;
this query does not import user operators, probe ports or create scan resources.

Named templates print their exact YAML to stdout; `--json` adds artifact kind
and prerequisites. Some examples depend on the files copied by `project create`;
Sampling fragments provide a full base-template query and instructions to
replace its `Sampling` block, then validate the composed card. A full card needs
`Scan`, `EnvReqs` and at least one `Calculators` or `Operas` section. Agents should edit the user's
file within the requested scope, use each validation issue's `manual_argv` to
look up the correction, then validate again. Preserve scientific choices and
existing files. Repeated explicit YAML keys in tasks and selected default files
are rejected with the field and both source locations; legitimate YAML merge
overrides remain allowed. Unrecognized manual subtopics fail with parent
navigation instead of returning a misleading successful parent lookup.

`validate` does not execute calculators or start Redis/Workers. A successful
report does not prove runtime or physics correctness. `check` and `run` execute
configured code and create runtime resources; bare `Jarvis TASK.yaml` also runs
the task. Existing user authorization remains valid, but YAML assistance alone
does not authorize execution, dependency installation or process termination.
CLI onboarding and dependency imports may initialize user caches even during
queries. JSON command output uses stdout; advisory diagnostics use stderr.

For a missing capability, produce an email draft with:

```bash
Jarvis agent request --title "Required feature" --description "Goal and constraints" --json
```

The recipient is `Jarvis.HEP.AI@mail.grokbot.com`, with the exact marker
`Resquest@Jarvis-HEP` (the spelling follows Jarvis-PLOT). Include installed version,
OS/Python/installation method, minimal YAML, commands and diagnostics, actual and
expected behavior, and scientific/compatibility constraints. Report directly in
a maintainer conversation; otherwise provide the draft for the user to submit.
The CLI does not send mail or write the draft to a file.
Use `--yaml`, `--commands`, `--diagnostics`, `--actual`, `--expected`,
`--constraints` and `--installation-method` to fill the matching draft sections.
These inputs are literal text; paths are not read and commands are not executed.
`missing_fields` identifies remaining blanks; `draft_ready` only means all required
text sections are filled, not that their contents were verified or mail was sent.

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
are not. See [SAMPLE storage](https://github.com/Pengxuan-Zhu-Phys/Jarvis-Books/blob/main/Docs/Jarvis-HEP/specs/sample-storage.md) for details.

## Documentation

- [Online documentation](https://pengxuan-zhu-phys.github.io/Jarvis-Docs/)
- `Jarvis man` — built-in manual for every task-card field, with examples
- [Installation guide](https://github.com/Pengxuan-Zhu-Phys/Jarvis-Books/blob/main/Docs/Jarvis-HEP/INSTALL.md) — installation, Redis setup, full command list, and
  sharing projects
- [SAMPLE storage](https://github.com/Pengxuan-Zhu-Phys/Jarvis-Books/blob/main/Docs/Jarvis-HEP/specs/sample-storage.md) — keeping or discarding per-point files
- [CHANGELOG](https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP/blob/master/CHANGELOG.md) — changes in each release
- [Developer documentation](https://github.com/Pengxuan-Zhu-Phys/Jarvis-Books/blob/main/Docs/Jarvis-HEP/README.md) — design, specifications and roadmap
- [Project template](https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP/blob/master/jarvishep2/project_template/README.md) — what
  `Jarvis project create` puts in a new project

### Download the PDF manual

```bash
Jarvis download-manual                              # latest stable release, current directory
Jarvis download-manual --release v2.1.0 -o manual.pdf # specific release
Jarvis download-manual -o ~/Downloads --force         # replace an existing copy
Jarvis download-manual --json                        # machine-readable result
```

The command downloads the `Jarvis-HEP-Manual.pdf` attachment from a public
[Jarvis-HEP GitHub Release](https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP/releases).
The default `latest` follows GitHub's latest stable release; use `--release TAG`
to select a specific release independently of your installed package version.
A plain version such as `2.1.0` is accepted as `v2.1.0`.

`--output` (`-o`) accepts a file path or an existing directory. Existing files
require `--force`; a failed or interrupted download leaves the original file
intact. Downloads have a 256 MiB size limit and a 30-second network operation
timeout, configurable with `--timeout`. The downloaded file receives basic PDF
header, trailer and transfer-length checks before it is saved. If the release
or attachment has not been published, the command reports this clearly.
Downloading the manual does not require Redis or a GitHub login.

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

Developer documentation, including the change rules, architecture and
roadmap, is in [Jarvis-Books](https://github.com/Pengxuan-Zhu-Phys/Jarvis-Books/blob/main/Docs/Jarvis-HEP/README.md).

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

Jarvis-HEP is released under the [MIT License](https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP/blob/master/LICENSE). Bundled third-party
code keeps its own license (for example
[Dynesty](https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP/blob/master/jarvishep2/sampling/Source/Dynesty/LICENSE)).

## Releasing

`jhrel` is a maintainer-local script, excluded from Git and distribution archives.

Record changes in `CHANGELOG.md`, then run from this checkout (replace the
example version with the intended release):

```bash
./jhrel 2.1.1 --dry
./jhrel 2.1.1
```

`--dry` performs read-only source and release-note checks. The normal command
requires a clean Git checkout, builds and checks the wheel and sdist, uploads
to PyPI, refreshes the local version with `python3 -m pip install -e .`, and
verifies the installed version. Only then does it commit the version, push the
`vX.Y.Z` tag, and create a GitHub release with changelog notes. Configure Twine
and `gh` authentication before publishing.
