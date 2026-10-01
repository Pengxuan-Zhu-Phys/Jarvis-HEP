# Jarvis-HEP V2 — Installation

How to install and use Jarvis-HEP V2. It is installed from PyPI as
`Jarvis-HEP`, imported in Python as `jarvishep2`, and run with the `Jarvis`
command. For the fields you can use in a task YAML, run `Jarvis man` or see the
[online documentation](https://pengxuan-zhu-phys.github.io/Jarvis-Docs/).

## Requirements

- **Python ≥ 3.10** (Workers use the `spawn` multiprocessing context; macOS and Linux supported)
- **Redis server** for distributed runs. V2 connects internally to the local
  `127.0.0.1:6379` service; it is not configured in task YAML. Tests do *not* need one — the
  suite runs on `fakeredis`.

## Install

From PyPI:

```bash
python3 -m pip install Jarvis-HEP
```

From a source checkout:

```bash
git clone https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP.git
cd Jarvis-HEP

# runtime (recommended)
python3 -m pip install -e .

# development tests
python3 -m pip install "pytest>=7.0" "fakeredis>=2.0" "colorlog>=6.0"
```

Core runtime depends on:

- **`Jarvis-HEP-Portal`** — calculator I/O formats (JSON, CSV, TSV, DAT, Wolfram, …)
- **`Jarvis-Operas`** — operator registry and qualified expression functions (e.g. `helper.eggbox2d`)
- **`JarvisPLOT`** (product name: **Jarvis-PLOT**) — YAML-driven plotting and flowchart rendering

These three packages are core dependencies. Their minimum versions track the
current PyPI releases so a fresh `pip install Jarvis-HEP` resolves the latest
compatible release available from the index.

If you are developing those packages too, clone them next to `Jarvis-HEP` and
install all of them from source:

```bash
python3 -m pip install -e ../Jarvis-Portal
python3 -m pip install -e ../Jarvis-Operas
python3 -m pip install -e .
python3 -m pip install "pytest>=7.0" "fakeredis>=2.0" "colorlog>=6.0"
```

All Jarvis runtime dependencies are installed by default; Jarvis does not
publish optional dependency groups. Test-only tools remain separate from the
runtime package.

## Task-card rules

A task card is the YAML file that describes one scan. All `Jarvis` commands
read it the same way:

- Required: `Scan.name`, `Sampling.Method`, and an `EnvReqs` section.
- At least one of `Calculators` (external programs) or `Operas` (Python
  functions) must be present.
- Optional: `LibDeps` (external libraries your calculators need).
- Settings for the chosen sampling method go under `Sampling.Bounds`, written
  in lowercase with underscores (e.g. `point_number`).
- Runtime settings (number of workers, checkpoint interval, …) go under
  `EnvReqs.V2`.
- The other blocks allowed under `EnvReqs` are checks run before the scan
  starts: `Python`, `CERN_ROOT`, `Check_default_dependencies`, and `OS`.
  `OS` lists the operating systems the scan may run on; each item has a
  `name` (as reported by Python's `platform.system()`, e.g. `Linux` or
  `Darwin`) and a minimum `version` such as `>=5.4`. An empty list allows any
  system. Any other key under `EnvReqs` is an error.
- `Jarvis check TASK.yaml` runs your calculators on a few points as a quick
  test, one at a time. To test specific points, list them in a CSV file and
  set its path in `EnvReqs.V2.check_modules.data`; otherwise a few points are
  drawn from your sampler.

Run `Jarvis man` before writing or editing a card, and `Jarvis validate` before
running it.

**Commands.** Jarvis-HEP V2 installs a single command, `Jarvis`. The V1-era
`Jarvis2` command no longer exists.

## Redis

```bash
# macOS
brew install redis && brew services start redis
# or containerized
docker run -d --name jarvis-redis -p 6379:6379 redis:7
# sanity check
redis-cli ping        # → PONG
```

### First-command Redis check

After the first `Jarvis` command, V2 checks for a Redis-compatible server
executable (`redis-server`, `redis6-server`, or `valkey-server`). If none is
installed, Jarvis prints an OS-specific command after the normal command
output. For example:

```bash
# macOS (Homebrew's `redis` formula provides `redis-server`)
brew install redis

# Ubuntu / Debian
sudo apt install -y redis-server
```

Linux hints also cover the common package families below; unknown Linux
distributions fall back to whichever package manager is detected on `PATH`:

| Distribution family | Suggested command |
| --- | --- |
| Fedora | `sudo dnf install -y valkey` |
| RHEL / CentOS / Rocky / Alma | `sudo dnf install -y redis` |
| Amazon Linux 2023 | `sudo dnf install -y valkey` |
| Amazon Linux 2 | `sudo amazon-linux-extras install redis6 -y` |
| Arch / Manjaro / Artix | `sudo pacman -S valkey` |
| openSUSE / SLES | `sudo zypper install -y redis` |
| Alpine | `sudo apk add redis` |
| Gentoo | `sudo emerge --ask dev-db/redis` |
| Void | `sudo xbps-install -S redis` |
| NixOS | `nix profile install nixpkgs#redis` |
| Solus / Mageia | `sudo eopkg install redis` / `sudo urpmi redis` |

Package names and service commands can vary by release. Fedora and Arch now
ship Valkey in their maintained repositories; Valkey is Redis-protocol
compatible and Jarvis recognizes its `valkey-server` binary. Jarvis only
prints the suggestion and never executes it.

This is a one-time advisory and never runs the package manager automatically.
The completion marker is `~/.jarvis/redis-install-check-v1`; remove that marker
if you intentionally want Jarvis to show the onboarding check again.

> Redis connection details are intentionally internal to V2. Ensure the local service above is
> running before launching a scan; V2 fails early with a focused connection error if it is not.

## Verify

```bash
python3 -m pytest -q          # full suite; fakeredis opens a local test socket
python3 -m pytest -q tests/test_d0_integration.py tests/test_worker_mvp.py   # quick subset
```

`tests/test_adaptive_bridson.py` is excluded from the default suite because it
contains long-running process/resume integration tests. Run it explicitly when
working on AdaptiveBridson:

```bash
python3 -m pytest -q tests/test_adaptive_bridson.py
```

`tests/test_ensemble_samplers.py` is also excluded from the default suite. Its
feedback-loop tests wait for durable Archiver acknowledgements; the mock-only
fixture has no Archiver and therefore pays a five-second barrier timeout per
generation. Run it explicitly when working on EnsembleMCMC/DEMCMC/PT:

```bash
python3 -m pytest -q tests/test_ensemble_samplers.py
```

The slow distributed acceptance gates in `tests/test_distributed_acceptance.py`
are also skipped by default. They launch multiple Workers/Archiver processes
and include machine-relative throughput thresholds, so run them explicitly
when validating distributed performance:

```bash
python3 -m pytest -q tests/test_distributed_acceptance.py
```

`tests/test_distributed_resume.py`, `tests/test_mcmc_sampler.py`, and
`tests/test_worker_pool.py` are also excluded from the default suite. They run
interruption/checkpoint-resume, multi-process sampler, and Worker calculator-
pool scenarios, so execute the relevant file manually when changing those
subsystems:

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

## Quickstart

Save as `quickstart.yaml` (any directory; outputs land in `<project-root>/outputs/<scan-name>/`):

```yaml
Scan:
  name: quickstart
EnvReqs:
  V2:
    workers: 2
    batch_size: 256
Sampling:
  Method: Random
  Bounds:
    point_number: 20
    seed: 7
  Variables:
    - name: x
      distribution: {type: Flat, parameters: {min: 0.0, max: 1.0}}
    - name: y
      distribution: {type: Flat, parameters: {min: 0.0, max: 1.0}}
  LogLikelihood:
    - name: LogL_Z
      expression: z
Operas:
  Modules:
    - name: TrivialEggbox
      operator: jarvishep2.testing.eggbox.eggbox2d_numpy
      call_mode: call
      input:
        - {name: x, expression: x}
        - {name: y, expression: y}
      output:
        - {name: z, entry: z}
```

Run and inspect:

```bash
Jarvis -v                            # logo + authors + package version
Jarvis --version
Jarvis --refs                        # framework and bundled-sampler citations

Jarvis run quickstart.yaml           # preferred
Jarvis run quickstart.yaml --console-level INFO   # show progress messages (default: WARNING)
Jarvis run quickstart.yaml --silence # no console log; files under logs/<scan>/ still written
Jarvis quickstart.yaml               # legacy alias → run
Jarvis check path/to/check.yaml      # quick test of calculators on a few points
Jarvis validate path/to/task.yaml    # check the card for mistakes without running
Jarvis man                           # YAML writing manuals (keys, paths, examples)
Jarvis man calculator.execution.output --type JSON
Jarvis man operas                    # Operas.Modules list YAML shape
Jarvis man yaml.Calculators.Modules.execution
Jarvis man yaml.EnvReqs.V2
# List-valued fields use their field name alone; man reports type `list`.
Jarvis monitor                       # one read-only status snapshot
Jarvis plot path/to/scene.yaml       # render Jarvis-PLOT scene (core dependency)
Jarvis portal man                    # Portal format *runtime* manuals (same as jportal man)
Jarvis portal man slha               # format runtime manual
Jarvis portal path/to/io.yaml        # run Portal IO YAML (same as jportal file)
Jarvis operas list                   # Operas operator catalog
Jarvis operas info helper.eggbox2d   # operator signature / return shape

# while a scan is running (another terminal), or after a hard kill:
Jarvis ps                            # running Jarvis* / Jarvis-Redis* processes
Jarvis kill                          # list + confirm [y/N], then terminate
Jarvis kill --yes                    # non-interactive kill

# legacy run/plot aliases (not YAML validation interfaces):
Jarvis quickstart.yaml --resume
Jarvis path/to/plot.yaml --plot      # deprecated warning → prefer `Jarvis plot`
```

## Project tools (scaffold, catalog, public / restricted packs)

Packing, encrypting and decrypting projects are all done with `Jarvis project`
commands; you never need to run `openssl` yourself. The list of official
example projects is kept in the
[Jarvis-Examples](https://github.com/Pengxuan-Zhu-Phys/Jarvis-Examples)
repository.

### Local project

```bash
Jarvis project create MyScan
cd MyScan
Jarvis run bin/quickstart_bridson_operas.yaml

Jarvis project pack . --share          # plain tarball
Jarvis project pack . --repro
Jarvis project pack . --full
Jarvis project pack . --man            # write pack manifest only
```

### Official library (GitHub JSON catalog — no PyPI package)

Default index:

```text
https://raw.githubusercontent.com/Pengxuan-Zhu-Phys/Jarvis-Examples/main/catalog/official_project_library.json
```

```bash
# Browse projects; columns include Access (public|restricted) and Key (no|required)
Jarvis project browse

Jarvis project info Eggbox
Jarvis project fetch Eggbox            # public — no key
```

### Restricted (encrypted) projects — fetch

```bash
# See Key: required in the browse output
Jarvis project browse

# Type the key without it being shown or saved in your shell history
read -rs JARVIS_PROJECT_FETCH_KEY && export JARVIS_PROJECT_FETCH_KEY

# Decrypt + unpack
Jarvis project fetch SecretName
```

You can also pass `--key 'YOUR_KEY'`, but then the key is part of the command
line: other users on the same machine can see it with `ps` while the command
runs, and it is saved in your shell history. On shared machines, use the
environment variable above.

Backend: OpenSSL-compatible AES-256-CBC (PBKDF2). Jarvis uses system `openssl` if
available, otherwise optional `pip install cryptography`. **You still only call
`Jarvis project fetch`.**

### Restricted projects — maintainers (encrypt)

```bash
# Same key variable as for fetch
read -rs JARVIS_PROJECT_FETCH_KEY && export JARVIS_PROJECT_FETCH_KEY

# Pack then encrypt → *.tar.gz.jenc
Jarvis project pack MyPrivate --repro --encrypt

# Or encrypt an existing archive
Jarvis project encrypt MyPrivate_repro_….tar.gz
```

Upload the `.jenc`, register it in `Jarvis-Examples/catalog/official_project_library.json`
with `access: restricted` and `requires_key: true`. Do **not** publish the plaintext tree
to the public Examples repo. Share the key out-of-band.

Optional overrides:

```bash
export JARVIS_OFFICIAL_LIBRARY_INDEX_URL=file:///path/to/catalog.json   # local test
export JARVIS_OFFICIAL_LIBRARY_TIMEOUT_SEC=30
```

Stop a running scan with **Ctrl+C** (SIGINT). Jarvis will shut down Workers, Archiver,
and any managed `Jarvis-Redis:<scan>` it started. Prefer not to use **Ctrl+Z** — suspend
leaves the job half-alive; Ctrl+Z is ignored during a scan for that reason.

Exit codes: **0** success (also version / empty `ps`); **1** any sample failed (including partial),
kill incomplete, or runtime error; **2** usage/config; **130** interrupted.

`--plot` currently expects a **JarvisPLOT YAML**, not the scan task YAML (scan-driven plot
scenes are partial). Portal / Operas / project tools are native `Jarvis` subcommands.

Outputs under `outputs/<scan>/` (example scan name depends on the YAML):

- `DATABASE/samples.hdf5` — one row per evaluated point: input parameters,
  calculator/operator outputs, `LogL`, and file paths if any
- `DATABASE/samples.csv` — the same table as CSV, written automatically when the
  scan finishes; run `Jarvis convert TASK.yaml` to regenerate it (for example
  after an interrupted scan)
- `SAMPLE/<bucket>/<uuid>/` or `SAMPLE/<bucket>.tar.gz` — files kept for each
  point (outputs marked `save: true`, per-point logs)
- `run_summary.json` / `.csv` / `.txt` — how many points succeeded or failed,
  and how fast the scan ran
- The screen output and main log end with a **`[Scan Performance]`** block
  (`samples / sec`, `avg sample (sec)`, …)
- Logs are under `logs/<scan>/`, one file per process: `core.log`,
  `worker-00.log`…, `archiver.log`. The log for a single point is in
  `SAMPLE/.../Sample_running.log`.
- Screen messages: default level **WARNING**; use `--console-level INFO` for more
  detail, or `--silence` / `-s` to turn screen output off. Log files always keep
  full detail.
- Resume data: `<project-root>/checkpoints/<scan>/<sampler>/state.pkl`

## Uninstall

```bash
python3 -m pip uninstall Jarvis-HEP
```
