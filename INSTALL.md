# Jarvis-HEP V2 — Installation

How to install and use Jarvis-HEP V2. It is installed from PyPI as
`Jarvis-HEP`, imported in Python as `jarvishep2`, and run with the `Jarvis`
command. For the fields you can use in a task YAML, run `Jarvis man` or see the
[online documentation](https://pengxuan-zhu-phys.github.io/Jarvis-Docs/).

## Requirements

- **Python ≥ 3.10** (Workers use the `spawn` multiprocessing context; macOS and Linux supported)
- **Redis server program** (`redis-server`, `redis6-server`, or the compatible
  `valkey-server`) installed on the machine; Redis must be 6.2 or newer. Jarvis
  starts and stops its own Redis for each scan, so you do not run a Redis service
  yourself. Many unit tests use `fakeredis`; real runtime/integration tests need
  a Redis-compatible server program.
- **Windows:** use WSL, and install and run both Jarvis-HEP and Redis inside the
  Linux environment. The scan runtime is not supported on native Windows.

## Install

HEP 2.0.17's wheel is public and has passed a fresh complete installation from
PyPI, `pip check`, SDK absence and installed agent/Redis checks. The source
archive upload is still being verified; see the [release record](docs/release-2.0.17.md).

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
- **`Jarvis-HEP-Sampler`** — plug-in samplers for `Jarvis sampler` and `Sampling.Method` (any release newer than 1)
- **`JarvisPLOT`** (product name: **Jarvis-PLOT**) — YAML-driven plotting and flowchart rendering

These packages are core dependencies. Their minimum versions track the
current PyPI releases so a fresh `pip install Jarvis-HEP` resolves the latest
compatible release available from the index. Jarvis-HEP-Sampler must be
newer than 1. `Jarvis sampler` is the same command as `jsampler`.

If you are developing those packages too, clone them next to `Jarvis-HEP` and
install all of them from source:

```bash
python3 -m pip install -e ../Jarvis-Portal
python3 -m pip install -e ../Jarvis-Operas
python3 -m pip install -e ../Jarvis-HEP-Sampler
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

Jarvis uses Redis to pass work between its processes. You only need to
**install** the Redis server program; a global Redis service is not required:

```bash
# macOS
brew install redis            # no `brew services start redis`
# Ubuntu 24.04 / Debian 12 (older releases: use `Jarvis redis install`)
sudo apt update
sudo apt install -y redis-server
```

For every scan, Jarvis starts its own Redis server on the local machine (it
shows up as `Jarvis-Redis:<scan>` in `Jarvis ps`) and stops it when the scan
ends or you press Ctrl+C.

Redis settings are configured in YAML under `EnvReqs.V2.redis` (`host`, `port`,
and `db`; defaults `127.0.0.1`, `6379`, and `0`). For project cards, the port
usually comes from the project's `deps/environment_default.yaml`. This corrects
the older installation guide's statement that Redis was not configured in task
YAML. If something else already uses that port, for example a Redis service you
started yourself:

- for task cards that load the project's default environment
  (`EnvReqs.Check_default_dependencies.default_yaml_path`, as in the cards made
  by `Jarvis project create`), Jarvis picks the next free port, writes it into
  `deps/environment_default.yaml`, and prints a warning;
- for other task cards, the run stops with
  `Redis port 127.0.0.1:6379 is already used by another program`. Stop the other Redis
  (`brew services stop redis`, or `sudo systemctl stop redis-server` on Linux)
  or make the card load the project defaults.

Some Linux packages start a Redis service automatically when installed; stop and
disable it the same way if you see that error.

### Redis environment commands

```bash
Jarvis redis                    # same as `Jarvis redis info`
Jarvis redis info               # OS/distro, executable, version and setup instructions
Jarvis redis check              # fresh PATH lookup and bounded --version check
Jarvis redis install            # installation instructions only
Jarvis redis start              # startup instructions only
Jarvis redis check --json        # --json is available for every topic
```

Use the same command-group help layout as the rest of the CLI:

```bash
Jarvis redis --help              # list the four commands and shared options
Jarvis redis info --help         # overview help
Jarvis redis check --help        # environment-check help and exit behavior
Jarvis redis install --help      # installation-guidance help
Jarvis redis start --help        # startup-guidance help
```

Each command has its own description and options, including `--json`. The JSON
flag works before or after the command: `Jarvis redis --json check` and
`Jarvis redis check --json` are equivalent. Bare `Jarvis redis` and
`Jarvis redis --json` still select `info`. Agents can discover the same exact
child commands through `Jarvis agent --command "redis check" --json` or
`Jarvis agent search "redis check" --json`.

`info` and `check` inspect the executable with a bounded `--version` subprocess.
`check` returns 0 for a usable compatible executable and 1 when it is absent,
cannot execute, reports an unknown/old version (Redis < 6.2), or runs on native
Windows. An environment check does not contact Redis, test a port, verify a
running service, or start a scan. `install` and `start` only print instructions,
and do not execute the server or package manager. All topics support one JSON
report on stdout with `--json`.

For normal scans, use `Jarvis run TASK.yaml`: Jarvis handles startup and cleanup.
`Jarvis redis start` explains this and prints an optional foreground debugging
command, using the detected server binary, equivalent to:

```bash
redis-server --bind 127.0.0.1 --port 6380 --save "" --appendonly no
```

This manual instance binds only to loopback, uses a separate example port, and
disables persistence. It is optional; stop it with Ctrl+C. The guide does not
check whether port 6380 is free. The printed command uses `valkey-server` or
`redis6-server` where appropriate.

### First-command Redis check

After the first `Jarvis` command following installation, V2 checks for a Redis-compatible server
executable (`redis-server`, `redis6-server`, or `valkey-server`). If none is
installed, Jarvis prints OS-specific installation and startup guidance to
stderr after the normal command output, keeping JSON stdout parseable. This
PATH-only check does not run a version subprocess. It runs on first CLI use,
not during `pip install`. Use `Jarvis redis check` to recheck at any time.
For example:

```bash
# macOS (Homebrew's `redis` formula provides `redis-server`)
brew install redis

# Ubuntu 24.04 / Debian 12 (older releases: use `Jarvis redis install`)
sudo apt update
sudo apt install -y redis-server
```

Linux hints use the concrete distribution ID first, then its ID_LIKE family;
unknown distributions fall back to the package manager detected on `PATH`.
The installation guidance includes `apt update`/`apt-get update` or `apk update`
where appropriate, followed by the package command below:

| Distribution family | Installation guidance |
| --- | --- |
| Ubuntu 20.04 / 22.04, Debian 11 | `Jarvis redis install` prints the signed Redis APT repository steps, then `sudo apt-get install -y redis` |
| Older Ubuntu / Debian releases without a verified recipe | Documentation and upgrade guidance; no unverified package command |
| Modern Debian / Ubuntu, e.g. Debian 12 / Ubuntu 24.04, and derivatives | `sudo apt install -y redis-server`; verify the installed version |
| Fedora | `sudo dnf install -y valkey` |
| RHEL / CentOS / Rocky / Alma family, version 10+ | `sudo dnf install -y valkey` |
| RHEL family, version 8 | `sudo dnf module install -y redis:6` with updated AppStream repositories |
| RHEL family, earlier or unknown version | `sudo dnf install -y redis` (or `yum` if `dnf` is absent) |
| Amazon Linux 2023 | `sudo dnf install -y valkey` |
| Amazon Linux 2 | `sudo amazon-linux-extras install redis6 -y` |
| Arch / Manjaro / Artix | `sudo pacman -S valkey` |
| openSUSE / SLES | `sudo zypper install -y redis` |
| Alpine 3.20+ or unknown version | `sudo apk add valkey` |
| Alpine before 3.20 or edge | `sudo apk add redis` |
| Gentoo | `sudo emerge --ask dev-db/redis` |
| Void | `sudo xbps-install -S redis` |
| NixOS | `nix profile install nixpkgs#redis` |
| Solus / Mageia | `sudo eopkg install redis` / `sudo urpmi redis` |

Ubuntu 20.04/22.04 and Debian 11 have default Redis packages below Jarvis's
minimum. Their printed recipe installs prerequisites, writes the official Redis
signing key and signed APT source for the known distribution codename, refreshes
APT, and installs `redis`. It changes system key/source files only when you run
the commands yourself. See the [official Redis APT installation guide](https://redis.io/docs/latest/operate/oss_and_stack/install/install-stack/apt/).
Other older Ubuntu/Debian releases without a verified repository recipe receive
documentation and upgrade guidance instead of an ordinary package command that
would install an incompatible version. A derivative's VERSION_ID is not used
as the parent distribution's repository release.

RHEL 8 requires the `redis:6` stream and current AppStream repositories:
[RHEL 8.7 rebased that stream to Redis 6.2.7](https://docs.redhat.com/en/documentation/red_hat_enterprise_linux/8/html/8.7_release_notes/new-features).
Repositories pinned to earlier releases can still provide Redis 6.0. Update the
enabled repositories and rerun `Jarvis redis check`; if `redis:5` is already
enabled or installed, follow the distribution's stream-switching instructions
first. Older or unknown RHEL-family versions retain a package-manager hint,
with no promise that the available package meets the minimum.

Package availability and versions vary by OS release and enabled repositories;
run `Jarvis redis check` after installation to verify the local executable.
Valkey is Redis-protocol compatible and Jarvis recognizes `valkey-server`.
The version-specific changes are documented by [RHEL 10](https://docs.redhat.com/en/documentation/red_hat_enterprise_linux/10/html/considerations_in_adopting_rhel_10/infrastructure-services),
[Alpine 3.20](https://pkgs.alpinelinux.org/package/v3.20/main/x86_64/valkey),
[Amazon Linux 2023](https://docs.aws.amazon.com/linux/al2023/ug/redis6-to-valkey-al2023.html),
and [Arch Linux](https://archlinux.org/news/valkey-to-replace-redis-in-the-extra-repository/).
Jarvis only prints these suggestions and never executes them.

This is a one-time advisory and never runs the package manager automatically.
The completion marker is `~/.jarvis/redis-install-check-v1`; remove that marker
if you intentionally want Jarvis to show the onboarding check again.

If no Redis-compatible server program is installed, `Jarvis run` fails early
with an installation hint. These instructions only print commands; users run
them themselves. On native Windows, the guide recommends
`wsl --install -d Ubuntu` in Administrator PowerShell, restarting when requested,
then installing and running Jarvis-HEP and Redis inside Ubuntu. See the
[Microsoft WSL installation guide](https://learn.microsoft.com/en-us/windows/wsl/install).

## Verify

```bash
python3 -m pytest -q          # many tests use fakeredis; runtime/integration tests need Redis
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
