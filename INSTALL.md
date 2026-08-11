# Jarvis-HEP V2 — Installation

Quick install guide for the V2 distributed runtime (PyPI distribution
`Jarvis-HEP`, Python package `jarvishep2`, CLI `Jarvis`).
For the architecture see
`~/Jarvis-Workshop/Jarvis-Books/Jarvis-HEP V2/DESIGN_2.0_DISTRIBUTED.md`;
for the task-YAML schema see
`~/Jarvis-Workshop/Jarvis-Books/Jarvis-HEP V2/YAML_REFERENCE_2.0.md`.

## Requirements

- **Python ≥ 3.10** (Workers use the `spawn` multiprocessing context; macOS and Linux supported)
- **Redis server** for distributed runs. V2 connects internally to the local
  `127.0.0.1:6379` service; it is not configured in task YAML. Tests do *not* need one — the
  suite runs on `fakeredis`.

## Install

From PyPI:

```bash
python3 -m pip install 'Jarvis-HEP[distributed]'
```

If you previously installed the short-lived `jarvishep2` distribution, remove
it first. It contains the same Python import package as `Jarvis-HEP` V2 and the
two distributions should not be installed together:

```bash
python3 -m pip uninstall jarvishep2
python3 -m pip install --upgrade 'Jarvis-HEP[distributed]'
```

From a source checkout:

```bash
cd ~/Jarvis-Workshop/Jarvis-HEP

# runtime + real-Redis extras (recommended)
python3 -m pip install -e '.[distributed]'

# development (tests, fakeredis, colorlog)
python3 -m pip install -e '.[distributed,dev]'
```

Core runtime depends on:

- **`Jarvis-HEP-Portal`** — calculator I/O formats (JSON, CSV, TSV, DAT, Wolfram, …)
- **`Jarvis-Operas`** — operator registry and qualified expression functions (e.g. `helper.eggbox2d`)

Local editable installs of those packages are fine during development:

```bash
python3 -m pip install -e ../Jarvis-Portal
python3 -m pip install -e ../Jarvis-Operas
python3 -m pip install -e '.[distributed,dev]'
```

Extras:

- `distributed` = `redis`, `msgpack`, `aiofiles`
- `plot` = `JarvisPLOT` (YAML-driven plotting; CLI bridge)
- `dev` = `pytest`, `fakeredis`, `colorlog`
- `operas` = **deprecated no-op alias** (Operas is core since D12.0; kept for old install scripts)

```bash
python3 -m pip install -e '.[distributed,plot,dev]'
```

## Task-card contract

Use the same contract for every task card and every CLI path:

- Required: `Scan.name`, `Sampling.Method`, and `EnvReqs`.
- Conditional: at least one of `Calculators` or `Operas`.
- Optional: `LibDeps`.
- Method-specific sampler settings are under `Sampling.Bounds` and use lower
  `snake_case` keys.
- Runtime settings are under `EnvReqs.V2`. `EnvReqs.OS` is the supported
  V1-compatible host preflight list; each item uses `name` and `version`
  (`>=X.Y`), while an empty list imposes no OS restriction. Other undeclared
  siblings are rejected.
- `Jarvis check TASK.yaml` is the only check entry point. The only fixed-point
  CSV setting is `EnvReqs.V2.check_modules.data`.
- Outputs go directly to `SAMPLE/`; task YAML has no `cleanup.strategy` or
  `archiver.handoff` interface.

`Jarvis man`, `Jarvis validate`, `Jarvis run`, and `Jarvis check` use this same
vocabulary. Use `Jarvis man` before writing or editing a card.

**CLI ownership.** The product command is **`Jarvis`** (this package). V1
(`Jarvis-HEP` / `jarvishep`) is retired from the CLI: it no longer installs a
console script. The old **`Jarvis2`** command is removed. Do not import
`jarvishep` from V2 code.

## Redis

```bash
# macOS
brew install redis && brew services start redis
# or containerized
docker run -d --name jarvis-redis -p 6379:6379 redis:7
# sanity check
redis-cli ping        # → PONG
```

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
Jarvis run quickstart.yaml --console-level WARNING
Jarvis run quickstart.yaml --silence # no console log; files under logs/<scan>/ still written
Jarvis quickstart.yaml               # legacy alias → run
Jarvis check path/to/check.yaml      # fixed-point calculator smoke (only check entry point)
Jarvis validate path/to/task.yaml    # schema + contract gate
Jarvis man                           # YAML writing manuals (keys, paths, examples)
Jarvis man calculator.execution.output --type JSON
Jarvis man operas                    # Operas.Modules list YAML shape
Jarvis man yaml.Calculators.Modules.execution
Jarvis man yaml.EnvReqs.V2
# List-valued fields use their field name alone; man reports type `list`.
Jarvis monitor                       # one read-only status snapshot
Jarvis plot path/to/scene.yaml       # render JarvisPLOT scene (plot extra)
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

All packing and **encrypt/decrypt** go through the CLI. End users do **not** run `openssl`
by hand. Design notes:
`Jarvis-Books/Jarvis-HEP V2/components/project_tools.md` and
`Jarvis-Examples/catalog/README.md`.

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
# List projects; columns include Access (public|restricted) and Key (no|required)
Jarvis project list

Jarvis project info Eggbox
Jarvis project fetch Eggbox            # public — no key
```

### Restricted (encrypted) projects — fetch

```bash
# See Key: required in the list
Jarvis project list

# Decrypt + unpack (preferred)
Jarvis project fetch SecretName --key 'YOUR_KEY'

# Or set the key once
export JARVIS_PROJECT_FETCH_KEY='YOUR_KEY'
Jarvis project fetch SecretName
```

Backend: OpenSSL-compatible AES-256-CBC (PBKDF2). Jarvis uses system `openssl` if
available, otherwise optional `pip install cryptography`. **You still only call
`Jarvis project fetch`.**

### Restricted projects — maintainers (encrypt)

```bash
# Pack then encrypt → *.tar.gz.jenc
Jarvis project pack MyPrivate --repro --encrypt --key 'YOUR_KEY'

# Or encrypt an existing archive
Jarvis project encrypt MyPrivate_repro_….tar.gz --key 'YOUR_KEY'
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

- `DATABASE/samples.hdf5` — full **observables** JSON rows (params + module outputs + `LogL` + optional paths)
- `DATABASE/samples.csv` — full-column CSV export of the same records (post-run plot emit)
- `SAMPLE/<bucket>/<uuid>/` or `SAMPLE/<bucket>.tar.gz` — per-sample artifacts (`save: true` files, logs)
- `run_summary.json` / `.csv` / `.txt` — end-of-run counters + throughput
- Console / main log ends with a **`[Scan Performance]`** block (`samples / sec`, `avg sample (sec)`, …)
- Logs under `logs/<scan>/` by component: `core.log`, `worker-00.log`…, `archiver.log`
  (style from packaged `jarvishep2/card/logging.yaml`; per-sample detail in `SAMPLE/.../Sample_running.log`)
- Console: default **INFO** on screen; `Jarvis run … --console-level WARNING`; `--silence` / `-s` turns screen off
- Checkpoints: `<project-root>/checkpoints/<scan>/<sampler>/state.pkl`

## Uninstall

```bash
python3 -m pip uninstall Jarvis-HEP
```
