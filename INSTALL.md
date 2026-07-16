# Jarvis-HEP V2 — Installation

Quick install guide for the V2 distributed runtime (`jarvishep2`, CLI `Jarvis2`).
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

```bash
cd ~/Jarvis-Workshop/Jarvis-HEP-v2

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

**Side-by-side with V1.** V2 uses a distinct distribution (`jarvishep2`) and CLI entry
(`Jarvis2`), so it installs alongside the frozen V1 line (`Jarvis-HEP` / `Jarvis` 1.7.4)
in the same environment without clobbering it. Do not import `jarvishep` from V2 code.

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

## Quickstart

Save as `quickstart.yaml` (any directory; outputs land in `<project-root>/outputs/<scan-name>/`):

```yaml
project_name: quickstart
Scan:
  name: quickstart
EnvReqs:
  V2:
    workers: 2
    batch_size: 256
Sampling:
  Method: Random
  "Point number": 20
  Seed: 7
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
Jarvis2 -v                            # logo + authors + package version
Jarvis2 --version

Jarvis2 run quickstart.yaml           # preferred
Jarvis2 run quickstart.yaml --console-level WARNING
Jarvis2 run quickstart.yaml --silence # no console log; files under logs/<scan>/ still written
Jarvis2 quickstart.yaml               # legacy alias → run
Jarvis2 check path/to/check.yaml      # fixed-point calculator smoke
Jarvis2 monitor                       # one read-only status snapshot
Jarvis2 plot path/to/scene.yaml       # render JarvisPLOT scene (plot extra)
Jarvis2 portal man                    # list Portal formats (same as jportal man)
Jarvis2 portal man slha               # format manual
Jarvis2 portal path/to/io.yaml        # run Portal IO YAML (same as jportal file)
Jarvis2 operas list                   # list Operas operators

# while a scan is running (another terminal), or after a hard kill:
Jarvis2 ps                            # running Jarvis2* / Jarvis-Redis* processes
Jarvis2 kill                          # list + confirm [y/N], then terminate
Jarvis2 kill --yes                    # non-interactive kill

# legacy (still work):
Jarvis2 quickstart.yaml --resume
Jarvis2 path/to/plot.yaml --plot      # deprecated warning → prefer `Jarvis2 plot`
```

## Project tools (scaffold, catalog, public / restricted packs)

All packing and **encrypt/decrypt** go through the CLI. End users do **not** run `openssl`
by hand. Design notes:
`Jarvis-Books/Jarvis-HEP V2/components/project_tools.md` and
`Jarvis-Examples/catalog/README.md`.

### Local project

```bash
Jarvis2 project create MyScan
cd MyScan
Jarvis2 run bin/quickstart_bridson_operas.yaml

Jarvis2 project pack . --share          # plain tarball
Jarvis2 project pack . --repro
Jarvis2 project pack . --full
Jarvis2 project pack . --man            # write pack manifest only
```

### Official library (GitHub JSON catalog — no PyPI package)

Default index:

```text
https://raw.githubusercontent.com/Pengxuan-Zhu-Phys/Jarvis-Examples/main/catalog/official_project_library.json
```

```bash
# List projects; columns include Access (public|restricted) and Key (no|required)
Jarvis2 project list

Jarvis2 project info Eggbox
Jarvis2 project fetch Eggbox            # public — no key
```

### Restricted (encrypted) projects — fetch

```bash
# See Key: required in the list
Jarvis2 project list

# Decrypt + unpack (preferred)
Jarvis2 project fetch SecretName --key 'YOUR_KEY'

# Or set the key once
export JARVIS_PROJECT_FETCH_KEY='YOUR_KEY'
Jarvis2 project fetch SecretName
```

Backend: OpenSSL-compatible AES-256-CBC (PBKDF2). Jarvis2 uses system `openssl` if
available, otherwise optional `pip install cryptography`. **You still only call
`Jarvis2 project fetch`.**

### Restricted projects — maintainers (encrypt)

```bash
# Pack then encrypt → *.tar.gz.jenc
Jarvis2 project pack MyPrivate --repro --encrypt --key 'YOUR_KEY'

# Or encrypt an existing archive
Jarvis2 project encrypt MyPrivate_repro_….tar.gz --key 'YOUR_KEY'
```

Upload the `.jenc`, register it in `Jarvis-Examples/catalog/official_project_library.json`
with `access: restricted` and `requires_key: true`. Do **not** publish the plaintext tree
to the public Examples repo. Share the key out-of-band.

Optional overrides:

```bash
export JARVIS_OFFICIAL_LIBRARY_INDEX_URL=file:///path/to/catalog.json   # local test
export JARVIS_OFFICIAL_LIBRARY_TIMEOUT_SEC=30
```

Stop a running scan with **Ctrl+C** (SIGINT). Jarvis2 will shut down Workers, Archiver,
and any managed `Jarvis-Redis:<scan>` it started. Prefer not to use **Ctrl+Z** — suspend
leaves the job half-alive; Ctrl+Z is ignored during a scan for that reason.

Exit codes: **0** success (also version / empty `ps`); **1** any sample failed (including partial),
kill incomplete, or runtime error; **2** usage/config; **130** interrupted.

`--plot` currently expects a **JarvisPLOT YAML**, not the scan task YAML (scan-driven plot
scenes are partial). Portal / Operas / project tools are native `Jarvis2` subcommands.

Outputs under `outputs/<scan>/` (example scan name depends on the YAML):

- `DATABASE/samples.hdf5` — full **observables** JSON rows (params + module outputs + `LogL` + optional paths)
- `DATABASE/samples.csv` — full-column CSV export of the same records (post-run plot emit)
- `SAMPLE/<bucket>/<uuid>/` or `SAMPLE/<bucket>.tar.gz` — per-sample artifacts (`save: true` files, logs)
- `run_summary.json` / `.csv` / `.txt` — end-of-run counters + throughput
- Console / main log ends with a **`[Scan Performance]`** block (`samples / sec`, `avg sample (sec)`, …)
- Logs under `logs/<scan>/` by component: `core.log`, `worker-00.log`…, `archiver.log`
  (style from packaged `jarvishep2/card/logging.yaml`; per-sample detail in `SAMPLE/.../Sample_running.log`)
- Console: default **INFO** on screen; `Jarvis2 run … --console-level WARNING`; `--silence` / `-s` turns screen off
- Checkpoints: `<project-root>/checkpoints/<scan>/<sampler>/state.pkl`

## Uninstall

```bash
python3 -m pip uninstall jarvishep2
```
