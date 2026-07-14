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
Jarvis2 run quickstart.yaml           # preferred
Jarvis2 quickstart.yaml               # legacy alias → run
Jarvis2 check path/to/check.yaml      # fixed-point calculator smoke
Jarvis2 monitor                       # one read-only status snapshot
Jarvis2 plot path/to/scene.yaml       # render JarvisPLOT scene (plot extra)
Jarvis2 portal man                    # list Portal formats (same as jportal man)
Jarvis2 portal man slha               # format manual
Jarvis2 portal path/to/io.yaml        # run Portal IO YAML (same as jportal file)
Jarvis2 operas list                   # list Operas operators

# legacy (still work):
Jarvis2 quickstart.yaml --resume
Jarvis2 path/to/plot.yaml --plot      # deprecated warning → prefer `Jarvis2 plot`
```

Stop a running scan with **Ctrl+C** (SIGINT). Jarvis2 will shut down Workers, Archiver,
and any managed `Jarvis-Redis:<scan>` it started. Prefer not to use **Ctrl+Z** — suspend
leaves the job half-alive; Ctrl+Z is ignored during a scan for that reason.

Exit codes: **0** all samples succeeded; **1** any sample failed (including partial) or runtime
error; **2** usage/config; **130** interrupted.

`--plot` currently expects a **JarvisPLOT YAML**, not the scan task YAML. Scan-to-plot generation
and post-run plotting are planned follow-ups. Portal and Operas discovery currently use their
package CLIs (`jportal`, `jopera`) until native `Jarvis2 portal …` / `Jarvis2 operas …`
subcommands land.

Outputs under `outputs/quickstart/`:

- `DATABASE/samples.hdf5` — archived records (params, observables, `LogL`)
- `SAMPLE/<uuid>/` — per-sample artifacts (empty for opera-only successes)
- `run_summary.json` / `.csv` / `.txt` — end-of-run summary
- top-level process logs in `./logs/`; checkpoints in `<project-root>/checkpoints/<scan>/<sampler>/state.pkl`

## Uninstall

```bash
python3 -m pip uninstall jarvishep2
```
