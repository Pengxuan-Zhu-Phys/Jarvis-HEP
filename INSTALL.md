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

Core runtime depends on **`Jarvis-HEP-Portal`** for calculator input/output formats
(JSON, CSV, TSV, DAT, Wolfram, …). Local editable install of Portal is fine during development:

```bash
python3 -m pip install -e ../Jarvis-Portal
python3 -m pip install -e '.[distributed,dev]'
```

Extras:

- `distributed` = `redis`, `msgpack`, `aiofiles`
- `operas` = `Jarvis-Operas` (catalog operators such as `helper.eggbox2d`)
- `plot` = `JarvisPLOT` (YAML-driven plotting; CLI bridge)
- `dev` = `pytest`, `fakeredis`, `colorlog`

```bash
python3 -m pip install -e '.[distributed,operas,plot,dev]'
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
Jarvis2 quickstart.yaml               # run the scan
Jarvis2 --monitor                     # one read-only status snapshot (separate terminal)
Jarvis2 quickstart.yaml --resume      # resume from checkpoint without the 30 s prompt
Jarvis2 path/to/plot.yaml --plot      # render an existing JarvisPLOT scene (plot extra)
```

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
