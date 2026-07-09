# Jarvis-HEP V2 — Installation

Quick install guide for the V2 distributed runtime (`jarvishep2`, CLI `Jarvis2`).
For the architecture see `~/Jarvis-Workshop/Jarvis-Books/Docs/DESIGN_2.0_DISTRIBUTED.md`;
for the task-YAML schema see `Jarvis-Books/Docs/YAML_REFERENCE_2.0.md`.

## Requirements

- **Python ≥ 3.10** (Workers use the `spawn` multiprocessing context; macOS and Linux supported)
- **Redis server** for distributed runs (`Runtime.mode: redis`). Tests do *not* need one — the
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

Extras: `distributed` = `redis`, `msgpack`, `aiofiles`; `dev` = `pytest`, `fakeredis`, `colorlog`.

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

> **Always set `Runtime.redis` explicitly in the task YAML.** If the block is omitted, the
> control process silently falls back to an in-process fakeredis while spawned Workers connect
> to `localhost:6379` — the run will hang (see `YAML_REFERENCE_2.0.md`, Appendix A.1).

## Verify

```bash
python3 -m pytest -q          # full suite, no Redis server needed (~4 min, 207 tests)
python3 -m pytest -q tests/test_d0_integration.py tests/test_worker_mvp.py   # quick subset
```

## Quickstart

Save as `quickstart.yaml` (any directory; outputs land in `<project-root>/outputs/<scan-name>/`):

```yaml
project_name: quickstart
Scan:
  name: quickstart
Runtime:
  mode: redis
  workers: 2
  redis:
    host: 127.0.0.1
    port: 6379
    db: 0
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
```

Outputs under `outputs/quickstart/`:

- `DATABASE/samples.hdf5` — archived records (params, observables, `LogL`)
- `SAMPLE/<uuid>/` — per-sample artifacts (empty for opera-only successes)
- `run_summary.json` / `.csv` / `.txt` — end-of-run summary
- top-level process logs in `./logs/`; checkpoints in `<project-root>/checkpoints/<scan>/<sampler>/state.pkl`

## Uninstall

```bash
python3 -m pip uninstall jarvishep2
```
