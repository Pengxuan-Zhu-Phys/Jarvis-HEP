# Jarvis Monitor TUI

Status: **accepted** (2026-09-10)
Audience: implementers of `Jarvis monitor`
Non-goal: this document does not describe the in-scan `TaskFactory._MonitorLoop`. That loop stays an internal factory snapshot; the TUI never starts it.

## 1. Product

`Jarvis monitor` is a **separate-process, read-only Textual app**. It looks like btop on the overview and like Jarvis-Agent on the sub-pages: one dense home screen, seven numbered pages, a right-hand detail pane, Jarvis blue/gold chrome.

It must not change the scan. No Redis writes, no signals, no log tailing, no pause/resume, no kill. Stopping a scan remains `Jarvis kill`.

Locked decisions:

| Decision | Choice |
| --- | --- |
| Engine | Textual, lazy-imported from the CLI only |
| Safety | Observer only |
| Launch | `Jarvis monitor` always opens the TUI |
| Visual | Dense overview, spacious sub-pages |
| Pages | 7: Overview, Workers, Factory, Sampler, Calculators, Samples, Host |
| Logs | v1 does not tail files |
| Drill-down | Same-page right detail pane |

## 2. What exists today (and what we keep)

Current `Jarvis monitor` lists scans or prints one plaintext snapshot via `format_monitor_view`. That snapshot is not the product.

Keep:

- `list_active_scans` / sticky `R1` refs from `process_cleanup` — TUI picker and attach target.
- `RedisQueue.snapshot_raw` / `SnapshotReader` / `MonitorView` — the read projection. Extend it; do not bypass it with ad-hoc `HGETALL` in widgets.
- `run_monitor` plaintext path — becomes `Jarvis monitor --once` (and `--json`) for scripts.
- Factory `_MonitorLoop` — untouched. It lives inside the scan process at high frequency and is not a TUI data source.

Throw away as UI: `format_monitor_view` as the interactive experience. `EnvReqs.V2.monitor.hz` is reserved YAML for the factory loop, not the TUI refresh rate.

## 3. Isolation contract

The TUI is a **viewer attached to a live scan**, never a collaborator.

### 3.1 Process

- Own OS process, started by `Jarvis monitor`.
- Own Redis client, using the existing control socket timeout (`CONTROL_SOCKET_TIMEOUT_SEC`).
- Worker / process identity comes from **OS inventory** (`list_jarvis_processes` / `list_active_scans`), never from `SCAN`/`KEYS` of `hep:worker:*` or `hep:proc:*`.
- Never construct a live `TaskFactory` just to start `_MonitorLoop`. Attach is Redis + OS + optional runtime-metadata JSON.
- Never claim `hep:control:lock`.
- Never import Textual (or the monitor package) from Worker / Archiver / Core runtime paths.

### 3.2 Redis: allowed vs forbidden

Allowed, pipelined, on the control client:

- `GET`, `HGETALL`, `LLEN`, `EXISTS`, `TTL`, `PING`, `INFO`
- `HGETALL` of known keys: `hep:proc:core`, `hep:proc:archiver`, `hep:proc:redis`, `hep:worker:status:{id}`, `hep:proc:worker:{id}`, `hep:proc:children:{id}`, `hep:calculator:status`, `hep:sample:stats`, `hep:sample:bucket:meta`, `hep:runtime:metadata_path`
- `LLEN` of known lists: `hep:task_queue`, `hep:archive_queue`, `hep:feedback`, `hep:feedback:chain:{id}` when chain ids are already known from metadata
- `GET` of `hep:control:lock` and its `TTL` (existence + age only)

Forbidden:

- Any write: `SET` / `HSET` / `INCR*` / `RPUSH` / `LPUSH` / `DEL` / `EXPIRE` / `BLPOP` / `BLMOVE` / `EVAL`
- `SCAN` / `KEYS` of worker, result, inflight, or proc namespaces
- `LRANGE` of task / archive / feedback queues (payloads are work, not telemetry)
- `GET hep:results:{uuid}` in v1 (result blobs are large)
- Bumping `hep:{kind}:op_count`

Read-only is a **client convention plus tests**, not a Redis ACL. Tests must wrap the client and fail the run if a write method is called.

### 3.3 Host / filesystem

Allowed:

- `psutil` on **the attached scan's process group only** (control, workers, file-operation children, archiver, managed redis) plus host aggregates (`cpu_percent`, `virtual_memory`, `swap`, `getloadavg`)
- `stat` of the runtime-metadata JSON named by `hep:runtime:metadata_path`, and reading that small JSON

Forbidden in v1:

- Walking every OS process like btop
- Tailing `logs/<scan>/*.log`
- Reading `DATABASE/samples.hdf5` or `SAMPLE/` trees on the refresh path

Host page therefore shows: machine totals + this scan's tree. It is not a system process monitor.

### 3.4 Poll budget

| Stream | Default | Notes |
| --- | --- | --- |
| Redis boards + queue lengths + sample/calc hashes | 2 Hz | One non-transactional pipeline per tick |
| `psutil` for known PIDs + host CPU/MEM | 2 Hz | Same tick as Redis |
| Redis `INFO` | 0.5 Hz | Memory, clients, ops/sec |
| Runtime metadata JSON | 0.2 Hz or on attach | Sampler name, worker count, chain ids, calculator names |

User-adjustable refresh: `[` slower, `]` faster. Allowed values: `0.5`, `1`, `2`, `4` Hz. Cap at 4 Hz. Never 120 Hz.

A tick that times out or errors must **keep the last good frame** and show a stale badge. It must not retry in a tight loop.

## 4. CLI

```text
Jarvis monitor                 # TUI splash (logo + Jarvis ps table)
Jarvis monitor R1              # skip chooser, attach that scan
Jarvis monitor --once [R1]     # no TUI; print the list or one snapshot
Jarvis monitor --once --json   # JSON list / snapshot, for scripts
```

TUI attach rules:

| Invocation | Behaviour |
| --- | --- |
| `Jarvis monitor` | Logo splash. Task list is `list_active_scans()` (same as `Jarvis ps`). `j`/`k` + `enter` attach. Empty list stays on splash with a start hint. |
| `Jarvis monitor R1` | Same splash underneath; session screen on top. `esc` returns to the chooser. Unknown ref stays on splash with the error. |
| not a TTY | Fall back to the old table / snapshot so scripts and tests do not hang. |

If Textual is not installed, TUI mode prints an install hint (`pip install 'Jarvis-HEP[monitor]'` or whatever extra we ship) and exits non-zero. `--once` / `--json` must still work without Textual.

`textual` is an **optional extra**, not a Worker/runtime dependency. Lazy import only inside `dispatch_monitor`.

## 5. Information architecture

Seven pages, numbered like btop. Header and footer are global.

```text
┌ Jarvis Monitor  <scan>  ● <mode>  <elapsed>  <REF>  redis:<host:port> ────── <hz> ┐
│ 1 Overview  2 Workers  3 Factory  4 Sampler  5 Calculators  6 Samples  7 Host    │
├──────────────────────────────────────────────────────────────────────────────────┤
│ page body                                                                        │
├──────────────────────────────────────────────────────────────────────────────────┤
│ Enter: attach  |  J/K: select  |  R: refresh  |  Q: quit                         │
└──────────────────────────────────────────────────────────────────────────────────┘
```

Header (always): scan name, `scan_mode`, elapsed since `started_at`, sticky REF, Redis endpoint, live/stale, refresh Hz.

Footer (always): the frozen key-hint component in `jarvishep2/monitor/hints.py`. Format is locked: bold gold `Key`, dim `: meaning`, groups joined by ` | `. Bindings change per page; the chrome does not. `?` still opens a help overlay.

### 5.1 Overview — dense, btop-like

No row selection. The whole page is a dashboard.

Panels:

- **Samples**: `completed` / `running` / `failed` from `hep:sample:stats`; bar if a total is known (`workers_total` is not a total; prefer sampler bound / archived prefix + remaining queue when honest). Never invent a percentage.
- **Queues**: `task_queue_length`, `archive_queue_length`, `feedback_queue_length` as bars + integers.
- **Workers**: alive/total, stale heartbeat count (age > watchdog `stale_sec` or board TTL).
- **Calculators**: one line per known calculator, `free` / `busy` from `hep:calculator:status`.
- **Sparklines**: samples/min, task-queue depth, host CPU. History is **local to the TUI process** (ring buffer, ~120 points). Not stored in Redis.
- **Liveness row**: core / archiver / redis / control-lock dots from proc boards + `GET hep:control:lock`.
- **Compact host**: one CPU bar, one MEM bar. Full graphs live on page 7.

### 5.2 Workers — table + right detail

Left table, one row per worker id from OS inventory (not from SCAN):

| Column | Source |
| --- | --- |
| id | worker id |
| pid | proc board / OS |
| status | `hep:worker:status:{id}` / proc board |
| uuid | `current_uuid` (short in the table, full in detail) |
| hb | heartbeat age from `ts` |
| cpu% / rss | `psutil` on that pid |
| calc | `held_calc_n` |
| alive | OS `is_running` |

Right pane for the highlighted row: full uuid, `file_operation_pid`, children pgids, heartbeat interval, host, held calc count. No log tail. No result payload.

Empty right pane when nothing is selected. Esc clears selection.

### 5.3 Factory

The factory page is the **control plane**: Core, Archiver, Redis, lease. Not a second worker table.

- **Core** (`hep:proc:core`): pid, host, scan_name, run_id, started_at, workers_total, scan_mode, lease_owner, lease_ts, archiver_alive / archiver_pid, redis_alive
- **Archiver** (`hep:proc:archiver`): pid, status, records_written, last_bucket_packed, db_path, scan_name, host
- **Redis** (`hep:proc:redis`): status, port, pid if started-by-us, last_pong_ts, host, title
- **Lock**: `GET`/`TTL` of `hep:control:lock` — owner string and remaining TTL only

Right pane: the selected role's full board key/value list (stringified). Still read-only.

### 5.4 Sampler

- Method, run_id, bounds summary from runtime-metadata JSON (not from guessing Redis)
- Queue depths: task / archive / feedback
- Per-chain feedback `LLEN` only when chain ids are listed in metadata. Do not `SCAN hep:feedback:chain:*`
- Archived contiguous prefix (`hep:archived-prefix:{scan}`) if present
- Checkpoint path mtime if metadata names a file (stat only)

Right pane: selected queue or chain — length, key name. Never the queued payloads.

### 5.5 Calculators

Left: one row per calculator known from task metadata + `hep:calculator:status` fields (`{name}:free`, `{name}:busy`).

Right: free/busy counts; optional `HGETALL calc:busy:{name}` **only for the selected calculator** (bounded by pack pool size, typically ≤ worker count). Show pack id → owner worker. Do not dump every pool every tick.

### 5.6 Samples

Honest limitation: Redis does **not** expose a cheap full UUID catalogue. This page is counters + currently running, not a sample browser.

- Totals: running / completed / failed
- Bucket meta (`hep:sample:bucket:meta`) if present
- Currently running: one row per worker whose `current_uuid` is non-empty
- Archived prefix progress when available

Right pane: selected running uuid, worker id, pid, status. v1 does not fetch `hep:results:{uuid}`.

### 5.7 Host

- Host: per-core CPU bars, memory, swap, load
- Redis `INFO`: `used_memory_human`, `connected_clients`, `instantaneous_ops_per_sec`
- Table: **this scan's process group only** — role, pid, pgid, cpu%, rss, num fds, cmdline short

Right pane: selected process — full cmdline, create time, ppid, num threads.

Not in v1: GPU, disk IO of SAMPLE/, network, other users' processes.

## 6. Interaction

| Key | Action |
| --- | --- |
| `1`…`7` | Jump to page |
| `tab` / `shift+tab` | Next / previous page |
| `j` `k` or arrows | Move table highlight |
| `enter` | Pin the current row in the right pane (highlight can move; pin stays until Esc) |
| `esc` | Clear pin and highlight |
| `[` `]` | Slower / faster refresh |
| `r` | Force one tick now |
| `?` | Help overlay |
| `q` or `ctrl+c` | Quit TUI (scan keeps running) |
| mouse | Click tabs, click rows (Textual default) |

No command composer. This is not Jarvis-Agent's workbench. No slash commands, no kill confirmations, no pause. Scan picking is the splash table (`j`/`k`/`enter`), matching `Jarvis ps`; `Jarvis monitor R1` skips the chooser.

Vim keys `h`/`l` are **not** bound to pages (they conflict with table scrolling). Pages are digits and Tab.

## 7. Visual language

Steal chrome from Jarvis-Agent / `jarvishep2.versioning`, not the transcript layout.

- Palette: logo colours `#2f7fd8` `#73b8f4` `#f6d33f` `#134a8d`, plus status greens/reds/ambers
- Overview: tight panels, sparklines, bars, liveness dots, little whitespace
- Pages 2–7: one table, one vertical rule, one detail pane, comfortable row height
- Empty / stale / disconnected states are first-class (dim copy + next action), not traceback dumps
- Splash: Jarvis-Agent home/hero, not a one-shot pulse. Left is the 8×8 ⬤ logo-monitor (same `card/logo` pattern, colour from `versioning.ICON_COLORS`). Right is the `Jarvis -v` wordmark, tagline, authors, version, then Monitor copy (scan / redis / workers / sampler). `enter` attaches (picker if several scans). Skip the animation if `stdout` is not a TTY or `NO_COLOR` is set; still show the settled frame. Plain-text canvases live in `docs/TUI/*/0-splash.txt`.

Theme lives in a Textual CSS file next to the app, not inline in Python.

## 8. Module layout (when we implement)

Proposed, not created in this change:

```text
jarvishep2/monitor/
  __init__.py          # public: run_tui, extra name
  app.py               # Textual App, header/footer, page switch
  collector.py         # tick loop, Redis+psutil, MonitorFrame dataclass
  theme.tcss
  screens/
    picker.py
    empty.py
    overview.py
    workers.py
    factory.py
    sampler.py
    calculators.py
    samples.py
    host.py
```

Rules:

- Widgets render `MonitorFrame`. They do not talk to Redis.
- `collector.py` is the only Redis/psutil reader. It extends `SnapshotReader` rather than duplicating key names.
- CLI (`dispatch_monitor`) chooses TUI vs `--once` vs `--json`. Runtime packages do not import `jarvishep2.monitor`.

`pyproject.toml`: optional extra, e.g. `monitor = ["textual>=0.86"]`. Core scan extras stay as they are.

## 9. Data frame

One frozen dataclass per tick, roughly:

```text
MonitorFrame
  timestamp, hz, stale
  scan: name, ref, run_id, mode, started_at, redis_endpoint
  queues: task, archive, feedback, chains[{id, length}]
  samples: running, completed, failed, bucket_meta, archived_prefix
  workers: [{id, pid, status, uuid, hb_age, cpu, rss, held_calc_n, alive, children}]
  factory: core{}, archiver{}, redis{}, lock{owner, ttl}
  calculators: [{name, free, busy, busy_packs?}]
  host: cpu_per_core[], mem, swap, load, redis_info{}, procs[{role, pid, ...}]
  sparks: samples_per_min[], task_queue[], cpu[]
```

`MonitorView` today is a subset. Extend it (or map it into `MonitorFrame`) instead of growing a second ad-hoc dict in the App.

## 10. Testing

Must hold before we call it done:

- Collector issues **zero Redis writes** (reuse the write-guard pattern in `tests/test_dashboard_reader.py` / `test_monitor_snapshot.py`).
- Collector never `SCAN`/`KEYS` worker or result namespaces.
- Worker rows come from an injected OS inventory, not from Redis key discovery.
- `--once` / `--json` still work without Textual.
- TUI unit tests drive the collector + frame projection with fakeredis; they do not require a real terminal.
- Optional: one Textual snapshot/pilot test for Overview CSS if we add `textual[dev]` later. Not a gate.

## 11. Non-goals (v1)

- Pause, resume, kill, scale workers
- Log tails, HDF5 browsers, SAMPLE/ directory walkers
- Full UUID history / result observables
- Multi-scan in one TUI (picker attaches to one; quit and re-run to switch)
- Using `TaskFactory.get_monitor_snapshot()` from inside the scan as the TUI feed
- System-wide process list
- GPU / disk / network dashboards
- Writing refresh Hz into the task YAML

## 12. Implementation order (next, not this change)

1. Extra + `jarvishep2.monitor` package + `MonitorFrame` collector on fakeredis.
2. CLI: TUI default, `--once`/`--json` preserved, picker/empty/attach.
3. Overview (dense) with local sparklines.
4. Workers and Host (tables + right pane + psutil).
5. Factory, Sampler, Calculators, Samples.
6. Isolation tests, then visual polish (palette, splash, help overlay).

Do not ship a page that reads payloads to look busy. Empty-but-honest beats a SCAN.
