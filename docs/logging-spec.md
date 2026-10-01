# Jarvis-HEP logging specification

This document is the rule book for log output in Jarvis-HEP V2, with a focus
on samplers. It has two parts:

1. **The Jarvis log format.** The V1 layout is the Jarvis look and is kept
   exactly. Nothing in this document changes it.
2. **What a sampler must log.** Which events every sampler records, at which
   level, with which wording, so that every sampler's log reads the same way.

Every new sampler, including samplers written with an AI assistant or by
users as plug-ins, must follow it. A contract test enforces the required
events (see [Enforcement](#enforcement)).

## 1. The Jarvis log format (V1 contract)

### 1.1 Record layout

Every record is rendered as:

```text

·•· <module label> 
	-> <MM-DD HH:mm:ss.SSS> - [<LEVEL>] >>> 
<message>
```

- A blank line, then the bullet `·•·`, a space, the module label, and a
  trailing space.
- A tab, `-> `, the timestamp with milliseconds, ` - [LEVEL] >>> `.
- The message on its own line(s).
- DATABASE / `samples.hdf5` records (module `Jarvis-HEP.DataRecorder`) use the
  bullet `Ϡ` instead of `·•·`.
- Records marked `raw` (calculator screen output) are written as-is, without
  the header.

This is the V1 layout (`jarvishep/core.py`, `custom_format`). In V2 it is
defined once in `jarvishep2/card/logging.yaml` and rendered by
`JarvisContextFormatter`. Tests pin it. **Do not change the layout, the
bullets, or the timestamp format, and do not build headers by hand.**

On the console the module label, timestamp, and level are colored. Log files
must not contain color codes. (Today the logo banner at the top of
`core.log` still does; see §4.)

### 1.2 Colors

The module label is colored by component. **Jarvis yellow and Jarvis blue are
the Jarvis colors (from the logo) and have top priority**: they belong to the
control process and the samplers, and no other component may use them.

| Component | Label | Color |
| --- | --- | --- |
| Control process | `Jarvis-HEP` | Jarvis yellow `#f6d33f` |
| Samplers | `Jarvis-HEP.Sampler.*` | Jarvis blue `#2f7fd8` |
| Factory | `Jarvis-HEP.Factory` | magenta `#d670d6` |
| Archiver | `Jarvis-HEP.Archiver` | green `#35c98a` |
| DATABASE writer | `Jarvis-HEP.DataRecorder` | green `#35c98a` (told apart by the `Ϡ` bullet) |
| Workers and samples | `Jarvis-HEP.Worker.NN`, `Sample@…` | lavender `#a78bfa` |

The timestamp stays green and the level keeps its level color, as in V1.
Colors are set in `card/logging.yaml` under `process.module_colors`; a new
component gets its own entry there, never a reuse of yellow or blue.

### 1.3 How to get a logger

Always log through the Jarvis logger. Never use `print()` or a bare
`logging.getLogger()` for anything a user should see.

```python
from jarvishep2.logging import get_jarvis_logger

self._logger = get_jarvis_logger("sampler.random")   # label: Jarvis-HEP.Sampler.Random
```

Module labels always start with `Jarvis-HEP` and use dots only:

| Component | Label | File |
| --- | --- | --- |
| Control process | `Jarvis-HEP` | `logs/<scan>/core.log` |
| Factory | `Jarvis-HEP.Factory` | `factory.log` |
| Sampler | `Jarvis-HEP.Sampler.<Method>` (sub-parts: `.Inner`, `.Pool`, …) | `sampler.log` |
| Archiver | `Jarvis-HEP.Archiver` | `archiver.log` |
| DATABASE writer | `Jarvis-HEP.DataRecorder` (bullet `Ϡ`) | `datarecorder.log` |
| Worker | `Jarvis-HEP.Worker.<NN>` | `worker-NN.log` |
| One sample | `Sample@<uuid>` … | `SAMPLE/<bucket>/<uuid>/Sample_running.log` |

The logger name decides the file (`sampler.*` → `sampler.log`); a
`Jarvis-HEP.Sampler…` label also routes there. A sampler that logs under any
other name ends up in `core.log` instead of `sampler.log`.

### 1.4 Levels

Jarvis uses levels the way V1 did. The screen shows WARNING and above by
default; files keep everything from DEBUG up.

| Level | Meaning in Jarvis | Shown on screen by default |
| --- | --- | --- |
| `ERROR` | Something failed: the scan, a stage, or a component. | yes |
| `WARNING` | A milestone the user should see while the scan runs (start, ready, every whole percent of progress, finish, files written), **or** a real problem the run recovered from. | yes |
| `INFO` | Detail for later reading: settings tables, every ‰ of progress, per-generation numbers, checkpoint saves. | no (file only) |
| `DEBUG` | Internals useful only when debugging Jarvis itself. | no (file only) |

Rules:

- A real problem logged at WARNING must say so in words, e.g. start with
  `... was adjusted`, `... ignored`, `... retrying`, so it cannot be mistaken
  for a milestone.
- Never log an exception at a level below WARNING.
- An unexpected exception is logged with its traceback: use
  `logger.exception(...)` (or `exc_info=True`) inside the `except` block. The
  traceback goes to the file; keep the first line of the message short and
  readable.

### 1.5 Message wording

Keep the V1 phrasing, with V1's spelling mistakes corrected (`initializaing`
→ `initializing`, `submited` → `submitted`). A user who knows one sampler's
log can read any other.

| Purpose | Pattern (V1 wording) | Example |
| --- | --- | --- |
| Start | `Initializing the <Method> Sampling` | `Initializing the Bridson Sampling` |
| Ready to submit | `WorkerFactory is ready for <Method> sampler` | `WorkerFactory is ready for Grid sampler` |
| Progress | `<‰>‰ of <done>/<total> <what> in <HH:MM:SS.mmm>` | `250‰ of 1000/4000 samples submitted in 00:00:12.481` |
| One-line facts | `<Subject> -> <key> -> <value> \| <key> -> <value>` | `DRAM Chain -> 3 \| accept rate -> 0.27 \| stage-2 accepts -> 41` |
| Table block | `<Title> ->` followed by a two-column table from `format_two_column_log` | `Random Sampler Summary ->` + table |
| Result | `<Method> Sampler obtains <N> samples in <T>` | `Grid Sampler obtains 4000 samples in 00:01:12.004` |
| File written | `<What> saved to <path>` | `Results saved to /…/DATABASE/dynesty_result.csv` |
| Failure | `<Method> Sampler meets error when <doing what> -> <reason>` | `CSV Sampler meets error when reading data/points.csv -> missing column x` |

Formatting helpers (do not reinvent them):

- Durations: `jarvishep2.log_kv.format_duration` → `HH:MM:SS.mmm`.
- Progress: `jarvishep2.log_kv.PermilleProgress`. It logs every ‰ at INFO and
  every whole percent at WARNING, as V1 did.
- Tables: `jarvishep2.log_kv.format_two_column_log(title, rows)`. Keep keys
  short, in plain words, with units in the value (`12.4 s`, `0.025`).

## 2. What every sampler must log

A run has the stages below. Rows marked **required** must appear in every
sampler's `sampler.log`; the contract test checks them.

| # | Event | Level | Required | Content |
| --- | --- | --- | --- | --- |
| S1 | Start | WARNING | **required** | `Initializing the <Method> Sampling` |
| S2 | Settings | INFO | **required** | `<Method> Sampler Settings ->` table: method, seed, number of variables, variable names, key `Bounds` values, resume (yes/no) |
| S3 | Ready | WARNING | **required** | `WorkerFactory is ready for <Method> sampler` |
| S4 | Progress | INFO / WARNING | **required** for samplers with a known total | `PermilleProgress` lines; samplers without a fixed total log S5 instead |
| S5 | Step | INFO | **required** for feedback samplers | one line per generation / barrier / batch with the quantities that drive the algorithm (see §3) |
| S6 | Checkpoint | INFO | **required** when checkpoints are on | each save and each load: `Checkpoint saved to <path> (<reason>)`, `Checkpoint loaded from <path>` |
| S7 | Adjustment | WARNING | when it happens | a setting that was changed, ignored, or clamped, with old → new value and why |
| S8 | Failure | ERROR | when it happens | V1 failure wording, logged with `logger.exception` so the traceback reaches the file |
| S9 | Result | WARNING | **required** | `<Method> Sampler obtains <N> samples in <T>` |
| S10 | Summary | WARNING | **required** | `<Method> Sampler Summary ->` table: stop reason, submitted, completed, failed, elapsed, plus the method fields in §3 |
| S11 | Files | WARNING | when it happens | `<What> saved to <path>` for every result file the sampler writes |

The **stop reason** in S10 is one short phrase, e.g. `point_number reached`,
`converged (dlogz < 0.5)`, `max_generations reached`, `interrupted`,
`all chains finished`.

S1, S6, S9 and S10 are emitted by the sampler base class, so a sampler only
supplies the values (see [Enforcement](#enforcement)). A sampler writes S5
and its method fields itself.

## 3. Method-specific content

| Family | S5 Step line | Extra fields in the S10 Summary |
| --- | --- | --- |
| Random, Grid | — | points generated, points rejected by `selection` |
| CSV | — | file, rows read, rows skipped (and why) |
| Bridson | — | radius, points generated |
| AdaptiveBridson | per generation: radius change, band, cores, points (already logged) | convergence, generations, cores, final points |
| MCMC family (MCMC, ToyMCMC, AMMCMC, DRAM, EnsembleMCMC, DEMCMC, PTMCMC, PTEnsemble) | at a fixed interval: steps done per chain, acceptance rate per chain; AMMCMC/DRAM: covariance updates; DRAM: stage-2 attempts and accepts; PT: swap attempts and accepts per pair of temperatures | proposed, accepted, acceptance rate, rejected outside the prior, failed evaluations, diagnostics files |
| Nested (Dynesty, MultiNest) | at a fixed interval: iteration, ncall, current logZ ± error, dlogz | nlive, iterations, ncall, efficiency, logZ ± error, result files |

"At a fixed interval" means time-based (default every 60 s), not per step, so
long runs do not flood `sampler.log`.

## 4. Status today

From a review of the V2 sources (2026-10). ✓ = present, partial, ✗ = missing.

| Sampler | S1 Start | S2 Settings | S4/S5 Progress | S6 Checkpoint | S9/S10 Result & summary |
| --- | --- | --- | --- | --- | --- |
| Random | ✗ | ✗ | ✓ | ✗ | ✗ |
| Grid | ✗ | partial | ✓ | ✗ | ✗ |
| CSV | ✗ | ✗ | ✗ | ✗ | ✗ (no log lines at all) |
| Bridson | ✓ | ✗ | ✗ | ✗ | ✓ result, ✗ summary |
| AdaptiveBridson | ✓ (own wording, not the V1 phrase) | ✓ | ✓ | ✗ | ✓ |
| MCMC family | ✓ | ✓ | partial (no per-chain acceptance) | ✗ | ✓ |
| Dynesty / MultiNest | partial | partial | from the dynesty engine | partial | ✓ |

Gaps that affect every sampler:

- The seed is never logged, so a log cannot tell how to reproduce a run.
- Periodic checkpoint saves are not logged; only the save on interrupt is.
- No traceback is ever written: no code uses `logger.exception` or
  `exc_info=True`.
- The logo banner is written to `core.log` with its terminal color codes.

## Enforcement

- **Base class.** The sampler base class emits S1, S2, S6, S9 and S10. A
  sampler provides two small hooks: one returning its settings rows (S2) and
  one returning its summary rows and stop reason (S10). The wording and
  tables then come out identical for every sampler.
- **Contract test.** For every method in the sampler catalog, a short run
  checks that `sampler.log` contains S1, S2, S3, S9 and S10 in the Jarvis
  layout. A sampler that skips them, built-in or plug-in, fails the test.

## Decisions

Settled by the maintainer (2026-10):

1. **Sampler label**: `Jarvis-HEP.Sampler.<Method>` (V2), not V1's
   `Jarvis-HEP.<Method>`.
2. **Colors**: per component, with Jarvis yellow for the control process and
   Jarvis blue for samplers (§1.2).
3. **V1 spelling mistakes are corrected**; the V1 wording is kept otherwise.
