# Jarvis-HEP V2 (`jarvishep2`)

Independent Python package for the Jarvis-HEP 2.0 distributed runtime (Redis + spawn Workers).

**Where we are (2026-07-16):** D0–D7 runtime, D11 UI, and **D12.0–D12.6 + D12.8** are in
(including `Jarvis2 project …`, Examples catalog, restricted pack crypto via CLI). Recent
as-built: **FileOperation** SAMPLE save, component logs under `logs/<scan>/`, file-driven
log style (`card/logging.yaml`), **[Scan Performance]** summary, CLI **`-v` / `ps` / `kill`**.
**Agent bridge (D8) is parked.** Authoritative ledger:
`Jarvis-Books/Jarvis-HEP V2/V2_DISTRIBUTED_PLAN.md`.

- Installation + CLI: [INSTALL.md](INSTALL.md) (includes `Jarvis2 -v`, `ps`, `kill`, logging flags)
- Project catalog & **encrypt/decrypt usage**:
  - [INSTALL.md § Project tools](INSTALL.md#project-tools-scaffold-catalog-public--restricted-packs)
  - Design: `Jarvis-Books/Jarvis-HEP V2/components/project_tools.md`
  - Catalog source: `Jarvis-Examples/catalog/README.md`
- Design docs: `~/Jarvis-Workshop/Jarvis-Books/Jarvis-HEP V2/`
- Task-YAML reference: `Jarvis-Books/Jarvis-HEP V2/YAML_REFERENCE_2.0.md`
- Portal IO + FileOperation: `Jarvis-Books/Jarvis-HEP V2/DESIGN_PORTAL_IO_2.0.md`,
  `components/io_system.md`
- Run summary / monitor: `components/monitor.md`

**Plugin packages** (upgrade them to extend HEP without a HEP release):

| Package | HEP bridge | Install |
|---------|------------|--------|
| Jarvis-HEP-Portal | `io_portal.py` calculator **formats** (variable R/W) | **core dependency** |
| Jarvis-Operas | `operas.py` operator registry + expression functions | **core dependency** (D12.0) |
| JarvisPLOT | `plot_bridge.py` / `Jarvis2 plot` | optional extra `[plot]` |

**Runtime defaults worth knowing:**

| Area | Default |
|------|---------|
| Archiver | `mode: process`, `pack_buckets: true`, log `logs/<scan>/archiver.log` |
| FileOperation | per-Worker process; YAML `save: true` → SAMPLE (not Portal) |
| Cleanup / handoff | `direct` (no `staging/` hop; optional `mv_to_staging`) |
| SAMPLE layout | buckets of 200 → `SAMPLE/000001/<uuid>/` then tar after archive |
| DATABASE | `samples.hdf5` full observables JSON rows; `samples.csv` full-column export |
| Process titles | `Jarvis2:<scan>`, `Jarvis2-Worker-N`, `Jarvis2-Archiver`, `Jarvis-Redis:<scan>` |
| Component logs | `logs/<scan>/core.log`, `factory.log`, `sampler.log`, `archiver.log`, `datarecorder.log`, `worker-NN.log` |
| CLI version | `Jarvis2 -v` / `--version` → logo + Author + Version |
| Process inspect | `Jarvis2 ps` / `Jarvis2 kill` (interactive confirm; `--yes` for scripts) |
| Console | default INFO; `--console-level`; `--silence` / `-s` |
| AdaptiveBridson tuning | public knobs: `outer_half_width`, `min_radius`; core/stop width = `outer_half_width / 8` |
| Official catalog | GitHub JSON in Jarvis-Examples (no PyPI catalog package) |
| Restricted packs | `Jarvis2 project fetch NAME --key …` / `pack --encrypt --key …` |

V1 (`jarvishep`) is frozen and must not be imported from this package.

## AdaptiveBridson minimal configuration

Besides the required target expression/value, AdaptiveBridson exposes two
algorithm-tuning controls:

```yaml
Sampling:
  Method: AdaptiveBridson
  AdaptiveBridson:
    target_expression: Omega_h2
    target_value: 0.12
    outer_half_width: 0.02
    min_radius: 0.002
```

The sampler derives `core_half_width` and the `t_max - t_min` stop threshold as
`outer_half_width / 8`. Radius shrink, coverage, bridging, endpoint exploration,
and safety budgets use internal defaults. Historical advanced keys remain
readable so existing scan cards and checkpoints can resume, but new cards should
use the two-knob interface.

At the default `INFO` log level, AdaptiveBridson reports only its start, `r_g`
scale transitions, and final summary. Per-fill band diagnostics, endpoint/root
search, bridging, and densify timing are available at `DEBUG`; abnormal safety
stops and output failures remain `WARNING`.

DataRecorder reports DATABASE progress every 500 rows (plus the forced final
count). Archiver keeps process start/final summaries at `INFO`; per-bucket pack,
loop, drain, and process-exit details are available at `DEBUG`.
