# Jarvis-HEP V2 (`jarvishep2`)

Independent Python package for the Jarvis-HEP 2.0 distributed runtime (Redis + spawn Workers).

**Where we are (2026-07-16):** D0–D7 runtime, D11 UI, and **D12.0–D12.6 + D12.8** are in
(including `Jarvis2 project …`, Examples catalog, restricted pack crypto via CLI). Recent
as-built: **FileOperation** owns SAMPLE save/copy/delete (Portal keeps format R/W only),
Archiver has its own logger, full `DATABASE/samples.csv` export, end-of-scan **[Scan Performance]**
log (`samples_per_sec`, `avg_sample_sec`, …). **Agent bridge (D8) is parked.** Authoritative
ledger: `Jarvis-Books/Jarvis-HEP V2/V2_DISTRIBUTED_PLAN.md`.

- Installation + CLI: [INSTALL.md](INSTALL.md)
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
| Archiver | `mode: process`, `pack_buckets: true`, log `logs/<scan>/jarvis_archiver.log` |
| FileOperation | per-Worker process; YAML `save: true` → SAMPLE (not Portal) |
| Cleanup / handoff | `direct` (no `staging/` hop; optional `mv_to_staging`) |
| SAMPLE layout | buckets of 200 → `SAMPLE/000001/<uuid>/` then tar after archive |
| DATABASE | `samples.hdf5` full observables JSON rows; `samples.csv` full-column export |
| Process titles | `Jarvis2:<scan>`, `Jarvis2-Worker-N`, `Jarvis2-Archiver`, `Jarvis-Redis:<scan>` |
| Component logs | `logs/<scan>/core.log`, `worker-NN.log`, `archiver.log` |
| Official catalog | GitHub JSON in Jarvis-Examples (no PyPI catalog package) |
| Restricted packs | `Jarvis2 project fetch NAME --key …` / `pack --encrypt --key …` |

V1 (`jarvishep`) is frozen and must not be imported from this package.
