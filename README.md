# Jarvis-HEP V2 (`jarvishep2`)

Independent Python package for the Jarvis-HEP 2.0 distributed runtime (Redis + spawn Workers).

**Where we are (2026-07-14, `jarvis2` / `64d7486`):** D0–D7 runtime core + D10 AdaptiveLevelSet
core + D12.0–D12.2 + **D12.8** (SAMPLE buckets, process Archiver, direct handoff, process titles,
managed Redis, clean SIGINT) are in. **Agent bridge (D8) is parked.** Next focus: **D11**
UI/integration and remaining **D12** (flowchart / project tools / catalog). Authoritative ledger:
`Jarvis-Books/Jarvis-HEP V2/V2_DISTRIBUTED_PLAN.md`.

- Installation: [INSTALL.md](INSTALL.md)
- Design docs: `~/Jarvis-Workshop/Jarvis-Books/Jarvis-HEP V2/`
- Execution plan (open work): `Jarvis-Books/Jarvis-HEP V2/V2_DISTRIBUTED_PLAN.md`
- Prototype closeout + D12 theme: `Jarvis-Books/Jarvis-HEP V2/PROTOTYPE_CLOSEOUT_REVIEW_2026-07-14.md`
- Task-YAML reference: `Jarvis-Books/Jarvis-HEP V2/YAML_REFERENCE_2.0.md`
- UI/integration review: `Jarvis-Books/Jarvis-HEP V2/USER_INTERFACE_INTEGRATION_REVIEW_2026-07-13.md`
- Portal IO design: `Jarvis-Books/Jarvis-HEP V2/DESIGN_PORTAL_IO_2.0.md`
- Agent bridge: `Jarvis-Books/Jarvis-HEP V2/DESIGN_AGENT_BRIDGE_2.0.md`

**Plugin packages** (upgrade them to extend HEP without a HEP release):

| Package | HEP bridge | Install |
|---------|------------|--------|
| Jarvis-HEP-Portal | `io_portal.py` calculator formats | **core dependency** |
| Jarvis-Operas | `operas.py` operator registry + expression functions | **core dependency** (D12.0) |
| JarvisPLOT | `plot_bridge.py` / `Jarvis2 --plot` | optional extra `[plot]` |

Current integration depth:

- **Portal**: production calculator I/O uses the Portal registry. The installed bridge exposes
  CSV, DAT, JSON, TSV, and Wolfram; HEP-specific SLHA/xSLHA adapters still need Portal-side
  exposure and end-to-end fixtures.
- **Operas**: registry names such as `helper.eggbox2d` and dotted Python callables work. Strict
  `call_mode` validation, V1 `Operas.Functions`, signature filtering, and sample-logger parity
  are follow-up work.
- **PLOT**: `Jarvis2 <plot.yaml> --plot` renders an existing JarvisPLOT scene. It does not yet
  generate a plot scene from a scan YAML, render the workflow graph, or overlay D10
  `levelset.json` output.

**Runtime defaults worth knowing:**

| Area | Default |
|------|---------|
| Archiver | `mode: process`, `pack_buckets: true` |
| Cleanup / handoff | `direct` (no `staging/` hop; optional `mv_to_staging`) |
| SAMPLE layout | buckets of 200 → `SAMPLE/000001/<uuid>/` then `000001.tar.gz` after archive |
| Process titles | `Jarvis2:<scan>`, `Jarvis2-Worker-N`, `Jarvis2-Archiver`, `Jarvis-Redis:<scan>` |
| Stop | Ctrl+C clean shutdown; Ctrl+Z ignored during a scan |

The current `Jarvis2` CLI is a compatibility surface, not the final user interface. In
particular, `--plot` changes the positional YAML from a scan task into a plot scene, `--pid` is
reserved but unused, and discovery subcommands are not implemented. See the UI review before
adding new CLI behavior.

V1 (`jarvishep`) is frozen and must not be imported from this package.
