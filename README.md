# Jarvis-HEP V2 (`jarvishep2`)

Independent Python package for the Jarvis-HEP 2.0 distributed runtime (Redis + spawn Workers).

- Installation: [INSTALL.md](INSTALL.md)
- Design docs: `~/Jarvis-Workshop/Jarvis-Books/Jarvis-HEP V2/`
- Prototype closeout review + D12 plan (Calculator/UX parity): `Jarvis-Books/Jarvis-HEP V2/PROTOTYPE_CLOSEOUT_REVIEW_2026-07-14.md`
- Current UI/integration review: `Jarvis-Books/Jarvis-HEP V2/USER_INTERFACE_INTEGRATION_REVIEW_2026-07-13.md`
- Task-YAML reference: `Jarvis-Books/Jarvis-HEP V2/YAML_REFERENCE_2.0.md`
- Portal IO design: `Jarvis-Books/Jarvis-HEP V2/DESIGN_PORTAL_IO_2.0.md`
- Agent bridge (machine-readable API for Jarvis-Agent): `Jarvis-Books/Jarvis-HEP V2/DESIGN_AGENT_BRIDGE_2.0.md`

**Plugin packages** (upgrade them to extend HEP without a HEP release):

| Package | HEP bridge | Extra |
|---------|------------|--------|
| Jarvis-HEP-Portal | `io_portal.py` calculator formats | core dependency |
| Jarvis-Operas | `operas.py` operator registry | `[operas]` |
| JarvisPLOT | `plot_bridge.py` / `Jarvis2 --plot` | `[plot]` |

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

The current `Jarvis2` CLI is a compatibility surface, not the final user interface. In
particular, `--plot` changes the positional YAML from a scan task into a plot scene, `--pid` is
reserved but unused, and discovery subcommands are not implemented. See the review above before
adding new CLI behavior.

V1 (`jarvishep`) is frozen and must not be imported from this package.
