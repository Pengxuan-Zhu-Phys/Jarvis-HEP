# Jarvis-HEP V2 (`jarvishep2`)

Independent Python package for the Jarvis-HEP 2.0 distributed runtime (Redis + spawn Workers).

- Installation: [INSTALL.md](INSTALL.md)
- Design docs: `~/Jarvis-Workshop/Jarvis-Books/Jarvis-HEP V2/`
- Task-YAML reference: `Jarvis-Books/Jarvis-HEP V2/YAML_REFERENCE_2.0.md`
- Portal IO design: `Jarvis-Books/Jarvis-HEP V2/DESIGN_PORTAL_IO_2.0.md`
- Agent bridge (machine-readable API for Jarvis-Agent): `Jarvis-Books/Jarvis-HEP V2/DESIGN_AGENT_BRIDGE_2.0.md`

**Plugin packages** (upgrade them to extend HEP without a HEP release):

| Package | HEP bridge | Extra |
|---------|------------|--------|
| Jarvis-HEP-Portal | `io_portal.py` calculator formats | core dependency |
| Jarvis-Operas | `operas.py` operator registry | `[operas]` |
| JarvisPLOT | `plot_bridge.py` / `Jarvis2 --plot` | `[plot]` |

V1 (`jarvishep`) is frozen and must not be imported from this package.