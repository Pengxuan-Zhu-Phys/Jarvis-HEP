# Jarvis-HEP V2 (`jarvishep2`)

Independent Python package for the Jarvis-HEP 2.0 distributed runtime (Redis + spawn Workers).

- Installation: [INSTALL.md](INSTALL.md)
- Design docs: `~/Jarvis-Workshop/Jarvis-Books/Jarvis-HEP V2/`
- Task-YAML reference: `Jarvis-Books/Jarvis-HEP V2/YAML_REFERENCE_2.0.md`
- Portal IO design: `Jarvis-Books/Jarvis-HEP V2/DESIGN_PORTAL_IO_2.0.md`
- Agent bridge (machine-readable API for Jarvis-Agent): `Jarvis-Books/Jarvis-HEP V2/DESIGN_AGENT_BRIDGE_2.0.md`

Calculator file formats are provided by **Jarvis-HEP-Portal** (`jarvishep2/io_portal.py`).
Upgrading Portal can add new IO types without changing this package.

V1 (`jarvishep`) is frozen and must not be imported from this package.