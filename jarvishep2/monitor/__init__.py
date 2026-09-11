"""Read-only Jarvis Monitor TUI (entry splash first; pages later)."""

from __future__ import annotations

TEXTUAL_EXTRA = "Jarvis-HEP[monitor]"


def run_tui(
    *,
    scan_ref: str | None = None,
    scan_lister: object | None = None,
) -> int:
    """Open the Textual monitor. Lazy-imports Textual so Workers never load it."""
    from jarvishep2.monitor.app import MonitorApp
    from jarvishep2.run_outcome import EXIT_OK

    MonitorApp(scan_ref=scan_ref, scan_lister=scan_lister).run()
    return EXIT_OK


__all__ = ["TEXTUAL_EXTRA", "run_tui"]
