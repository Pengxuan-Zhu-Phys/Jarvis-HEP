#!/usr/bin/env python3
"""Truthful run result contract for CLI / run_summary / future Agent JSON (D11.1)."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Mapping


# Process exit codes (human CLI).
EXIT_OK = 0
EXIT_RUN_FAILED = 1
EXIT_USAGE = 2
EXIT_INTERRUPTED = 130


@dataclass
class RunOutcome:
    """Single source of truth for whether a scan "succeeded".

    Product policy (UI review 2026-07-13):
    - any sample failure → non-zero exit (even if some samples completed);
    - all-failed must not report success solely because archive rows were written;
    - interrupt → 130.
    """

    submitted: int = 0
    completed: int = 0
    failed: int = 0
    cancelled: int = 0
    archived: int = 0
    run_id: str = ""
    status: str = "unknown"  # success | partial_failure | failed | interrupted | error
    error_type: str | None = None
    error: str | None = None
    failed_module: str | None = None
    extras: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.submitted = max(0, int(self.submitted or 0))
        self.completed = max(0, int(self.completed or 0))
        self.failed = max(0, int(self.failed or 0))
        self.cancelled = max(0, int(self.cancelled or 0))
        self.archived = max(0, int(self.archived or 0))
        self.run_id = str(self.run_id or "")
        if self.status in {"unknown", "", None}:
            self.status = self._infer_status()

    def _infer_status(self) -> str:
        if self.error_type == "KeyboardInterrupt" or self.status == "interrupted":
            return "interrupted"
        if self.failed > 0 and self.completed <= 0:
            return "failed"
        if self.failed > 0 and self.completed > 0:
            return "partial_failure"
        if self.completed > 0 or (self.submitted == 0 and not self.error):
            return "success"
        if self.error:
            return "error"
        return "failed"

    @property
    def ok(self) -> bool:
        return self.status == "success"

    @property
    def exit_code(self) -> int:
        if self.status == "interrupted":
            return EXIT_INTERRUPTED
        if self.status in {"failed", "partial_failure", "error"}:
            return EXIT_RUN_FAILED
        if self.status == "success":
            return EXIT_OK
        return EXIT_RUN_FAILED

    def __int__(self) -> int:
        """Backward-compatible numeric view: submitted sample count."""
        return int(self.submitted)

    def __index__(self) -> int:
        return int(self.submitted)

    def __gt__(self, other: object) -> bool:
        if isinstance(other, (int, float)):
            return int(self.submitted) > other
        if isinstance(other, RunOutcome):
            return int(self.submitted) > int(other.submitted)
        return NotImplemented

    def __eq__(self, other: object) -> bool:
        if isinstance(other, RunOutcome):
            return self.to_dict() == other.to_dict()
        if isinstance(other, int):
            # Tests historically compared core.run() to an expected submit count.
            return int(self.submitted) == other
        return NotImplemented

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["ok"] = self.ok
        payload["exit_code"] = self.exit_code
        return payload

    @classmethod
    def from_counters(
        cls,
        *,
        submitted: int,
        completed: int,
        failed: int,
        archived: int = 0,
        cancelled: int = 0,
        run_id: str = "",
        interrupted: bool = False,
        error: str | None = None,
        error_type: str | None = None,
        failed_module: str | None = None,
        extras: Mapping[str, Any] | None = None,
    ) -> RunOutcome:
        status = "interrupted" if interrupted else "unknown"
        return cls(
            submitted=submitted,
            completed=completed,
            failed=failed,
            cancelled=cancelled,
            archived=archived,
            run_id=run_id,
            status=status,
            error=error,
            error_type=error_type,
            failed_module=failed_module,
            extras=dict(extras or {}),
        )

    @classmethod
    def from_exception(
        cls,
        exc: BaseException,
        *,
        run_id: str = "",
        submitted: int = 0,
        failed_module: str | None = None,
    ) -> RunOutcome:
        if isinstance(exc, KeyboardInterrupt):
            return cls(
                submitted=submitted,
                run_id=run_id,
                status="interrupted",
                error_type="KeyboardInterrupt",
                error=str(exc) or "interrupted",
                failed_module=failed_module,
            )
        return cls(
            submitted=submitted,
            run_id=run_id,
            status="error",
            error_type=type(exc).__name__,
            error=str(exc),
            failed_module=failed_module,
        )


__all__ = [
    "EXIT_INTERRUPTED",
    "EXIT_OK",
    "EXIT_RUN_FAILED",
    "EXIT_USAGE",
    "RunOutcome",
]
