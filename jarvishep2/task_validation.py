#!/usr/bin/env python3
"""Early task-YAML validation gate (D13.9).

Pure functions — no Redis, no workers, no ``sys.exit``.  Callers translate
:class:`ValidationReport` into logs / CLI exit codes.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Literal

from jarvishep2.contracts.common import is_check_modules_mode, sampling_block
from jarvishep2.contracts.methods import validate_method_sampling
from jarvishep2.contracts.operational import validate_operational_blocks
from jarvishep2.contracts.variables import validate_variables
from jarvishep2.diagnostic_guidance import guidance_for
from jarvishep2.distributor import Distributor
from jarvishep2.task_schema import validate_task_card_schema
from jarvishep2.task_config import get_check_modules_settings

Level = Literal["error", "warning"]


def _ascii_diagnostic(value: str | None) -> str | None:
    """Escape user-provided non-ASCII text before it reaches logs or tables."""
    if value is None:
        return None
    return value.encode("ascii", "backslashreplace").decode("ascii")


@dataclass(frozen=True)
class ValidationIssue:
    """One diagnostic produced by the config gate."""

    level: Level
    code: str
    path: str
    message: str
    hint: str | None = None
    suggestion: str | None = None
    example: str | None = None

    def format_line(self) -> str:
        code = _ascii_diagnostic(self.code) or self.code
        path = _ascii_diagnostic(self.path) or self.path
        message = _ascii_diagnostic(self.message) or self.message
        hint = _ascii_diagnostic(self.hint)
        suggestion = _ascii_diagnostic(self.suggestion)
        example = _ascii_diagnostic(self.example)
        base = f"  [{self.level}] {code}  {path}\n          {message}"
        if hint:
            base += f"\n          hint: {hint}"
        if suggestion:
            base += f"\n          suggestion: {suggestion}"
        if example:
            base += f"\n          example:\n{_indent_example(example)}"
        return base


@dataclass
class ValidationReport:
    """Collected validation issues for one task config."""

    issues: list[ValidationIssue] = field(default_factory=list)

    def add(self, item: ValidationIssue) -> None:
        self.issues.append(item)

    def extend(self, items: Sequence[ValidationIssue]) -> None:
        self.issues.extend(items)

    def errors(self) -> list[ValidationIssue]:
        return [i for i in self.issues if i.level == "error"]

    def warnings(self) -> list[ValidationIssue]:
        return [i for i in self.issues if i.level == "warning"]

    @property
    def ok(self) -> bool:
        """True when there are no errors (warnings are allowed)."""
        return not self.errors()

    def promote_warnings_to_errors(self) -> None:
        """Used by ``--strict``: every warning becomes an error."""
        self.issues = [
            ValidationIssue(
                level="error",
                code=i.code,
                path=i.path,
                message=i.message,
                hint=i.hint,
                suggestion=i.suggestion,
                example=i.example,
            )
            if i.level == "warning"
            else i
            for i in self.issues
        ]


class ConfigValidationError(ValueError):
    """Raised when the task card fails the validation gate."""

    def __init__(self, report: ValidationReport) -> None:
        self.report = report
        super().__init__(format_report(report))


def issue(
    level: Level,
    code: str,
    path: str,
    message: str,
    hint: str | None = None,
    suggestion: str | None = None,
    example: str | None = None,
) -> ValidationIssue:
    if suggestion is None:
        suggestion, automatic_example = guidance_for(code, path, message)
        if example is None:
            example = automatic_example
    return ValidationIssue(
        level=level, code=_ascii_diagnostic(code) or code,
        path=_ascii_diagnostic(path) or path,
        message=_ascii_diagnostic(message) or message,
        hint=_ascii_diagnostic(hint), suggestion=_ascii_diagnostic(suggestion),
        example=_ascii_diagnostic(example),
    )


def _indent_example(example: str) -> str:
    return "\n".join(f"            {line}" for line in example.splitlines())


def _ellipsize(value: str, width: int) -> str:
    """Keep a terminal-table cell readable without losing the detailed entry."""
    value = " ".join(value.split())
    return value if len(value) <= width else value[: width - 3] + "..."


def _format_issue_summary_table(issues: Sequence[ValidationIssue]) -> str:
    """Render a compact, plain-text table suitable for stderr and log files."""
    rows = [
        (
            str(index),
            _ascii_diagnostic(item.code) or item.code,
            _ellipsize(_ascii_diagnostic(item.path) or item.path, 38),
            _ellipsize(_ascii_diagnostic(item.message) or item.message, 72),
        )
        for index, item in enumerate(issues, start=1)
    ]
    headers = ("#", "Code", "YAML path", "Problem")
    widths = [
        max((len(headers[index]), *(len(row[index]) for row in rows)), default=len(headers[index]))
        for index in range(len(headers))
    ]

    def row(values: tuple[str, ...]) -> str:
        return "|" + "|".join(
            f" {value:<{width}} " for value, width in zip(values, widths)
        ) + "|"

    return "\n".join([
        "Error summary (reference only):",
        "+" + "+".join("-" * (width + 2) for width in widths) + "+",
        row(headers),
        "+" + "+".join("-" * (width + 2) for width in widths) + "+",
        *(row(values) for values in rows),
        "+" + "+".join("-" * (width + 2) for width in widths) + "+",
        "This table is only a quick reference. Please consult the Jarvis-HEP YAML settings documentation before editing the task card.",
    ])


def format_report(report: ValidationReport) -> str:
    n_err = len(report.errors())
    n_warn = len(report.warnings())
    header = f"Config validation failed ({n_err} error(s), {n_warn} warning(s)):"
    if report.ok and n_warn:
        header = f"Config validation warnings ({n_warn}):"
    elif report.ok:
        return "Config validation successful."
    lines = [header]
    errors = report.errors()
    if len(errors) >= 2:
        lines.extend(["", _format_issue_summary_table(errors), "", "Detailed diagnostics:"])
    for item in report.issues:
        lines.append(item.format_line())
    return "\n".join(lines)


def format_report_success(report: ValidationReport) -> str:
    n_warn = len(report.warnings())
    if n_warn:
        return format_report(report)
    return "Config validation successful. The task YAML meets the V2 contract."


def raise_if_errors(report: ValidationReport) -> None:
    if not report.ok:
        raise ConfigValidationError(report)


def _validate_dead_keys(config: Mapping[str, Any]) -> list[ValidationIssue]:
    """Known dead / ignored user-facing keys (warnings)."""
    issues: list[ValidationIssue] = []
    runtime = config.get("Runtime")
    if isinstance(runtime, Mapping) and "Subprocess" in runtime:
        issues.append(
            issue(
                "warning",
                "JV2-DEAD-001",
                "Runtime.Subprocess",
                "key is ignored in V2 (internal dead key)",
            )
        )

    return issues


def _yaml_path_part(path: str, part: Any) -> str:
    """Append one parsed YAML key or list position to a diagnostic path."""
    return f"{path}[{part}]" if isinstance(part, int) else f"{path}.{part}"


def _non_ascii_positions(value: str) -> list[int]:
    """Return one-based character positions of non-ASCII code points."""
    return [index for index, char in enumerate(value, start=1) if ord(char) > 0x7F]


def validate_task_card_encoding(config: Any) -> list[ValidationIssue]:
    """Reject non-ASCII keys and string values in the parsed task card.

    YAML comments are intentionally absent here: PyYAML removes them before the
    document reaches this gate, which keeps Chinese explanations legal in
    ``#`` comments while task-card data itself remains portable ASCII.
    """
    issues: list[ValidationIssue] = []

    def add_issue(path: str, value: str, *, location: str) -> None:
        positions = ", ".join(str(position) for position in _non_ascii_positions(value))
        issues.append(issue(
            "error", "JV2-ENC-001", path,
            f"non-ASCII character(s) at position(s) {positions} in {location} {value!r}",
            hint="Task-card keys and string values must be ASCII. Put Chinese text in a # comment, which is fully supported.",
            suggestion="Replace this text with English ASCII, or move Chinese explanatory text to a # comment.",
        ))

    def visit(value: Any, path: str) -> None:
        if isinstance(value, Mapping):
            for key, child in value.items():
                key_path = _yaml_path_part(path, key)
                if isinstance(key, str) and _non_ascii_positions(key):
                    add_issue(key_path, key, location="key")
                visit(child, key_path)
        elif isinstance(value, str):
            if _non_ascii_positions(value):
                add_issue(path, value, location="string value")
        elif isinstance(value, Sequence):
            for index, child in enumerate(value):
                visit(child, _yaml_path_part(path, index))

    visit(config, "$")
    return issues


def _non_ascii_key_names(config: Any) -> frozenset[str]:
    """Return raw user-authored keys whose schema unknown-key error is redundant."""
    names: set[str] = set()

    def visit(value: Any) -> None:
        if isinstance(value, Mapping):
            for key, child in value.items():
                if isinstance(key, str) and _non_ascii_positions(key):
                    names.add(key)
                visit(child)
        elif isinstance(value, Sequence) and not isinstance(value, str):
            for child in value:
                visit(child)

    visit(config)
    return frozenset(names)


def _validate_operas_call_mode(config: Mapping[str, Any]) -> list[ValidationIssue]:
    issues: list[ValidationIssue] = []
    operas = config.get("Operas")
    if not isinstance(operas, Mapping):
        return issues
    modules = operas.get("Modules")
    if not isinstance(modules, list):
        return issues
    allowed = frozenset({"call", "acall"})
    for index, mod in enumerate(modules):
        if not isinstance(mod, Mapping) or "call_mode" not in mod:
            continue
        raw = str(mod.get("call_mode") or "").strip().lower()
        if raw not in allowed:
            issues.append(
                issue(
                    "error",
                    "JV2-OPR-001",
                    f"Operas.Modules[{index}].call_mode",
                    f"{mod.get('call_mode')!r} invalid; expected one of: call, acall",
                )
            )
    return issues


def _finalize_report(report: ValidationReport, *, strict: bool) -> ValidationReport:
    """Apply the stable CLI ordering after every validation branch."""
    if strict:
        report.promote_warnings_to_errors()
    report.issues.sort(key=lambda item: (item.level != "error", item.path, item.code, item.message))
    return report


def validate_task_config(
    config: Mapping[str, Any],
    *,
    strict: bool = False,
    check_modules: bool | None = None,
) -> ValidationReport:
    """Validate a loaded task config (post-``load_task_yaml``).

    Parameters
    ----------
    config:
        Mapping produced by :func:`jarvishep2.task_config.load_task_yaml` (or
        equivalent in-memory structure).
    strict:
        Promote warnings to errors.
    check_modules:
        When True, treat as check-modules mode (ignore Sampling.Method).
        When None, detect from ``Sampling.mode``.
    """
    report = ValidationReport()
    if not isinstance(config, Mapping):
        report.add(
            issue(
                "error",
                "JV2-LOAD-001",
                "$",
                f"task config must be a mapping, got {type(config).__name__}",
            )
        )
        return _finalize_report(report, strict=strict)

    # Keep the raw task document for user-facing encoding diagnostics. The
    # normalized config contains injected copies such as ``scan_name`` and
    # ``task_result_dir`` that the user cannot edit in their YAML.
    raw_task_card = getattr(config, "raw_task_card", config)
    report.extend(validate_task_card_encoding(raw_task_card))
    non_ascii_keys = _non_ascii_key_names(raw_task_card)

    # Structural validation must also see the user-authored document, rather
    # than runtime-only fields injected by ``load_task_yaml``.
    report.extend(
        validate_task_card_schema(
            raw_task_card, suppress_additional_property_keys=non_ascii_keys,
        )
    )

    cm = (
        bool(check_modules)
        if check_modules is not None
        else is_check_modules_mode(config)
    )

    sampling = sampling_block(config)
    if sampling is None:
        if not cm:
            report.add(
                issue(
                    "error",
                    "JV2-MTH-001",
                    "Sampling",
                    "Sampling block is required for scan runs",
                )
            )
        # Operational checks still useful.
        report.extend(validate_operational_blocks(config))
        report.extend(_validate_dead_keys(config))
        report.extend(_validate_operas_call_mode(config))
        return _finalize_report(report, strict=strict)

    if not isinstance(sampling, Mapping):
        report.add(
            issue(
                "error",
                "JV2-MTH-001",
                "Sampling",
                f"expected a mapping, got {type(sampling).__name__}",
            )
        )
        report.extend(validate_operational_blocks(config))
        return _finalize_report(report, strict=strict)

    method_raw = sampling.get("Method")
    method = str(method_raw).strip() if method_raw is not None else ""

    if cm:
        # check-modules: CSV fixed points (optional) OR sampler-drawn smoke points.
        # Prefer validating Method when present so fallback sampling can run.
        data = str(sampling.get("data") or sampling.get("points_csv") or "").strip()
        settings = get_check_modules_settings(
            {"EnvReqs": config.get("EnvReqs"), "Sampling": sampling}
        )
        default_data = str(settings.get("data") or "").strip()

        if method:
            available_list = Distributor.available_methods()
            if method not in available_list:
                # Invalid Method is only fatal if there is also no CSV path configured;
                # a dedicated check card may leave a leftover Method and still use CSV.
                if not data and not default_data:
                    report.add(
                        issue(
                            "error",
                            "JV2-MTH-003",
                            "Sampling.Method",
                            f"{method!r} is not implemented; for check-modules either "
                            f"fix Method or set Sampling.data / EnvReqs.V2.check_modules.data. "
                            f"Available: {', '.join(available_list) or 'none'}",
                        )
                    )
                # else: CSV path may still work at runtime
            else:
                report.extend(validate_method_sampling(sampling, method=method))
        elif "Variables" in sampling:
            report.extend(
                validate_variables(
                    sampling,
                    method=None,
                    require_nonempty=False,
                )
            )

        if not method and not data and not default_data:
            report.add(
                issue(
                    "error",
                    "JV2-MTH-040",
                    "Sampling.data",
                    "check-modules needs either Sampling.Method (sampler smoke points) "
                    "or Sampling.data / EnvReqs.V2.check_modules.data (CSV fixed points)",
                    hint=(
                        "Add Sampling.Method for V1-like N-point smoke, or set "
                        "Sampling.data / EnvReqs.V2.check_modules.data to a CSV path."
                    ),
                )
            )
    else:
        if not method:
            available = ", ".join(Distributor.available_methods()) or "none"
            report.add(
                issue(
                    "error",
                    "JV2-MTH-002",
                    "Sampling.Method",
                    "Sampling.Method is required for scan runs "
                    f"(available: {available})",
                )
            )
        else:
            available_list = Distributor.available_methods()
            if method not in available_list:
                report.add(
                    issue(
                        "error",
                        "JV2-MTH-003",
                        "Sampling.Method",
                        f"{method!r} is not implemented in Jarvis-HEP V2. "
                        f"Available: {', '.join(available_list) or 'none'}",
                    )
                )
            else:
                report.extend(validate_method_sampling(sampling, method=method))

    report.extend(validate_operational_blocks(config))
    report.extend(_validate_dead_keys(config))
    report.extend(_validate_operas_call_mode(config))

    return _finalize_report(report, strict=strict)


def validate_or_raise(
    config: Mapping[str, Any],
    *,
    strict: bool = False,
    check_modules: bool | None = None,
) -> ValidationReport:
    """Validate and raise :class:`ConfigValidationError` on errors."""
    report = validate_task_config(
        config, strict=strict, check_modules=check_modules
    )
    raise_if_errors(report)
    return report


__all__ = [
    "ConfigValidationError",
    "ValidationIssue",
    "ValidationReport",
    "format_report",
    "format_report_success",
    "issue",
    "raise_if_errors",
    "validate_task_card_encoding",
    "validate_or_raise",
    "validate_task_config",
]
