#!/usr/bin/env python3
"""Early task-YAML validation gate (D13.9).

Pure functions — no Redis, no workers, no ``sys.exit``.  Callers translate
:class:`ValidationReport` into logs / CLI exit codes.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Literal

from jarvishep2.contracts.common import sampling_block
from jarvishep2.contracts.mapper import validate_mapper
from jarvishep2.contracts.methods import validate_method_sampling
from jarvishep2.contracts.operational import validate_operational_blocks
from jarvishep2.contracts.variables import validate_variables
from jarvishep2.calculator_modes import validate_calculator_modes
from jarvishep2.diagnostic_guidance import guidance_for
from jarvishep2.distributor import Distributor
from jarvishep2.task_schema import validate_task_card_schema
from jarvishep2.task_config import get_check_modules_settings

Level = Literal["error", "warning"]


def _ascii_diagnostic(value: str | None) -> str | None:
    """Escape text that may embed *user-provided* non-ASCII (paths / messages).

    Author-written guidance (``suggestion`` / ``example`` / ``hint``) must **not**
    go through this helper — otherwise punctuation such as ``…`` ``→`` ``≥``
    becomes ``\\u2026`` in CLI output (D23.14).  ``message`` still escapes because
    it often quotes user YAML values (see JV2-ENC-001).
    """
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
        # path/message may embed user YAML; suggestion/example/hint are author prose.
        path = _ascii_diagnostic(self.path) or self.path
        message = _ascii_diagnostic(self.message) or self.message
        base = f"  [{self.level}] {self.code}  {path}\n          {message}"
        if self.hint:
            base += f"\n          hint: {self.hint}"
        if self.suggestion:
            base += f"\n          suggestion: {self.suggestion}"
        if self.example:
            base += f"\n          example:\n{_indent_example(self.example)}"
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
    # D24.10: every diagnostic carries an executable Jarvis man command.
    if hint is None or (isinstance(hint, str) and "docs/" in hint):
        try:
            from jarvishep2.man_codes import man_command_for

            hint = man_command_for(code, path)
        except Exception:
            hint = hint or "Run: Jarvis man"
    return ValidationIssue(
        level=level,
        code=code,
        # Path/message may quote user YAML (escape). Guidance fields stay Unicode (D23.14).
        path=_ascii_diagnostic(path) or path,
        message=_ascii_diagnostic(message) or message,
        hint=hint,
        suggestion=suggestion,
        example=example,
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
            item.code,
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


def _operas_constant_full_names() -> frozenset[str]:
    """Return registered constant full names, or empty if Operas is unavailable."""
    try:
        from jarvis_operas import build_constant_dicts
    except ImportError:
        return frozenset()
    try:
        return frozenset(str(name) for name in build_constant_dicts())
    except Exception:
        # Validation must never fail open/closed on a transient registry error.
        return frozenset()


def _validate_operas_constant_operators(
    config: Mapping[str, Any],
) -> list[ValidationIssue]:
    """Reject Modules[].operator that names a Jarvis-Operas constant (D23.4).

    ``operator: pdg.mZ`` + ``call_mode: call`` can run at runtime (arity-0
    callable), but it is meaningless as a module step — write the constant in
    an expression instead.
    """
    issues: list[ValidationIssue] = []
    operas = config.get("Operas")
    if not isinstance(operas, Mapping):
        return issues
    modules = operas.get("Modules")
    if not isinstance(modules, list):
        return issues
    constants = _operas_constant_full_names()
    if not constants:
        return issues
    for index, mod in enumerate(modules):
        if not isinstance(mod, Mapping):
            continue
        operator = str(mod.get("operator") or "").strip()
        if not operator or operator not in constants:
            continue
        issues.append(
            issue(
                "error",
                "JV2-OPR-002",
                f"Operas.Modules[{index}].operator",
                f"{operator!r} is a Jarvis-Operas constant, not a module operator. "
                f"Write it directly in an expression (for example "
                f'expression: "{operator}") instead of Operas.Modules.',
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
        Internal flag set only by the ``Jarvis check`` command.
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
            raw_task_card,
            suppress_additional_property_keys=non_ascii_keys,
            check_modules=bool(check_modules),
        )
    )
    report.extend(
        issue(
            item.level, item.code, item.path, item.message,
            hint=item.hint, suggestion=item.suggestion,
        )
        for item in validate_calculator_modes(
            raw_task_card, resolved_config=config,
        )
    )

    cm = bool(check_modules)

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
        report.extend(_validate_operas_call_mode(config))
        report.extend(_validate_operas_constant_operators(config))
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
        # Jarvis check: CSV fixed points (optional) OR sampler-drawn smoke points.
        # Prefer validating Method when present so fallback sampling can run.
        settings = get_check_modules_settings(config)
        default_data = str(settings.get("data") or "").strip()

        if method:
            available_list = Distributor.available_methods()
            if method not in available_list:
                # Invalid Method is only fatal if there is also no CSV path configured;
                # a dedicated check card may leave a leftover Method and still use CSV.
                if not default_data:
                    report.add(
                        issue(
                            "error",
                            "JV2-MTH-003",
                            "Sampling.Method",
                            f"{method!r} is not implemented; for Jarvis check either "
                            f"fix Method or set EnvReqs.V2.check_modules.data. "
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

        if not method and not default_data:
            report.add(
                issue(
                    "error",
                    "JV2-MTH-040",
                    "EnvReqs.V2.check_modules.data",
                    "Jarvis check needs either Sampling.Method (sampler smoke points) "
                    "or EnvReqs.V2.check_modules.data (CSV fixed points)",
                    hint=(
                        "Add Sampling.Method for V1-like N-point smoke, or set "
                        "EnvReqs.V2.check_modules.data to a CSV path."
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
    report.extend(_validate_operas_call_mode(config))
    report.extend(_validate_operas_constant_operators(config))
    # D22: optional Sampling.Mapper (list of {name, expression}, closed namespace).
    report.extend(validate_mapper(config, method=method or None))

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
