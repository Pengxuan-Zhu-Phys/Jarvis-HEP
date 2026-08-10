#!/usr/bin/env python3
"""V1-compatible Python and CERN ROOT environment requirement checks."""

from __future__ import annotations

import platform
import re
import shlex
import subprocess
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from importlib import metadata as importlib_metadata
from typing import Any


def _version_tuple(value: Any) -> tuple[int, ...]:
    """Parse the numeric portion of Python/ROOT versions (e.g. ``6.30/04``)."""
    numbers = re.findall(r"\d+", str(value))
    return tuple(int(number) for number in numbers) or (0,)


def _meets_minimum(installed: Any, minimum: Any) -> bool:
    left = _version_tuple(installed)
    right = _version_tuple(minimum)
    width = max(len(left), len(right))
    return left + (0,) * (width - len(left)) >= right + (0,) * (width - len(right))


def _minimum_version(specification: Any) -> str | None:
    match = re.fullmatch(r"\s*>=\s*(\d+(?:\.\d+)*)\s*", str(specification or ""))
    return match.group(1) if match else None


@dataclass
class EnvironmentRequirementReport:
    """Outcome of checking the V1-compatible ``EnvReqs`` blocks."""

    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    summary: dict[str, Any] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return not self.errors

    def format(self) -> str:
        return "Environment requirements failed:\n" + "\n".join(
            f"  - {message}" for message in self.errors
        )


class EnvironmentRequirementError(RuntimeError):
    """Raised when declared OS, Python, or CERN ROOT requirements are not met."""

    def __init__(self, report: EnvironmentRequirementReport) -> None:
        self.report = report
        super().__init__(report.format())


def _run_command(command: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        shlex.split(command), capture_output=True, text=True, check=False
    )


def _check_python(requirement: Mapping[str, Any], report: EnvironmentRequirementReport) -> None:
    installed = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    report.summary["Python version"] = installed
    minimum = _minimum_version(requirement.get("version"))
    if requirement.get("version") is not None and minimum is None:
        report.warnings.append(
            f"Python has unsupported version requirement {requirement['version']!r}; skipping it"
        )
    elif minimum and not _meets_minimum(installed, minimum):
        report.errors.append(f"Python {installed} does not meet the requirement >= {minimum}")

    dependencies = requirement.get("Dependencies")
    if not isinstance(dependencies, Sequence) or isinstance(dependencies, (str, bytes)):
        return
    for item in dependencies:
        if not isinstance(item, Mapping):
            report.warnings.append("Python dependency entry is not a mapping; skipping it")
            continue
        name = str(item.get("name") or "").strip()
        if not name:
            report.warnings.append("Python dependency entry without name; skipping it")
            continue
        required = bool(item.get("required", False))
        minimum = _minimum_version(item.get("version"))
        if item.get("version") is not None and minimum is None:
            report.warnings.append(
                f"Python dependency {name} has unsupported version requirement {item['version']!r}; skipping it"
            )
            continue
        try:
            version = importlib_metadata.version(name)
        except importlib_metadata.PackageNotFoundError:
            report.summary[name] = None
            message = f"Python package {name}{'>=' + minimum if minimum else ''} is not installed"
            (report.errors if required else report.warnings).append(message)
            continue
        report.summary[name] = version
        if minimum and not _meets_minimum(version, minimum):
            message = f"Python package {name}=={version} does not meet the requirement >= {minimum}"
            (report.errors if required else report.warnings).append(message)


def _check_os(requirements: Sequence[Any], report: EnvironmentRequirementReport) -> None:
    """Check the V1 ``EnvReqs.OS`` list against the current platform.

    V1 matched ``platform.system()`` case-insensitively and compared the
    numeric portions of ``platform.release()`` with a ``>=`` requirement.
    Keep that behavior, but return diagnostics instead of terminating the
    process so callers can render one complete preflight report.
    """
    current_os = platform.system()
    current_version = platform.release()
    report.summary["OS"] = f"{current_os}-{current_version}"
    if not requirements:
        return

    matching: Mapping[str, Any] | None = None
    for item in requirements:
        if not isinstance(item, Mapping):
            report.warnings.append("OS requirement entry is not a mapping; skipping it")
            continue
        if str(item.get("name") or "").strip().casefold() == current_os.casefold():
            matching = item
            break

    if matching is None:
        report.errors.append(
            f"Operating system {current_os} is not listed in EnvReqs.OS"
        )
        return

    minimum = _minimum_version(matching.get("version"))
    if minimum is None:
        report.errors.append(
            f"OS requirement for {current_os} must use a supported >= version"
        )
        return
    if not _meets_minimum(current_version, minimum):
        report.errors.append(
            f"OS {current_os} release {current_version} does not meet the requirement >= {minimum}"
        )


def _check_root(requirement: Mapping[str, Any], report: EnvironmentRequirementReport) -> None:
    if not bool(requirement.get("required", False)):
        # A task may deliberately avoid a ROOT version check but still need its
        # configured prefix in a LibDeps build command (`@{ROOT path}`).
        path = requirement.get("path")
        if isinstance(path, str) and path.strip():
            report.summary["ROOT path"] = path.strip()
        return
    try:
        result = _run_command("root-config --version")
    except (FileNotFoundError, ValueError) as exc:
        report.summary["ROOT"] = False
        report.errors.append(f"CERN ROOT is not installed or root-config is unavailable: {exc}")
        return
    if result.returncode != 0 or not result.stdout.strip():
        report.summary["ROOT"] = False
        detail = result.stderr.strip() or f"exit {result.returncode}"
        report.errors.append(f"CERN ROOT is not installed or root-config is unavailable ({detail})")
        return

    root_version = result.stdout.strip()
    report.summary["ROOT"] = True
    report.summary["ROOT version"] = root_version
    minimum = _minimum_version(requirement.get("version"))
    if requirement.get("version") is not None and minimum is None:
        report.warnings.append(
            f"CERN_ROOT has unsupported version requirement {requirement['version']!r}; skipping it"
        )
    elif minimum and not _meets_minimum(root_version, minimum):
        report.errors.append(f"CERN ROOT {root_version} does not meet the requirement >= {minimum}")

    path_command = requirement.get("get_path_command")
    if isinstance(path_command, str) and path_command.strip():
        try:
            path_result = _run_command(path_command)
        except (FileNotFoundError, ValueError) as exc:
            report.errors.append(f"CERN ROOT path command failed: {exc}")
        else:
            if path_result.returncode == 0 and path_result.stdout.strip():
                report.summary["ROOT path"] = path_result.stdout.strip()
            else:
                report.errors.append("CERN ROOT path command failed")
    elif isinstance(requirement.get("path"), str):
        report.summary["ROOT path"] = requirement["path"]

    dependencies = requirement.get("Dependencies")
    if not isinstance(dependencies, Sequence) or isinstance(dependencies, (str, bytes)):
        return
    for item in dependencies:
        if not isinstance(item, Mapping):
            report.warnings.append("CERN ROOT dependency entry is not a mapping; skipping it")
            continue
        name = str(item.get("name") or "").strip() or "unnamed"
        command = item.get("check_command")
        required = bool(item.get("required", False))
        if not isinstance(command, str) or not command.strip():
            report.warnings.append(f"CERN ROOT dependency {name} has no check_command; skipping it")
            continue
        try:
            dependency_result = _run_command(command)
        except (FileNotFoundError, ValueError):
            dependency_result = None
        expected = str(item.get("expected_output") or "").strip()
        found = bool(
            dependency_result is not None
            and dependency_result.returncode == 0
            and dependency_result.stdout.strip() == expected
        )
        report.summary[f"ROOT-{name}"] = found
        if not found:
            message = f"CERN ROOT feature {name} does not meet the requirement"
            (report.errors if required else report.warnings).append(message)


def check_environment_requirements(config: Mapping[str, Any]) -> EnvironmentRequirementReport:
    """Check V1-compatible OS, Python, and CERN ROOT requirement blocks."""
    report = EnvironmentRequirementReport()
    envreqs = config.get("EnvReqs")
    if not isinstance(envreqs, Mapping):
        return report
    python = envreqs.get("Python")
    if isinstance(python, Mapping):
        _check_python(python, report)
    operating_system = envreqs.get("OS")
    if isinstance(operating_system, Sequence) and not isinstance(operating_system, (str, bytes)):
        _check_os(operating_system, report)
    root = envreqs.get("CERN_ROOT")
    if isinstance(root, Mapping):
        _check_root(root, report)
    return report


def raise_if_environment_requirements_fail(config: Mapping[str, Any]) -> EnvironmentRequirementReport:
    report = check_environment_requirements(config)
    if not report.ok:
        raise EnvironmentRequirementError(report)
    return report


__all__ = [
    "EnvironmentRequirementError",
    "EnvironmentRequirementReport",
    "check_environment_requirements",
    "raise_if_environment_requirements_fail",
]
