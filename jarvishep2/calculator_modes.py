"""Shared-PackID calculator multi-mode configuration (D20).

``modes`` is user-facing declaration sugar.  It expands to ordinary dotted
execution steps (``Prospino.ng``), while every mode of one parent shares the
same physical PackID namespace (``Prospino/001``).  Redis and Workers consume
the internal metadata below to serialize same-parent modes and switch a pack's
build only while that PackID is exclusively held.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from copy import deepcopy
from dataclasses import dataclass
from difflib import get_close_matches
from typing import Any


MODE_PARENT = "_mode_parent"
MODE_NAME = "_mode_name"
SHARED_MODE_PACK = "_shared_mode_pack"
PARENT_INSTALLATION = "_parent_installation"
MODE_INSTALLATION = "_mode_installation"

_MODE_INHERITED_FIELDS = frozenset({
    "path", "source", "clone_shadow", "env_setup", "timeout",
    "make_parallel", "selection", "required_modules", "symlink_name",
})
_MODE_PHYSICAL_FIELDS = frozenset({
    "source", "clone_shadow", "make_parallel", "symlink_name",
})


@dataclass(frozen=True)
class ModeIssue:
    code: str
    path: str
    message: str
    hint: str | None = None
    suggestion: str | None = None
    level: str = "error"


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes))


def _command_list(value: Any) -> list[Any]:
    if isinstance(value, str):
        return [value]
    return list(value) if _is_sequence(value) else []


def _module_children(modules: Sequence[Mapping[str, Any]]) -> dict[str, list[str]]:
    children: dict[str, list[str]] = {}
    for item in modules:
        name = str(item.get("name") or "").strip()
        if not name:
            continue
        modes = item.get("modes")
        if _is_sequence(modes):
            children[name] = [
                f"{name}.{str(mode.get('name') or '').strip()}"
                for mode in modes
                if isinstance(mode, Mapping) and str(mode.get("name") or "").strip()
            ]
        else:
            children[name] = [name]
    return children


def mode_info(module: Mapping[str, Any] | None) -> tuple[str, str] | None:
    """Return ``(parent, mode)`` for an expanded shared-pack mode module."""
    if not isinstance(module, Mapping) or not bool(module.get(SHARED_MODE_PACK)):
        return None
    parent = str(module.get(MODE_PARENT) or "").strip()
    mode = str(module.get(MODE_NAME) or "").strip()
    return (parent, mode) if parent and mode else None


def shared_mode_groups(modules: Sequence[Mapping[str, Any]] | None) -> dict[str, list[str]]:
    """Return parent -> declared mode names for expanded Worker modules."""
    grouped: dict[str, list[str]] = defaultdict(list)
    for module in modules or []:
        info = mode_info(module)
        if info is not None:
            parent, mode = info
            if mode not in grouped[parent]:
                grouped[parent].append(mode)
    return dict(grouped)


def _final_mode_config(parent: Mapping[str, Any], mode: Mapping[str, Any]) -> dict[str, Any]:
    parent_name = str(parent.get("name") or "").strip()
    mode_name = str(mode.get("name") or "").strip()
    child: dict[str, Any] = {
        key: deepcopy(value)
        for key, value in parent.items()
        if key not in {"name", "modes", "execution", "installation", "initialization"}
    }
    for key in _MODE_INHERITED_FIELDS:
        if key in parent:
            child[key] = deepcopy(parent[key])
    for key, value in mode.items():
        if key in {"name", "path", "installation", "initialization", *_MODE_PHYSICAL_FIELDS}:
            continue
        child[key] = deepcopy(value)
    # Parent path always anchors the shared physical @PackID runtime. A mode
    # may nevertheless select a different *execution* directory: explicit
    # execution.path wins, otherwise mode.path is its default.
    execution = mode.get("execution")
    if isinstance(execution, Mapping):
        execution = deepcopy(dict(execution))
        if "path" not in execution and mode.get("path") not in (None, ""):
            execution["path"] = deepcopy(mode["path"])
        child["execution"] = execution
    child["name"] = f"{parent_name}.{mode_name}"
    child[MODE_PARENT] = parent_name
    child[MODE_NAME] = mode_name
    child[SHARED_MODE_PACK] = True
    # Keep the public combined representation for consumers which only need a
    # complete declaration. RuntimePreparer uses the two internal lists to run
    # the shared base only for a zero/forced build, not for every mode switch.
    parent_installation = _command_list(parent.get("installation"))
    mode_installation = _command_list(mode.get("installation"))
    child[PARENT_INSTALLATION] = deepcopy(parent_installation)
    child[MODE_INSTALLATION] = deepcopy(mode_installation)
    child["installation"] = parent_installation + mode_installation
    child["initialization"] = _command_list(parent.get("initialization")) + _command_list(mode.get("initialization"))
    return child


def expand_calculator_modes(modules: Sequence[Mapping[str, Any]] | None) -> list[dict[str, Any]]:
    """Expand mode declarations while retaining one PackID namespace per parent."""
    source = [dict(item) for item in (modules or []) if isinstance(item, Mapping)]
    expanded: list[dict[str, Any]] = []
    for parent in source:
        modes = parent.get("modes")
        if not _is_sequence(modes):
            expanded.append(deepcopy(parent))
            continue
        for mode in modes:
            if isinstance(mode, Mapping):
                expanded.append(_final_mode_config(parent, mode))

    children = _module_children(source)
    known = {str(item.get("name") or "").strip() for item in expanded}
    for child in expanded:
        own_name = str(child.get("name") or "").strip()
        raw_deps = child.get("required_modules")
        if not _is_sequence(raw_deps):
            continue
        dependencies: list[str] = []
        for dependency in raw_deps:
            name = str(dependency).strip()
            if not name:
                continue
            if name in children and children[name] != [name]:
                dependencies.extend(
                    dependency
                    for dependency in children[name]
                    if dependency != own_name
                )
            elif name in known or name in children:
                if name != own_name:
                    dependencies.append(name)
            else:
                dependencies.append(name)
        child["required_modules"] = list(dict.fromkeys(dependencies))
    return expanded


def _mode_path_issue(path: str, message: str, suggestion: str) -> ModeIssue:
    return ModeIssue("JV2-MOD-006", path, message, suggestion=suggestion)


def _iter_string_values(value: Any, path: str):
    """Yield every scalar string below a raw multi-mode declaration."""
    if isinstance(value, Mapping):
        for key, child in value.items():
            yield from _iter_string_values(child, f"{path}.{key}")
    elif _is_sequence(value):
        for index, child in enumerate(value):
            yield from _iter_string_values(child, f"{path}[{index}]")
    elif isinstance(value, str):
        yield path, value


def _unsupported_mode_tokens(parent: Mapping[str, Any], path: str) -> list[ModeIssue]:
    """Reject removed per-mode-directory tokens before any runtime work starts."""
    issues: list[ModeIssue] = []
    for value_path, text in _iter_string_values(parent, path):
        for token in ("@Mode", "${mode_dir}"):
            if token not in text:
                continue
            issues.append(ModeIssue(
                "JV2-MOD-006", value_path,
                f"{token} is not supported for multi-mode calculators",
                hint=(
                    "All modes share the parent physical path and its @PackID "
                    "slot; mode is not part of the directory name."
                ),
                suggestion=(
                    "Use the parent path such as runtime/Prospino/@PackID and "
                    "write the mode-specific command directly in that mode's block."
                ),
            ))
    return issues


def validate_calculator_modes(
    config: Mapping[str, Any],
    *,
    resolved_config: Mapping[str, Any] | None = None,
) -> list[ModeIssue]:
    """Validate shared-PackID multi-mode rules on the raw task card."""
    calculators = config.get("Calculators")
    if not isinstance(calculators, Mapping):
        return []
    raw_modules = calculators.get("Modules")
    if not _is_sequence(raw_modules):
        return []
    modules = [item for item in raw_modules if isinstance(item, Mapping)]
    issues: list[ModeIssue] = []
    children = _module_children(modules)
    known = {child for values in children.values() for child in values}
    parent_names: set[str] = set()

    for index, parent in enumerate(modules):
        parent_path = f"Calculators.Modules[{index}]"
        parent_name = str(parent.get("name") or "").strip()
        if parent_name in parent_names:
            issues.append(ModeIssue(
                "JV2-MOD-003", f"{parent_path}.name",
                f"duplicate calculator module name {parent_name!r}",
                suggestion="Give every calculator module a distinct name.",
            ))
        if parent_name:
            parent_names.add(parent_name)
        if parent_name and "." in parent_name:
            issues.append(ModeIssue(
                "JV2-MOD-001", f"{parent_path}.name",
                "calculator module names cannot contain '.' because '.' separates a module from its mode",
                suggestion="Rename the module without '.', then use modes[].name for the mode suffix.",
            ))
        modes = parent.get("modes")
        if not _is_sequence(modes):
            continue
        if "execution" in parent:
            issues.append(ModeIssue(
                "JV2-MOD-002", f"{parent_path}.execution",
                "a multi-mode parent must not define execution; each mode owns its execution block",
                suggestion="Move execution into every modes[] entry.",
            ))
        parent_runtime = str(parent.get("path") or "")
        issues.extend(_unsupported_mode_tokens(parent, parent_path))
        if not bool(parent.get("clone_shadow")):
            issues.append(_mode_path_issue(
                f"{parent_path}.clone_shadow",
                "multi-mode calculators require clone_shadow: true so an acquired PackID owns one mutable runtime directory",
                "Set clone_shadow: true and use a path containing @PackID.",
            ))
        if "@PackID" not in parent_runtime:
            issues.append(_mode_path_issue(
                f"{parent_path}.path",
                "multi-mode calculator path must contain @PackID for its shared physical pack slots",
                "Use a path such as runtime/Prospino/@PackID.",
            ))

        mode_names: set[str] = set()
        output_names: dict[str, str] = {}
        for mode_index, mode in enumerate(modes):
            mode_path = f"{parent_path}.modes[{mode_index}]"
            if not isinstance(mode, Mapping):
                continue
            mode_name = str(mode.get("name") or "").strip()
            if mode_name and "." in mode_name:
                issues.append(ModeIssue(
                    "JV2-MOD-001", f"{mode_path}.name",
                    "mode names cannot contain '.' because '.' separates a module from its mode",
                    suggestion="Use a simple mode name such as ng or ns.",
                ))
            if mode_name in mode_names:
                issues.append(ModeIssue(
                    "JV2-MOD-003", f"{mode_path}.name",
                    f"duplicate mode name {mode_name!r} in calculator {parent_name!r}",
                    suggestion="Give each mode of one calculator a distinct name.",
                ))
            mode_names.add(mode_name)
            for field in sorted(_MODE_PHYSICAL_FIELDS):
                if field in mode:
                    issues.append(_mode_path_issue(
                        f"{mode_path}.{field}",
                        f"a mode cannot override {field} because all modes share one physical PackID pool",
                        f"Move {field} to the parent calculator declaration.",
                    ))
            execution = mode.get("execution")
            if not isinstance(execution, Mapping):
                continue
            outputs = execution.get("output")
            if _is_sequence(outputs):
                for output_index, output in enumerate(outputs):
                    if not isinstance(output, Mapping):
                        continue
                    output_name = str(output.get("name") or "").strip()
                    if output_name in output_names:
                        issues.append(ModeIssue(
                            "JV2-MOD-004", f"{mode_path}.execution.output[{output_index}].name",
                            f"output name {output_name!r} duplicates mode {output_names[output_name]!r}",
                            suggestion="Give every mode output a distinct name.",
                        ))
                    elif output_name:
                        output_names[output_name] = mode_name

    for parent_index, parent in enumerate(modules):
        modes = parent.get("modes")
        candidates: list[tuple[str, Mapping[str, Any]]] = []
        if _is_sequence(modes):
            for mode_index, mode in enumerate(modes):
                if isinstance(mode, Mapping):
                    candidates.append((f"Calculators.Modules[{parent_index}].modes[{mode_index}]", mode))
        else:
            candidates.append((f"Calculators.Modules[{parent_index}]", parent))
        for base_path, candidate in candidates:
            deps = candidate.get("required_modules")
            if deps is None and _is_sequence(modes):
                deps = parent.get("required_modules")
            if not _is_sequence(deps):
                continue
            for dep_index, raw_dependency in enumerate(deps):
                dependency = str(raw_dependency).strip()
                # Bare dependencies are intentionally cross-block: production
                # cards name LibDeps, Operas, and the Parameters pseudo-module.
                # This validator owns only the dotted Parent.mode vocabulary.
                if "." not in dependency:
                    continue
                dependency_parent = dependency.split(".", 1)[0]
                parent_children = children.get(dependency_parent)
                if not parent_children or parent_children == [dependency_parent]:
                    continue
                if dependency in parent_children:
                    continue
                suggestions = get_close_matches(
                    dependency, sorted(parent_children), n=1, cutoff=0.55,
                )
                suffix = f" Did you mean {suggestions[0]!r}?" if suggestions else ""
                issues.append(ModeIssue(
                    "JV2-MOD-005", f"{base_path}.required_modules[{dep_index}]",
                    f"unknown calculator module or mode dependency {dependency!r}.{suffix}",
                    hint="Use Parent for all modes, or Parent.mode for one or more explicit modes.",
                    suggestion="Use an existing calculator module name, or an existing Module.mode name.",
                ))

    pools = calculators.get("Pools")
    if isinstance(pools, Mapping):
        physical_pool_names = {
            str(parent.get("name") or "").strip()
            for parent in modules
            if str(parent.get("name") or "").strip()
        }
        for raw_name in pools:
            name = str(raw_name).strip()
            if name in known and name not in children:
                issues.append(ModeIssue(
                    "JV2-MOD-007", f"Calculators.Pools.{name}",
                    f"pool {name!r} is a mode name; shared-pack modes use their parent calculator pool",
                    suggestion=f"Use Calculators.Pools.{name.split('.', 1)[0]} instead.",
                ))
            elif name not in physical_pool_names:
                suggestions = get_close_matches(name, sorted(physical_pool_names), n=1, cutoff=0.55)
                suffix = f" Did you mean {suggestions[0]!r}?" if suggestions else ""
                issues.append(ModeIssue(
                    "JV2-MOD-008", f"Calculators.Pools.{name}",
                    f"unknown physical calculator pool {name!r}.{suffix}",
                    suggestion="Use the parent calculator name for a multi-mode pool.",
                ))

    # Performance guidance, not a correctness gate. The resolved config is
    # used when EnvReqs defaults supplied the Worker count outside this card.
    runtime_source = resolved_config if isinstance(resolved_config, Mapping) else config
    runtime = runtime_source.get("Runtime")
    workers = 0
    if isinstance(runtime, Mapping):
        try:
            workers = max(0, int(runtime.get("workers") or 0))
        except (TypeError, ValueError):
            workers = 0
    if workers <= 0:
        envreqs = runtime_source.get("EnvReqs")
        v2 = envreqs.get("V2") if isinstance(envreqs, Mapping) else None
        if isinstance(v2, Mapping):
            try:
                workers = max(0, int(v2.get("workers") or v2.get("worker") or 0))
            except (TypeError, ValueError):
                workers = 0
    pool_map = pools if isinstance(pools, Mapping) else {}
    global_slots = calculators.get("make_parallel", 1)
    for parent_index, parent in enumerate(modules):
        parent_name = str(parent.get("name") or "").strip()
        mode_children = children.get(parent_name, [])
        if not parent_name or mode_children == [parent_name] or not mode_children:
            continue
        raw_slots = pool_map.get(
            parent_name, parent.get("make_parallel", global_slots),
        )
        try:
            slots = max(1, int(raw_slots or 1))
        except (TypeError, ValueError):
            continue
        mode_count = len(mode_children)
        recommended = mode_count + workers if workers > 0 else mode_count
        if slots >= recommended:
            continue
        if slots < mode_count:
            message = (
                f"shared pool {parent_name!r} has {slots} PackID slot(s) for "
                f"{mode_count} modes; mode rebuild thrashing is unavoidable"
            )
        else:
            message = (
                f"shared pool {parent_name!r} has {slots} PackID slot(s); "
                f"approximately {recommended} is recommended for {mode_count} "
                f"modes and {workers} Workers"
            )
        issues.append(ModeIssue(
            "JV2-MOD-009",
            f"Calculators.Modules[{parent_index}].modes",
            message,
            hint=(
                "Shared mode affinity cannot eliminate rebuilds when the pool is "
                "smaller than the active mode/Worker demand."
            ),
            suggestion=f"Set Calculators.Pools.{parent_name}: {recommended} when resources allow.",
            level="warning",
        ))
    return issues


def expand_calculator_modes_in_config(config: Mapping[str, Any]) -> dict[str, Any]:
    """Copy *config* with modes expanded while preserving the raw task card."""
    raw_task_card = getattr(config, "raw_task_card", None)
    result = deepcopy(config)
    if not isinstance(result, dict):
        result = deepcopy(dict(config))
    if raw_task_card is not None:
        try:
            result.raw_task_card = deepcopy(raw_task_card)
        except AttributeError:
            pass
    calculators = result.get("Calculators")
    if isinstance(calculators, Mapping):
        block = dict(calculators)
        if _is_sequence(block.get("Modules")):
            block["Modules"] = expand_calculator_modes(block["Modules"])
        result["Calculators"] = block
    return result


__all__ = [
    "MODE_NAME", "MODE_PARENT", "MODE_INSTALLATION", "PARENT_INSTALLATION",
    "SHARED_MODE_PACK", "ModeIssue",
    "expand_calculator_modes", "expand_calculator_modes_in_config", "mode_info",
    "shared_mode_groups", "validate_calculator_modes",
]
