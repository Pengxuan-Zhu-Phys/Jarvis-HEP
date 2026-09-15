"""Sampler-specific SAMPLES block shared by SIMU and live Monitor views."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any, Mapping

from jarvishep2.monitor.bars import bar_rail, visible_width


BUILTIN_SAMPLERS = (
    "Random",
    "Grid",
    "CSV",
    "Bridson",
    "AdaptiveBridson",
    "MCMC",
    "ToyMCMC",
    "AMMCMC",
    "DRAM",
    "EnsembleMCMC",
    "DEMCMC",
    "PTMCMC",
    "PTEnsemble",
    "Dynesty",
    "MultiNest",
)


@dataclass(frozen=True)
class SamplerProgress:
    kind: str = "none"
    label: str = ""
    current: float | int | None = None
    target: float | int | None = None
    eta: str = "—"


@dataclass(frozen=True)
class SamplerDisplay:
    method: str = "Unknown"
    family: str = ""
    state: str = "initial"
    progress: SamplerProgress = field(default_factory=SamplerProgress)
    metrics: Mapping[str, Any] = field(default_factory=dict)
    saved: int | None = None


@dataclass(frozen=True)
class SamplesBlockData:
    completed: int = 0
    running: int = 0
    failed: int = 0
    rate: str = "—"
    sampler: SamplerDisplay = field(default_factory=SamplerDisplay)


@dataclass(frozen=True)
class SamplesBlockRender:
    border_title: str
    body: str
    border_subtitle: str


def sample_border_title(outer_width: int, method: str) -> str:
    """Keep SAMPLES left and the canonical method on the top-right edge."""
    base = "SAMPLES"
    available = max(len(base), int(outer_width) - 6)
    method = str(method or "Unknown")
    max_method = max(1, available - len(base) - 3)
    if len(method) > max_method:
        method = method[: max(1, max_method - 1)] + "…"
    gap = max(1, available - len(base) - len(method) - 2)
    return f"{base} {'─' * gap} {method}"


def _clip_left(text: str, width: int) -> str:
    if width <= 0:
        return ""
    if len(text) <= width:
        return text
    if width == 1:
        return "…"
    return text[: width - 1] + "…"


def _edge_align(left: str, right: str, width: int) -> str:
    """Align plain text while preserving the complete right-hand value."""
    if width <= 0:
        return ""
    right = _clip_left(str(right), width)
    left_width = max(0, width - len(right) - (1 if right else 0))
    left = _clip_left(str(left), left_width if right else width)
    gap = max(0, width - len(left) - len(right))
    return left + " " * gap + right


def _edge_align_markup(left: str, right: str, width: int) -> str:
    """Align a row whose left side may contain Rich markup."""
    right = _clip_left(str(right), width)
    gap = max(0, width - visible_width(left) - len(right))
    return left + " " * gap + right


def _metric(display: SamplerDisplay, name: str, default: Any = None) -> Any:
    value = display.metrics.get(name, default)
    return default if value is None else value


def _integer(value: Any, default: str = "—") -> str:
    try:
        return str(int(value))
    except (TypeError, ValueError, OverflowError):
        return default


def _number(value: Any, *, digits: int = 3, default: str = "—") -> str:
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return default
    return f"{number:.{digits}g}"


def _percent(value: Any, default: str = "—") -> str:
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return default
    if 0.0 <= number <= 1.0:
        number *= 100.0
    return f"{number:0.0f}%"


def _saved(display: SamplerDisplay) -> str:
    return _integer(display.saved)


def _telemetry_waiting(display: SamplerDisplay) -> bool:
    return (
        str(display.state or "").lower()
        in {"initial", "initializing", "unavailable"}
        and display.progress.current is None
    )


def _common_row(width: int, data: SamplesBlockData) -> str:
    left = f"Acc {data.completed} · RUN {data.running} · FAIL {data.failed}"
    right = f"RATE {data.rate}"
    if len(left) + len(right) > width:
        left = f"Acc {data.completed} R{data.running} F{data.failed}"
    return _edge_align(left, right, width)


def _finite_rows(
    width: int,
    display: SamplerDisplay,
    *,
    label: str,
) -> tuple[str, str] | None:
    progress = display.progress
    try:
        current = max(0.0, float(progress.current))
    except (TypeError, ValueError, OverflowError):
        return None
    try:
        target = float(progress.target)
    except (TypeError, ValueError, OverflowError):
        target = 0.0
    if target <= 0:
        state = str(display.state or "running").upper()
        return (
            _edge_align(f"{label} {_integer(current)} / —", state, width),
            _edge_align("PROGRESS —", "TARGET UNAVAILABLE", width),
        )
    fraction = max(0.0, min(1.0, current / target))
    first = _edge_align(
        f"{label} {_integer(current)} / {_integer(target)}",
        f"{fraction * 100:0.0f}%",
        width,
    )
    eta = f"ETA {progress.eta or '—'}"
    rail_width = min(31, max(8, width - len(eta) - 1))
    second = _edge_align_markup(bar_rail(fraction, rail_width), eta, width)
    return first, second


def _waiting_rows(width: int, display: SamplerDisplay) -> tuple[str, str]:
    state = str(display.state or "initial").upper()
    return (
        _edge_align("SAMPLER TELEMETRY", "—", width),
        _edge_align("WAITING FOR CORE HEARTBEAT", state, width),
    )


def _mcmc_rows(width: int, display: SamplerDisplay, *, role: str, unit: str) -> tuple[str, str]:
    progress = display.progress
    try:
        current = max(0.0, float(progress.current))
        target = float(progress.target)
    except (TypeError, ValueError, OverflowError):
        return _waiting_rows(width, display)
    if target <= 0:
        return _waiting_rows(width, display)
    fraction = max(0.0, min(1.0, current / target))
    count = _integer(_metric(display, "chains", _metric(display, "walkers", 0)), "0")
    first = _edge_align(
        f"{role} {count} · {unit} FLOOR {_integer(current)} / {_integer(target)}",
        f"{fraction * 100:0.0f}%",
        width,
    )
    eta = f"ETA {progress.eta or '—'}"
    rail_width = min(31, max(8, width - len(eta) - 1))
    return first, _edge_align_markup(bar_rail(fraction, rail_width), eta, width)


def _nested_rows(width: int, display: SamplerDisplay) -> tuple[str, str]:
    state = str(display.state or "running").upper()
    nlive = _integer(_metric(display, "nlive"))
    niter = _integer(_metric(display, "niter"))
    ncall = _integer(_metric(display, "ncall"))
    eff = _percent(_metric(display, "efficiency"))
    tolerance = _number(_metric(display, "dlogz_target"), digits=2)
    return (
        _edge_align(f"LIVE {nlive} · ITER {niter}", state, width),
        _edge_align(f"CALLS {ncall} · EFF {eff}", f"DLOGZ≤ {tolerance}", width),
    )


def _subtitle(display: SamplerDisplay) -> str:
    method = display.method
    saved = _saved(display)
    if _telemetry_waiting(display):
        return f"METHOD KNOWN · SAVED {saved}"
    if method == "Random":
        return (
            f"ACCEPTED {_integer(_metric(display, 'accepted'))} · "
            f"SEED {_integer(_metric(display, 'seed'))} · SAVED {saved}"
        )
    if method == "Grid":
        return (
            f"SHAPE {_metric(display, 'shape', '—')} · "
            f"DIMS {_integer(_metric(display, 'dimensions'))} · SAVED {saved}"
        )
    if method == "CSV":
        return (
            f"ACCEPTED {_integer(_metric(display, 'accepted'))} · "
            f"{_metric(display, 'basename', 'CSV')} · SAVED {saved}"
        )
    if method == "Bridson":
        return (
            f"R {_number(_metric(display, 'radius'))} · "
            f"K {_integer(_metric(display, 'max_attempt'))} · SAVED {saved}"
        )
    if method == "AdaptiveBridson":
        return (
            f"R {_number(_metric(display, 'radius'))} · "
            f"CORE {_integer(_metric(display, 'core'))} · SAVED {saved}"
        )
    accept = _percent(_metric(display, "accept_rate"))
    if method == "MCMC":
        return f"ACCEPT {accept} · R-HAT {_number(_metric(display, 'rhat'), digits=3)} · SAVED {saved}"
    if method == "ToyMCMC":
        return f"ACCEPT {accept} · SCALE {_number(_metric(display, 'proposal_scale'))} · SAVED {saved}"
    if method == "AMMCMC":
        adapt = str(_metric(display, "adapt_state", "—")).upper()
        return f"ACCEPT {accept} · ADAPT {adapt} · SAVED {saved}"
    if method == "DRAM":
        retry = _percent(_metric(display, "retry_rate"))
        return f"ACCEPT {accept} · RETRY {retry} · SAVED {saved}"
    if method == "EnsembleMCMC":
        stretch = _number(_metric(display, "stretch_a"))
        return f"ACCEPT {accept} · STRETCH {stretch} · SAVED {saved}"
    if method == "DEMCMC":
        raw_gamma = _metric(display, "de_gamma")
        gamma = "AUTO" if raw_gamma in (None, 0, 0.0, "0", "0.0") else _number(raw_gamma)
        return f"ACCEPT {accept} · GAMMA {gamma} · SAVED {saved}"
    if method in {"PTMCMC", "PTEnsemble"}:
        swaps = (
            f"{_integer(_metric(display, 'swap_accepts'), '0')}/"
            f"{_integer(_metric(display, 'swap_attempts'), '0')}"
        )
        return f"ACCEPT {accept} · SWAP {swaps} · SAVED {saved}"
    if method in {"Dynesty", "MultiNest"}:
        family = "DYNAMIC" if method == "Dynesty" else "STATIC"
        return f"{family} · NLIVE {_integer(_metric(display, 'nlive'))} · SAVED {saved}"
    return f"CUSTOM METHOD · SAVED {saved}"


def render_samples_block(width: int, data: SamplesBlockData) -> SamplesBlockRender:
    """Render one of the 15 dedicated layouts into a fixed three-row body."""
    width = max(1, int(width))
    display = data.sampler
    method = str(display.method or "Unknown")
    progress = display.progress

    if _telemetry_waiting(display):
        rows = _waiting_rows(width, display)
    elif method == "Random":
        rows = _finite_rows(width, display, label="CANDIDATES")
    elif method == "Grid":
        rows = _finite_rows(width, display, label="CELLS")
    elif method == "CSV":
        rows = _finite_rows(width, display, label="ROWS")
    elif method == "Bridson":
        rows = _finite_rows(width, display, label="CLOUD CURSOR")
    elif method == "AdaptiveBridson":
        rows = _finite_rows(width, display, label="GENERATION")
        try:
            target = float(progress.target)
            fraction = max(0.0, min(1.0, float(progress.current) / target))
        except (TypeError, ValueError, ZeroDivisionError):
            target = 0.0
            fraction = 0.0
        if rows is not None and target > 0:
            state = str(display.state or "partial").upper()
            right = f"OPEN {_integer(_metric(display, 'open'), '0')} · {state}"
            rail_width = min(31, max(8, width - len(right) - 1))
            rows = (rows[0], _edge_align_markup(bar_rail(fraction, rail_width), right, width))
    elif method in {"MCMC", "ToyMCMC", "AMMCMC", "DRAM"}:
        rows = _mcmc_rows(width, display, role="CHAINS", unit="ITER")
    elif method in {"EnsembleMCMC", "DEMCMC"}:
        rows = _mcmc_rows(width, display, role="WALKERS", unit="GEN")
    elif method == "PTMCMC":
        rows = _mcmc_rows(width, display, role="REPLICAS", unit="ITER")
    elif method == "PTEnsemble":
        rows = _mcmc_rows(width, display, role="REPLICAS", unit="GEN")
    elif method in {"Dynesty", "MultiNest"}:
        rows = _nested_rows(width, display)
    else:
        state = str(display.state or "running").upper()
        rows = (
            _edge_align(f"COMPLETED {data.completed}", state, width),
            _edge_align("METHOD PROGRESS —", "TELEMETRY N/A", width),
        )

    if rows is None:
        rows = _waiting_rows(width, display)
    body = "\n".join((rows[0], rows[1], _common_row(width, data)))
    return SamplesBlockRender(
        border_title=sample_border_title(width + 4, method),
        body=body,
        border_subtitle=_subtitle(display),
    )


def sampler_display_from_sources(
    metadata: Mapping[str, Any] | None,
    proc_core: Mapping[str, Any] | None,
) -> SamplerDisplay:
    """Normalize one-shot metadata and bounded Core heartbeat telemetry."""
    sampler_meta = dict(metadata or {})
    method = str(sampler_meta.get("method") or "Unknown")
    family = str(sampler_meta.get("family") or "")
    config = sampler_meta.get("config")
    metrics = dict(config) if isinstance(config, Mapping) else {}
    if "dimensions" in sampler_meta:
        metrics.setdefault("dimensions", sampler_meta.get("dimensions"))
    if isinstance(metrics.get("shape"), (list, tuple)):
        metrics["shape"] = "×".join(str(item) for item in metrics["shape"])

    payload: Mapping[str, Any] | None = None
    raw = dict(proc_core or {}).get("sampler_status")
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8", errors="replace")
    if isinstance(raw, str) and len(raw.encode("utf-8")) <= 4096:
        try:
            decoded = json.loads(raw)
        except (TypeError, ValueError):
            decoded = None
        if isinstance(decoded, Mapping) and decoded.get("schema") == 1:
            decoded_method = str(decoded.get("method") or "")
            if method == "Unknown" and decoded_method:
                method = decoded_method
            if decoded_method == method:
                payload = decoded

    if payload is None:
        return SamplerDisplay(
            method=method,
            family=family,
            state="initial",
            metrics=metrics,
        )
    raw_progress = payload.get("progress")
    progress_map = dict(raw_progress) if isinstance(raw_progress, Mapping) else {}
    raw_metrics = payload.get("metrics")
    if isinstance(raw_metrics, Mapping):
        metrics.update(raw_metrics)

    target = progress_map.get("target")
    if target is None:
        if method == "Random":
            target = metrics.get("point_number")
        elif method == "Grid":
            target = metrics.get("total")
        elif method == "AdaptiveBridson":
            target = metrics.get("max_generations")
        elif method in {
            "MCMC",
            "ToyMCMC",
            "AMMCMC",
            "DRAM",
            "EnsembleMCMC",
            "DEMCMC",
            "PTMCMC",
            "PTEnsemble",
        }:
            target = metrics.get("num_iters")
    if method in {
        "MCMC",
        "ToyMCMC",
        "AMMCMC",
        "DRAM",
        "PTMCMC",
        "PTEnsemble",
    }:
        metrics.setdefault("chains", metrics.get("num_chains"))
    elif method in {"EnsembleMCMC", "DEMCMC"}:
        metrics.setdefault("walkers", metrics.get("num_chains"))

    return SamplerDisplay(
        method=method,
        family=family,
        state=str(payload.get("state") or "running"),
        progress=SamplerProgress(
            kind=str(progress_map.get("kind") or "none"),
            label=str(progress_map.get("label") or ""),
            current=progress_map.get("current"),
            target=target,
            eta=str(progress_map.get("eta") or "—"),
        ),
        metrics=metrics,
    )


__all__ = [
    "BUILTIN_SAMPLERS",
    "SamplerDisplay",
    "SamplerProgress",
    "SamplesBlockData",
    "SamplesBlockRender",
    "render_samples_block",
    "sampler_display_from_sources",
    "sample_border_title",
]
