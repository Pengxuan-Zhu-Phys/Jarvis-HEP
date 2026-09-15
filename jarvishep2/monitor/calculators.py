"""PackID activity constellation for the simulated CALCULATORS Overview block."""

from __future__ import annotations

import math
import random
from dataclasses import dataclass

from rich.cells import cell_len

from jarvishep2.monitor.bars import visible_width
from jarvishep2.monitor.styles import paint


@dataclass(frozen=True)
class CalculatorPoolData:
    """One physical calculator pool and its admitted busy PackIDs."""

    name: str
    slots: int
    busy_packs: tuple[int, ...] = ()


@dataclass(frozen=True)
class CalculatorsBlockData:
    """Pool snapshots; busy_packs carry positions, never Worker ownership."""

    pools: tuple[CalculatorPoolData, ...] = ()


@dataclass(frozen=True)
class CalculatorsBlockRender:
    border_title: str
    body: str
    border_subtitle: str


def _clip(text: str, width: int) -> str:
    if width <= 0:
        return ""
    if cell_len(text) <= width:
        return text
    if width == 1:
        return "…"
    chars: list[str] = []
    used = 0
    for char in text:
        size = cell_len(char)
        if used + size > width - 1:
            break
        chars.append(char)
        used += size
    return "".join(chars) + "…"


def _title(outer_width: int) -> str:
    base, right = "CALCULATORS", "ACTIVITY"
    available = max(len(base), int(outer_width) - 6)
    gap = max(1, available - len(base) - len(right) - 2)
    return f"{base} {'─' * gap} {right}"


def _star(pack_id: int, phase: float) -> str:
    """Fast, phase-shifted pulse that never changes the PackID's position."""
    bright = math.sin(phase * 4.0 + pack_id * 1.37) > 0.15
    if bright:
        return "[bold #e6e8eb]✦[/]"
    return "[#134a8d]✦[/]"


def _slot_positions(pool: CalculatorPoolData, width: int, slots: int) -> list[int]:
    """Return stable jittered positions that preserve PackID order."""
    if slots <= 0:
        return []
    if slots == 1:
        return [width // 2]
    step = (width - 1) / (slots - 1)
    rng = random.Random(f"{pool.name}:{slots}:{width}")
    positions = [0]
    for index in range(1, slots - 1):
        ideal = index * step
        jitter = rng.uniform(-0.34 * step, 0.34 * step)
        lower = positions[-1] + 1
        upper = round((index + 1) * step) - 1
        positions.append(max(lower, min(upper, round(ideal + jitter))))
    positions.append(width - 1)
    return positions


def _constellation(pool: CalculatorPoolData, width: int, phase: float) -> str:
    slots = max(0, int(pool.slots))
    busy = {pack for pack in pool.busy_packs if 1 <= int(pack) <= slots}
    if slots <= width:
        cells = [" "] * width
        positions = _slot_positions(pool, width, slots)
        for pack in range(1, slots + 1):
            # One physical PackID owns one stable, irregularly-spaced cell.
            position = positions[pack - 1]
            cells[position] = _star(pack, phase) if pack in busy else "[#134a8d]☆[/]"
        return "".join(cells)
    # A narrow terminal cannot truthfully give every PackID a fixed cell.
    # Use explicit counts rather than inventing a compressed ordering.
    return paint("dim", _clip(f"✦×{len(busy)}  ·×{slots - len(busy)}", width))


def _pool_rows(pool: CalculatorPoolData, width: int, phase: float) -> tuple[str, str]:
    slots = max(0, int(pool.slots))
    busy = len({pack for pack in pool.busy_packs if 1 <= int(pack) <= slots})
    state = "ACTIVE" if busy else "QUIET"
    state_width = len("ACTIVE")
    count_width = max(6, cell_len(f"{busy} / {slots}"))
    state_column = max(1, width - state_width - count_width - 1)
    label = _clip(pool.name, state_column)
    status_color = "#73b8f4" if busy else "#8d93a1"
    status_gap = " " * max(0, state_column - cell_len(label))
    status = (
        paint("panel-title", label)
        + status_gap
        + f"[bold {status_color}]{state:<{state_width}}[/]"
        + " "
        + f"[bold #e6e8eb]{f'{busy} / {slots}':>{count_width}}[/]"
    )
    stars = _constellation(pool, width, phase)
    return status, stars + " " * max(0, width - visible_width(stars))


def render_calculators_block(
    width: int,
    data: CalculatorsBlockData,
    *,
    outer_width: int | None = None,
    pulse_phase: float = 0.0,
) -> CalculatorsBlockRender:
    """Render a star per physical PackID, never a utilisation progress bar."""
    width = max(1, int(width))
    pools = tuple(data.pools)
    busy = sum(
        len({pack for pack in pool.busy_packs if 1 <= int(pack) <= max(0, int(pool.slots))})
        for pool in pools
    )
    slots = sum(max(0, int(pool.slots)) for pool in pools)
    summary = _clip(f"{busy} active  ·  {max(0, slots - busy)} quiet", width)
    rows = [
        paint("panel-title", "  FLEET")
        + " " * max(0, width - cell_len("  FLEET") - cell_len(summary))
        + f"[bold #e6e8eb]{summary}[/]"
    ]
    for pool in pools:
        rows.extend(_pool_rows(pool, width, pulse_phase))
    return CalculatorsBlockRender(
        border_title=_title(outer_width or width + 4),
        body="\n".join(rows),
        border_subtitle="✦ white pulse = occupied · ☆ quiet = free",
    )


__all__ = [
    "CalculatorPoolData",
    "CalculatorsBlockData",
    "CalculatorsBlockRender",
    "render_calculators_block",
]
