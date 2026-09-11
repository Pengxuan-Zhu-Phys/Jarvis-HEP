"""Shared workspace chrome: topbar, frozen folder tabs, frozen hint."""

from __future__ import annotations

from textual.app import ComposeResult
from textual.containers import Horizontal
from textual.widgets import Static

from jarvishep2.monitor.hints import OVERVIEW_KEYS, render_key_hint
from jarvishep2.monitor.styles import paint
from jarvishep2.monitor.topbar import render_topbar_clock, render_topbar_left

PAGES: tuple[tuple[str, str], ...] = (
    ("overview", "Overview"),
    ("workers", "Workers"),
    ("factory", "Factory"),
    ("sampler", "Sampler"),
    ("calculators", "Calculators"),
    ("samples", "Samples"),
    ("host", "Host"),
)

SHORT_LABEL = {
    "Overview": "Over",
    "Workers": "Wrk",
    "Factory": "Fac",
    "Sampler": "Spl",
    "Calculators": "Cal",
    "Samples": "Sam",
    "Host": "Host",
}


IDLE_GAP = 5
MIN_GAP_BEFORE_SHRINK = 2


def _abbrev_levels(full: str) -> list[str]:
    """Name tokens: Calculators → Cal. → Ca. → C. → ''."""
    stem = SHORT_LABEL[full]
    dotted = [stem[:n] + "." for n in range(len(stem), 0, -1)]
    if stem == full:
        dotted = [stem[:n] + "." for n in range(len(stem) - 1, 0, -1)]
    return [full, *dotted, ""]


def _idle_for_level(active_index: int, level: int) -> list[str | None]:
    texts: list[str | None] = []
    for j, (_slug, name) in enumerate(PAGES):
        if j == active_index:
            texts.append(None)
            continue
        forms = _abbrev_levels(name)
        token = forms[min(level, len(forms) - 1)]
        texts.append(f"{j + 1} {token}".rstrip() if token else f"{j + 1}")
    return texts


def _max_name_level() -> int:
    return max(len(_abbrev_levels(name)) - 1 for _slug, name in PAGES)


def folder_tab_lines(
    width: int,
    active: int,
    *,
    compact: bool | None = None,
) -> tuple[str, str, str, list[tuple[str, int, int]]]:
    """Folder tab on a rounded page frame.

    Active inner is `` 1 Overview `` (space after the bar, after the
    number, and after the name). Idle tabs: full name → ``Cal.`` while
    shrinking gaps down to 2 → then ``Ca.`` / ``C.`` / empty. Frozen.
    """
    del compact
    if width <= 0:
        width = 80
    count = len(PAGES)
    index = max(1, min(count, active)) - 1
    full = PAGES[index][1]
    inner = f" {index + 1} {full} "
    box_w = len(inner) + 2
    if box_w > width:
        inner = inner[: max(0, width - 2)].ljust(max(0, width - 2))
        box_w = width

    n_gaps = count - 1
    idle_indices = [j for j in range(count) if j != index]

    def _fits(texts: list[str | None], gap: int) -> bool:
        used = box_w + n_gaps * gap
        used += sum(len(text) for text in texts if text is not None)
        return used <= width

    idle = _idle_for_level(index, 0)
    gap_w = IDLE_GAP
    if not _fits(idle, IDLE_GAP):
        idle = _idle_for_level(index, 1)
        picked = False
        for gap_w in range(IDLE_GAP, MIN_GAP_BEFORE_SHRINK - 1, -1):
            if _fits(idle, gap_w):
                picked = True
                break
        if not picked:
            gap_w = MIN_GAP_BEFORE_SHRINK
            chosen = idle
            for level in range(2, _max_name_level() + 1):
                candidate = _idle_for_level(index, level)
                if _fits(candidate, MIN_GAP_BEFORE_SHRINK):
                    chosen = candidate
                    break
                chosen = candidate
            idle = chosen
            if not _fits(idle, MIN_GAP_BEFORE_SHRINK):
                idle_len = sum(len(text) for text in idle if text is not None)
                leftover = width - box_w - idle_len
                gap_w = max(0, leftover // n_gaps) if n_gaps else 0

    extra = max(
        0,
        width
        - box_w
        - n_gaps * gap_w
        - sum(len(text) for text in idle if text is not None),
    )
    pads = [0] * count
    if idle_indices and extra:
        base, more = divmod(extra, len(idle_indices))
        for offset, j in enumerate(idle_indices):
            pads[j] = base + (1 if offset < more else 0)

    gap = " " * gap_w
    gap_rule = "─" * gap_w
    top_parts: list[str] = []
    mid_parts: list[str] = []
    join_parts: list[str] = []
    hits: list[tuple[str, int, int]] = []
    column = 0
    for j, (slug, _name) in enumerate(PAGES):
        if j == index:
            top_parts.append("╭" + ("─" * len(inner)) + "╮")
            mid_parts.append("│" + inner + "│")
            join_parts.append("╯" + (" " * len(inner)) + "╰")
            hits.append((slug, column, column + box_w))
            column += box_w
        else:
            label = (idle[j] or "") + (" " * pads[j])
            slot = len(label)
            top_parts.append(" " * slot)
            mid_parts.append(label)
            join_parts.append("─" * slot)
            hits.append((slug, column, column + slot))
            column += slot
        if j != count - 1:
            top_parts.append(gap)
            mid_parts.append(gap)
            join_parts.append(gap_rule)
            column += gap_w

    top = _clip("".join(top_parts), width)
    mid = _clip("".join(mid_parts), width)
    join = list(_clip("".join(join_parts), width))
    if index != 0:
        join[0] = "╭"
    else:
        join[0] = "│"
    if index != count - 1:
        join[-1] = "╮"
    else:
        join[-1] = "│"
    return top, mid, "".join(join), hits


def _clip(text: str, width: int) -> str:
    if len(text) > width:
        return text[:width]
    return text + (" " * (width - len(text)))


class MonitorTopbar(Horizontal):
    """Git / path / clock. Lives in every workspace view."""

    def __init__(self) -> None:
        super().__init__(id="topbar")

    def compose(self) -> ComposeResult:
        yield Static(id="topbar-left")
        yield Static(id="topbar-clock")

    def on_mount(self) -> None:
        self.paint()
        self.set_interval(1.0, self.paint)

    def paint(self) -> None:
        self.query_one("#topbar-left", Static).update(
            paint("topbar", render_topbar_left())
        )
        self.query_one("#topbar-clock", Static).update(
            paint("topbar", render_topbar_clock())
        )


class TabBar(Static):
    """Frozen folder-tab strip. See docs/TUI/STYLES.txt (Tabs)."""

    def __init__(self) -> None:
        super().__init__(id="tab-bar")
        self._slug = "overview"
        self._hits: list[tuple[str, int, int]] = []

    def on_mount(self) -> None:
        self.set_active(self._slug)

    def on_resize(self) -> None:
        self.set_active(self._slug)

    def set_active(self, slug: str) -> None:
        self._slug = slug
        width = self.size.width or 80
        active = 1
        for index, (page_slug, _) in enumerate(PAGES, start=1):
            if page_slug == slug:
                active = index
                break
        top, mid, join, hits = folder_tab_lines(width, active)
        self._hits = hits
        gold = "tab-active"
        dim = "tab-idle"
        # Colour the active label gold, idle names dim; box drawing stays.
        _, full = PAGES[active - 1]
        needle = f" {active} {full} "
        if needle in mid:
            mid = mid.replace(
                needle,
                paint(gold, needle),
                1,
            )
        for index, (_, idle_full) in enumerate(PAGES, start=1):
            if index == active:
                continue
            tokens = []
            for form in _abbrev_levels(idle_full):
                tokens.append(f"{index} {form}".rstrip() if form else str(index))
            for token in sorted(set(tokens), key=len, reverse=True):
                if token and token in mid:
                    mid = mid.replace(token, paint(dim, token), 1)
                    break
        mid = mid.replace("│", paint("frame", "│"))
        self.update(
            f"{paint('frame', top)}\n{mid}\n{paint('frame', join)}"
        )

    def on_click(self, event) -> None:  # noqa: ANN001
        x = int(getattr(event, "x", -1))
        switch = getattr(self.screen, "show_page", None)
        if not callable(switch) or x < 0:
            return
        for slug, start, end in self._hits:
            if start <= x < end:
                switch(slug)
                return


class MonitorHint(Static):
    def __init__(self) -> None:
        super().__init__(render_key_hint(OVERVIEW_KEYS), id="hint")

    def set_keys(self, pairs: tuple[tuple[str, str], ...]) -> None:
        self.update(render_key_hint(pairs))


__all__ = [
    "PAGES",
    "MonitorHint",
    "MonitorTopbar",
    "TabBar",
    "folder_tab_lines",
]
