#!/usr/bin/env python3
"""Mint exact-WxH TUI canvases into docs/TUI/.

WARNING: this overwrites every canvas from fake data. After you start drawing
by hand, stop running it. Use _check.py to verify dimensions instead.
"""

from __future__ import annotations

import re
from pathlib import Path

from jarvishep2.monitor.chrome import folder_tab_lines

ROOT = Path(__file__).resolve().parent
LOGO_PATH = ROOT.parent.parent / "jarvishep2" / "card" / "logo"
VERSION = "2.0.10"
ICON_TEMPLATE_RE = re.compile(r"^([BWY.]{8})( {2})(.*)$")

# ---------------------------------------------------------------------------
# Fake scan (stable across every canvas)
# ---------------------------------------------------------------------------
SCAN = "iDM_Vector_V1"
MODE = "running"
ELAPSED = "00:42:18"
REF = "R1"
REDIS = "127.0.0.1:6379"
HZ = "2.0 Hz"
RUN_ID = "run-20260910-120104"
METHOD = "AdaptiveBridson"
HOST = "nersc-login-02"
STARTED = "2026-09-10 12:01:04 UTC"
WORKERS_ALIVE = 47
WORKERS_TOTAL = 48
COMPLETED = 1842
RUNNING_N = 12
FAILED = 7
TASK_Q = 24
ARCH_Q = 6
FB_Q = 0
CPU_PCT = 41.0
MEM_USED_G = 18.4
MEM_TOTAL_G = 64.0
LOAD = (12.41, 10.02, 8.77)

CALCS = [
    ("SoftSUSY", 12, 4),
    ("micrOMEGAs", 8, 0),
    ("HiggsBounds", 10, 2),
    ("HiggsSignals", 11, 1),
    ("SModelS", 16, 0),
    ("HiggsTools", 8, 2),
    ("FeynHiggs", 4, 0),
    ("SPheno", 6, 1),
]

SPARK_SAMPLES = [
    3, 4, 5, 6, 8, 10, 9, 11, 12, 14, 13, 12, 15, 16, 14, 13,
    15, 17, 18, 16, 15, 14, 16, 18, 19, 17, 16, 18, 20, 19, 18, 17,
    16, 18, 19, 21, 20, 19, 18, 17, 19, 20, 18, 17, 16, 18, 19, 20,
]
SPARK_QUEUE = [
    8, 9, 12, 14, 18, 22, 20, 16, 14, 11, 9, 10, 13, 17, 24, 22,
    19, 15, 12, 10, 8, 9, 12, 16, 20, 24, 21, 18, 14, 12, 10, 11,
    13, 16, 19, 22, 24, 20, 16, 13, 11, 9, 10, 14, 18, 21, 24, 22,
]
SPARK_CPU = [
    22, 28, 35, 40, 44, 38, 33, 41, 48, 52, 47, 43, 39, 36, 42, 49,
    55, 51, 46, 41, 37, 34, 39, 45, 50, 47, 41, 38, 36, 40, 44, 41,
    39, 43, 48, 46, 42, 40, 38, 41, 45, 49, 44, 40, 37, 41, 43, 41,
]

PAGES = [
    (1, "Overview", "Over"),
    (2, "Workers", "Wrk"),
    (3, "Factory", "Fac"),
    (4, "Sampler", "Spl"),
    (5, "Calculators", "Cal"),
    (6, "Samples", "Sam"),
    (7, "Host", "Host"),
]

SPARK_CHARS = "▁▂▃▄▅▆▇█"


def fit(text: str, n: int, align: str = "left") -> str:
    if n <= 0:
        return ""
    if len(text) > n:
        return text[:n]
    pad = n - len(text)
    if align == "right":
        return (" " * pad) + text
    if align == "center":
        a, b = pad // 2, pad - pad // 2
        return (" " * a) + text + (" " * b)
    return text + (" " * pad)


def bar(frac: float, width: int, fill: str = "█", empty: str = "░") -> str:
    if width <= 0:
        return ""
    frac = max(0.0, min(1.0, frac))
    k = int(round(frac * width))
    k = min(width, max(0, k))
    return fill * k + empty * (width - k)


def spark(values: list[int], width: int) -> str:
    if width <= 0:
        return ""
    seq = (values * ((width // len(values)) + 1))[:width]
    lo, hi = min(seq), max(seq)
    span = max(1, hi - lo)
    out = []
    for v in seq:
        idx = int(round((v - lo) / span * (len(SPARK_CHARS) - 1)))
        out.append(SPARK_CHARS[idx])
    return "".join(out)


def hline(width: int, left: str, mid: str, right: str, fill: str = "─") -> str:
    return left + (fill * (width - 2)) + right


def status_inner(inner: int, compact: bool) -> str:
    if compact:
        left, right = f" Jarvis  {SCAN}  ● run  {ELAPSED}  {REF}", HZ.replace(" ", "")
    else:
        left = f" Jarvis Monitor  {SCAN}  ● {MODE}  {ELAPSED}  {REF}"
        right = f"{REDIS}  {HZ}"
    space = inner - len(left) - len(right)
    if space < 3:
        left = f" Jarvis  {SCAN}  ● {MODE}  {REF}"
        right = HZ
        space = inner - len(left) - len(right)
    if space < 3:
        left = f" {SCAN}  ● {MODE}"
        right = HZ
        space = inner - len(left) - len(right)
    if space < 1:
        return fit(f"{SCAN} {MODE} {HZ}", inner)
    dash = max(1, space - 2)
    return left + " " + ("─" * dash) + " " + right


def tabs_inner(inner: int, current: int, compact: bool) -> str:
    parts: list[str] = []
    for num, full, short in PAGES:
        label = short if compact else full
        token = f"{num} {label}"
        parts.append(f"[{token}]" if num == current else f" {token} ")
    return fit(" " + " ".join(parts), inner)


def footer_inner(inner: int, compact: bool) -> str:
    if compact:
        text = "1-7 pages   tab   j/k   enter pin   esc   [ ] hz   q   ?"
    else:
        text = (
            "1-7 pages   tab next   j/k rows   enter pin   esc clear"
            "   [ ] hz   r refresh   q quit   ?"
        )
    return fit(text, inner)


def _notch_folder_tab(line: str, left: int, right: int) -> str:
    """Open the active tab into the first body rule with rounded ╯ ╰."""
    chars = list(line)
    inner_left = 1 + left
    inner_right = 1 + right
    if left > 0 and 0 <= inner_left < len(chars):
        chars[inner_left] = "╯"
        start = inner_left + 1
    else:
        start = inner_left + 1
    for index in range(start, min(inner_right, len(chars) - 1)):
        chars[index] = " "
    if 0 <= inner_right < len(chars):
        chars[inner_right] = "╰"
    return "".join(chars)


def chrome(width: int, height: int, current: int, body: list[str], compact: bool) -> list[str]:
    """`body` is already `height-6` lines, each exactly `width` chars (includes side bars)."""
    inner = width - 2
    tab_top, tab_mid, _, hits = folder_tab_lines(inner, current, compact=compact)
    lines = [
        "┌" + status_inner(inner, compact) + "┐",
        "│" + tab_top + "│",
        "│" + tab_mid + "│",
    ]
    if len(body) != height - 6:
        raise RuntimeError(f"body rows {len(body)} != {height - 6}")
    body = list(body)
    if body and hits:
        from jarvishep2.monitor.chrome import PAGES as TAB_PAGES

        slug = TAB_PAGES[current - 1][0]
        left = right_ex = 0
        for hit_slug, hit_left, hit_right in hits:
            if hit_slug == slug:
                left, right_ex = hit_left, hit_right
                break
        body[0] = _notch_folder_tab(body[0], left, right_ex - 1)
    lines.extend(body)
    lines.append("├" + ("─" * inner) + "┤")
    lines.append("│" + footer_inner(inner, compact) + "│")
    lines.append("└" + ("─" * inner) + "┘")
    for i, line in enumerate(lines):
        if len(line) != width:
            raise RuntimeError(f"chrome line {i} width {len(line)} != {width}: {line!r}")
    if len(lines) != height:
        raise RuntimeError(f"page height {len(lines)} != {height}")
    return lines


def inner_row(width: int, text: str, align: str = "left") -> str:
    return "│" + fit(text, width - 2, align) + "│"


def sep(width: int, kind: str = "mid") -> str:
    inner = width - 2
    if kind == "top":
        return "├" + ("─" * inner) + "┤"
    if kind == "bot":
        return "├" + ("─" * inner) + "┤"
    return "├" + ("─" * inner) + "┤"


def split_sep(width: int, left_w: int, joints: str = "├┬┤") -> str:
    right_w = width - 2 - left_w - 1
    a, b, c = joints
    return a + ("─" * left_w) + b + ("─" * right_w) + c


def fake_workers() -> list[dict[str, str]]:
    rows = []
    uuids = [
        "a3f2c91e-7b14-4d2a-9c88-e0b1d44a1101",
        "91bc4e12-0aa1-4f77-b3de-2c91f08c4e12",
        "c019d33a-55e0-41b2-a901-77ab12cc90de",
        "e77a01b4-9c22-4a10-8f01-bb4455667788",
        "11f0aa92-cd31-4e90-a712-90ff00aa11bb",
        "b2b2b2b2-0001-4000-8000-000000000007",
        "d4e5f607-89ab-4cde-a012-3456789abcde",
        "00aa11bb-22cc-43dd-8eef-ff0011223344",
    ]
    pids = [44100 + i for i in range(WORKERS_TOTAL)]
    for i in range(WORKERS_TOTAL):
        running = i % 5 != 4
        stale = i == 17
        idle = not running and not stale
        status = "stale" if stale else ("run" if running else "idle")
        uuid = uuids[i % len(uuids)] if running else "—"
        if running:
            uuid = uuid[:-2] + f"{i:02d}"
        hb = 9.4 if stale else (0.12 + (i % 7) * 0.07)
        cpu = 3.0 if idle else (12.0 if stale else 55.0 + (i * 3) % 40)
        rss = 88 if idle else (140 if stale else 360 + (i * 13) % 220)
        held = 0 if not running else (1 + i % 3)
        rows.append(
            {
                "id": f"{i:02d}",
                "pid": str(pids[i]),
                "status": status,
                "uuid": uuid,
                "uuid_short": uuid[:8] if uuid != "—" else "—",
                "hb": f"{hb:0.1f}s",
                "cpu": f"{cpu:0.0f}%",
                "rss": f"{rss}M",
                "calc": str(held),
                "fo_pid": str(pids[i] + 80) if running else "—",
                "children": f"{pids[i] + 80},{pids[i] + 81}" if running else "—",
            }
        )
    return rows


WORKERS = fake_workers()
SELECTED = 7


def kv_block(pairs: list[tuple[str, str]], width: int, height: int, title: str) -> list[str]:
    lines = [fit(title, width), fit("", width)]
    key_w = min(14, max(10, width // 3))
    for key, value in pairs:
        if len(lines) >= height:
            break
        room = max(1, width - key_w - 1)
        text = str(value)
        lines.append(fit(key, key_w) + " " + fit(text[:room], room))
        rest = text[room:]
        while rest and len(lines) < height:
            lines.append(fit("", key_w) + " " + fit(rest[:room], room))
            rest = rest[room:]
    while len(lines) < height:
        lines.append(fit("", width))
    return lines[:height]


def table_row(cols: list[tuple[str, int]]) -> str:
    return " ".join(fit(text, w) for text, w in cols)


# ---------------------------------------------------------------------------
# Overview
# ---------------------------------------------------------------------------
def overview_body(width: int, height: int, compact: bool) -> list[str]:
    """height is body rows (page-4). Each line length == width."""
    inner = width - 2
    four_col = width >= 150
    if four_col:
        return overview_four_col(width, height, compact)
    left_w = inner // 2
    right_w = inner - left_w - 1
    seps = 4
    live_h = 3 if height >= 22 else 2
    spark_h = 5 if height >= 28 else (4 if height >= 22 else 3)
    remaining = height - live_h - spark_h - seps
    if remaining < 8:
        spark_h = max(3, spark_h - (8 - remaining))
        remaining = height - live_h - spark_h - seps
    top_h = max(4, remaining // 2)
    bot_h = max(4, remaining - top_h)
    extra = height - (seps + top_h + bot_h + spark_h + live_h)
    if extra > 0:
        spark_h += extra

    samples = panel_samples(left_w, top_h, compact)
    queues = panel_queues(right_w, top_h, compact)
    workers = panel_workers_sum(left_w, bot_h, compact)
    calcs = panel_calcs(right_w, bot_h, compact)
    sparks = panel_sparks(inner, spark_h, compact)
    live = panel_live(inner, live_h, compact)

    lines: list[str] = [split_sep(width, left_w, "├┬┤")]
    for a, b in zip(samples, queues):
        lines.append("│" + a + "│" + b + "│")
    lines.append(split_sep(width, left_w, "├┼┤"))
    for a, b in zip(workers, calcs):
        lines.append("│" + a + "│" + b + "│")
    lines.append(split_sep(width, left_w, "├┴┤"))
    for row in sparks:
        lines.append("│" + row + "│")
    lines.append(sep(width))
    for row in live:
        lines.append("│" + row + "│")
    # pad / trim to exact body height
    while len(lines) < height:
        lines.append(inner_row(width, ""))
    return lines[:height]


def overview_four_col(width: int, height: int, compact: bool) -> list[str]:
    inner = width - 2
    # 4 columns + 3 joints
    usable = inner - 3
    w1 = w2 = w3 = usable // 4
    w4 = usable - 3 * w1
    widths = [w1, w2, w3, w4]
    if height <= 44:
        top_h, spark_h, extra_h, live_h = 12, 8, 12, 4
    else:
        top_h, spark_h, extra_h, live_h = 14, 10, 16, 4
    used = 1 + top_h + 1 + spark_h + 1 + extra_h + 1 + live_h
    extra = height - used
    if extra > 0:
        extra_h += extra
    elif extra < 0:
        extra_h = max(6, extra_h + extra)

    p1 = panel_samples(w1, top_h, compact)
    p2 = panel_queues(w2, top_h, compact)
    p3 = panel_workers_sum(w3, top_h, compact)
    p4 = panel_calcs(w4, top_h, compact)
    sparks = panel_sparks(inner, spark_h, compact)
    extra_panel = panel_running_strip(inner, extra_h, compact)
    live = panel_live(inner, live_h, compact)

    def grid_sep(joints: str) -> str:
        a, m, z = joints[0], joints[1], joints[2]
        return a + m.join("─" * w for w in widths) + z

    lines = [grid_sep("├┬┤")]
    for i in range(top_h):
        lines.append("│" + "│".join([p1[i], p2[i], p3[i], p4[i]]) + "│")
    lines.append(grid_sep("├┴┤")[:width].ljust(width) if False else "├" + "┴".join("─" * w for w in widths) + "┤")
    for row in sparks:
        lines.append("│" + row + "│")
    lines.append(sep(width))
    for row in extra_panel:
        lines.append("│" + row + "│")
    lines.append(sep(width))
    for row in live:
        lines.append("│" + row + "│")
    while len(lines) < height:
        lines.append(inner_row(width, ""))
    return lines[:height]


def titled(width: int, title: str, rows: list[str], height: int) -> list[str]:
    lines = [fit(title, width)]
    lines.extend(fit(r, width) for r in rows)
    while len(lines) < height:
        lines.append(fit("", width))
    return lines[:height]


def panel_samples(w: int, h: int, compact: bool) -> list[str]:
    target = 3000
    done_frac = COMPLETED / target
    bw = max(6, w - 20)
    rows = [
        f"target {target:>6}",
        f"done   {COMPLETED:>6}  " + bar(done_frac, bw) + f"  {done_frac * 100:0.0f}%",
        f"run    {RUNNING_N:>6}",
        f"fail   {FAILED:>6}",
        f"remain {target - COMPLETED - RUNNING_N:>6}",
        "rate    12.4 /min",
        "eta     1:31:40",
        "avg     3.82 s",
    ]
    return titled(w, "SAMPLES", rows, h)


def panel_queues(w: int, h: int, compact: bool) -> list[str]:
    cap = 32
    rows = [
        f"task      {TASK_Q:>4}  " + bar(TASK_Q / cap, max(6, w - 18)),
        f"archive   {ARCH_Q:>4}  " + bar(ARCH_Q / cap, max(6, w - 18)),
        f"feedback  {FB_Q:>4}  " + bar(0.02, max(6, w - 18)),
        "",
        "chain-0      0",
        "chain-1      0",
        "chain-2      0",
    ]
    return titled(w, "QUEUES", rows, h)


def panel_workers_sum(w: int, h: int, compact: bool) -> list[str]:
    rows = [
        f"alive     {WORKERS_ALIVE}/{WORKERS_TOTAL}",
        f"stale     1",
        f"busy      {RUNNING_N}",
        f"idle      {WORKERS_ALIVE - RUNNING_N}",
        "",
        "cpu  " + bar(CPU_PCT / 100.0, max(6, w - 16)) + f"  {CPU_PCT:0.0f}%",
        "mem  " + bar(MEM_USED_G / MEM_TOTAL_G, max(6, w - 16)) + f"  {MEM_USED_G:0.1f}G",
    ]
    return titled(w, "WORKERS", rows, h)


def panel_calcs(w: int, h: int, compact: bool) -> list[str]:
    rows = []
    for name, free, busy in CALCS:
        total = free + busy
        rows.append(
            fit(name, 12) + f" {busy:>2}/{total:<2}  " + bar(busy / max(1, total), max(4, w - 22))
        )
    return titled(w, "CALCULATORS", rows, h)


def panel_sparks(w: int, h: int, compact: bool) -> list[str]:
    if compact or w < 40:
        labels = [("smp", 3), ("que", 3), ("cpu", 3)]
    else:
        labels = [("samples/min", 12), ("task queue", 12), ("host cpu", 12)]
    series = [SPARK_SAMPLES, SPARK_QUEUE, SPARK_CPU]
    rows: list[str] = []
    if h >= 4:
        rows.append("SPARKS  (local TUI ring, not Redis)")
    for (label, label_w), data in zip(labels, series):
        plot_w = max(8, w - label_w - 1)
        rows.append(fit(label, label_w) + " " + spark(data, plot_w))
    if h > 6:
        rows.append("")
        rows.append("each tick 2.0 Hz   history ~120 pts   never written back")
    out = [fit(r, w) for r in rows]
    while len(out) < h:
        out.append(fit("", w))
    return out[:h]


def panel_live(w: int, h: int, compact: bool) -> list[str]:
    dots = "core ●   archiver ●   redis ●   lock ●"
    rest = f"started {STARTED}   method {METHOD}   run {RUN_ID}"
    if compact:
        rest = f"{METHOD}   {RUN_ID}"
    rows = [dots, rest, f"host {HOST}   workers {WORKERS_ALIVE}/{WORKERS_TOTAL}   stale 1"]
    return [fit(r, w) for r in rows][:h] + [fit("", w)] * max(0, h - min(h, len(rows)))


def panel_running_strip(w: int, h: int, compact: bool) -> list[str]:
    header = fit("RUNNING SAMPLES", w)
    cols = "  wrk   pid      uuid                                  hb     cpu   calc"
    lines = [header, fit(cols, w)]
    shown = 0
    for row in WORKERS:
        if row["status"] != "run":
            continue
        line = (
            f"  {row['id']}   {row['pid']:<8} {fit(row['uuid'], 36)}  "
            f"{row['hb']:>5}  {row['cpu']:>4}    {row['calc']}"
        )
        lines.append(fit(line, w))
        shown += 1
        if shown >= h - 3:
            break
    while len(lines) < h:
        lines.append(fit("", w))
    return lines[:h]


# ---------------------------------------------------------------------------
# Master-detail pages
# ---------------------------------------------------------------------------
def master_detail(
    width: int,
    height: int,
    left_lines: list[str],
    right_lines: list[str],
    left_frac: float = 0.64,
) -> list[str]:
    inner = width - 2
    left_w = int(inner * left_frac)
    if left_w < 28:
        left_w = max(22, inner // 2)
    if inner - left_w - 1 < 18:
        left_w = inner - 1 - 18
    right_w = inner - left_w - 1
    body_h = height - 2  # split top sep + we already have tabs; last body is not footer
    # first line split cap, then rows, no bottom split (footer chrome handles)
    # height here is body height. Use 1 sep + (height-1) rows.
    rows_h = height - 1
    L = [fit(x, left_w) for x in left_lines]
    R = [fit(x, right_w) for x in right_lines]
    while len(L) < rows_h:
        L.append(fit("", left_w))
    while len(R) < rows_h:
        R.append(fit("", right_w))
    lines = [split_sep(width, left_w, "├┬┤")]
    for i in range(rows_h):
        lines.append("│" + L[i] + "│" + R[i] + "│")
    return lines[:height]


def workers_page_body(width: int, height: int, compact: bool) -> list[str]:
    inner = width - 2
    left_frac = 0.62 if width >= 120 else 0.58
    left_w = int(inner * left_frac)
    if inner - left_w - 1 < 20:
        left_w = inner - 21
    # table
    if compact or width < 100:
        spec = [("SEL", 3), ("PID", 6), ("STAT", 5), ("UUID", 8), ("HB", 5), ("CPU", 4)]
        header_spec = [(" ID", 3), ("PID", 6), ("STAT", 5), ("UUID", 8), ("HB", 5), ("CPU", 4)]
    elif width < 150:
        spec = [
            ("SEL", 3),
            ("PID", 6),
            ("STAT", 5),
            ("UUID", 10),
            ("HB", 5),
            ("CPU", 4),
            ("RSS", 5),
            ("CALC", 4),
        ]
        header_spec = [
            (" ID", 3),
            ("PID", 6),
            ("STAT", 5),
            ("UUID", 10),
            ("HB", 5),
            ("CPU", 4),
            ("RSS", 5),
            ("CALC", 4),
        ]
    else:
        spec = [
            ("SEL", 3),
            ("PID", 7),
            ("STAT", 5),
            ("UUID", 12),
            ("HB", 6),
            ("CPU", 5),
            ("RSS", 6),
            ("CALC", 4),
            ("FO", 6),
        ]
        header_spec = [
            (" ID", 3),
            ("PID", 7),
            ("STAT", 5),
            ("UUID", 12),
            ("HB", 6),
            ("CPU", 5),
            ("RSS", 6),
            ("CALC", 4),
            ("FO", 6),
        ]
    header = table_row([(h, w) for h, w in header_spec])
    left = [f"WORKERS  {WORKERS_ALIVE}/{WORKERS_TOTAL} alive   1 stale", header]
    rows_h = height - 1 - 2  # sep + title + header, remaining in left after master_detail first sep
    # master_detail uses height-1 data rows; first two are title+header
    max_rows = max(3, height - 3)
    start = 0
    # keep selected visible
    if SELECTED >= max_rows - 2:
        start = SELECTED - (max_rows - 3)
    view = WORKERS[start : start + max_rows]
    for row in view:
        idx = int(row["id"])
        cols = []
        mapping = {
            "SEL": (">" if idx == SELECTED else " ") + row["id"],
            "PID": row["pid"],
            "STAT": row["status"],
            "UUID": row["uuid_short"],
            "HB": row["hb"],
            "CPU": row["cpu"],
            "RSS": row["rss"],
            "CALC": row["calc"],
            "FO": row["fo_pid"],
        }
        for name, w in spec:
            cols.append((mapping[name], w))
        left.append(table_row(cols))
    sel = WORKERS[SELECTED]
    right = kv_block(
        [
            ("id", sel["id"]),
            ("pid", sel["pid"]),
            ("status", sel["status"]),
            ("uuid", sel["uuid"]),
            ("heartbeat", sel["hb"]),
            ("cpu", sel["cpu"]),
            ("rss", sel["rss"]),
            ("held_calc", sel["calc"]),
            ("fo_pid", sel["fo_pid"]),
            ("children", sel["children"]),
            ("host", HOST),
            ("board", f"hep:proc:worker:{sel['id']}"),
        ],
        80,
        40,
        f"worker-{sel['id']}",
    )
    return master_detail(width, height, left, right, left_frac=left_frac)


def factory_page_body(width: int, height: int, compact: bool) -> list[str]:
    left = [
        "FACTORY  control plane",
        "",
        "ROLE       PID      STAT      HOST",
        ">core      44001    running   " + HOST,
        " archiver  44012    running   " + HOST,
        " redis     44000    running   " + HOST,
        " lock      —        held      hep:control:lock",
        "",
        "scan_mode     running",
        "workers_total 48",
        "archiver_alive 1",
        "redis_alive    1",
        "lease_ttl      96 s",
    ]
    right = kv_block(
        [
            ("role", "core"),
            ("pid", "44001"),
            ("host", HOST),
            ("scan_name", SCAN),
            ("run_id", RUN_ID),
            ("started_at", STARTED),
            ("scan_mode", MODE),
            ("workers_total", str(WORKERS_TOTAL)),
            ("lease_owner", f"core:{44001}"),
            ("lease_ts", "12:43:21"),
            ("archiver_pid", "44012"),
            ("archiver_alive", "1"),
            ("redis_alive", "1"),
            ("board", "hep:proc:core"),
        ],
        80,
        40,
        "hep:proc:core",
    )
    return master_detail(width, height, left, right, left_frac=0.58)


def sampler_page_body(width: int, height: int, compact: bool) -> list[str]:
    left = [
        f"SAMPLER  {METHOD}",
        "",
        "QUEUE              LEN",
        ">hep:task_queue      24",
        " hep:archive_queue    6",
        " hep:feedback         0",
        " hep:feedback:chain:0 0",
        " hep:feedback:chain:1 0",
        " hep:feedback:chain:2 0",
        "",
        "bounds.point_number   3000",
        "bounds.seed           7",
        "archived_prefix       1842",
        "checkpoint            checkpoints/iDM_Vector_V1/AdaptiveBridson/state.pkl",
        "checkpoint mtime      12:42:01",
    ]
    right = kv_block(
        [
            ("key", "hep:task_queue"),
            ("llen", "24"),
            ("role", "pending samples"),
            ("payloads", "not listed (v1)"),
            ("method", METHOD),
            ("run_id", RUN_ID),
            ("feedback", "0"),
            ("chains", "3"),
        ],
        80,
        40,
        "hep:task_queue",
    )
    return master_detail(width, height, left, right, left_frac=0.62)


def calculators_page_body(width: int, height: int, compact: bool) -> list[str]:
    left = ["CALCULATORS  pack pools", "", "NAME           BUSY  FREE  TOTAL  UTIL"]
    for i, (name, free, busy) in enumerate(CALCS):
        total = free + busy
        util = f"{busy / total * 100:0.0f}%"
        mark = ">" if i == 0 else " "
        left.append(
            f"{mark}{fit(name, 14)} {busy:>4}  {free:>4}  {total:>5}  {util:>4}  "
            + bar(busy / total, 12)
        )
    right = [
        "SoftSUSY",
        "",
        fit("pack   owner", 40),
        fit("001    worker-03", 40),
        fit("002    worker-07", 40),
        fit("003    worker-12", 40),
        fit("004    worker-22", 40),
        fit("005-16 free", 40),
        "",
        "hash   hep:calculator:status",
        "busy   calc:busy:SoftSUSY",
        "free   calc:free:SoftSUSY",
        "mode   exclusive",
    ]
    return master_detail(width, height, left, right, left_frac=0.62)


def samples_page_body(width: int, height: int, compact: bool) -> list[str]:
    left = [
        "SAMPLES  counters + currently running",
        "Redis has no cheap UUID catalogue. This is not a browser.",
        "",
        f"completed  {COMPLETED}    running  {RUNNING_N}    failed  {FAILED}",
        "bucket     4 open    11 sealed    last pack SAMPLE/00011.tar",
        "archived_prefix  1842",
        "",
        "WRK  UUID                                  STAT   PID",
    ]
    n = 0
    sel_uuid = None
    for row in WORKERS:
        if row["status"] != "run":
            continue
        mark = ">" if n == 1 else " "
        if n == 1:
            sel_uuid = row
        left.append(
            f"{mark}{row['id']}  {fit(row['uuid'], 36)}  run   {row['pid']}"
        )
        n += 1
        if n >= 16:
            break
    sel = sel_uuid or WORKERS[SELECTED]
    right = kv_block(
        [
            ("uuid", sel["uuid"]),
            ("worker", sel["id"]),
            ("pid", sel["pid"]),
            ("status", "running"),
            ("hb", sel["hb"]),
            ("held_calc", sel["calc"]),
            ("result", "not fetched (v1)"),
            ("logs", "not tailed (v1)"),
        ],
        80,
        40,
        "running sample",
    )
    return master_detail(width, height, left, right, left_frac=0.66)


def host_page_body(width: int, height: int, compact: bool) -> list[str]:
    cores = [62, 48, 71, 33, 12, 9, 44, 51, 80, 77, 15, 18, 40, 42, 8, 11]
    core_lines = []
    per = 4 if width >= 120 else 2
    chunk = []
    for i, c in enumerate(cores):
        chunk.append(f"{i:02d} " + bar(c / 100.0, 8) + f" {c:2d}%")
        if len(chunk) == per:
            core_lines.append("  ".join(chunk))
            chunk = []
    if chunk:
        core_lines.append("  ".join(chunk))
    left = [
        f"HOST  {HOST}   load {LOAD[0]:0.2f} {LOAD[1]:0.2f} {LOAD[2]:0.2f}",
        "cpu  " + bar(CPU_PCT / 100.0, 28) + f"  {CPU_PCT:0.0f}%",
        "mem  " + bar(MEM_USED_G / MEM_TOTAL_G, 28) + f"  {MEM_USED_G:0.1f}/{MEM_TOTAL_G:0.0f} G",
        "swp  " + bar(0.02, 28) + "  0.3/8 G",
        "",
        "REDIS INFO",
        "used_memory              42.1M",
        "connected_clients        8",
        "instantaneous_ops/sec    310",
        "",
        "ROLE       PID     CPU   RSS    FDS  CMDLINE",
        ">core      44001   12%   210M    32  Jarvis run bin/iDM.yaml",
        " worker-00 44100   88%   412M    41  Jarvis-Worker-00",
        " worker-07 44121   91%   508M    44  Jarvis-Worker-07",
        " archiver  44012    7%   180M    28  Jarvis-Archiver",
        " redis     44000    2%    48M    12  redis-server *:6379",
        " fo-07     44210   40%    66M    18  file-operation",
    ]
    left[5:5] = core_lines[:4]
    right = kv_block(
        [
            ("role", "core"),
            ("pid", "44001"),
            ("ppid", "1"),
            ("cpu", "12.4%"),
            ("rss", "210 MiB"),
            ("threads", "18"),
            ("fds", "32"),
            ("created", STARTED),
            ("cmdline", "Jarvis run bin/iDM_Vector_V1.yaml"),
            ("scope", "this scan only"),
        ],
        80,
        40,
        "process 44001",
    )
    return master_detail(width, height, left, right, left_frac=0.68)


def load_hep_logo() -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Same parse as jarvishep2.versioning / Jarvis-Agent branding."""
    pattern: list[str] = []
    banners: list[str] = []
    text = LOGO_PATH.read_text(encoding="utf-8")
    for line in text.splitlines():
        matched = ICON_TEMPLATE_RE.match(line)
        if not matched:
            continue
        pattern.append(matched.group(1))
        rest = matched.group(3)
        if "Version:" in rest:
            rest = f"{rest.split('Version:', 1)[0]}Version:  {VERSION}"
        if "Jarvis-HEP" in rest and "V2" not in rest:
            rest = rest.replace("Jarvis-HEP", "Jarvis-HEP V2", 1)
        banners.append(rest)
    return tuple(pattern[:8]), tuple(banners[:8])


def logo_monitor_row(pattern_row: str) -> str:
    """Agent #logo-monitor: eight ⬤ cells, 16 columns (⬤ + space)."""
    cells = pattern_row[:8].ljust(8, ".")
    return "".join("⬤ " for _ in cells)


def _splash_topbar(inner: int) -> str:
    left = " ⎇ main · ~/Jarvis-Workshop/Jarvis-Examples/Eggbox "
    right = "· 8:53 PM "
    fill = inner - len(left) - len(right)
    if fill < 1:
        left = " ⎇ main · "
        fill = max(1, inner - len(left) - len(right))
    return fit(left + (" " * fill) + right, inner)


def splash_page(width: int, height: int, compact: bool) -> list[str]:
    """Live SplashScreen: round hero (logo|banner+tag), full-width ps table, hint."""
    inner = width - 2
    pattern, banners = load_hep_logo()
    icon_w = 16
    gap = 2
    hint = "Enter: attach  |  J/K: select  |  R: refresh  |  Q: quit"
    side = 1 if compact else 2
    hero_w = inner - 2 * side
    if hero_w < 40:
        side = 1
        hero_w = inner - 2
    hero_inner = hero_w - 2
    inset = 1
    content_w = max(20, hero_inner - 2 * inset)
    panel_w = max(0, content_w - icon_w - gap)
    gutter = icon_w + gap

    hero_body: list[str] = []
    for i in range(8):
        icon = logo_monitor_row(pattern[i] if i < len(pattern) else "........")
        banner = banners[i] if i < len(banners) else ""
        hero_body.append(fit(icon + (" " * gap) + fit(banner, panel_w), content_w))
    hero_body.append(fit("", content_w))
    hero_body.append(fit((" " * gutter) + "Jarvis Monitor", content_w))
    hero_body.append(fit((" " * gutter) + "read-only runtime dashboard", content_w))

    def wrap_hero() -> list[str]:
        out = [
            fit((" " * side) + "╭" + ("─" * (hero_w - 2)) + "╮", inner),
            fit((" " * side) + "│" + fit("", hero_inner) + "│", inner),
        ]
        for row in hero_body:
            out.append(
                fit(
                    (" " * side) + "│" + fit((" " * inset) + row, hero_inner) + "│",
                    inner,
                )
            )
        if not compact:
            out.append(fit((" " * side) + "│" + fit("", hero_inner) + "│", inner))
        out.append(fit((" " * side) + "╰" + ("─" * (hero_w - 2)) + "╯", inner))
        return out

    def chooser_line(ref: str, scan: str, control: str, procs: str) -> str:
        pad_left = 1
        pad_right = 1
        usable = inner - pad_left - pad_right
        right = f"{control:>10}  {procs:>5}"
        left = f"{ref:<3}  "
        scan_w = max(1, usable - len(left) - 2 - len(right))
        return fit(
            (" " * pad_left) + left + fit(scan, scan_w) + "  " + right + (" " * pad_right),
            inner,
        )

    rows = [
        chooser_line(">R1", "iDM_Vector_V1", "44001", "52"),
        chooser_line(" R2", "eggbox-bridson", "51020", "6"),
    ]
    chooser: list[str] = [
        fit(" running scans", inner),
        chooser_line("REF", "SCAN", "CONTROL", "PROCS"),
        *rows,
    ]

    topbar = _splash_topbar(inner)
    top_pad = 0
    # outer top + top_pad + hero + gap + chooser + notice + hint + outer bot
    hero = wrap_hero()
    reserved = 2 + 1 + top_pad + len(hero) + 1 + 2 + 1 + 1
    extra = height - reserved - len(rows)
    if extra > 0:
        chooser.extend(fit("", inner) for _ in range(extra))
    elif extra < 0:
        # Drop blank gap before chooser first, then hero inner pads.
        reserved_tight = 2 + len(hero) + 2 + 1 + 1
        extra = height - reserved_tight - len(rows)
        top_pad = 0
        if extra > 0:
            chooser.extend(fit("", inner) for _ in range(extra))

    lines = [hline(width, "┌", "─", "┐")]
    lines.append("│" + topbar + "│")
    for _ in range(top_pad):
        lines.append("│" + fit("", inner) + "│")
    for row in hero:
        lines.append("│" + row + "│")
    lines.append("│" + fit("", inner) + "│")
    for row in chooser:
        lines.append("│" + row + "│")
    lines.append("│" + fit("", inner) + "│")
    lines.append("│" + fit(" " + hint, inner) + "│")
    lines.append(hline(width, "└", "─", "┘"))
    # Trim / pad to exact height without eating the hint.
    if len(lines) > height:
        overflow = len(lines) - height
        # Remove blank rows just above the hint (notice spacer).
        hint_idx = len(lines) - 2
        i = hint_idx - 1
        while overflow > 0 and i > 0:
            if lines[i] == "│" + fit("", inner) + "│":
                del lines[i]
                overflow -= 1
                hint_idx -= 1
            i -= 1
        lines = lines[:height]
        lines[-1] = hline(width, "└", "─", "┘")
    while len(lines) < height:
        lines.insert(-2, "│" + fit("", inner) + "│")
    return [fit(line, width) for line in lines[:height]]


def picker_page(width: int, height: int, compact: bool) -> list[str]:
    inner = width - 2
    content = [
        "Jarvis Monitor  pick a live scan",
        "",
        " REF   SCAN                 CONTROL    PROCS    REDIS",
        ">R1    iDM_Vector_V1        44001         52    127.0.0.1:6379",
        " R2    eggbox-bridson       51020          6    127.0.0.1:6379",
        "",
        "Enter: attach  |  J/K: select  |  Q: quit",
    ]
    lines = [hline(width, "┌", "─", "┐")]
    pad_top = max(1, (height - 2 - len(content)) // 3)
    for _ in range(pad_top):
        lines.append("│" + fit("", inner) + "│")
    for row in content:
        lines.append("│" + fit(row, inner) + "│")
    while len(lines) < height - 1:
        lines.append("│" + fit("", inner) + "│")
    lines.append(hline(width, "└", "─", "┘"))
    return [fit(l, width) for l in lines[:height]]


def empty_page(width: int, height: int, compact: bool) -> list[str]:
    inner = width - 2
    msg = [
        "No running Jarvis scan.",
        "",
        "Redis is reachable at 127.0.0.1:6379.",
        "Start one, then reopen the monitor:",
        "",
        "  Jarvis run bin/quickstart_bridson_operas.yaml",
        "  Jarvis monitor",
        "",
        "q quit",
    ]
    lines = [hline(width, "┌", "─", "┐")]
    pad_top = max(1, (height - 2 - len(msg)) // 2)
    for _ in range(pad_top):
        lines.append("│" + fit("", inner) + "│")
    for row in msg:
        lines.append("│" + fit(row, inner, "center" if not row.startswith("  ") else "left") + "│")
    while len(lines) < height - 1:
        lines.append("│" + fit("", inner) + "│")
    lines.append(hline(width, "└", "─", "┘"))
    return [fit(l, width) for l in lines[:height]]


def help_overlay(width: int, height: int, compact: bool) -> list[str]:
    """Help box drawn on top of a dim overview-like frame."""
    base = chrome(width, height, 1, overview_body(width, height - 6, compact), compact)
    box_w = min(64, width - 8)
    box_h = 12
    left = (width - box_w) // 2
    top = (height - box_h) // 2
    help_lines = [
        "┌" + fit(" keys", box_w - 2, "left").replace("keys", " keys ", 1) + "┐",
    ]
    # rebuild first line properly
    help_lines = ["┌" + fit(" KEYS ", box_w - 2, "center").replace(" ", "─") + "┐"]
    inner_h = box_h - 2
    content = [
        "1-7        jump to page",
        "tab        next page",
        "j / k      move highlight",
        "enter      pin row in the right pane",
        "esc        clear pin",
        "[ / ]      slower / faster refresh",
        "r          force one tick",
        "q          quit (scan keeps running)",
        "",
        "observer only. no kill, no pause, no log tail.",
    ]
    while len(content) < inner_h:
        content.append("")
    for row in content[:inner_h]:
        help_lines.append("│" + fit(" " + row, box_w - 2) + "│")
    help_lines.append("└" + ("─" * (box_w - 2)) + "┘")
    out = [list(line) for line in base]
    for i, hline_s in enumerate(help_lines):
        y = top + i
        for x, ch in enumerate(hline_s):
            out[y][left + x] = ch
    return ["".join(row) for row in out]


def write_page(path: Path, lines: list[str], width: int, height: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if len(lines) != height:
        raise RuntimeError(f"{path}: height {len(lines)} != {height}")
    for i, line in enumerate(lines):
        if len(line) != width:
            raise RuntimeError(f"{path}:{i+1}: width {len(line)} != {width} {line!r}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def page(width: int, height: int, current: int, body_fn, compact: bool) -> list[str]:
    return chrome(width, height, current, body_fn(width, height - 6, compact), compact)


BENCHMARKS = [
    ("min-80x24", 80, 24, True),
    ("narrow-80x48", 80, 48, True),
    ("ref-120x36", 120, 36, False),
    ("wide-160x48", 160, 48, False),
    ("ultra-220x64", 220, 64, False),
]


def main() -> None:
    # All 7 pages + picker/empty at min and ref.
    full_pages = {
        "1-overview": (1, overview_body),
        "2-workers": (2, workers_page_body),
        "3-factory": (3, factory_page_body),
        "4-sampler": (4, sampler_page_body),
        "5-calculators": (5, calculators_page_body),
        "6-samples": (6, samples_page_body),
        "7-host": (7, host_page_body),
    }
    archetypes = {"1-overview", "2-workers"}

    for folder, w, h, compact in BENCHMARKS:
        dest = ROOT / folder
        names = full_pages if folder in {"min-80x24", "ref-120x36"} else {k: full_pages[k] for k in archetypes}
        for name, (idx, fn) in names.items():
            write_page(dest / f"{name}.txt", page(w, h, idx, fn, compact), w, h)
        write_page(dest / "0-splash.txt", splash_page(w, h, compact), w, h)
        if folder in {"min-80x24", "ref-120x36"}:
            write_page(dest / "0-picker.txt", picker_page(w, h, compact), w, h)
            write_page(dest / "0-empty.txt", empty_page(w, h, compact), w, h)
        if folder == "ref-120x36":
            write_page(dest / "help.txt", help_overlay(w, h, compact), w, h)

    print("minted:")
    for p in sorted(ROOT.rglob("*.txt")):
        if p.name in {"BENCHMARKS.txt"}:
            continue
        lines = p.read_text(encoding="utf-8").splitlines()
        print(f"  {p.relative_to(ROOT)}  {len(lines[0]) if lines else 0}x{len(lines)}")


if __name__ == "__main__":
    main()
