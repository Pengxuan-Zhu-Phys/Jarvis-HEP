#!/usr/bin/env python3
"""Monitor TUI entry: splash chooser, branding, CLI --once."""

from __future__ import annotations

import asyncio
import io
import json
import math
import os
import re
import time
import unittest
from contextlib import redirect_stdout
from unittest import mock

from jarvishep2.client import build_parser, dispatch_monitor, main
from jarvishep2.monitor.branding import (
    load_branding,
    render_banner_markup,
    render_logo_monitor_frame,
)
from jarvishep2.monitor.hints import SPLASH_KEYS, render_key_hint, render_key_hint_plain
from jarvishep2.monitor.bars import bar_rail, eighths, visible_width
from jarvishep2.monitor.overview import (
    _balanced_spark_width,
    _fit_spark_series,
    _health_row,
    _spark_column_widths,
    _spark_area_rows,
    _sample_border_title,
    HealthItem,
    PacmanGame,
    render_adaptive_bridson_samples,
    render_health_block,
    render_queue_block,
    render_resource_row,
    render_status_block,
    render_spark_block,
)
from jarvishep2.monitor.metrics import format_sample_rate
from jarvishep2.monitor.queues import (
    QueueBlockData,
    QueueDirection,
    QueueDirectionWindow,
)
from jarvishep2.monitor.calculators import (
    CalculatorPoolData,
    CalculatorsBlockData,
    render_calculators_block,
)
from jarvishep2.monitor.workers import WorkerBlockData, render_workers_block
from jarvishep2.monitor.samples import (
    BUILTIN_SAMPLERS,
    SamplesBlockData,
    render_samples_block,
)
from jarvishep2.monitor.chrome import folder_tab_lines
from jarvishep2.monitor.collector import Collector
from jarvishep2.monitor.simu import SimuEngine
from jarvishep2.monitor.scans import (
    ScanChoice,
    choice_from_scan,
    list_scan_choices,
    resolve_choice,
    simulated_choice,
)
from jarvishep2.redis_queue import MONITOR_WANT, make_fakeredis_queue
from jarvishep2.monitor.topbar import (
    clock_label,
    compact_path,
    render_topbar_clock,
    render_topbar_left,
)
from jarvishep2.process_cleanup import JarvisProcess, JarvisScan
from jarvishep2.run_outcome import EXIT_OK, EXIT_USAGE


def _scan(reference: str, name: str, control_pid: int, extra: int = 2) -> JarvisScan:
    processes = [JarvisProcess(control_pid, f"Jarvis:{name}")]
    for index in range(extra):
        processes.append(
            JarvisProcess(
                control_pid + 1 + index,
                f"Jarvis-Worker-{index:02d}:{name}",
            )
        )
    return JarvisScan(
        reference=reference,
        name=name,
        processes=tuple(processes),
    )


class ResourceRowTests(unittest.TestCase):
    def test_cdot_is_on_the_center_column(self) -> None:
        row = render_resource_row(76, 41.0, 18.4, 64.0)
        self.assertEqual(visible_width(row), 76)
        plain = re.sub(r"\[/?[^\]]*\]", "", row)
        self.assertEqual(plain[76 // 2], "·")
        cpu_len = plain.index("┤") - plain.index("├")
        mem_len = plain.rindex("┤") - plain.rindex("├")
        self.assertGreater(cpu_len, mem_len)


class ProgressBarTests(unittest.TestCase):
    def test_rail_gradient_and_visible_width(self) -> None:
        self.assertEqual(len(eighths(0.5, 10, empty="─")), 10)
        rail = bar_rail(0.63, 20)
        shown = visible_width(rail)
        self.assertEqual(shown, 20)
        self.assertTrue(shown >= 2)
        self.assertIn("├", rail)
        self.assertIn("┤", rail)
        self.assertIn("#134a8d", rail)
        self.assertIn("#73b8f4", bar_rail(1.0, 12))


class MetricFormattingTests(unittest.TestCase):
    def test_sample_rate_uses_the_monitor_wide_adaptive_units(self) -> None:
        self.assertEqual(format_sample_rate(1.0), "1.0 / sec")
        self.assertEqual(format_sample_rate(0.5), "30.0 / min")
        self.assertEqual(format_sample_rate(1.0 / 120.0), "30.0 / hour")
        self.assertEqual(format_sample_rate(None), "—")


class StatusBlockTests(unittest.TestCase):
    def test_status_has_two_edge_aligned_rows_and_no_blank_row(self) -> None:
        block = render_status_block(
            78,
            "simu-iDM_Vector_V1",
            "RUNNING",
            "00:42:18",
            "simu-iDM_Vector_V1",
            "127.0.0.1:6379",
        )
        lines = block.splitlines()
        self.assertEqual(len(lines), 2)
        self.assertTrue(all(visible_width(line) == 78 for line in lines))
        self.assertTrue(lines[0].startswith("simu-iDM_Vector_V1"))
        self.assertTrue(lines[0].rstrip().endswith("00:42:18"))
        self.assertEqual(lines[0].count("●[/]"), 3)
        self.assertIn(" running", lines[0])
        self.assertTrue(lines[1].startswith("REF  simu-iDM_Vector_V1"))
        self.assertTrue(lines[1].rstrip().endswith("Redis  127.0.0.1:6379"))

    def test_running_status_dot_uses_a_breathing_color(self) -> None:
        white = render_status_block(78, "scan", "running", "00:00:01", "simu-iDM_Vector_V1", "redis", breathe_phase=0.0)
        yellow = render_status_block(78, "scan", "running", "00:00:01", "simu-iDM_Vector_V1", "redis", breathe_phase=math.pi)
        self.assertEqual(white.count("●[/]"), 3)
        self.assertEqual(yellow.count("●[/]"), 3)
        self.assertIn("[#e6e8eb]●[/]", white)
        self.assertIn("[#f6d33f]●[/]", yellow)


class HealthBlockTests(unittest.TestCase):
    def test_health_has_four_component_rows_and_healthy_title(self) -> None:
        from jarvishep2.monitor.simu import SimuEngine

        rendered = render_health_block(52, SimuEngine().empty_frame())
        lines = rendered.body.splitlines()
        self.assertEqual(len(lines), 4)
        self.assertTrue(all(visible_width(line) == 52 for line in lines))
        self.assertTrue(lines[0].startswith("[bold #f6d33f]💡[/] CORE"))
        self.assertIn("FACTORY", lines[1])
        self.assertIn("WORKERS", lines[2])
        self.assertIn("ARCHIVER", lines[3])
        self.assertIn("[#35c98a]●[/]", rendered.body)
        self.assertIn("HEALTHY · 4/4", rendered.border_title)
        self.assertEqual(rendered.border_subtitle, "No active warnings")
        narrow_lines = render_health_block(31, SimuEngine().empty_frame()).body.splitlines()
        self.assertTrue(all(visible_width(line) == 31 for line in narrow_lines))
        self.assertTrue(narrow_lines[3].endswith("queue 6"))

    def test_health_reports_a_stale_worker_as_the_active_warning(self) -> None:
        from dataclasses import replace
        from jarvishep2.monitor.simu import SimuEngine

        frame = replace(SimuEngine().empty_frame(), stale=1, workers_alive=47)
        rendered = render_health_block(52, frame)
        self.assertIn("DEGRADED · 3/4", rendered.border_title)
        self.assertIn("WORKERS", rendered.body)
        self.assertIn("1 stale", rendered.body)
        self.assertIn("[#f6d33f]●[/]", rendered.body)
        self.assertEqual(rendered.border_subtitle, "WARN · workers 1 stale")

    def test_help_trigger_stays_in_one_column_when_state_changes(self) -> None:
        marker_positions = []
        for state in ("HEALTHY", "DEGRADED", "CRITICAL"):
            row = _health_row(HealthItem("CORE", state, "dynamic fact", "fact"), 52)
            marker_prefix = row.split("💡", 1)[0]
            marker_positions.append(visible_width(marker_prefix))
        self.assertEqual(marker_positions, [0, 0, 0])


class QueueBlockTests(unittest.TestCase):
    def test_queue_block_shows_independent_core_paths_without_bars_or_total(self) -> None:
        rendered = render_queue_block(
            52,
            QueueBlockData(
                task_depth=24,
                archive_depth=6,
                feedback_mode="inactive",
                task_direction=QueueDirection("draining", -3.0),
                archive_direction=QueueDirection("rising", 1.0),
            ),
            outer_width=56,
            breathe_phase=math.pi,
        )
        lines = rendered.body.splitlines()
        self.assertEqual(len(lines), 6)
        self.assertTrue(all(visible_width(line) == 52 for line in lines))
        self.assertTrue(lines[0].startswith("[#35c98a]●[/]"))
        self.assertIn("TASK     → WORKERS", lines[0])
        self.assertTrue(re.sub(r"\[/?[^\]]*\]", "", lines[0]).endswith("24  ↓ 3 / min"))
        self.assertIn("ARCHIVE  → ARCHIVER", lines[2])
        self.assertIn("FEEDBACK → SAMPLER", lines[4])
        self.assertTrue(lines[4].endswith("—  inactive"))
        self.assertNotIn("pending", rendered.body.lower())
        self.assertNotIn("├", rendered.body)
        self.assertTrue(rendered.border_title.endswith("FLOW"))
        self.assertEqual(rendered.border_subtitle, "depths are separate paths")

    def test_queue_direction_window_suppresses_one_item_jitter(self) -> None:
        window = QueueDirectionWindow()
        self.assertEqual(
            window.observe("task", 24, observed_at=0.0).state,
            "warming",
        )
        self.assertEqual(
            window.observe("task", 25, observed_at=8.0).state,
            "flat",
        )
        direction = window.observe("task", 20, observed_at=8.0)
        self.assertEqual(direction.state, "draining")
        self.assertEqual(direction.label(), "↓ 30 / min")

    def test_simu_feedback_mode_matches_the_selected_sampler(self) -> None:
        from jarvishep2.monitor.simu import SimuEngine

        engine = SimuEngine()
        self.assertEqual(engine.empty_frame().queue_block.feedback_mode, "inactive")
        engine.set_sampler_method("MCMC")
        sharded = engine.tick(6)
        self.assertEqual(sharded.queue_block.feedback_mode, "sharded")
        self.assertEqual(len(sharded.queue_block.feedback_shards), 4)
        self.assertEqual(
            sharded.queue_block.feedback_depth,
            sum(sharded.queue_block.feedback_shards.values()),
        )
        engine.set_sampler_method("PTMCMC")
        shared = engine.tick(6)
        self.assertEqual(shared.queue_block.feedback_mode, "shared")
        self.assertEqual(shared.queue_block.feedback_shards, {})
        self.assertIsNotNone(shared.queue_block.feedback_depth)

    def test_live_collector_projects_only_metadata_known_feedback_shards(self) -> None:
        queue = make_fakeredis_queue()
        queue.connect()
        assert queue.r is not None
        queue.r.rpush("hep:feedback:chain:0", "known")
        queue.r.rpush("hep:feedback:chain:7", "not-admitted")
        collector = Collector(
            queue,
            sampler_metadata={"method": "MCMC", "config": {"num_chains": 2}},
            monotonic=lambda: 1.0,
        )
        block = collector.tick("overview").queue_block
        self.assertEqual(block.feedback_mode, "sharded")
        self.assertEqual(block.feedback_shards, {"0": 1, "1": 0})
        self.assertEqual(block.feedback_depth, 1)


class WorkersBlockTests(unittest.TestCase):
    def test_workers_canvas_has_no_padding_row(self) -> None:
        rendered = render_workers_block(52, SimuEngine().empty_frame().worker_block)
        self.assertEqual(len(rendered.body.splitlines()), 6)

    def test_workers_block_makes_file_operator_a_first_class_fleet_signal(self) -> None:
        rendered = render_workers_block(
            52,
            WorkerBlockData(
                expected=48,
                present=48,
                busy=12,
                idle=36,
                assigned=12,
                file_ops_expected=48,
                file_ops_alive=47,
                file_ops_missing=1,
                heartbeat_recent=48,
            ),
            outer_width=56,
            breathe_phase=math.pi,
        )
        lines = rendered.body.splitlines()
        self.assertEqual(len(lines), 6)
        self.assertTrue(all(visible_width(line) == 52 for line in lines))
        self.assertIn("FILE OPS", lines[4])
        self.assertIn("47 / 48 alive", lines[4])
        self.assertTrue(lines[4].startswith("[#"))
        self.assertNotIn("[#35c98a]●[/]", lines[4])
        self.assertTrue(rendered.border_title.endswith("FLEET"))

    def test_collector_checks_only_advertised_file_operator_pids(self) -> None:
        queue = make_fakeredis_queue()
        queue.connect()
        queue.publish_proc_board("core", role="core", workers_total=1)
        queue.heartbeat(
            "3",
            status="busy",
            pid=110,
            ts=1.0,
            board_ttl_sec=30,
            current_sample="sample-3",
            file_operation_pid=220,
            file_operation_pgid=220,
            file_operation_mode="process",
        )
        checked: list[int] = []

        def _alive(pid: int) -> bool:
            checked.append(pid)
            return pid == 110

        block = Collector(queue, owner_ids=["3"], pid_alive=_alive).tick("overview").worker_block
        self.assertEqual(block.expected, 1)
        self.assertEqual(block.present, 1)
        self.assertEqual(block.file_ops_expected, 1)
        self.assertEqual(block.file_ops_alive, 0)
        self.assertEqual(block.file_ops_missing, 1)
        self.assertEqual(checked, [110, 220])


class CalculatorsBlockTests(unittest.TestCase):
    def test_constellation_keeps_busy_pack_positions_without_progress_rails(self) -> None:
        rendered = render_calculators_block(
            52,
            CalculatorsBlockData(
                pools=(CalculatorPoolData("SoftSUSY", 8, (2, 6)),),
            ),
            outer_width=56,
            pulse_phase=0.0,
        )
        lines = rendered.body.splitlines()
        plain = re.sub(r"\[/?[^\]]*\]", "", lines[2])
        self.assertEqual(len(lines), 3)
        self.assertEqual(visible_width(lines[1]), 52)
        self.assertEqual(visible_width(lines[2]), 52)
        self.assertIn("SoftSUSY", lines[1])
        self.assertIn("ACTIVE", lines[1])
        self.assertEqual(plain.count("✦"), 2)
        self.assertEqual(plain.count("☆"), 6)
        self.assertEqual(plain.count("✦") + plain.count("☆"), 8)
        self.assertTrue(re.sub(r"\[/?[^\]]*\]", "", lines[1]).endswith("2 / 8"))
        self.assertNotIn("█", rendered.body)
        self.assertNotIn("░", rendered.body)
        self.assertEqual(rendered.border_subtitle, "✦ white pulse = occupied · ☆ quiet = free")

    def test_simu_busy_packs_are_not_synthesized_as_a_left_aligned_bar(self) -> None:
        frame = SimuEngine().empty_frame()
        pools = {pool.name: pool for pool in frame.calculator_block.pools}
        softsusy = pools["SoftSUSY"]
        self.assertNotEqual(softsusy.busy_packs, tuple(range(1, len(softsusy.busy_packs) + 1)))

    def test_simu_gives_every_pool_a_status_line_and_a_star_line(self) -> None:
        frame = SimuEngine().empty_frame()
        rendered = render_calculators_block(52, frame.calculator_block)
        lines = rendered.body.splitlines()
        self.assertEqual(len(lines), 1 + 2 * len(frame.calculator_block.pools))
        for index, pool in enumerate(frame.calculator_block.pools):
            status = re.sub(r"\[/?[^\]]*\]", "", lines[1 + 2 * index])
            stars = re.sub(r"\[/?[^\]]*\]", "", lines[2 + 2 * index])
            self.assertIn(pool.name, status)
            self.assertRegex(status, r"(?:ACTIVE|QUIET)\s+\d+ / \d+$")
            self.assertNotIn(pool.name, stars)
            self.assertEqual(visible_width(lines[2 + 2 * index]), 52)
        states = [
            re.sub(r"\[/?[^\]]*\]", "", lines[1 + 2 * index])
            for index in range(len(frame.calculator_block.pools))
        ]
        state_columns = [
            line.index("ACTIVE") if "ACTIVE" in line else line.index("QUIET")
            for line in states
        ]
        self.assertEqual(len(set(state_columns)), 1)

    def test_star_positions_are_irregular_but_stable_while_they_pulse(self) -> None:
        frame = SimuEngine().empty_frame()
        first = render_calculators_block(52, frame.calculator_block, pulse_phase=0.0).body
        second = render_calculators_block(52, frame.calculator_block, pulse_phase=0.65).body

        def star_positions(body: str) -> tuple[int, ...]:
            plain = re.sub(r"\[/?[^\]]*\]", "", body).splitlines()[2]
            return tuple(index for index, char in enumerate(plain) if char in "✦☆")

        positions = star_positions(first)
        self.assertEqual(positions, star_positions(second))
        self.assertEqual(len(positions), 16)
        gaps = {right - left for left, right in zip(positions, positions[1:])}
        self.assertGreater(len(gaps), 1)


class SamplesBlockTests(unittest.TestCase):
    def test_adaptive_bridson_uses_three_dense_edge_aligned_rows(self) -> None:
        from jarvishep2.monitor.simu import SimuEngine

        frame = SimuEngine().empty_frame()
        block = render_adaptive_bridson_samples(52, frame)
        lines = block.splitlines()
        self.assertEqual(len(lines), 3)
        self.assertTrue(all(visible_width(line) == 52 for line in lines))
        self.assertTrue(lines[0].startswith("GENERATION 8 / 25"))
        self.assertTrue(lines[0].endswith("32%"))
        self.assertIn("OPEN 4 · PARTIAL", lines[1])
        self.assertTrue(lines[2].startswith("Acc 1842 · RUN 12 · FAIL 7"))
        self.assertTrue(lines[2].endswith("RATE —"))
        self.assertNotIn("remain", block.lower())
        self.assertNotIn("avg", block.lower())

    def test_sampler_name_occupies_the_right_side_of_the_top_edge(self) -> None:
        title = _sample_border_title(56, "AdaptiveBridson")
        self.assertEqual(len(title), 50)
        self.assertTrue(title.startswith("SAMPLES ─"))
        self.assertTrue(title.endswith("AdaptiveBridson"))

    def test_every_builtin_sampler_has_a_fixed_dedicated_canvas(self) -> None:
        from jarvishep2.monitor.simu import SimuEngine

        engine = SimuEngine()
        for method in BUILTIN_SAMPLERS:
            with self.subTest(method=method):
                engine.set_sampler_method(method)
                frame = engine.empty_frame()
                rendered = render_samples_block(
                    52,
                    SamplesBlockData(
                        completed=frame.done,
                        running=frame.running,
                        failed=frame.failed,
                        rate=frame.rate,
                        sampler=frame.sampler,
                    ),
                )
                lines = rendered.body.splitlines()
                self.assertEqual(len(lines), 3)
                self.assertTrue(all(visible_width(line) == 52 for line in lines))
                self.assertTrue(lines[2].startswith("Acc "))
                self.assertTrue(rendered.border_title.endswith(method))
                self.assertTrue(rendered.border_subtitle)

class SparkBlockTests(unittest.TestCase):
    def test_spark_columns_use_the_available_width_with_one_right_padding(self) -> None:
        self.assertEqual(_balanced_spark_width(69, 2), 69)
        self.assertEqual(_spark_column_widths(69, 2), [32, 32])
        self.assertEqual(_spark_column_widths(68, 2), [31, 32])

    def test_display_order_keeps_newest_value_at_the_right_edge(self) -> None:
        self.assertEqual(_fit_spark_series((4.0, 3.0, 2.0, 1.0), 4), [1.0, 2.0, 3.0, 4.0])
        with self.assertRaises(ValueError):
            _fit_spark_series((4.0, 3.0), 4)

    def test_area_spark_fills_only_from_the_bottom(self) -> None:
        rows = _spark_area_rows((0, 1, 2, 3, 4), 5)
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(len(row) == 5 for row in rows))
        for column in range(5):
            bottom_to_top = [row[column] != " " for row in reversed(rows)]
            seen_empty = False
            for occupied in bottom_to_top:
                if not occupied:
                    seen_empty = True
                elif seen_empty:
                    self.fail(f"spark column {column} has a filled gap")

    def test_spark_block_has_two_columns_and_four_plot_rows(self) -> None:
        block = render_spark_block(
            72,
            (
                ("SAMPLES", "12.4 /min", tuple(range(33)), 20.0),
                ("TASK QUEUE", "24", tuple(range(34)), 32.0),
            ),
        )
        lines = block.splitlines()
        self.assertEqual(len(lines), 5)
        self.assertTrue(all(visible_width(line) == 72 for line in lines))
        cells = lines[0].split("  ·  ")
        self.assertTrue(cells[0].startswith("SAMPLES"))
        self.assertTrue(cells[0].rstrip().endswith("12.4 /min"))
        self.assertTrue(cells[1].lstrip().startswith("TASK QUEUE"))
        self.assertTrue(cells[1].rstrip().endswith("24"))
        self.assertNotIn("HOST CPU", lines[0])
        self.assertIn("[#73b8f4]", lines[1])
        self.assertIn("[#134a8d]", lines[4])
        self.assertEqual(lines[1].count("[#73b8f4]"), 2)
        self.assertEqual(lines[4].count("[#134a8d]"), 2)


class PacmanGameTests(unittest.TestCase):
    def test_each_row_starts_with_one_or_two_breakable_blocks(self) -> None:
        game = PacmanGame()
        game._ensure_board(12, 3)

        self.assertTrue(all(1 <= len(blocks) <= 2 for blocks in game._blocks))

    def test_pacman_breaks_blocks_but_ghosts_turn_around(self) -> None:
        game = PacmanGame()
        game._ensure_board(12, 3)
        game._blocks = [{0}, set(), set()]

        game._advance_game()

        self.assertNotIn(0, game._blocks[0])
        ghost = game._ghosts[0]
        ghost.row, ghost.column, ghost.direction = 0, 10, 1
        game._blocks[1].add(0)
        game._move_ghosts(2, 0, 12)

        self.assertEqual((ghost.row, ghost.column, ghost.direction), (0, 10, -1))

    def test_caught_pacman_does_not_pause_pellet_recovery(self) -> None:
        game = PacmanGame()
        game._ensure_board(12, 2)
        game._pellets[0][1] = False
        game._restore_queue.append((0, 1))
        game._caught_ticks = 1

        game._advance_game()

        self.assertTrue(game._pellets[0][1])
        self.assertFalse(game._restore_queue)

    def test_pellet_recovery_spawns_one_or_two_new_blocks_in_one_row(self) -> None:
        game = PacmanGame()
        game._ensure_board(12, 3)
        game._blocks = [set(), set(), set()]
        game._pellets[0][0] = False

        game._queue_row_restore(0)

        self.assertIn((0, 0), game._restore_queue)
        block_count = sum(len(blocks) for blocks in game._blocks)
        self.assertIn(block_count, (1, 2))

    def test_ghosts_loop_between_rows_without_turning_at_edges(self) -> None:
        game = PacmanGame()
        game._ensure_board(12, 3)
        game._blocks = [set(), set(), set()]
        ghost = game._ghosts[0]
        ghost.row, ghost.column, ghost.direction = 0, 10, 1

        game._move_ghosts(2, 0, 12)
        self.assertEqual((ghost.row, ghost.column, ghost.direction), (1, 0, 1))

        ghost.row, ghost.column, ghost.direction = 0, 0, -1
        game._move_ghosts(2, 0, 12)
        self.assertEqual((ghost.row, ghost.column, ghost.direction), (2, 10, -1))

    def test_giant_spawns_as_a_random_large_dot(self) -> None:
        game = PacmanGame()
        game._ensure_board(12, 2)
        game._giant_spawn_ticks = 0
        initial_position = game._path[game._path_index]

        game._advance_game()

        self.assertIsNotNone(game._giant)
        assert game._giant is not None
        self.assertNotEqual(game._giant, initial_position)

    def test_caught_pacman_holds_position_for_one_second(self) -> None:
        game = PacmanGame()
        game._ensure_board(12, 2)
        row, column = game._path[game._path_index]
        ghost = game._ghosts[0]
        ghost.row = row
        ghost.column = column
        ghost.direction = 1

        self.assertTrue(game._handle_collisions(row, column))
        self.assertEqual(game._caught_ticks, game._CAUGHT_TICKS)
        self.assertEqual(ghost.direction, 1)
        path_index = game._path_index
        game._advance_game()

        self.assertEqual(game._path_index, path_index)
        self.assertEqual(game._caught_ticks, game._CAUGHT_TICKS - 1)


class FolderTabTests(unittest.TestCase):
    def test_release_tab_has_no_number_and_joins_the_right_frame(self) -> None:
        _top, mid, join, hits = folder_tab_lines(
            78,
            1,
            pages=(("overview", "Overview"),),
        )
        self.assertIn("│     Overview     │", mid)
        self.assertNotIn("1 Overview", mid)
        self.assertEqual(hits, [("overview", 0, len("│     Overview     │"))])
        self.assertEqual(join, "│" + " " * 18 + "╰" + "─" * 57 + "╮")

    def test_active_stays_in_order_and_idles_share_space(self) -> None:
        _top, mid, join, hits = folder_tab_lines(78, 1)
        self.assertIn("│ 1 Overview │", mid)
        self.assertEqual(hits[0][0], "overview")
        self.assertEqual(hits[0][1], 0)
        self.assertEqual(join[0], "│")
        self.assertEqual(join[-1], "╮")

        _top, mid, join, hits = folder_tab_lines(78, 2)
        self.assertIn("│ 2 Workers │", mid)
        self.assertEqual(hits[1][0], "workers")
        self.assertGreater(hits[1][1], 0)
        self.assertLess(hits[0][1], hits[1][1])
        self.assertEqual(join[0], "╭")
        self.assertIn("╯", join)
        self.assertIn("╰", join)
        self.assertNotIn("┌", _top)
        self.assertIn("╭", _top)

    def test_idle_names_shorten_only_when_needed_then_gaps(self) -> None:
        _top, mid, _join, _hits = folder_tab_lines(120, 1)
        self.assertIn("Calculators", mid)
        self.assertNotIn("Cal.", mid)
        without_active = mid.split("│", 2)[-1]
        self.assertIn("     ", without_active)

        _top, mid, _join, _hits = folder_tab_lines(78, 1)
        self.assertIn("Cal.", mid)
        self.assertNotIn("Calculators", mid)
        self.assertNotIn("Cal.6", mid)

        _top, mid, _join, _hits = folder_tab_lines(61, 1)
        self.assertNotIn("Calculators", mid)
        self.assertNotIn("Cal.", mid)
        self.assertIn("Ca.", mid)


class SimulatedScanTests(unittest.TestCase):
    def test_fixture_is_marked_simulated(self) -> None:
        choice = simulated_choice()
        self.assertTrue(choice.simulated)
        self.assertEqual(choice.reference, "SIM")
        self.assertTrue(choice.name.startswith("simu-"))

    def test_engine_ticks_change_the_frame(self) -> None:
        from jarvishep2.monitor.simu import SimuEngine

        engine = SimuEngine()
        first = engine.tick()
        for _ in range(6):
            second = engine.tick()
        self.assertGreaterEqual(second.done, first.done)
        self.assertNotEqual(first.spark_samples, second.spark_samples)

    def test_spark_history_inserts_newest_at_right_edge_and_drops_oldest(self) -> None:
        from jarvishep2.monitor.simu import SimuEngine

        engine = SimuEngine()
        self.assertEqual(engine.empty_frame().spark_samples, ())
        first = engine.tick(4)
        second = engine.tick(4)
        self.assertEqual(len(first.spark_samples), 4)
        self.assertEqual(len(second.spark_samples), 4)
        self.assertIsNone(first.spark_samples[-1])
        self.assertAlmostEqual(second.spark_samples[1], first.spark_samples[0])

    def test_spark_history_width_follows_the_display_column(self) -> None:
        from jarvishep2.monitor.simu import SimuEngine

        engine = SimuEngine()
        engine.tick(6)
        self.assertEqual(len(engine.current_frame().spark_samples), 6)
        engine.set_history_width(3)
        resized = engine.current_frame()
        self.assertEqual(len(resized.spark_samples), 3)
        engine.set_history_width(8)
        expanded = engine.current_frame()
        self.assertEqual(len(expanded.spark_samples), 8)
        self.assertEqual(expanded.spark_samples[:3], resized.spark_samples)
        self.assertTrue(all(value is None for value in expanded.spark_samples[3:]))

    def test_spark_histories_can_use_the_odd_width_remainder_independently(self) -> None:
        from jarvishep2.monitor.simu import SimuEngine

        engine = SimuEngine()
        frame = engine.tick((3, 4))
        self.assertEqual(len(frame.spark_samples), 3)
        self.assertEqual(len(frame.spark_queue), 4)


class ScanChoiceTests(unittest.TestCase):
    def test_projects_ps_fields(self) -> None:
        scan = _scan("R1", "iDM_Vector_V1", 44001, extra=3)
        choice = choice_from_scan(scan)
        self.assertEqual(choice.reference, "R1")
        self.assertEqual(choice.name, "iDM_Vector_V1")
        self.assertEqual(choice.control_pid, 44001)
        self.assertEqual(choice.process_count, 4)

    def test_list_and_resolve(self) -> None:
        scans = [
            _scan("R1", "iDM_Vector_V1", 44001),
            _scan("R2", "eggbox-bridson", 51020),
        ]
        choices = list_scan_choices(lambda: scans)
        self.assertEqual([choice.reference for choice in choices], ["R1", "R2"])
        picked = resolve_choice("R2", choices, scans=scans)
        self.assertEqual(picked.name, "eggbox-bridson")
        picked = resolve_choice("44001", choices, scans=scans)
        self.assertEqual(picked.reference, "R1")


class KeyHintTests(unittest.TestCase):
    def test_markup_is_key_colon_meaning_with_bars(self) -> None:
        markup = render_key_hint(SPLASH_KEYS)
        self.assertIn("[bold #f6d33f]Enter[/]", markup)
        self.assertIn("[#8d93a1]: attach[/]", markup)
        self.assertIn("|", markup)
        self.assertIn("J/K", markup)
        self.assertIn("Q", markup)
        self.assertEqual(
            render_key_hint_plain(SPLASH_KEYS),
            "Enter: attach  |  J/K: select  |  R: refresh  |  Q: quit",
        )
        plain = render_key_hint_plain(SPLASH_KEYS)
        self.assertNotIn("S:", plain)
        self.assertNotIn("simu", plain.lower())


class TopbarTests(unittest.TestCase):
    def test_clock_and_left_piece(self) -> None:
        from datetime import datetime

        self.assertEqual(clock_label(datetime(2026, 9, 10, 20, 53)), "8:53 PM")
        left = render_topbar_left(branch="main", path="~/Jarvis-HEP")
        self.assertIn("⎇ main", left)
        self.assertIn("~/Jarvis-HEP", left)
        self.assertIn(" · ", left)
        self.assertFalse(left.endswith("·"))
        bare = render_topbar_left(branch="", path="~/Jarvis-HEP")
        self.assertEqual(bare, "~/Jarvis-HEP")
        self.assertNotIn("⎇", bare)
        self.assertNotIn("·", bare)
        self.assertEqual(
            render_topbar_clock(datetime(2026, 9, 10, 20, 53)),
            "· 8:53 PM",
        )
        self.assertEqual(compact_path("abcdefghij", 7), "...ghij")


class BrandingTests(unittest.TestCase):
    def test_logo_and_banner_come_from_hep_card(self) -> None:
        branding = load_branding()
        self.assertEqual(len(branding.logo_pattern), 8)
        self.assertTrue(any("Version:" in line for line in branding.banner_lines))
        markup = render_banner_markup(branding.banner_lines)
        self.assertIn("2.0.", markup)
        frame = render_logo_monitor_frame(40, branding.logo_pattern, animate=False)
        lines = frame.splitlines()
        self.assertEqual(len(lines), 8)
        self.assertTrue(all("⬤" in line for line in lines))


class CliMonitorOnceTests(unittest.TestCase):
    def test_once_without_ref_prints_scan_table(self) -> None:
        scans = [_scan("R1", "iDM_Vector_V1", 44001)]
        args = build_parser().parse_args(["monitor", "--once"])
        buf = io.StringIO()
        with mock.patch(
            "jarvishep2.process_cleanup.list_active_scans",
            return_value=scans,
        ):
            with redirect_stdout(buf):
                code = dispatch_monitor(args)
        self.assertEqual(code, EXIT_OK)
        self.assertIn("iDM_Vector_V1", buf.getvalue())
        self.assertIn("R1", buf.getvalue())

    def test_once_json_emits_rows(self) -> None:
        scans = [_scan("R1", "iDM_Vector_V1", 44001)]
        args = build_parser().parse_args(["monitor", "--once", "--json"])
        buf = io.StringIO()
        with mock.patch(
            "jarvishep2.process_cleanup.list_active_scans",
            return_value=scans,
        ):
            with redirect_stdout(buf):
                code = dispatch_monitor(args)
        self.assertEqual(code, EXIT_OK)
        payload = json.loads(buf.getvalue())
        self.assertEqual(payload[0]["reference"], "R1")
        self.assertEqual(payload[0]["control_pid"], 44001)

    def test_unknown_ref_with_once_is_usage(self) -> None:
        args = build_parser().parse_args(["monitor", "R9", "--once"])
        with mock.patch(
            "jarvishep2.process_cleanup.list_active_scans",
            return_value=[],
        ):
            code = dispatch_monitor(args)
        self.assertEqual(code, EXIT_USAGE)

    def test_non_tty_monitor_lists_scans(self) -> None:
        scans = [_scan("R1", "iDM_Vector_V1", 44001)]
        buf = io.StringIO()
        with mock.patch(
            "jarvishep2.process_cleanup.list_active_scans",
            return_value=scans,
        ):
            with mock.patch("sys.stdout.isatty", return_value=False):
                with redirect_stdout(buf):
                    code = main(["monitor"])
        self.assertEqual(code, EXIT_OK)
        self.assertIn("iDM_Vector_V1", buf.getvalue())


class MonitorAppTests(unittest.TestCase):
    def test_splash_lists_ps_rows_and_quit(self) -> None:
        try:
            from jarvishep2.monitor.app import MonitorApp
        except ImportError:
            self.skipTest("textual is not installed")

        scans = [
            _scan("R1", "iDM_Vector_V1", 44001),
            _scan("R2", "eggbox-bridson", 51020),
        ]

        async def _run() -> None:
            app = MonitorApp(scan_lister=lambda: scans)
            async with app.run_test() as pilot:
                await pilot.pause()
                from textual.widgets import DataTable, Static

                logo = app.screen.query_one("#logo-monitor", Static)
                self.assertIn("⬤", str(logo.render()))
                table = app.screen.query_one("#scans", DataTable)
                self.assertEqual(table.row_count, 2)
                hint = app.screen.query_one("#hint", Static)
                hint_text = str(hint.render())
                self.assertIn("Enter", hint_text)
                self.assertIn("attach", hint_text)
                self.assertNotIn("simu", hint_text.lower())
                self.assertNotRegex(hint_text, r"\bS:")
                top = str(app.screen.query_one("#topbar-left", Static).render())
                self.assertIn("⎇", top)
                clock = str(app.screen.query_one("#topbar-clock", Static).render())
                self.assertRegex(clock, r"· \d{1,2}:\d{2} (AM|PM)")
                await pilot.press("q")

        asyncio.run(_run())

    def test_escape_does_not_quit_splash(self) -> None:
        try:
            from jarvishep2.monitor.app import MonitorApp
        except ImportError:
            self.skipTest("textual is not installed")

        scans = [_scan("R1", "iDM_Vector_V1", 44001)]

        async def _run() -> None:
            app = MonitorApp(scan_lister=lambda: scans)
            async with app.run_test() as pilot:
                await pilot.pause()
                from textual.widgets import DataTable, Static

                await pilot.press("escape")
                await pilot.pause()
                self.assertFalse(app._exit)
                app.screen.query_one("#logo-monitor", Static)
                table = app.screen.query_one("#scans", DataTable)
                scan_col = next(
                    col for col in table.ordered_columns if col.key.value == "scan"
                )
                self.assertFalse(scan_col.auto_width)
                self.assertGreaterEqual(scan_col.width, 8)
                await pilot.press("q")

        asyncio.run(_run())

    def test_empty_scan_list_shows_dash_not_start_hint(self) -> None:
        try:
            from jarvishep2.monitor.app import MonitorApp
        except ImportError:
            self.skipTest("textual is not installed")

        async def _run() -> None:
            app = MonitorApp(scan_lister=lambda: [])
            async with app.run_test() as pilot:
                await pilot.pause()
                from textual.widgets import DataTable, Static

                table = app.screen.query_one("#scans", DataTable)
                self.assertGreaterEqual(table.row_count, 1)
                scan_cell = str(table.get_row_at(0)[1])
                self.assertEqual(scan_cell.strip(), "-")
                notice = str(app.screen.query_one("#notice", Static).render())
                self.assertNotIn("Start one", notice)
                await pilot.press("q")

        asyncio.run(_run())

    def test_direct_ref_opens_workspace(self) -> None:
        try:
            from jarvishep2.monitor.app import MonitorApp
            from jarvishep2.monitor.overview_workspace import OverviewReleaseScreen
        except ImportError:
            self.skipTest("textual is not installed")

        scans = [_scan("R1", "iDM_Vector_V1", 44001)]

        async def _run() -> None:
            app = MonitorApp(scan_ref="R1", scan_lister=lambda: scans)
            async with app.run_test() as pilot:
                await pilot.pause()
                self.assertIsInstance(app.screen, OverviewReleaseScreen)
                self.assertEqual(app.screen.choice.reference, "R1")
                self.assertEqual(app.screen.choice.name, "iDM_Vector_V1")
                app.screen.query_one("#tab-bar")
                await pilot.press("q")

        asyncio.run(_run())

    def test_live_attach_uses_overview_only_release_workspace(self) -> None:
        try:
            from jarvishep2.monitor.app import MonitorApp
            from jarvishep2.monitor.overview_workspace import OverviewReleaseScreen
            from textual.widgets import ContentSwitcher
        except ImportError:
            self.skipTest("textual is not installed")

        queue = make_fakeredis_queue()
        queue.register_calc_pool("DemoCalc", 4)
        queue._acquire_calc("DemoCalc", timeout=1, worker_id="0")
        queue.publish_proc_board(
            "core",
            role="core",
            pid=44001,
            host="test-host",
            status="running",
            scan_mode="running",
            ts=time.time(),
            started_at=time.time() - 65,
            workers_total=1,
            sampler_status=json.dumps(
                {
                    "schema": 1,
                    "method": "Random",
                    "state": "running",
                    "progress": {"kind": "finite", "current": 4, "target": 10},
                    "metrics": {"accepted": 4, "seed": 7},
                }
            ),
        )
        queue.publish_proc_board("archiver", role="archiver", status="running")
        queue.publish_proc_board("redis", role="redis", status="running")
        queue.publish_proc_board(
            "worker",
            owner_id="0",
            pid=44004,
            status="busy",
            current_uuid="sample-4",
            board_ttl_sec=30,
            file_operation_mode="inline",
        )
        queue.r.hset(
            "hep:sample:stats",
            mapping={"completed": 4, "running": 1, "failed": 0},
        )
        collector = Collector(
            queue,
            pid=9,
            owner_ids=["0"],
            pool_names=["DemoCalc"],
            slots={"DemoCalc": 4},
            process_inventory=[{"role": "worker", "pid": 44004}],
            pid_alive=lambda _pid: True,
            host_snapshotter=lambda: {
                "available": True,
                "cpu_percent": 25.0,
                "memory_used": 8 * 1024**3,
                "memory_total": 32 * 1024**3,
                "processes": [{"fds": 12, "fd_limit": 1024}],
            },
            sampler_metadata={
                "method": "Random",
                "family": "simple",
                "dimensions": 2,
                "config": {"point_number": 10, "seed": 7},
            },
        )
        choice = ScanChoice(
            reference="R1",
            name="live-scan",
            control_pid=os.getpid(),
            process_count=1,
            pids=(os.getpid(),),
            simulated=False,
        )

        async def _run() -> None:
            app = MonitorApp(scan_lister=lambda: [])
            async with app.run_test() as pilot:
                await pilot.pause()
                app.attach_scan(choice, collector=collector)
                await pilot.pause()
                self.assertIsInstance(app.screen, OverviewReleaseScreen)
                self.assertEqual(
                    app.screen.query_one("#pages", ContentSwitcher).current,
                    "overview",
                )
                self.assertEqual(len(app.screen.query("#live-overview-body")), 0)
                status = str(app.screen.query_one("#ov-status").render())
                self.assertIn("live-scan", status)
                self.assertIn("REF  live-scan", status)
                self.assertIn("Redis  localhost:6379", status)
                resources = app.screen.query_one("#ov-resources")
                self.assertIn("25%", str(resources.render()))
                self.assertIn("8.0/32G", str(resources.render()))
                self.assertIn("12/1024", str(resources.border_subtitle))
                samples = app.screen.query_one("#ov-samples")
                self.assertTrue(str(samples.border_title).endswith("Random"))
                self.assertIn("CANDIDATES 4 / 10", str(samples.render()))
                queues = app.screen.query_one("#ov-queues")
                self.assertIn("TASK", str(queues.render()))
                workers = app.screen.query_one("#ov-workers")
                self.assertIn("1 / 1 present", str(workers.render()))
                calculators = app.screen.query_one("#ov-calcs")
                self.assertIn("DemoCalc", str(calculators.render()))
                self.assertIn("1 / 4", str(calculators.render()))
                health = app.screen.query_one("#ov-health")
                self.assertIn("HEALTHY", str(health.border_title))
                for hidden_id in (
                    "#workers",
                    "#factory",
                    "#sampler",
                    "#calculators",
                    "#samples",
                    "#host",
                ):
                    self.assertEqual(len(app.screen.query(hidden_id)), 0)
                await pilot.press("2", "7", "right")
                await pilot.pause()
                self.assertEqual(
                    app.screen.query_one("#pages", ContentSwitcher).current,
                    "overview",
                )
                await pilot.press("q")

        asyncio.run(_run())
        self.assertEqual(int(queue.r.exists(MONITOR_WANT) or 0), 0)

    def test_s_attaches_simulated_scan_without_hint(self) -> None:
        try:
            from jarvishep2.monitor.app import MonitorApp
        except ImportError:
            self.skipTest("textual is not installed")

        async def _run() -> None:
            app = MonitorApp(scan_lister=lambda: [])
            async with app.run_test() as pilot:
                await pilot.pause()
                from textual.widgets import Static

                hint = str(app.screen.query_one("#hint", Static).render())
                self.assertNotIn("simu", hint.lower())
                await pilot.press("s")
                await pilot.pause()
                self.assertEqual(
                    app.screen.__class__.__name__,
                    "OverviewReleaseScreen",
                )
                status = app.screen.query_one("#ov-status")
                self.assertEqual(status.border_title, "STATUS")
                self.assertIn("simu-iDM_Vector_V1", str(status.render()))
                self.assertIn("REF  simu-iDM_Vector_V1", str(status.render()))
                samples = app.screen.query_one("#ov-samples")
                self.assertTrue(str(samples.border_title).startswith("SAMPLES"))
                self.assertTrue(str(samples.border_title).endswith("AdaptiveBridson"))
                sample_lines = str(samples.render()).splitlines()
                self.assertEqual(len(sample_lines), 3)
                self.assertIn("GENERATION", sample_lines[0])
                self.assertIn("OPEN", sample_lines[1])
                self.assertIn("RATE", sample_lines[2])
                self.assertTrue(sample_lines[2].startswith("Acc "))
                self.assertNotIn("avg", str(samples.render()).lower())
                self.assertIn("CORE", str(samples.border_subtitle))
                self.assertIn("SAVED", str(samples.border_subtitle))
                await pilot.press("m")
                await pilot.pause()
                self.assertTrue(str(samples.border_title).endswith("MCMC"))
                resources = app.screen.query_one("#ov-resources")
                self.assertEqual(resources.border_title, "RESOURCES")
                self.assertIn("CPU", str(resources.render()))
                self.assertIn("MEM", str(resources.render()))
                self.assertIn("/", resources.border_subtitle)
                self.assertIn("├", str(resources.render()))
                health = app.screen.query_one("#ov-health")
                self.assertTrue(str(health.border_title).startswith("HEALTH"))
                self.assertTrue(str(app.screen.query_one("#ov-health-core").render()).startswith("💡 CORE"))
                self.assertTrue(str(app.screen.query_one("#ov-health-factory").render()).startswith("💡 FACTORY"))
                self.assertTrue(str(app.screen.query_one("#ov-health-workers").render()).startswith("💡 WORKERS"))
                self.assertTrue(str(app.screen.query_one("#ov-health-archiver").render()).startswith("💡 ARCHIVER"))
                self.assertEqual(health.border_subtitle, "No active warnings")
                queues = app.screen.query_one("#ov-queues")
                queue_lines = str(queues.render()).splitlines()
                self.assertTrue(str(queues.border_title).endswith("FLOW"))
                self.assertEqual(len(queue_lines), 6)
                self.assertIn("TASK", queue_lines[0])
                self.assertIn("WORKERS", queue_lines[0])
                self.assertIn("ARCHIVE", queue_lines[2])
                self.assertIn("ARCHIVER", queue_lines[2])
                self.assertIn("FEEDBACK", queue_lines[4])
                self.assertIn("SAMPLER", queue_lines[4])
                self.assertNotIn("pending", str(queues.render()).lower())
                self.assertEqual(queues.region.height, 8)
                calculators = app.screen.query_one("#ov-calcs")
                workers = app.screen.query_one("#ov-workers")
                self.assertEqual(calculators.region.width, samples.region.width)
                self.assertEqual(
                    calculators.content_region.width,
                    samples.content_region.width,
                )
                self.assertEqual(health.region.width, queues.region.width)
                self.assertEqual(queues.region.width, workers.region.width)
                self.assertEqual(calculators.region.x, samples.region.x)
                self.assertEqual(health.region.x, queues.region.x)
                self.assertEqual(queues.region.x, workers.region.x)
                self.assertLess(samples.region.x, health.region.x)
                sparks = app.screen.query_one("#ov-sparks")
                overview = app.screen.query_one("#overview")
                self.assertEqual(
                    overview.spark_history_widths(),
                    (
                        len(overview.frame.spark_samples),
                        len(overview.frame.spark_queue),
                    ),
                )
                self.assertEqual(sparks.styles.padding.bottom, 0)
                self.assertEqual(sparks.styles.padding.right, 1)
                self.assertEqual(sparks.styles.border_subtitle_align, "right")
                self.assertEqual(sparks.styles.border_subtitle_color.hex.lower(), "#8d93a1")
                self.assertIn("Peak : samples", sparks.border_subtitle)
                self.assertIn("queue", sparks.border_subtitle)
                app.screen.query_one("#topbar")
                app.screen.query_one("#tab-bar")
                app.screen.query_one("#hint")
                from textual.widgets import ContentSwitcher

                tab_text = str(app.screen.query_one("#tab-bar").render())
                hint_text = str(app.screen.query_one("#hint").render())
                self.assertIn("Overview", tab_text)
                self.assertNotIn("1 Overview", tab_text)
                for hidden_name in (
                    "Workers",
                    "Factory",
                    "Sampler",
                    "Calculators",
                    "Samples",
                    "Host",
                ):
                    self.assertNotIn(hidden_name, tab_text)
                self.assertIn("R: refresh", hint_text)
                self.assertNotIn("1-7", hint_text)
                self.assertNotIn("page", hint_text.lower())
                self.assertEqual(
                    app.screen.query_one("#pages", ContentSwitcher).current,
                    "overview",
                )
                for hidden_id in (
                    "#workers",
                    "#factory",
                    "#sampler",
                    "#calculators",
                    "#samples",
                    "#host",
                ):
                    self.assertEqual(len(app.screen.query(hidden_id)), 0)
                await pilot.press("tab", "left", "right", "2", "7")
                await pilot.pause()
                self.assertEqual(
                    app.screen.query_one("#pages", ContentSwitcher).current,
                    "overview",
                )
                await pilot.press("q")

        asyncio.run(_run())

    def test_simu_workspace_never_holds_a_collector(self) -> None:
        try:
            from jarvishep2.monitor.workspace import WorkspaceScreen
        except ImportError:
            self.skipTest("textual is not installed")
        queue = make_fakeredis_queue()
        screen = WorkspaceScreen(simulated_choice(), collector=Collector(queue, pid=1))
        self.assertIsNone(screen.collector)
        self.assertIsNotNone(screen.engine)

    def test_health_info_trigger_has_hover_hint_and_centered_detail_modal(self) -> None:
        try:
            from jarvishep2.monitor.app import MonitorApp
        except ImportError:
            self.skipTest("textual is not installed")

        async def _run() -> None:
            app = MonitorApp(scan_lister=lambda: [])
            async with app.run_test(size=(120, 64)) as pilot:
                await pilot.pause()
                await pilot.press("s")
                await pilot.pause()
                row = app.screen.query_one("#ov-health-core")
                # The marker lives at the leading edge, regardless of refreshes.
                await pilot.hover("#ov-health-core", offset=(0, 0))
                await pilot.pause()
                self.assertIn("health-info-hover", row.classes)
                self.assertEqual(row.tooltip, "💡 Click to explain CORE")
                # A neighboring blank cell is deliberately not a click target.
                await pilot.click("#ov-health-core", offset=(19, 0))
                await pilot.pause()
                self.assertEqual(
                    app.screen.__class__.__name__,
                    "OverviewReleaseScreen",
                )
                await pilot.click("#ov-health-core", offset=(0, 0))
                await pilot.pause()
                body = str(app.screen.query_one("#health-detail-body").render())
                self.assertIn("Jarvis2Core", body)
                self.assertIn("HEALTHY MEANS", body)
                self.assertIn("NEEDS ATTENTION WHEN", body)
                await pilot.press("escape")
                await pilot.pause()
                self.assertEqual(
                    app.screen.__class__.__name__,
                    "OverviewReleaseScreen",
                )
                await pilot.press("q")

        asyncio.run(_run())

    def test_health_rows_reflow_after_terminal_resize(self) -> None:
        try:
            from jarvishep2.monitor.app import MonitorApp
        except ImportError:
            self.skipTest("textual is not installed")

        async def _run() -> None:
            app = MonitorApp(scan_lister=lambda: [])
            async with app.run_test(size=(120, 64)) as pilot:
                await pilot.pause()
                await pilot.press("s")
                await pilot.pause()
                await pilot.resize_terminal(80, 64)
                await pilot.pause()
                for row_id, suffixes in (
                    ("#ov-health-core", ("hb 0.3s", "0.3s")),
                    ("#ov-health-factory", ("dispatching", "run")),
                    ("#ov-health-workers", ("48/48 alive", "48/48")),
                    ("#ov-health-archiver", ("queue",)),
                ):
                    row = app.screen.query_one(row_id)
                    rendered = str(row.render())
                    self.assertEqual(visible_width(rendered), row.content_region.width)
                    if suffixes == ("queue",):
                        self.assertIn("queue ", rendered)
                    else:
                        self.assertTrue(rendered.rstrip().endswith(suffixes))
                await pilot.press("q")

        asyncio.run(_run())

    def test_lower_blocks_use_two_columns_without_blank_rows(self) -> None:
        try:
            from jarvishep2.monitor.app import MonitorApp
        except ImportError:
            self.skipTest("textual is not installed")

        async def _run() -> None:
            app = MonitorApp(scan_lister=lambda: [])
            async with app.run_test(size=(120, 64)) as pilot:
                await pilot.pause()
                await pilot.press("s")
                await pilot.pause()
                sparks = app.screen.query_one("#ov-sparks")
                samples = app.screen.query_one("#ov-samples")
                health = app.screen.query_one("#ov-health")
                self.assertEqual(samples.region.y, sparks.region.y + sparks.region.height)
                self.assertEqual(health.region.y, samples.region.y)
                self.assertGreater(health.region.x, samples.region.x)
                calculators = app.screen.query_one("#ov-calcs")
                self.assertEqual(calculators.region.x, samples.region.x)
                self.assertEqual(
                    calculators.region.y,
                    samples.region.y + samples.region.height,
                )
                queues = app.screen.query_one("#ov-queues")
                self.assertEqual(queues.region.x, health.region.x)
                self.assertEqual(queues.region.y, health.region.y + health.region.height)
                workers = app.screen.query_one("#ov-workers")
                self.assertEqual(workers.region.x, health.region.x)
                self.assertEqual(
                    workers.region.y,
                    queues.region.y + queues.region.height,
                )
                pacman = app.screen.query_one("#ov-pacman")
                lower_right = app.screen.query_one("#ov-lower-right")
                self.assertEqual(pacman.region.x, health.region.x)
                self.assertEqual(
                    pacman.region.y,
                    workers.region.y + workers.region.height,
                )
                self.assertEqual(
                    pacman.region.y + pacman.region.height,
                    lower_right.region.y + lower_right.region.height,
                )
                pacman_render = str(pacman.render())
                self.assertIn("·", pacman_render)
                self.assertIn("👻", pacman_render)
                self.assertTrue(
                    any(glyph in pacman_render for glyph in ("ᗧ", "◯"))
                )
                await pilot.press("q")

        asyncio.run(_run())

    def test_live_workspace_writes_want_only_on_calculators_and_samples(self) -> None:
        try:
            from jarvishep2.monitor.app import MonitorApp
            from jarvishep2.monitor.workspace import WorkspaceScreen
        except ImportError:
            self.skipTest("textual is not installed")

        queue = make_fakeredis_queue()
        collector = Collector(queue, pid=9)
        choice = ScanChoice(
            reference="R1",
            name="live-scan",
            control_pid=44001,
            process_count=1,
            pids=(44001,),
            simulated=False,
        )

        async def _run() -> None:
            app = MonitorApp(scan_lister=lambda: [])
            async with app.run_test() as pilot:
                await pilot.pause()
                app.push_screen(WorkspaceScreen(choice, collector=collector))
                await pilot.pause()
                self.assertEqual(int(queue.r.exists(MONITOR_WANT) or 0), 0)
                await pilot.press("5")
                await pilot.pause()
                self.assertEqual(queue.r.hget(MONITOR_WANT, "ch:calc"), "1")
                await pilot.press("1")
                await pilot.pause()
                self.assertEqual(int(queue.r.exists(MONITOR_WANT) or 0), 0)
                await pilot.press("6")
                await pilot.pause()
                self.assertEqual(queue.r.hget(MONITOR_WANT, "ch:sample"), "1")
                await pilot.press("q")

        asyncio.run(_run())
        self.assertEqual(int(queue.r.exists(MONITOR_WANT) or 0), 0)

    def test_live_inflight_pages_render_sidecars_without_task_payload_reads(self) -> None:
        try:
            from jarvishep2.monitor.app import MonitorApp
            from jarvishep2.monitor.workspace import WorkspaceScreen
            from textual.widgets import Static
        except ImportError:
            self.skipTest("textual is not installed")

        queue = make_fakeredis_queue()
        queue.register_calc_pool("DemoCalc", 1)
        queue.publish_proc_board(
            "worker",
            owner_id="3",
            pid=44003,
            status="busy",
            current_sample="sample-uuid-3",
            held_calc_n=1,
        )
        queue.r.hset("hep:sample:stats", mapping={"completed": 4, "running": 1, "failed": 0})
        collector = Collector(
            queue,
            owner_ids=["3"],
            pid=9,
            pool_names=["DemoCalc"],
            slots={"DemoCalc": 1},
        )
        choice = ScanChoice(
            reference="R1",
            name="live-scan",
            control_pid=44001,
            process_count=1,
            pids=(44001,),
            simulated=False,
        )

        async def _run() -> None:
            app = MonitorApp(scan_lister=lambda: [])
            async with app.run_test() as pilot:
                await pilot.pause()
                app.push_screen(WorkspaceScreen(choice, collector=collector))
                await pilot.pause()
                await pilot.press("5")
                await pilot.pause()
                queue._acquire_calc("DemoCalc", timeout=1, worker_id="3")
                await pilot.press("r")
                await pilot.pause()
                owners = str(app.screen.query_one("#calculators-detail", Static).render())
                self.assertIn("worker-03", owners)
                self.assertIn("hep:monitor:calc:busy:DemoCalc", owners)

                await pilot.press("6")
                await pilot.pause()
                queue.publish_sample_overlay(
                    "3", {"uuid": "sample-uuid-3", "step": "DemoCalc", "t0": 12}
                )
                await pilot.press("r")
                await pilot.pause()
                rows = str(app.screen.query_one("#samples-list", Static).render())
                detail = str(app.screen.query_one("#samples-detail", Static).render())
                self.assertIn("sample-uuid-3", rows)
                self.assertIn("DemoCalc", detail)
                self.assertIn("held_calc      1", detail)
                await pilot.press("q")

        asyncio.run(_run())
        self.assertEqual(int(queue.r.exists(MONITOR_WANT) or 0), 0)

    def test_live_standard_pages_render_collector_snapshot(self) -> None:
        try:
            from jarvishep2.monitor.app import MonitorApp
            from jarvishep2.monitor.workspace import WorkspaceScreen
            from textual.widgets import Static
        except ImportError:
            self.skipTest("textual is not installed")

        queue = make_fakeredis_queue()
        queue.publish_proc_board(
            "core",
            role="core",
            pid=44001,
            status="running",
            sampler_status=json.dumps(
                {
                    "schema": 1,
                    "method": "Random",
                    "state": "running",
                    "progress": {"kind": "finite", "current": 23, "target": 100},
                    "metrics": {"accepted": 19, "seed": 7},
                }
            ),
        )
        queue.publish_proc_board("archiver", role="archiver", pid=44002, status="idle")
        queue.publish_proc_board("redis", role="redis", pid=44003, status="ready")
        queue.publish_proc_board(
            "worker",
            owner_id="3",
            pid=44004,
            status="busy",
            current_sample="active-sample",
            held_calc_n=1,
        )
        queue.push_task({"uuid": "queued-sample", "u_coords": [0.0]})
        collector = Collector(
            queue,
            owner_ids=["3"],
            pid=9,
            host_snapshotter=lambda: {
                "available": True,
                "cpu_percent": 25.0,
                "memory_used": 8 * 1024**3,
                "memory_total": 32 * 1024**3,
                "swap_used": 0,
                "swap_total": 8 * 1024**3,
                "load": (1.0, 0.8, 0.6),
                "processes": [
                    {
                        "role": "worker",
                        "pid": 44004,
                        "alive": True,
                        "cpu_percent": 1.0,
                        "rss": 256 * 1024**2,
                        "fds": 8,
                        "ppid": 44001,
                        "threads": 2,
                        "cmdline": "Jarvis-Worker-03:live-scan",
                    }
                ],
            },
            sampler_metadata={
                "method": "Random",
                "family": "simple",
                "dimensions": 2,
                "config": {"point_number": 100, "seed": 7},
            },
        )
        choice = ScanChoice(
            reference="R1",
            name="live-scan",
            control_pid=44001,
            process_count=1,
            pids=(44001,),
            simulated=False,
        )

        async def _run() -> None:
            app = MonitorApp(scan_lister=lambda: [])
            async with app.run_test() as pilot:
                await pilot.pause()
                app.push_screen(WorkspaceScreen(choice, collector=collector))
                await pilot.pause()
                await pilot.press("r")
                await pilot.pause()
                overview = str(app.screen.query_one("#live-overview-body", Static).render())
                self.assertIn("live scan snapshot", overview)
                self.assertIn("CONTROL PLANE", overview)
                self.assertNotIn("simu-iDM", overview)
                samples_block = app.screen.query_one("#ov-samples")
                self.assertTrue(str(samples_block.border_title).endswith("Random"))
                self.assertIn("CANDIDATES 23 / 100", str(samples_block.render()))
                self.assertIn("Acc ", str(samples_block.render()))
                self.assertIn("SEED 7", str(samples_block.border_subtitle))
                queues_block = app.screen.query_one("#ov-queues")
                self.assertTrue(str(queues_block.border_title).endswith("FLOW"))
                self.assertIn("TASK", str(queues_block.render()))
                self.assertIn("FEEDBACK", str(queues_block.render()))
                self.assertNotIn("pending", str(queues_block.render()).lower())
                await pilot.press("2")
                await pilot.pause()
                workers = str(app.screen.query_one("#workers-list", Static).render())
                self.assertIn("active-sample", workers)

                await pilot.press("3")
                await pilot.pause()
                factory = str(app.screen.query_one("#factory-list", Static).render())
                self.assertIn("running", factory)
                self.assertIn("ready", factory)

                await pilot.press("4")
                await pilot.pause()
                sampler = str(app.screen.query_one("#sampler-list", Static).render())
                self.assertIn("task", sampler)
                self.assertIn("1", sampler)
                await pilot.press("7")
                await pilot.pause()
                host = str(app.screen.query_one("#host-list", Static).render())
                detail = str(app.screen.query_one("#host-detail", Static).render())
                self.assertIn("25.0%", host)
                self.assertIn("Jarvis-Worker-03", detail)
                await pilot.press("q")

        asyncio.run(_run())
        self.assertEqual(int(queue.r.exists(MONITOR_WANT) or 0), 0)


if __name__ == "__main__":
    unittest.main()
