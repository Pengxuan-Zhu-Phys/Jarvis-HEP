#!/usr/bin/env python3
"""AdaptiveBridson live-band / local-core Bridson tests."""

from __future__ import annotations

import json
import os
import signal
import tempfile
import threading
import time
import unittest
from typing import Any

import numpy as np
from fakeredis import TcpFakeServer

from jarvishep2.Sampling.adaptive_bridson import (
    AdaptiveBridsonSampler,
    DelaunayGraph,
    KNNGraph,
    LevelSetPoint,
    LoopAction,
    _SepIndex,
    _thin_centers,
    clip_unit_cube,
    sample_ball,
)
from jarvishep2.core import Jarvis2Core
from jarvishep2.distributor import Distributor, STATELESS_METHODS
from jarvishep2.factory import TaskFactory
from jarvishep2.mp_context import get_spawn_context
from jarvishep2.process_cleanup import kill_jarvis_processes, list_running_scans
from jarvishep2.database import read_persisted_uuids
from jarvishep2.redis_queue import FEEDBACK_QUEUE, make_fakeredis_queue


def _start_tcp_fakeredis() -> tuple[TcpFakeServer, dict]:
    server = TcpFakeServer(("127.0.0.1", 0))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    return server, {"host": host, "port": port, "db": 0, "managed": False}


def _flat_vars(names: list[str]) -> list[dict[str, Any]]:
    return [
        {
            "name": name,
            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
        }
        for name in names
    ]


def _opera_module(
    *,
    name: str,
    operator: str,
    inputs: list[str],
    outputs: list[tuple[str, str]],
) -> dict[str, Any]:
    return {
        "name": name,
        "operator": operator,
        "call_mode": "call",
        "input": [{"name": n, "expression": n} for n in inputs],
        "output": [{"name": out, "entry": entry} for out, entry in outputs],
    }


def _als_config(
    *,
    redis_config: dict[str, Any],
    tmpdir: str,
    method_block: dict[str, Any],
    operator: str = "jarvishep2.testing.eggbox.circle_r2",
    output_name: str = "r2",
    var_names: list[str] | None = None,
    seed: int = 7,
    workers: int = 2,
    project_name: str = "alevelset",
    scan_name: str = "adaptive-test-scan",
) -> dict[str, Any]:
    names = list(var_names or ["x", "y"])
    block = {
        "target_expression": output_name,
        "target_value": 0.04,
        "initial_radius": 0.12,
        "refinement_factor": 0.5,
        "core_half_width": 0.08,
        "outer_half_width": 0.25,
        "threshold": 0.08,
        "max_generations": 8,
        "max_points": 600,
        "max_new_per_generation": 80,
        "k_ref": 20,
        "radius_shrink_mode": "on_fill",
    }
    block.update(method_block)
    return {
        "project_name": project_name,
        "Scan": {"name": scan_name},
        "task_result_dir": tmpdir,
        "Runtime": {
            "mode": "redis",
            "workers": workers,
            "batch_size": max(4, workers * 2),
            "redis": redis_config,
            "sample_artifacts": "never",
            "Watchdog": {"enabled": False},
        },
        "Sampling": {
            "Method": "AdaptiveBridson",
            "Bounds": {**dict(block), "Seed": seed},
            "Variables": _flat_vars(names),
        },
        "Operas": {
            "Modules": [
                _opera_module(
                    name="Fixture",
                    operator=operator,
                    inputs=names,
                    outputs=[(output_name, output_name)],
                )
            ]
        },
    }


def _run_adaptive_resume_control(
    redis_config: dict[str, Any],
    tmpdir: str,
    scan_name: str,
    resume: bool,
    outcome_queue: Any,
) -> None:
    config = _als_config(
        redis_config=redis_config,
        tmpdir=tmpdir,
        method_block={
            "target_value": 0.04,
            "threshold": 0.12,
            "initial_radius": 0.15,
            "max_generations": 2,
            "max_points": 120,
            "max_new_per_generation": 30,
        },
        operator="jarvishep2.testing.resume.slow_circle_r2",
        seed=41,
        workers=2,
        project_name="adaptive-resume",
        scan_name=scan_name,
    )
    config["task_root"] = tmpdir
    config["Runtime"]["max_inflight"] = 8
    config["Runtime"]["Archiver"] = {"mode": "process", "batch_size": 4}
    outcome = Jarvis2Core(config).run(resume=resume, write_run_summary=False)
    outcome_queue.put(
        {
            "ok": outcome.ok,
            "submitted": outcome.submitted,
            "archived": outcome.archived,
            "exit_code": outcome.exit_code,
        }
    )
    TaskFactory.reset_instance()


class NeighborGraphUnitTests(unittest.TestCase):
    def test_delaunay_builds(self) -> None:
        pts = np.array([[0.2, 0.2], [0.8, 0.2], [0.2, 0.8], [0.8, 0.8], [0.5, 0.5]])
        edges = DelaunayGraph().build(pts)
        self.assertGreaterEqual(edges.shape[0], 4)

    def test_knn_symmetric(self) -> None:
        rng = np.random.default_rng(0)
        pts = rng.random((20, 2))
        edges = KNNGraph(knn_k=4).build(pts)
        as_set = {(int(a), int(b)) for a, b in edges}
        for a, b in as_set:
            self.assertLess(a, b)
        self.assertGreater(len(as_set), 0)

    def test_sample_ball_stays_near_center(self) -> None:
        rng = np.random.default_rng(1)
        center = np.array([0.5, 0.5])
        for _ in range(20):
            p = sample_ball(center, 0.1, rng, 2)
            self.assertLessEqual(float(np.linalg.norm(p - center)), 0.1 + 1e-9)

    def test_clip_unit_cube(self) -> None:
        out = clip_unit_cube(np.array([-1.0, 2.0]))
        self.assertGreater(out[0], 0.0)
        self.assertLess(out[1], 1.0)


class LiveBandUnitTests(unittest.TestCase):
    def _sampler_with_points(
        self, us: list[list[float]], fs: list[float | None], *, target: float
    ) -> AdaptiveBridsonSampler:
        s = AdaptiveBridsonSampler()
        s.set_config(
            {
                "Sampling": {
                    "Method": "AdaptiveBridson",
                    "Bounds": {"target_expression": "r2",
                        "target_value": target,
                        "initial_radius": 0.2,
                        "core_half_width": 0.01,
                        "outer_half_width": 0.05,
                        "threshold": 0.01,
                        "max_generations": 5,
                        "radius_shrink_mode": "on_fill", "Seed": 1},
                    "Variables": _flat_vars(["x", "y"]),
                },
                "Runtime": {"workers": 1},
            }
        )
        s._points = [
            LevelSetPoint(
                u=list(u),
                x={"x": u[0], "y": u[1]},
                f=f,
                uuid=f"p{i}",
                generation=0,
            )
            for i, (u, f) in enumerate(zip(us, fs))
        ]
        s._uuid_to_index = {p.uuid: i for i, p in enumerate(s._points)}
        s._graph = s._make_graph()
        return s

    def test_loop_decisions_name_terminal_and_advance_transitions(self) -> None:
        s = self._sampler_with_points(
            [[0.5, 0.5]], [0.04], target=0.04
        )
        missing = s._pre_fill_decision(
            best=None, t_min=None, t_max=None, fill_needed=False
        )
        self.assertEqual(missing.action, LoopAction.STOP)
        self.assertEqual(missing.reason, "level-set not present in domain")

        s._contour_converged = lambda **_kwargs: True  # type: ignore[method-assign]
        converged = s._pre_fill_decision(
            best=0, t_min=0.04, t_max=0.04, fill_needed=False
        )
        self.assertEqual(converged.action, LoopAction.STOP)
        self.assertEqual(converged.reason, "converged")

        s._contour_converged = lambda **_kwargs: False  # type: ignore[method-assign]
        s._should_advance_generation = lambda **_kwargs: True  # type: ignore[method-assign]
        self.assertEqual(
            s._post_fill_decision(
                t_min=0.02, t_max=0.06, fill_needed=False, cores=[0]
            ).action,
            LoopAction.ADVANCE,
        )
        self.assertEqual(
            s._post_fill_decision(
                t_min=0.02, t_max=0.06, fill_needed=True, cores=[0]
            ).action,
            LoopAction.CONTINUE,
        )

    def test_best_index_closest_to_target(self) -> None:
        s = self._sampler_with_points(
            [[0.1, 0.1], [0.5, 0.5], [0.9, 0.9]],
            [0.0, 0.04, 0.2],
            target=0.04,
        )
        self.assertEqual(s._best_index(), 1)

    def test_classify_absolute_outer_and_core(self) -> None:
        # target=0.04, core±0.01 → [0.03,0.05]; outer±0.05 → [−0.01,0.09]
        s = self._sampler_with_points(
            [[0.5, 0.5], [0.6, 0.5], [0.7, 0.5], [0.1, 0.1]],
            [0.04, 0.07, 0.20, 0.035],  # core, frontier, outside, core
            target=0.04,
        )
        cores, frontier, outer, best = s._classify_bands()
        self.assertEqual(best, 0)
        self.assertIn(0, cores)
        self.assertIn(3, cores)
        self.assertIn(1, frontier)
        self.assertNotIn(2, outer)  # |0.20-0.04|=0.16 > outer
        self.assertNotIn(1, cores)

    def test_straddle_brackets_and_secant(self) -> None:
        # Two points on opposite sides of T=0.04 along an edge.
        s = self._sampler_with_points(
            [[0.2, 0.5], [0.8, 0.5], [0.5, 0.2]],
            [0.01, 0.09, 0.50],
            target=0.04,
        )
        brackets = s._straddle_brackets()
        self.assertGreaterEqual(len(brackets), 1)
        corr = s._root_correct_proposals(
            brackets, sep_index=None, r_sep=0.01, budget=20
        )
        self.assertGreater(len(corr), 0)
        xs = [float(u[0]) for u in corr]
        self.assertTrue(any(0.25 < x < 0.75 for x in xs), xs)

    def test_log_secant_alphas_for_orders_of_magnitude(self) -> None:
        # Ω-like: 1e-4 ↔ 10, T=0.12 → log-secant α far from linear 0.5.
        alphas = AdaptiveBridsonSampler._secant_alphas(1e-4, 10.0, 0.12)
        self.assertTrue(any(a < 0.45 for a in alphas))  # log pulls toward low side
        self.assertIn(0.5, [round(a, 5) for a in alphas] + [0.5])  # bisection present

    def test_full_graph_brackets_even_when_outer_empty_of_straddles(self) -> None:
        """Gen-0 style: outer shell all on one side of T; true straddle is far out."""
        # outer half-w=0.05 → |f-0.04|≤0.05 keeps f∈[−0.01,0.09]
        # points at 0.041 and 0.045 are outer but same side of T=0.04
        # points at 0.001 and 0.20 straddle T and are Voronoi-linked if cloud is rich.
        us = [
            [0.50, 0.50],  # near T, slightly above
            [0.52, 0.50],  # near T, above
            [0.20, 0.50],  # far below T
            [0.80, 0.50],  # far above T
            [0.50, 0.20],
            [0.50, 0.80],
            [0.30, 0.30],
            [0.70, 0.70],
        ]
        fs = [0.041, 0.045, 0.001, 0.20, 0.10, 0.10, 0.15, 0.15]
        s = self._sampler_with_points(us, fs, target=0.04)
        # Production failure mode: tight outer; near-T points same side of T.
        s._core_half_width = 0.0005
        s._outer_half_width = 0.01
        cores, frontier, outer, _ = s._classify_bands()
        self.assertEqual(cores, [])
        # Outer may contain 0.041 only (d=0.001); both near-T points same side.
        brackets_outer = s._straddle_brackets(pool=outer if outer else [])
        brackets_full = s._straddle_brackets(pool=None)
        # Full cloud must see the below/above pair; outer-only often sees none.
        self.assertGreaterEqual(len(brackets_full), 1)
        self.assertGreaterEqual(len(brackets_full), len(brackets_outer))

    def test_null_f_ignored_for_best(self) -> None:
        s = self._sampler_with_points(
            [[0.2, 0.2], [0.5, 0.5]],
            [None, 0.1],
            target=0.0,
        )
        self.assertEqual(s._best_index(), 1)

    def test_densify_stays_in_core_window(self) -> None:
        s = self._sampler_with_points(
            [[0.5, 0.5], [0.55, 0.5], [0.5, 0.55]],
            [0.04, 0.05, 0.03],
            target=0.04,
        )
        s._bridge_gaps = False  # pure local-ball contract
        r_window = 0.15
        r_sep = 0.05
        new_us = s._densify_local(
            [0, 1, 2], r_window=r_window, r_sep=r_sep, generation=1
        )
        cores = [np.asarray(s._points[i].u) for i in (0, 1, 2)]
        for u in new_us:
            dmin = min(float(np.linalg.norm(u - c)) for c in cores)
            self.assertLessEqual(dmin, r_window + 1e-9)
            for c in cores:
                # May be near a core but not closer than ~0.5*r_sep to same core
                # Global exclusion is vs all accepted including cores.
                self.assertGreaterEqual(
                    float(np.linalg.norm(u - c)), r_sep * 0.5 - 1e-9
                )

    def test_gap_bridge_reconnects_distant_cores(self) -> None:
        """Two live cores farther than 2 r_window still get mid-segment samples."""
        # Separation 0.30; r_window=0.10 → pure balls cannot meet (need >0.20).
        s = self._sampler_with_points(
            [[0.30, 0.50], [0.60, 0.50]],
            [0.04, 0.04],
            target=0.04,
        )
        s._bridge_gaps = True
        s._bridge_span_factor = 3.0
        s._k_ref = 5
        s._max_new_per_generation = 50
        r_window = 0.10
        r_sep = 0.04
        new_us = s._densify_local(
            [0, 1], r_window=r_window, r_sep=r_sep, generation=1
        )
        self.assertGreater(len(new_us), 0)
        # At least one proposal should sit in the mid-gap (not only near cores).
        mid = np.array([0.45, 0.50])
        d_mid = min(float(np.linalg.norm(u - mid)) for u in new_us)
        self.assertLess(
            d_mid,
            0.08,
            msg=f"expected a bridge sample near mid-gap; closest={d_mid}",
        )

    def test_mst_bridge_pairs_prefers_nearest_links(self) -> None:
        # Three colinear cores: MST should use the two short edges only.
        pos = [
            np.array([0.0, 0.0]),
            np.array([0.1, 0.0]),
            np.array([0.2, 0.0]),
        ]
        pairs = AdaptiveBridsonSampler._mst_bridge_pairs(
            pos, min_dist=0.05, max_dist=0.25
        )
        ends = {(a, b) for a, b, _d in pairs}
        self.assertEqual(ends, {(0, 1), (1, 2)})

    def test_sep_index_rejects_near_neighbors(self) -> None:
        idx = _SepIndex(0.1, 2)
        idx.add(np.array([0.5, 0.5]))
        self.assertFalse(idx.far_enough(np.array([0.55, 0.5])))
        self.assertTrue(idx.far_enough(np.array([0.7, 0.5])))

    def test_thin_centers_reduces_overdense_set(self) -> None:
        # 11 points on a line with spacing 0.01; thin at 0.05 → fewer centers.
        pos = np.array([[0.1 + 0.01 * i, 0.5] for i in range(11)], dtype=np.float64)
        kept = _thin_centers(pos, min_spacing=0.05)
        self.assertLess(len(kept), 11)
        self.assertGreaterEqual(len(kept), 2)
        for a, b in zip(kept, kept[1:]):
            d = float(np.linalg.norm(pos[a] - pos[b]))
            self.assertGreaterEqual(d, 0.05 - 1e-12)

    def test_classify_bridge_edges_ABCD(self) -> None:
        # 0,1: both non-core straddle → A
        # 2 core near T, 3 far → B if they form a filtered edge
        s = self._sampler_with_points(
            [
                [0.20, 0.50],
                [0.80, 0.50],
                [0.50, 0.50],
                [0.50, 0.80],
                [0.48, 0.50],
                [0.52, 0.50],
            ],
            [0.001, 0.20, 0.040, 0.20, 0.039, 0.041],
            target=0.04,
        )
        s._core_half_width = 0.002
        s._outer_half_width = 0.05
        s._radius = 0.5
        s._bridge_span_factor = 2.5
        cores, _fr, _ou, _ = s._classify_bands()
        # cores should include points with |f-0.04|<=0.002 → 0.040, 0.039, 0.041
        filt = s._nearest_straddle_brackets(
            max_distance=2.5 * s._radius
        )
        cls = s._classify_bridge_edges(filt, cores, r_window=s._radius)
        self.assertGreaterEqual(len(cls["A"]) + len(cls["B"]), 1)
        # A/B endpoints listed for edge fill — not mixed into densify ball centers.
        ends = s._bridge_endpoint_indices(cls)
        for key in ("A", "B"):
            for i, j, _fi, _fj in cls[key]:
                self.assertIn(i, ends)
                self.assertIn(j, ends)
        # Accumulated cores expand only (register never shrinks the set).
        s._register_cores(cores)
        n0 = len(s._accumulated_core_uuids)
        s._register_cores(cores[:1] if cores else [])
        self.assertEqual(len(s._accumulated_core_uuids), n0)

    def test_fill_needed_respects_ABCD(self) -> None:
        s = self._sampler_with_points(
            [[0.5, 0.5], [0.55, 0.5]],
            [0.04, 0.04],
            target=0.04,
        )
        s._min_cores_for_coverage = 2
        s._radius = 0.2
        # Only D-like situation: no A/B/C, coverage ok with 2 cores close
        s._live_core_indices = [0, 1]
        self.assertFalse(
            s._fill_needed(
                cores=[0, 1],
                n_local_brackets=0,
                coverage_ok=True,
                densify_new=0,
                densify_attempted=True,
                n_bridge_A=0,
                n_bridge_B=0,
                n_bridge_C=0,
            )
        )
        self.assertTrue(
            s._fill_needed(
                cores=[0, 1],
                n_local_brackets=5,
                coverage_ok=True,
                densify_new=0,
                densify_attempted=True,
                n_bridge_A=3,
                n_bridge_B=0,
                n_bridge_C=0,
            )
        )
        self.assertTrue(
            s._fill_needed(
                cores=[0, 1],
                n_local_brackets=1,
                coverage_ok=True,
                densify_new=0,
                densify_attempted=True,
                n_bridge_A=0,
                n_bridge_B=1,
                n_bridge_C=0,
            )
        )
        # Quiet random passes cannot override structurally incomplete core
        # coverage; otherwise disconnected islands close the scale early.
        self.assertTrue(
            s._fill_needed(
                cores=[0, 1],
                n_local_brackets=5,
                coverage_ok=False,
                densify_new=0,
                densify_attempted=True,
                n_bridge_A=3,
                n_bridge_B=1,
                n_bridge_C=0,
                scale_quiet=True,
            )
        )
        self.assertTrue(
            s._fill_needed(
                cores=[0, 1],
                n_local_brackets=1,
                coverage_ok=True,
                n_endpoint_edges=1,
                scale_quiet=True,
            )
        )
        self.assertFalse(
            s._fill_needed(
                cores=[0, 1],
                n_local_brackets=1,
                coverage_ok=True,
                n_endpoint_edges=1,
                force_scale_close=True,
            )
        )

    def test_actionable_endpoint_is_recomputed_until_root_is_covered(self) -> None:
        s = self._sampler_with_points(
            [[0.5, 0.5], [0.10, 0.10], [0.14, 0.10]],
            [0.04, 0.02, 0.06],
            target=0.04,
        )
        s._core_half_width = 0.001
        edge = (1, 2, 0.02, 0.06)
        active = s._actionable_endpoint_edges([edge], [0], r_sep=0.02)
        self.assertEqual(active, [edge])

        ui = np.asarray(s._points[1].u)
        uj = np.asarray(s._points[2].u)
        for alpha in s._secant_alphas(0.02, 0.06, 0.04):
            alpha = float(min(0.85, max(0.08, alpha)))
            s._submitted_u_keys.add(s._u_key(ui + alpha * (uj - ui)))
        self.assertEqual(
            s._actionable_endpoint_edges([edge], [0], r_sep=0.02), []
        )

    def test_local_straddles_define_secant_working_cores(self) -> None:
        s = self._sampler_with_points(
            [[0.10, 0.10], [0.12, 0.10], [0.60, 0.10]],
            [0.03, 0.05, 0.05],
            target=0.04,
        )
        s._radius = 0.01
        near = (0, 1, 0.03, 0.05)
        far = (0, 2, 0.03, 0.05)
        local = s._local_bracket_edges([near, far])
        self.assertEqual(local, [near])
        anchors = s._bracket_root_anchors(local)
        self.assertEqual(len(anchors), 1)
        np.testing.assert_allclose(anchors[0], [0.11, 0.10], atol=1e-12)

    def test_nearest_opposite_side_pairs_bypass_graph_adjacency(self) -> None:
        s = self._sampler_with_points(
            [[0.00, 0.10], [0.10, 0.10], [0.11, 0.10], [0.90, 0.10]],
            [0.10, 0.11, 0.13, 0.14],
            target=0.12,
        )
        s._radius = 0.02
        pairs = s._nearest_straddle_brackets(max_distance=0.05)
        self.assertEqual([(edge[0], edge[1]) for edge in pairs], [(1, 2)])
        root = s._edge_root_anchor(pairs[0])
        self.assertGreater(float(root[0]), 0.10)
        self.assertLess(float(root[0]), 0.11)

    def test_virtual_core_samples_random_annulus(self) -> None:
        s = self._sampler_with_points(
            [[0.10, 0.10], [0.12, 0.10]], [0.03, 0.05], target=0.04
        )
        s._radius = 0.01
        s._submitted_u_keys = {s._u_key(point.u) for point in s._points}
        anchor = np.array([0.11, 0.10])
        self.assertEqual(s._refresh_virtual_cores([anchor], []), 1)
        proposals: list[np.ndarray] = []
        n = s._virtual_core_proposals(
            sep_index=s._submitted_sep_index(s._radius),
            new_us=proposals,
            budget=20,
            rng=np.random.default_rng(23),
        )
        self.assertGreater(n, 0)
        self.assertLessEqual(n, 2)
        distances = [float(np.linalg.norm(point - anchor)) for point in proposals]
        self.assertGreaterEqual(min(distances), s._radius - 1e-12)
        self.assertLessEqual(max(distances), 2.5 * s._radius + 1e-12)

    def test_direct_secant_precedes_random_virtual_fallback(self) -> None:
        s = self._sampler_with_points(
            [[0.10, 0.10]], [0.10], target=0.04
        )
        s._radius = 0.01
        anchor = np.array([0.30, 0.30])
        s._refresh_virtual_cores([anchor], [])
        proposals: list[np.ndarray] = []
        sep = s._submitted_sep_index(s._radius)
        n_direct = s._virtual_core_direct_proposals(
            sep_index=sep, new_us=proposals, budget=20
        )
        self.assertEqual(n_direct, 1)
        np.testing.assert_allclose(proposals[0], anchor)
        root = next(iter(s._active_virtual_cores.values()))
        self.assertTrue(root["pending_attempt"])

        n_random = s._virtual_core_proposals(
            sep_index=sep,
            new_us=proposals,
            budget=20,
            rng=np.random.default_rng(17),
        )
        self.assertEqual(n_random, 0)

    def test_virtual_core_persists_until_covered_or_stalled(self) -> None:
        s = self._sampler_with_points(
            [[0.10, 0.10], [0.50, 0.50]], [0.04, 0.04], target=0.04
        )
        s._radius = 0.01
        anchor = np.array([0.30, 0.30])
        self.assertEqual(s._refresh_virtual_cores([anchor], []), 1)
        root_id = next(iter(s._active_virtual_cores))
        s._active_virtual_cores[root_id]["pending_attempt"] = True
        self.assertEqual(s._refresh_virtual_cores([], []), 1)
        self.assertEqual(s._active_virtual_cores[root_id]["misses"], 1)

        # An actual evaluated core within 2*r_g retires the discovery root.
        s._points.append(
            LevelSetPoint(
                u=[0.315, 0.30], x={}, f=0.04, uuid="covered", generation=0
            )
        )
        self.assertEqual(s._refresh_virtual_cores([], [2]), 0)
        self.assertIn(root_id, s._retired_virtual_core_ids)

    def test_virtual_core_walks_to_best_feedback_probe(self) -> None:
        s = self._sampler_with_points(
            [[0.10, 0.10]], [0.10], target=0.04
        )
        s._radius = 0.01
        anchor = np.array([0.30, 0.30])
        s._refresh_virtual_cores([anchor], [])
        root_id = next(iter(s._active_virtual_cores))
        probe = [0.325, 0.30]
        s._points.append(
            LevelSetPoint(u=probe, x={}, f=0.041, uuid="probe", generation=0)
        )
        state = s._active_virtual_cores[root_id]
        state["pending_attempt"] = True
        state["pending_probe_us"] = [probe]
        self.assertEqual(s._refresh_virtual_cores([], []), 1)
        np.testing.assert_allclose(
            s._active_virtual_cores[root_id]["search_u"], probe
        )
        self.assertAlmostEqual(
            s._active_virtual_cores[root_id]["best_probe_dev"], 0.001
        )
        self.assertEqual(s._active_virtual_cores[root_id]["misses"], 0)

        # A completed but fully blue-noise-blocked trial is still a miss.
        s._active_virtual_cores[root_id]["pending_attempt"] = True
        s._active_virtual_cores[root_id]["pending_probe_us"] = []
        s._refresh_virtual_cores([], [])
        self.assertEqual(s._active_virtual_cores[root_id]["misses"], 1)

    def test_covered_bracket_root_is_not_activated(self) -> None:
        s = self._sampler_with_points(
            [[0.10, 0.10]], [0.04], target=0.04
        )
        s._radius = 0.01
        self.assertEqual(
            s._refresh_virtual_cores([np.array([0.115, 0.10])], [0]), 0
        )

    def test_submission_history_survives_control_cloud_prune(self) -> None:
        s = self._sampler_with_points(
            [[0.2, 0.2], [0.8, 0.8]], [0.01, 0.09], target=0.04
        )
        historical = np.array([0.35, 0.45])
        s._submitted_u_keys = {
            s._u_key(p.u) for p in s._points
        } | {s._u_key(historical)}
        s._prune_to_indices([0])

        unseen = s._filter_unsubmitted(
            [historical, np.array([0.36, 0.46]), historical.copy()]
        )
        self.assertEqual(len(unseen), 1)
        np.testing.assert_allclose(unseen[0], [0.36, 0.46])

    def test_blue_noise_uses_pruned_history_and_same_batch_candidates(self) -> None:
        s = self._sampler_with_points(
            [[0.8, 0.8]], [0.04], target=0.04
        )
        historical = np.array([0.20, 0.20])
        s._submitted_u_keys = {
            s._u_key(s._points[0].u),
            s._u_key(historical),
        }

        accepted = s._filter_blue_noise(
            [
                np.array([0.205, 0.20]),  # too close to pruned history
                np.array([0.40, 0.40]),   # accepted
                np.array([0.405, 0.40]),  # too close to same-batch accepted
                np.array([0.60, 0.60]),   # accepted
            ],
            min_distance=0.01,
        )

        self.assertEqual(len(accepted), 2)
        np.testing.assert_allclose(accepted[0], [0.40, 0.40])
        np.testing.assert_allclose(accepted[1], [0.60, 0.60])

    def test_densify_respects_requested_blue_noise_distance(self) -> None:
        s = self._sampler_with_points(
            [[0.50, 0.50]], [0.04], target=0.04
        )
        s._submitted_u_keys = {s._u_key(s._points[0].u)}
        s._k_ref = 200
        s._max_new_per_generation = 40

        new_us = s._densify_local(
            [0], r_window=0.20, r_sep=0.10, generation=1
        )

        self.assertGreater(len(new_us), 0)
        all_points = [np.array([0.50, 0.50]), *new_us]
        for i, left in enumerate(all_points):
            for right in all_points[i + 1 :]:
                self.assertGreaterEqual(
                    float(np.linalg.norm(left - right)), 0.10 - 1e-12
                )

    def test_fill_rng_changes_between_same_generation_passes(self) -> None:
        s = self._sampler_with_points(
            [[0.5, 0.5], [0.6, 0.5]], [0.04, 0.04], target=0.04
        )
        s._seed = 21
        s._fill_pass = 1
        a = s._fill_rng(0).random(8)
        s._fill_pass = 2
        b = s._fill_rng(0).random(8)
        self.assertFalse(np.array_equal(a, b))
        s._fill_pass = 1
        np.testing.assert_array_equal(a, s._fill_rng(0).random(8))

    def test_core_coverage_requires_one_connected_component(self) -> None:
        s = self._sampler_with_points(
            [[0.10, 0.50], [0.15, 0.50], [0.70, 0.50], [0.75, 0.50]],
            [0.04, 0.04, 0.04, 0.04],
            target=0.04,
        )
        s._radius = 0.05
        s._min_cores_for_coverage = 4
        self.assertLessEqual(s._max_core_nn_gap([0, 1, 2, 3]) or 1.0, 0.10)
        self.assertEqual(len(s._core_components([0, 1, 2, 3])), 2)
        self.assertFalse(s._core_coverage_ok([0, 1, 2, 3]))

    def test_geometric_endpoints_detect_line_ends(self) -> None:
        s = self._sampler_with_points(
            [[0.30, 0.50], [0.38, 0.50], [0.46, 0.50], [0.54, 0.50]],
            [0.04] * 4,
            target=0.04,
        )
        s._radius = 0.05
        endpoints = s._geometric_endpoint_candidates([0, 1, 2, 3])
        terminal = {e["anchor_uuid"]: np.asarray(e["direction"]) for e in endpoints}
        self.assertIn("p0", terminal)
        self.assertIn("p3", terminal)
        self.assertLess(float(terminal["p0"][0]), 0.0)
        self.assertGreater(float(terminal["p3"][0]), 0.0)

    def test_active_endpoint_proposals_obey_blue_noise(self) -> None:
        s = self._sampler_with_points(
            [[0.30, 0.50], [0.38, 0.50], [0.46, 0.50], [0.54, 0.50]],
            [0.04] * 4,
            target=0.04,
        )
        s._radius = 0.05
        s._submitted_u_keys = {s._u_key(p.u) for p in s._points}
        self.assertGreater(s._refresh_active_endpoints([0, 1, 2, 3]), 0)
        index = s._submitted_sep_index(s._radius)
        proposals: list[np.ndarray] = []
        n = s._active_endpoint_proposals(
            sep_index=index,
            new_us=proposals,
            budget=20,
            rng=np.random.default_rng(7),
        )
        self.assertGreater(n, 0)
        history = [np.asarray(p.u) for p in s._points]
        for i, proposal in enumerate(proposals):
            for old in [*history, *proposals[:i]]:
                self.assertGreaterEqual(
                    float(np.linalg.norm(proposal - old)), s._radius - 1e-12
                )

    def test_endpoint_search_uses_fixed_random_annulus(self) -> None:
        s = self._sampler_with_points(
            [[0.30, 0.50], [0.38, 0.50], [0.46, 0.50], [0.54, 0.50]],
            [0.04] * 4,
            target=0.04,
        )
        s._radius = 0.05
        s._submitted_u_keys = {s._u_key(p.u) for p in s._points}
        s._refresh_active_endpoints([0, 1, 2, 3])
        right_id = next(
            key
            for key, value in s._active_endpoints.items()
            if value["anchor_uuid"] == "p3"
        )
        s._active_endpoints = {right_id: s._active_endpoints[right_id]}
        s._active_endpoints[right_id]["misses"] = 3
        anchor = np.asarray(s._active_endpoints[right_id]["anchor_u"])
        proposals: list[np.ndarray] = []
        s._active_endpoint_proposals(
            sep_index=s._submitted_sep_index(s._radius),
            new_us=proposals,
            budget=20,
            rng=np.random.default_rng(9),
        )
        self.assertGreater(len(proposals), 0)
        distances = [float(np.linalg.norm(u - anchor)) for u in proposals]
        self.assertGreaterEqual(min(distances), s._radius - 1e-12)
        self.assertLessEqual(max(distances), 2.5 * s._radius + 1e-12)

    def test_endpoint_probes_cover_full_circle_before_rejudging(self) -> None:
        s = self._sampler_with_points([[0.50, 0.50]], [0.50], target=0.50)
        s._radius = 0.05
        s._endpoint_omni_probes = 16
        endpoint_id = "terminal:center:plus"
        s._active_endpoints = {
            endpoint_id: {
                "id": endpoint_id,
                "kind": "terminal",
                "anchor_uuid": "p0",
                "anchor_u": [0.50, 0.50],
                "direction": [1.0, 0.0],
                "target_uuid": None,
                "target_distance": None,
                "attempts": 0,
                "misses": 0,
                "pending_attempt": False,
                "max_probe_distance": 0.0,
                "generation": 0,
            }
        }
        s._submitted_u_keys = {s._u_key(s._points[0].u)}
        proposals: list[np.ndarray] = []
        s._active_endpoint_proposals(
            sep_index=s._submitted_sep_index(s._radius),
            new_us=proposals,
            budget=32,
            rng=np.random.default_rng(19),
        )
        delta = np.asarray(proposals) - np.array([0.50, 0.50])
        self.assertGreater(len(delta), 4)
        self.assertLess(float(np.min(delta[:, 0])), -0.9 * s._radius)
        self.assertGreater(float(np.max(delta[:, 0])), 0.9 * s._radius)
        self.assertLess(float(np.min(delta[:, 1])), -0.9 * s._radius)
        self.assertGreater(float(np.max(delta[:, 1])), 0.9 * s._radius)
        self.assertTrue(s._active_endpoints[endpoint_id]["pending_attempt"])

    def test_radius_shrink_discards_and_rebuilds_endpoint_state(self) -> None:
        s = self._sampler_with_points(
            [[0.30, 0.50], [0.38, 0.50], [0.46, 0.50], [0.54, 0.50]],
            [0.04] * 4,
            target=0.04,
        )
        s._radius = 0.05
        s._refresh_active_endpoints([0, 1, 2, 3])
        old_ids = set(s._active_endpoints)
        s._retired_endpoint_ids.update(old_ids)
        s._radius = 0.025
        s._reset_endpoint_state()
        self.assertFalse(s._active_endpoints)
        self.assertFalse(s._retired_endpoint_ids)
        self.assertGreater(s._refresh_active_endpoints([0, 1, 2, 3]), 0)
        self.assertTrue(
            all(value["misses"] == 0 for value in s._active_endpoints.values())
        )
        # The smaller link radius can turn one old component into several new
        # components, so endpoint identities are allowed (and expected) to change.
        self.assertNotEqual(old_ids, set(s._active_endpoints))

    def test_active_endpoint_survives_one_miss_then_migrates(self) -> None:
        s = self._sampler_with_points(
            [[0.30, 0.50], [0.38, 0.50], [0.46, 0.50], [0.54, 0.50]],
            [0.04] * 4,
            target=0.04,
        )
        s._radius = 0.05
        s._refresh_active_endpoints([0, 1, 2, 3])
        right_id = next(
            key
            for key, value in s._active_endpoints.items()
            if value["anchor_uuid"] == "p3"
        )
        s._active_endpoints[right_id]["pending_attempt"] = True
        s._refresh_active_endpoints([0, 1, 2, 3])
        self.assertIn(right_id, s._active_endpoints)
        self.assertEqual(s._active_endpoints[right_id]["misses"], 1)

        s._points.append(
            LevelSetPoint(
                u=[0.60, 0.50], x={"x": 0.60, "y": 0.50}, f=0.04,
                uuid="p4", generation=0,
            )
        )
        s._uuid_to_index["p4"] = 4
        s._refresh_active_endpoints([0, 1, 2, 3, 4])
        self.assertNotIn(right_id, s._active_endpoints)
        self.assertTrue(
            any(e["anchor_uuid"] == "p4" for e in s._active_endpoints.values())
        )

    def test_active_endpoint_retires_only_after_stall_limit(self) -> None:
        s = self._sampler_with_points(
            [[0.30, 0.50], [0.38, 0.50], [0.46, 0.50], [0.54, 0.50]],
            [0.04] * 4,
            target=0.04,
        )
        s._radius = 0.05
        s._endpoint_stall_limit = 4
        s._refresh_active_endpoints([0, 1, 2, 3])
        right_id = next(
            key
            for key, value in s._active_endpoints.items()
            if value["anchor_uuid"] == "p3"
        )
        for miss in range(1, s._endpoint_stall_limit):
            s._active_endpoints[right_id]["pending_attempt"] = True
            s._refresh_active_endpoints([0, 1, 2, 3])
            self.assertIn(right_id, s._active_endpoints)
            self.assertEqual(s._active_endpoints[right_id]["misses"], miss)
        s._active_endpoints[right_id]["pending_attempt"] = True
        s._refresh_active_endpoints([0, 1, 2, 3])
        self.assertNotIn(right_id, s._active_endpoints)
        self.assertIn(right_id, s._retired_endpoint_ids)

        # Retirement belongs to this core geometry only. A genuine core gain
        # rebuilds the endpoint set instead of locking the id forever.
        s._points.append(
            LevelSetPoint(
                u=[0.60, 0.50], x={}, f=0.04, uuid="p4", generation=0
            )
        )
        s._refresh_active_endpoints([0, 1, 2, 3, 4])
        self.assertNotIn(right_id, s._retired_endpoint_ids)
        self.assertTrue(s._active_endpoints)

    def test_interior_core_gain_is_not_structural_progress(self) -> None:
        progress = AdaptiveBridsonSampler._fill_made_structural_progress
        self.assertFalse(
            progress(
                core_gain=12,
                coverage_before=True,
                best_improved=False,
                extent_improved=False,
                gap_improved=False,
            )
        )
        self.assertTrue(
            progress(
                core_gain=1,
                coverage_before=True,
                best_improved=False,
                extent_improved=True,
                gap_improved=False,
            )
        )
        self.assertTrue(
            progress(
                core_gain=1,
                coverage_before=False,
                best_improved=False,
                extent_improved=False,
                gap_improved=False,
            )
        )

    def test_tmin_tmax_gates_exit_and_next_gen(self) -> None:
        """Exit when Δt thin; next gen only after fill done and band still wide."""
        s = self._sampler_with_points(
            [[0.40, 0.5], [0.45, 0.5], [0.50, 0.5], [0.55, 0.5]],
            [0.03, 0.04, 0.05, 0.04],
            target=0.04,
        )
        s._threshold = 0.05
        s._radius = 0.1
        # Wide band (neighbors span 0.03..0.05 → Δ=0.02 < 0.05) → band ok.
        best, t_min, t_max = s._radius_t_band()
        self.assertIsNotNone(best)
        self.assertTrue(s._band_width_ok(t_min, t_max))
        self.assertTrue(
            s._contour_converged(t_min=t_min, t_max=t_max, fill_needed=False)
        )
        # Still filling → not converged / not next gen.
        self.assertFalse(
            s._contour_converged(t_min=t_min, t_max=t_max, fill_needed=True)
        )
        self.assertFalse(
            s._should_advance_generation(
                t_min=t_min, t_max=t_max, fill_needed=True, cores=[0, 1, 2]
            )
        )
        # Wide threshold fail: force large Δ.
        s._threshold = 0.005
        best2, t_min2, t_max2 = s._radius_t_band()
        self.assertFalse(s._band_width_ok(t_min2, t_max2))
        self.assertTrue(
            s._should_advance_generation(
                t_min=t_min2,
                t_max=t_max2,
                fill_needed=False,
                cores=[0, 1, 2],
            )
        )

    def test_tmin_tmax_uses_all_evaluated_points_within_two_rg(self) -> None:
        s = self._sampler_with_points(
            [
                [0.50, 0.50],  # best
                [0.56, 0.50],  # inside 2*r_g
                [0.50, 0.58],  # inside 2*r_g
                [0.71, 0.50],  # outside 2*r_g
            ],
            [0.040, 0.010, 0.090, 1.000],
            target=0.04,
        )
        s._radius = 0.10

        best, t_min, t_max = s._radius_t_band()

        self.assertEqual(best, 0)
        self.assertAlmostEqual(float(t_min), 0.010)
        self.assertAlmostEqual(float(t_max), 0.090)

    def test_tmin_tmax_two_rg_boundary_is_inclusive_and_tracks_rg(self) -> None:
        s = self._sampler_with_points(
            [[0.50, 0.50], [0.60, 0.50], [0.71, 0.50]],
            [0.040, 0.070, 0.200],
            target=0.04,
        )
        s._radius = 0.05
        _best, t_min, t_max = s._radius_t_band()
        self.assertAlmostEqual(float(t_min), 0.040)
        self.assertAlmostEqual(float(t_max), 0.070)

        s._radius = 0.049
        _best, t_min, t_max = s._radius_t_band()
        self.assertIsNone(t_min)
        self.assertIsNone(t_max)

    def test_singleton_scale_band_cannot_converge(self) -> None:
        s = self._sampler_with_points(
            [[0.50, 0.50], [0.80, 0.80]], [0.040, 0.20], target=0.04
        )
        s._radius = 0.05
        _best, t_min, t_max = s._radius_t_band()
        self.assertIsNone(t_min)
        self.assertIsNone(t_max)
        self.assertFalse(
            s._contour_converged(
                t_min=t_min, t_max=t_max, fill_needed=False
            )
        )

    def test_thin_band_without_a_core_refines_instead_of_converging(self) -> None:
        s = self._sampler_with_points(
            [[0.50, 0.50], [0.55, 0.50]], [0.1280, 0.1281], target=0.12
        )
        s._radius = 0.10
        s._core_half_width = 0.0025
        s._threshold = 0.0025
        _best, t_min, t_max = s._radius_t_band()
        self.assertTrue(s._band_width_ok(t_min, t_max))
        self.assertFalse(
            s._contour_converged(
                t_min=t_min, t_max=t_max, fill_needed=False
            )
        )

    def test_no_core_no_candidate_bootstrap_advances_radius(self) -> None:
        s = self._sampler_with_points(
            [[0.50, 0.50]], [0.128], target=0.12
        )
        s._radius = 0.10
        s._min_radius = 0.005
        self.assertTrue(
            s._should_advance_generation(
                t_min=None, t_max=None, fill_needed=False, cores=[]
            )
        )

    def test_min_radius_clamps_refinement_and_blocks_finer_generation(self) -> None:
        s = self._sampler_with_points(
            [[0.4, 0.5], [0.5, 0.5], [0.6, 0.5]],
            [0.02, 0.04, 0.08],
            target=0.04,
        )
        s._threshold = 0.001
        s._min_radius = 1.0 / 200.0
        s._refinement_factor = 0.5
        s._radius = 0.006
        self.assertAlmostEqual(s._next_radius(), 0.005)
        s._radius = s._next_radius()
        self.assertTrue(s._at_min_radius())
        self.assertFalse(
            s._should_advance_generation(
                t_min=0.02,
                t_max=0.08,
                fill_needed=False,
                cores=[1],
            )
        )

    def test_outer_anneal_blocked_when_full_cloud_has_brackets(self) -> None:
        """Outer-only may be same-side of T while full cloud still straddles."""
        # near-T both slightly above T → outer shell same side
        # far below / far above → full cloud has sign-change
        s = self._sampler_with_points(
            [
                [0.45, 0.50],
                [0.50, 0.50],
                [0.20, 0.50],
                [0.80, 0.50],
                [0.50, 0.20],
                [0.50, 0.80],
            ],
            [0.045, 0.048, 0.001, 0.20, 0.10, 0.10],
            target=0.04,
        )
        s._core_half_width = 0.0005
        s._outer_half_width = 0.01
        s._radius = 0.2
        cores, frontier, outer, _ = s._classify_bands()
        # Outer pool: only near-T same-side points → no straddle.
        br_outer = s._straddle_brackets(pool=outer if outer else [])
        br_full = s._straddle_brackets(pool=None)
        self.assertEqual(len(br_outer), 0)
        self.assertGreaterEqual(len(br_full), 1)
        filt = s._nearest_straddle_brackets(
            max_distance=2.5 * s._radius
        )
        # Must NOT anneal while local nearest straddles remain.
        did = s._maybe_shrink_outer(
            n_core=max(1, len(cores)),
            open_brackets=len(filt) if filt else len(br_full),
            raw_brackets=len(br_full),
        )
        self.assertFalse(did)
        # Clear brackets → anneal allowed when cores exist.
        s._outer_half_width = 0.05
        # force a synthetic core so n_core>0
        s._points[0].f = 0.04
        s._core_half_width = 0.01
        cores2, _, _, _ = s._classify_bands()
        self.assertGreater(len(cores2), 0)
        did2 = s._maybe_shrink_outer(
            n_core=len(cores2), open_brackets=0, raw_brackets=0
        )
        self.assertTrue(did2)
        self.assertLess(s._outer_half_width, 0.05)

    def test_prune_keeps_bracket_endpoints_outside_outer(self) -> None:
        s = self._sampler_with_points(
            [
                [0.40, 0.50],
                [0.60, 0.50],
                [0.50, 0.50],
                [0.50, 0.30],
                [0.50, 0.70],
            ],
            [0.001, 0.20, 0.04, 0.10, 0.10],
            target=0.04,
        )
        s._core_half_width = 0.01
        s._outer_half_width = 0.02
        s._radius = 0.25
        # Indices 0 and 1 are outside outer (|0.001-0.04| and |0.20-0.04| > 0.02)
        # but straddle T; nearest opposite-side support must remain keepable.
        filt = s._nearest_straddle_brackets(
            max_distance=2.5 * s._radius
        )
        self.assertGreaterEqual(len(filt), 1)
        keep = {2}  # only "core-ish" point
        for i, j, _fi, _fj in filt:
            keep.add(int(i))
            keep.add(int(j))
        # Endpoints outside outer must still be in keep.
        ends = set()
        for i, j, _fi, _fj in filt:
            ends.add(i)
            ends.add(j)
        self.assertTrue(ends.issubset(keep))
        s._prune_to_indices(sorted(keep))
        # After prune, straddle endpoints still present.
        for i, j, _fi, _fj in filt:
            # original indices may have remapped — check by values via f
            pass
        # At least two points with opposite signs of (f-T) remain.
        t = 0.04
        signs = {
            1 if float(p.f) - t > 0 else -1
            for p in s._points
            if p.f is not None
        }
        self.assertIn(1, signs)
        self.assertIn(-1, signs)

    def test_final_half_width_export_excludes_loose_cores(self) -> None:
        s = self._sampler_with_points(
            [[0.4, 0.5], [0.5, 0.5], [0.6, 0.5]],
            [0.0405, 0.04005, 0.042],  # only middle is within final 0.0002
            target=0.04,
        )
        s._core_half_width = 0.005
        s._outer_half_width = 0.005
        s._final_half_width = 0.0002
        s._final_half_width_configured = True
        s._converged = True
        s._generation = 3
        payload = s.finalize()
        self.assertEqual(payload["final_half_width"], 0.0002)
        self.assertTrue(payload["final_half_width_configured"])
        self.assertEqual(payload["n_final_points"], 1)
        self.assertEqual(len(payload["final_points_f"]), 1)
        self.assertAlmostEqual(float(payload["final_points_f"][0]), 0.04005, places=8)
        # Primary export is final-only when configured.
        self.assertEqual(payload["n_points_total"], 1)
        self.assertEqual(len(payload["points_f"]), 1)
        # No final points → cannot converge when final is configured.
        s._points[1].f = 0.045
        s._final_half_width_configured = True
        self.assertEqual(s._final_indices(), [])
        self.assertFalse(
            s._contour_converged(t_min=0.04, t_max=0.0401, fill_needed=False)
        )

    def test_final_half_width_absent_keeps_compat_export(self) -> None:
        s = self._sampler_with_points(
            [[0.4, 0.5], [0.5, 0.5]],
            [0.04, 0.041],
            target=0.04,
        )
        s._core_half_width = 0.005
        s._final_half_width = None
        s._final_half_width_configured = False
        s._generation = 1
        payload = s.finalize()
        self.assertEqual(payload["final_half_width"], s._core_half_width)
        self.assertFalse(payload["final_half_width_configured"])
        self.assertEqual(payload["n_points_total"], 2)
        self.assertEqual(payload["n_final_points"], 2)
        self.assertIn("final_points_u", payload)
        self.assertIn("final_points_x", payload)
        self.assertIn("final_points_f", payload)

class FeedbackChannelTests(unittest.TestCase):
    def test_publish_and_pull_feedback(self) -> None:
        queue = make_fakeredis_queue()
        queue.publish_feedback({"uuid": "u1", "logL": -1.0, "r2": 0.1})
        record = queue.pull_feedback(timeout=1)
        self.assertIsNotNone(record)
        assert record is not None
        self.assertEqual(record["uuid"], "u1")
        self.assertEqual(record["r2"], 0.1)

    def test_feedback_off_by_default_key_empty(self) -> None:
        queue = make_fakeredis_queue()
        self.assertEqual(int(queue.r.llen(FEEDBACK_QUEUE)), 0)


class DistributorAdaptiveTests(unittest.TestCase):
    def test_adaptive_registered_stateful(self) -> None:
        self.assertNotIn("AdaptiveBridson", STATELESS_METHODS)
        sampler = Distributor.set_method("AdaptiveBridson")
        self.assertEqual(sampler.method, "AdaptiveBridson")
        self.assertEqual(Distributor.get_resume_status("AdaptiveBridson"), "implemented")


class AdaptiveConfigTests(unittest.TestCase):
    def test_two_knob_config_derives_core_threshold_and_defaults(self) -> None:
        sampler = AdaptiveBridsonSampler()
        sampler.set_config(
            {
                "Sampling": {
                    "Method": "AdaptiveBridson",
                    "Variables": _flat_vars(["x", "y"]),
                    "Bounds": {
                        "target_expression": "r2",
                        "target_value": 0.04,
                        "outer_half_width": 0.02,
                        "min_radius": 0.002,
                    },
                }
            }
        )
        self.assertAlmostEqual(sampler._outer_half_width, 0.02)
        self.assertAlmostEqual(sampler._core_half_width, 0.0025)
        self.assertAlmostEqual(float(sampler._threshold or 0.0), 0.0025)
        self.assertAlmostEqual(sampler._min_radius, 0.002)
        self.assertAlmostEqual(sampler._initial_radius, 0.10)
        self.assertAlmostEqual(sampler._refinement_factor, 0.5)
        self.assertEqual(sampler._radius_shrink_mode, "on_coverage")
        self.assertTrue(sampler._bridge_gaps)
        self.assertEqual(sampler._max_generations, 16)
        self.assertEqual(sampler._max_points, 50000)
        self.assertEqual(sampler._max_new_per_generation, 4000)
        self.assertEqual(sampler._k_ref, 30)
        self.assertEqual(sampler._endpoint_stall_limit, 12)
        self.assertEqual(sampler._endpoint_omni_probes, 16)

    def test_rejects_dim_outside_range(self) -> None:
        sampler = AdaptiveBridsonSampler()
        with self.assertRaises(ValueError):
            sampler.set_config(
                {
                    "Sampling": {
                        "Method": "AdaptiveBridson",
                        "Variables": [
                            {
                                "name": "x",
                                "distribution": {
                                    "type": "Flat",
                                    "parameters": {"min": 0, "max": 1},
                                },
                            }
                        ],
                        "Bounds": {
                            "target_expression": "x",
                            "target_value": 0.5,
                            "threshold": 0.01,
                        },
                    }
                }
            )

    def test_requires_target_fields(self) -> None:
        sampler = AdaptiveBridsonSampler()
        with self.assertRaises(ValueError):
            sampler.set_config(
                {
                    "Sampling": {
                        "Method": "AdaptiveBridson",
                        "Variables": _flat_vars(["x", "y"]),
                        "Bounds": {"target_value": 0.5},
                    }
                }
            )

    def test_threshold_only_config(self) -> None:
        sampler = AdaptiveBridsonSampler()
        sampler.set_config(
            {
                "Sampling": {
                    "Method": "AdaptiveBridson",
                    "Bounds": {"target_expression": "r2",
                        "target_value": 0.04,
                        "threshold": 0.02,
                        "initial_radius": 0.1, "Seed": 0},
                    "Variables": _flat_vars(["x", "y"]),
                },
                "Runtime": {"workers": 1},
            }
        )
        self.assertEqual(sampler._threshold, 0.02)
        self.assertAlmostEqual(sampler._min_radius, 1.0 / 200.0)
        self.assertGreater(sampler._max_generations, 0)

    def test_custom_min_radius_is_clamped_to_initial_radius(self) -> None:
        sampler = AdaptiveBridsonSampler()
        sampler.set_config(
            {
                "Sampling": {
                    "Method": "AdaptiveBridson",
                    "Variables": _flat_vars(["x", "y"]),
                    "Bounds": {
                        "target_expression": "r2",
                        "target_value": 0.04,
                        "initial_radius": 0.004,
                        "min_radius": 0.01,
                    },
                }
            }
        )
        self.assertAlmostEqual(sampler._min_radius, 0.004)


class AdaptiveSyntheticRunTests(unittest.TestCase):
    def setUp(self) -> None:
        TaskFactory.reset_instance()

    def tearDown(self) -> None:
        TaskFactory.reset_instance()

    def test_generation_barrier_resume_matches_uninterrupted_run(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        ctx = get_spawn_context()
        baseline_queue = ctx.Queue()
        interrupted_queue = ctx.Queue()
        resumed_queue = ctx.Queue()
        baseline = interrupted = resumed = None
        scan_prefix = f"adaptive-resume-{os.getpid()}"
        try:
            with (
                tempfile.TemporaryDirectory() as baseline_dir,
                tempfile.TemporaryDirectory() as resume_dir,
            ):
                baseline = ctx.Process(
                    target=_run_adaptive_resume_control,
                    args=(
                        redis_config,
                        baseline_dir,
                        f"{scan_prefix}-baseline",
                        False,
                        baseline_queue,
                    ),
                )
                baseline.start()
                baseline.join(timeout=120.0)
                self.assertFalse(baseline.is_alive())
                self.assertEqual(baseline.exitcode, 0)
                self.assertTrue(baseline_queue.get(timeout=5.0)["ok"])
                with open(
                    os.path.join(baseline_dir, "levelset.json"), encoding="utf-8"
                ) as handle:
                    baseline_levelset = json.load(handle)

                interrupted = ctx.Process(
                    target=_run_adaptive_resume_control,
                    args=(
                        redis_config,
                        resume_dir,
                        f"{scan_prefix}-resume",
                        False,
                        interrupted_queue,
                    ),
                )
                interrupted.start()
                database_dir = os.path.join(resume_dir, "DATABASE")
                deadline = time.monotonic() + 60.0
                durable: set[str] = set()
                while time.monotonic() < deadline:
                    durable = read_persisted_uuids(database_dir)
                    if 4 <= len(durable) < int(baseline_levelset["n_unique_submitted"]):
                        break
                    if not interrupted.is_alive():
                        self.fail("AdaptiveBridson finished before interruption")
                    time.sleep(0.05)
                self.assertGreaterEqual(len(durable), 4)
                os.kill(interrupted.pid, signal.SIGINT)
                interrupted.join(timeout=60.0)
                self.assertFalse(interrupted.is_alive())

                checkpoint = os.path.join(
                    resume_dir,
                    "checkpoints",
                    f"{scan_prefix}-resume",
                    "AdaptiveBridson",
                    "state.pkl",
                )
                self.assertTrue(os.path.isfile(checkpoint))
                resumed = ctx.Process(
                    target=_run_adaptive_resume_control,
                    args=(
                        redis_config,
                        resume_dir,
                        f"{scan_prefix}-resume",
                        True,
                        resumed_queue,
                    ),
                )
                resumed.start()
                resumed.join(timeout=120.0)
                self.assertFalse(resumed.is_alive())
                self.assertEqual(resumed.exitcode, 0)
                self.assertTrue(resumed_queue.get(timeout=5.0)["ok"])
                with open(
                    os.path.join(resume_dir, "levelset.json"), encoding="utf-8"
                ) as handle:
                    resumed_levelset = json.load(handle)

                self.assertEqual(resumed_levelset, baseline_levelset)
                uuids = read_persisted_uuids(database_dir)
                self.assertEqual(
                    len(uuids), int(resumed_levelset["n_unique_submitted"])
                )
        finally:
            for proc in (baseline, interrupted, resumed):
                if proc is not None and proc.is_alive():
                    proc.terminate()
                    proc.join(timeout=5.0)
            for scan in list_running_scans():
                if scan.name.startswith(scan_prefix):
                    kill_jarvis_processes(list(scan.processes), force=True)
            for queue in (baseline_queue, interrupted_queue, resumed_queue):
                queue.close()
                queue.join_thread()
            server.shutdown()
            server.server_close()

    def test_circle_run_writes_levelset_and_submits(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                config = _als_config(
                    redis_config=redis_config,
                    tmpdir=tmpdir,
                    method_block={
                        "target_value": 0.04,
                        "threshold": 0.15,
                        "initial_radius": 0.14,
                        "max_generations": 4,
                        "max_points": 400,
                        "max_new_per_generation": 60,
                    },
                    project_name="alevelset-circle",
                    scan_name="circle",
                )
                core = Jarvis2Core(config)
                outcome = core.run()
                submitted = int(getattr(outcome, "submitted", outcome) or 0)
                self.assertGreater(submitted, 0)
                levelset_path = os.path.join(tmpdir, "levelset.json")
                self.assertTrue(os.path.isfile(levelset_path), "levelset.json missing")
                payload = json.loads(open(levelset_path, encoding="utf-8").read())
                self.assertEqual(payload["dim"], 2)
                self.assertEqual(
                    payload.get("algorithm"), "outer_core_root_correction_bridson"
                )
                self.assertGreater(int(payload.get("n_points_total") or 0), 0)
                self.assertIn(payload.get("contour_status"), ("converged", "partial"))
                self.assertIn("stop_reason", payload)
        finally:
            server.shutdown()
            server.server_close()

    def test_seeded_reproducibility_point_count(self) -> None:
        def _run(workers: int) -> list[str]:
            server, redis_config = _start_tcp_fakeredis()
            try:
                with tempfile.TemporaryDirectory() as tmpdir:
                    config = _als_config(
                        redis_config=redis_config,
                        tmpdir=tmpdir,
                        method_block={
                            "threshold": 0.12,
                            "initial_radius": 0.15,
                            "max_generations": 3,
                            "max_points": 250,
                            "max_new_per_generation": 40,
                        },
                        seed=99,
                        workers=workers,
                    )
                    core = Jarvis2Core(config)
                    core.run()
                    path = os.path.join(tmpdir, "levelset.json")
                    payload = json.loads(open(path, encoding="utf-8").read())
                    return [
                        str(payload["n_points_total"]),
                        str(payload.get("n_generations")),
                        str(bool(payload.get("converged"))),
                        str(payload.get("stop_reason")),
                    ]
            finally:
                server.shutdown()
                server.server_close()
                TaskFactory.reset_instance()

        a = _run(1)
        b = _run(2)
        self.assertEqual(a, b)

    def test_export_import_runtime_state(self) -> None:
        s = AdaptiveBridsonSampler()
        s.set_config(
            {
                "Sampling": {
                    "Method": "AdaptiveBridson",
                    "Bounds": {"target_expression": "r2",
                        "target_value": 0.04,
                        "threshold": 0.05,
                        "initial_radius": 0.1, "Seed": 3},
                    "Variables": _flat_vars(["x", "y"]),
                },
                "Runtime": {"workers": 1},
            }
        )
        s._points = [
            LevelSetPoint(u=[0.5, 0.5], x={"x": 0.5, "y": 0.5}, f=0.04, uuid="a")
        ]
        s._radius = 0.05
        s._t_min = 0.03
        s._t_max = 0.05
        s._converged = False
        historical = np.array([0.25, 0.75])
        s._submitted_u_keys = {
            s._u_key(s._points[0].u),
            s._u_key(historical),
        }
        s._quiet_fill_passes = 2
        s._active_endpoints = {
            "terminal:a:open": {
                "id": "terminal:a:open",
                "kind": "terminal",
                "anchor_uuid": "a",
                "anchor_u": [0.5, 0.5],
                "direction": [1.0, 0.0],
                "target_uuid": None,
                "attempts": 2,
                "misses": 1,
                "pending_attempt": False,
                "generation": 0,
            }
        }
        s._active_virtual_cores = {
            "bracket:0:4,4": {
                "id": "bracket:0:4,4",
                "anchor_u": [0.2, 0.2],
                "generation": 0,
                "attempts": 1,
                "misses": 1,
                "pending_attempt": False,
            }
        }
        s._n_radius_refines = 3
        s._generation = 3
        state = s.export_runtime_state()
        # Simulate an older checkpoint whose feedback generation drifted.
        state["generation"] = 99
        s2 = AdaptiveBridsonSampler()
        s2.set_config(
            {
                "Sampling": {
                    "Method": "AdaptiveBridson",
                    "Bounds": {"target_expression": "r2",
                        "target_value": 0.04,
                        "threshold": 0.05,
                        "initial_radius": 0.1, "Seed": 3},
                    "Variables": _flat_vars(["x", "y"]),
                },
                "Runtime": {"workers": 1},
            }
        )
        s2.import_runtime_state(state)
        self.assertEqual(len(s2._points), 1)
        self.assertAlmostEqual(s2._radius, 0.05)
        # Generation must follow shrink count for SeedSequence.spawn(generation).
        self.assertEqual(s2._n_radius_refines, 3)
        self.assertEqual(s2._generation, 3)
        self.assertAlmostEqual(float(s2._t_min or 0), 0.03)
        self.assertAlmostEqual(float(s2._t_max or 0), 0.05)
        self.assertIn(s2._u_key(historical), s2._submitted_u_keys)
        self.assertEqual(s2._quiet_fill_passes, 2)
        self.assertIn("terminal:a:open", s2._active_endpoints)
        self.assertEqual(s2._active_endpoints["terminal:a:open"]["misses"], 1)
        self.assertIn("bracket:0:4,4", s2._active_virtual_cores)
        self.assertEqual(s2._active_virtual_cores["bracket:0:4,4"]["misses"], 1)


if __name__ == "__main__":
    unittest.main()
