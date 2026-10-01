"""Sampler logging contract (docs/logging-spec.md §2).

Every method in the sampler catalog must produce the required lifecycle
records (S1 start, S2 settings, S3 ready, S9 result, S10 summary) under its
own ``Jarvis-HEP.Sampler.<Method>`` label, in the Jarvis (V1) layout, and its
log hooks must work. A new sampler, built-in or plug-in, that breaks this
fails here.
"""

from __future__ import annotations

import logging
import os
import re
import tempfile
import unittest
from types import SimpleNamespace
from typing import Any

from jarvishep2.distributor import Distributor
from jarvishep2.logging.toplevel import JarvisContextFormatter, get_jarvis_logger
from jarvishep2.sampling.lifecycle_log import SamplerLifecycleLog

_STATELESS_BOUNDS: dict[str, dict[str, Any]] = {
    "Random": {"point_number": 10},
    "Grid": {},
    "Bridson": {"radius": 0.3, "max_attempt": 30},
    "CSV": {},
    "AdaptiveBridson": {
        "target_expression": "z",
        "target_value": 1.0,
        "outer_half_width": 0.2,
        "min_radius": 0.05,
    },
    "Dynesty": {"nlive": 20},
    "MultiNest": {"nlive": 20},
}
_MCMC_BOUNDS = {"num_chains": 4, "num_iters": 10, "proposal_scale": 0.1}

_V1_HEADER = re.compile(
    r"^\n·•· (?P<module>\S+) \n\t-> \d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{3} - "
    r"\[(?P<level>[A-Z]+)\] >>> \n"
)


def _config(method: str, tmp: str) -> dict[str, Any]:
    variables = []
    for name in ("x", "y"):
        parameters: dict[str, Any] = {"min": 0.0, "max": 1.0}
        if method == "Grid":
            parameters["num"] = 3
        variables.append(
            {"name": name, "distribution": {"type": "Flat", "parameters": parameters}}
        )
    bounds = dict(_STATELESS_BOUNDS.get(method, _MCMC_BOUNDS))
    bounds.setdefault("seed", 7)
    sampling: dict[str, Any] = {
        "Method": method,
        "Bounds": bounds,
        "Variables": variables,
        "LogLikelihood": [{"name": "LogL_Z", "expression": "LogGauss(z, 1, 1)"}],
    }
    if method == "CSV":
        path = os.path.join(tmp, "points.csv")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("x,y\n0.1,0.2\n0.3,0.4\n")
        sampling["Variables"] = []
        sampling["CSV"] = {"path": path}
        bounds["path"] = path
    return {
        "project_name": "contract",
        "Scan": {"name": f"{method.lower()}-contract"},
        "task_root": tmp,
        "task_result_dir": tmp,
        "Runtime": {"mode": "redis", "workers": 2, "batch_size": 2},
        "Sampling": sampling,
        "Operas": {
            "Modules": [
                {
                    "name": "EggBox",
                    "operator": "jarvishep2.testing.eggbox.eggbox2d_numpy",
                    "call_mode": "call",
                    "input": [
                        {"name": "x", "expression": "x"},
                        {"name": "y", "expression": "y"},
                    ],
                    "output": [{"name": "z", "entry": "z"}],
                }
            ]
        },
    }


class _Capture(logging.Handler):
    def __init__(self) -> None:
        super().__init__(level=logging.DEBUG)
        self.records: list[logging.LogRecord] = []
        self.setFormatter(JarvisContextFormatter(colorize=False))

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


class SamplerLoggingContractTests(unittest.TestCase):
    def _run_lifecycle(self, method: str) -> tuple[_Capture, Any]:
        sampler = Distributor.set_method(method)
        tmp = tempfile.mkdtemp(prefix="jarvis-log-contract-")
        sampler.set_config(_config(method, tmp))
        capture = _Capture()
        base = logging.getLogger("jarvis_hep")
        base.addHandler(capture)
        old_level = base.level
        base.setLevel(logging.DEBUG)
        try:
            log = SamplerLifecycleLog(sampler)
            log.start(resume=False, workers=2)
            log.ready()
            log.checkpoint_saved("/tmp/state.pkl", "heartbeat")
            log.result(
                SimpleNamespace(status="success", submitted=4, completed=3, failed=1)
            )
        finally:
            base.removeHandler(capture)
            base.setLevel(old_level)
        return capture, sampler

    def test_every_catalog_method_emits_the_required_records(self) -> None:
        methods = Distributor.available_methods()
        self.assertGreaterEqual(len(methods), 15)
        for method in methods:
            with self.subTest(method=method):
                capture, _ = self._run_lifecycle(method)
                rendered = [capture.format(record) for record in capture.records]
                messages = [record.getMessage() for record in capture.records]
                levels = {record.getMessage().split("\n")[0]: record.levelname for record in capture.records}

                # No hook may fail.
                self.assertFalse(
                    [m for m in messages if "could not report" in m], messages
                )
                # Every record is labelled Jarvis-HEP.Sampler.<Method> in the V1 layout.
                for text in rendered:
                    match = _V1_HEADER.match(text)
                    self.assertIsNotNone(match, text)
                    self.assertEqual(match.group("module"), f"Jarvis-HEP.Sampler.{method}")

                start = f"Initializing the {method} Sampling"
                ready = f"WorkerFactory is ready for {method} sampler"
                result = f"{method} Sampler obtains 4 samples in "
                self.assertEqual(levels.get(start), "WARNING")
                self.assertEqual(levels.get(ready), "WARNING")
                self.assertTrue(
                    any(m.startswith(result) for m in messages if m in levels), messages
                )
                settings = [m for m in messages if m.startswith(f"{method} Sampler Settings ->")]
                summary = [m for m in messages if m.startswith(f"{method} Sampler Summary ->")]
                self.assertEqual(len(settings), 1, messages)
                self.assertEqual(len(summary), 1, messages)
                for key in ("method", "seed", "variables", "resume"):
                    self.assertIn(key, settings[0])
                for key in ("stop reason", "submitted", "completed", "failed", "elapsed"):
                    self.assertIn(key, summary[0])
                self.assertIn("Checkpoint saved to /tmp/state.pkl (heartbeat)", messages)

    def test_seed_is_recorded_in_settings(self) -> None:
        capture, _ = self._run_lifecycle("Random")
        settings = next(
            r.getMessage() for r in capture.records
            if r.getMessage().startswith("Random Sampler Settings ->")
        )
        self.assertRegex(settings, r"seed\s+\S*\s*7")

    def test_result_is_logged_once(self) -> None:
        sampler = Distributor.set_method("Random")
        with tempfile.TemporaryDirectory() as tmp:
            sampler.set_config(_config("Random", tmp))
            log = SamplerLifecycleLog(sampler)
            capture = _Capture()
            logger = logging.getLogger("jarvis_hep")
            logger.addHandler(capture)
            try:
                outcome = SimpleNamespace(status="interrupted", submitted=1, completed=1, failed=0)
                log.result(outcome)
                log.result(outcome)
            finally:
                logger.removeHandler(capture)
        results = [r for r in capture.records if "Sampler obtains" in r.getMessage()]
        self.assertEqual(len(results), 1)
        summary = next(r.getMessage() for r in capture.records if "Summary ->" in r.getMessage())
        self.assertIn("interrupted (continue with --resume)", summary)

    def test_failure_carries_traceback_to_files_only(self) -> None:
        sampler = Distributor.set_method("Random")
        log = SamplerLifecycleLog(sampler)
        capture = _Capture()
        logger = logging.getLogger("jarvis_hep")
        logger.addHandler(capture)
        try:
            try:
                raise ValueError("bad point")
            except ValueError as exc:
                log.failure("proposing points", exc)
        finally:
            logger.removeHandler(capture)
        record = capture.records[-1]
        self.assertEqual(record.levelname, "ERROR")
        self.assertEqual(
            record.getMessage(),
            "Random Sampler meets error when proposing points -> bad point",
        )
        file_text = JarvisContextFormatter(colorize=False).format(record)
        self.assertIn("Traceback (most recent call last)", file_text)
        record.exc_text = None
        screen_text = JarvisContextFormatter(colorize=True, show_traceback=False).format(record)
        self.assertNotIn("Traceback", screen_text)


class SamplerHookTests(unittest.TestCase):
    def test_mcmc_family_summary_rows(self) -> None:
        for method in ("MCMC", "ToyMCMC", "AMMCMC", "DRAM", "PTMCMC", "EnsembleMCMC", "DEMCMC", "PTEnsemble"):
            with self.subTest(method=method), tempfile.TemporaryDirectory() as tmp:
                sampler = Distributor.set_method(method)
                sampler.set_config(_config(method, tmp))
                rows = dict(sampler.log_summary_rows())
                for key in ("proposed", "accepted", "acceptance rate",
                            "rejected outside the prior", "failed evaluations"):
                    self.assertIn(key, rows)
                if method == "DRAM":
                    self.assertIn("stage-2 accepted", rows)
                if method in ("AMMCMC", "DRAM"):
                    self.assertIn("covariance updates", rows)
                if method.startswith("PT"):
                    self.assertIn("swaps accepted", rows)

    def test_mcmc_step_line_is_time_based(self) -> None:
        sampler = Distributor.set_method("ToyMCMC")
        with tempfile.TemporaryDirectory() as tmp:
            sampler.set_config(_config("ToyMCMC", tmp))
            lines: list[str] = []
            sampler._logger = SimpleNamespace(  # type: ignore[assignment]
                info=lambda msg, *a, **k: lines.append(msg % a),
                warning=lambda msg, *a, **k: None,
                debug=lambda msg, *a, **k: None,
            )
            sampler._ensure_registry()
            sampler._maybe_log_step()  # arms the timer, logs nothing
            sampler._maybe_log_step()  # within the interval: nothing
            self.assertEqual(lines, [])
            sampler._maybe_log_step(force=True)
        self.assertEqual(len(lines), 1)
        self.assertTrue(lines[0].startswith("ToyMCMC Chains -> steps -> "), lines)
        self.assertIn("accept rate -> 0: ", lines[0])

    def test_sampler_logger_label(self) -> None:
        logger = get_jarvis_logger("sampler.adaptive_bridson")
        record = logging.LogRecord("jarvis_hep.sampler.adaptive_bridson", logging.INFO, "", 0, "x", None, None)
        record.__dict__.update(getattr(logger, "extra", {}) or {})
        text = JarvisContextFormatter(colorize=False).format(record)
        self.assertIn("Jarvis-HEP.Sampler.AdaptiveBridson", text)


if __name__ == "__main__":
    unittest.main()
