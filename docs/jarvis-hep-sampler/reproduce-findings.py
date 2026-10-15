"""Review evidence only; no production fixes. Run from the repository root.

    python docs/jarvis-hep-sampler/reproduce-findings.py

Exercises the reviewed code using temporary files and in-memory fakeredis.
The JSON describes observed bugs; exit 0 means the evidence script completed,
not that the reviewed implementation passed a correctness test.
"""

from __future__ import annotations

import json
import logging
import sys
import tempfile
import threading
import warnings
from pathlib import Path
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from jarvishep2.sampling.Source.Diver import defaults, run
from jarvishep2.sampling.multinest_sampler import MultiNestSampler
from jarvishep2.redis_queue import make_fakeredis_queue


def constant(value):
    def objective(params, calls, quit, valid, context):
        return value, calls + 1, False

    return objective


def options(path, **kwargs):
    result = defaults([-1.0], [1.0])
    result.update(seed=7, verbose=0, maxgen=5, convsteps=2, path=str(path),
                  disableIO=True, outputRaw=False, outputSam=False)
    result.update(kwargs)
    return result


def reproduce():
    evidence = {}
    with tempfile.TemporaryDirectory(prefix="jhep-review-") as tmp:
        root = Path(tmp)
        try:
            run(constant(0.0), options(root / "zero"))
        except ZeroDivisionError as exc:
            evidence["diver_zero_mean"] = {"error": type(exc).__name__, "message": str(exc)}
        else:
            evidence["diver_zero_mean"] = {"error": None}

        cfg = options(root / "converged", disableIO=False, convsteps=1)
        first = run(constant(1.0), cfg)
        cfg["resume"] = True
        second = run(constant(1.0), cfg)
        evidence["diver_completed_resume"] = {
            "first": {"generation": first.generations, "calls": first.function_calls,
                      "converged": first.converged},
            "resumed": {"generation": second.generations, "calls": second.function_calls,
                        "converged": second.converged},
        }

        cfg = options(root / "atomic", disableIO=False, maxgen=1, convsteps=20)
        run(constant(1.0), cfg)
        checkpoint = root / "atomic.diver.npz"
        original = checkpoint.read_bytes()
        cfg.update(resume=True, maxgen=2)
        try:
            with patch("numpy.lib.format.write_array", side_effect=OSError("simulated disk failure")):
                run(constant(1.0), cfg)
        except OSError:
            pass
        current = checkpoint.read_bytes()
        try:
            run(constant(1.0), cfg)
        except Exception as exc:
            recovery_error = type(exc).__name__
        else:
            recovery_error = None
        evidence["diver_nonatomic_save"] = {
            "before_bytes": len(original), "after_bytes": len(current),
            "previous_checkpoint_preserved": original == current,
            "recovery_error": recovery_error,
        }

        queue = make_fakeredis_queue()
        stop = threading.Event()

        def worker():
            while not stop.is_set():
                task = queue.pull_task(timeout=1)
                if task is not None:
                    u = np.asarray(task["u_coords"], dtype=float)
                    queue.publish_feedback({"uuid": task["uuid"],
                                            "logL": -float(np.sum((u - 0.5) ** 2)) * 1000})

        thread = threading.Thread(target=worker, daemon=True)
        thread.start()
        messages = []

        class Capture(logging.Handler):
            def emit(self, record):
                messages.append(record.getMessage())

        capture = Capture()
        base_logger = logging.getLogger("jarvis_hep")
        base_logger.addHandler(capture)
        try:
            sampler = MultiNestSampler()
            sampler.set_config({
                "Sampling": {
                    "Method": "MultiNest",
                    "Variables": [{"name": "x", "distribution": {
                        "type": "Flat", "parameters": {"min": 0, "max": 1}}}],
                    "Bounds": {"seed": 7, "nlive": 10, "dlogz": 1e-8,
                               "run_nested": {"maxiter": 1, "print_progress": False}},
                },
                "Runtime": {"workers": 1, "batch_size": 1},
                "EnvReqs": {"V2": {"checkpoint": {"enabled": False}}},
                "task_root": tmp, "task_result_dir": tmp, "scan_name": "review",
            })
            sampler.set_redis(queue)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                sampler.run_adaptive(generation_timeout=30)
            evidence["nested_budget_stop"] = {
                "engine_budget_warning": [m for m in
                                          [*(str(w.message) for w in caught), *messages]
                                          if "maxiter/maxcall" in m],
                "finished": sampler._finished,
                "reported_reason": sampler.log_stop_reason(),
            }
        finally:
            base_logger.removeHandler(capture)
            stop.set()
            thread.join(timeout=3)
    return evidence


if __name__ == "__main__":
    print(json.dumps(reproduce(), ensure_ascii=False, indent=2))
