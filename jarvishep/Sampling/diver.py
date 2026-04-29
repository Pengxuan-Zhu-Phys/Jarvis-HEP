#!/usr/bin/env python3

from __future__ import annotations

import concurrent.futures
import sys
from copy import deepcopy
from dataclasses import asdict
from typing import Any

import numpy as np

from jarvishep.log_kv import format_two_column_log
from jarvishep.Sampling.sampler import SamplingVirtial
from jarvishep.Sampling.Source.Diver.de import DEConfig, DifferentialEvolution
from jarvishep.sample import Sample


class Diver(SamplingVirtial):
    """Diver sampler integrated into Jarvis-HEP.

    Implementation notes:
    - Pure Python/NumPy DE backend (no Fortran dependency).
    - Internal optimization target is ``cost = -LogL`` (minimization), aligned
      with original Diver conventions.
    - Parameter vectors are evolved in unit-cube space and mapped into physical
      parameter distributions before workflow evaluation.
    """

    _BNDRY_MAP = {
        "not enforced": 0,
        "none": 0,
        "off": 0,
        "brick wall": 1,
        "brick": 1,
        "random re-initialization": 2,
        "random reinitialization": 2,
        "random": 2,
        "reflection": 3,
        "reflect": 3,
    }

    def __init__(self) -> None:
        super().__init__()
        self.method = "Diver"
        self.load_schema_file()
        self._D = 0
        self._run_cfg: dict[str, Any] = {}
        self._de_cfg: DEConfig | None = None
        self._de_result = None
        self._best_params: dict[str, float] | None = None
        self._de_state: dict[str, Any] | None = None

    def load_schema_file(self):
        self.schema = self.path.get("DiverSchema", self.path.get("RandomSchema"))

    def set_logger(self, logger) -> None:
        super().set_logger(logger)
        self.logger.warning("Sampling method initializaing ...")

    def set_factory(self, factory) -> None:
        self.factory = factory
        self.logger.warning("WorkerFactory is ready for Diver sampler")

    def supports_runtime_checkpointing(self) -> bool:
        return True

    def set_config(self, config_info) -> None:
        self.config = config_info
        self.set_bucket_alloc()
        self.init_generator()

    def set_likelihood(self, loglike):
        super().set_likelihood(loglike)

    def init_generator(self) -> None:
        self.load_variable()
        self._D = len(self.vars)

        self._run_cfg = self._extract_run_config()
        self._de_cfg = DEConfig(
            dim=self._D,
            pop_size=self._run_cfg["NP"],
            mode=self._run_cfg["mode"],
            max_gen=self._run_cfg["maxgen"],
            max_civ=self._run_cfg["maxciv"],
            F=self._run_cfg["F"],
            Cr=self._run_cfg["Cr"],
            lambda_=self._run_cfg["lambda"],
            current=self._run_cfg["current"],
            expon=self._run_cfg["expon"],
            bndry=self._run_cfg["bndry"],
            convthresh=self._run_cfg["convthresh"],
            convsteps=self._run_cfg["convsteps"],
            remove_duplicates=self._run_cfg["removeDuplicates"],
            discard_unfit_points=self._run_cfg["discard_unfit_points"],
            max_acceptable_value=self._run_cfg["max_acceptable_value"],
            seed=self._run_cfg["seed"],
        )

        self.total_core = int(self.config.get("Sampling", {}).get("MaxWorker", self.total_core))
        self.load_nuisance_sampler()
        if self._with_nuisance:
            self.logger.error("Diver sampler currently does not support Sampling.Nuisance.")
            sys.exit(2)

    def initialize(self):
        self.logger.warning("Initializing the Diver Sampling")
        self.logger.warning("Diver Sampling\nCopyright Elinore Roebber and Pat Scott")
        self.logger.warning(
            "If you use Diver Sampling method in your work, please cite\n"
            "\t1.) arXiv:1705.07959\n"
            "\t2.) arXiv:2101.04525"
        )

    def _export_sampler_state(self):
        if self._de_result is not None and hasattr(self._de_result, "_asdict"):
            de_result = dict(self._de_result._asdict())
        elif self._de_result is not None and isinstance(self._de_result, dict):
            de_result = dict(self._de_result)
        elif self._de_result is not None:
            de_result = {
                "best_vector": getattr(self._de_result, "best_vector", None),
                "best_cost": float(getattr(self._de_result, "best_cost", np.inf)),
                "best_loglike": float(getattr(self._de_result, "best_loglike", float("-inf"))),
                "best_civilization": int(getattr(self._de_result, "best_civilization", 0)),
                "best_generation": int(getattr(self._de_result, "best_generation", 0)),
                "evaluations": int(getattr(self._de_result, "evaluations", 0)),
                "acceptance_rate": float(getattr(self._de_result, "acceptance_rate", 0.0)),
                "population": getattr(self._de_result, "population", None),
                "costs": getattr(self._de_result, "costs", None),
            }
        else:
            de_result = None
        return {
            "D": int(self._D),
            "run_cfg": dict(self._run_cfg),
            "de_cfg": None if self._de_cfg is None else asdict(self._de_cfg),
            "de_result": de_result,
            "best_params": None if self._best_params is None else dict(self._best_params),
            "de_state": None if self._de_state is None else dict(self._de_state),
            "numpy_random_state": np.random.get_state(),
        }

    def _import_sampler_state(self, payload):
        self._D = int(payload.get("D", self._D or 0))
        self._run_cfg = dict(payload.get("run_cfg", self._run_cfg))
        de_cfg = payload.get("de_cfg")
        if de_cfg is not None:
            self._de_cfg = DEConfig(**de_cfg)
        best_params = payload.get("best_params")
        self._best_params = dict(best_params) if isinstance(best_params, dict) else None
        self._de_state = deepcopy(payload.get("de_state", self._de_state))
        np_state = payload.get("numpy_random_state")
        if np_state is not None:
            np.random.set_state(np_state)
        self._de_result = payload.get("de_result")

    def __iter__(self):
        return self

    def __next__(self):
        # Keep one-point-check mode compatible with existing Core.test_assembly_line.
        unit = np.random.random(self._D)
        return self.map_point_into_distribution(unit)

    def next_sample(self):
        return self.__next__()

    def map_point_into_distribution(self, row) -> dict[str, float]:
        mapped = {}
        eps = 1e-12
        for i, var in enumerate(self.vars):
            u = float(np.clip(row[i], eps, 1.0 - eps))
            value = var.map_standard_random_to_distribution(u)
            if isinstance(value, np.generic):
                value = value.item()
            mapped[var.name] = value
        return mapped

    def _extract_run_config(self) -> dict[str, Any]:
        sampling_cfg = self.config.get("Sampling", {}) or {}
        run = sampling_cfg.get("run", sampling_cfg.get("Bounds", {})) or {}

        mode = str(run.get("mode", "current"))
        NP = int(run.get("NP", max(10 * self._D, 40)))
        maxgen = int(run.get("maxgen", 300))
        maxciv = int(run.get("maxciv", 1))

        F_raw = run.get("F", 0.7)
        if isinstance(F_raw, (list, tuple, np.ndarray)):
            F = float(F_raw[0]) if len(F_raw) else 0.7
        else:
            F = float(F_raw)

        Cr = float(run.get("Cr", 0.9))
        lambda_ = float(run.get("lambda", 0.0))
        current = bool(run.get("current", False))
        expon = bool(run.get("expon", False))

        bndry_raw = run.get("bndry", "brick wall")
        if isinstance(bndry_raw, str):
            bndry = self._BNDRY_MAP.get(bndry_raw.strip().lower(), 1)
        else:
            bndry = int(bndry_raw)
            if bndry not in (0, 1, 2, 3):
                bndry = 1

        convthresh = float(run.get("convthresh", 1e-3))
        convsteps = int(run.get("convsteps", 10))
        remove_duplicates = bool(run.get("removeDuplicates", False))
        discard_unfit_points = bool(run.get("discard_unfit_points", False))
        max_acceptable_value = float(run.get("max_acceptable_value", np.inf))

        seed = run.get("seed", None)
        if seed is not None:
            seed = int(seed)

        return {
            "mode": mode,
            "NP": NP,
            "maxgen": maxgen,
            "maxciv": maxciv,
            "F": F,
            "Cr": Cr,
            "lambda": lambda_,
            "current": current,
            "expon": expon,
            "bndry": bndry,
            "convthresh": convthresh,
            "convsteps": convsteps,
            "removeDuplicates": remove_duplicates,
            "discard_unfit_points": discard_unfit_points,
            "max_acceptable_value": max_acceptable_value,
            "seed": seed,
        }

    def _evaluate_population_loglike(self, unit_population: np.ndarray) -> np.ndarray:
        """Evaluate a block of unit-cube vectors via Jarvis-HEP workflow."""
        npts = int(unit_population.shape[0])
        if npts == 0:
            return np.empty(0, dtype=float)

        values = np.full(npts, -np.inf, dtype=float)
        future_to_idx: dict[concurrent.futures.Future, int] = {}
        future_to_sample: dict[concurrent.futures.Future, Sample] = {}

        base_sample_cfg = self.info["sample"]

        for idx, row in enumerate(unit_population):
            params = self.map_point_into_distribution(row)
            sample = Sample(params)
            sample_cfg = self.build_sample_config(
                base_sample_cfg,
                save_dir=self._next_bucket_dir_for_sample(),
            )
            sample.set_config(sample_cfg)
            try:
                future = self.factory.submit_task(sample.info)
            except Exception:
                self._on_sample_completed(sample.info)
                sample.close()
                raise
            future_to_idx[future] = idx
            future_to_sample[future] = sample

        pending = set(future_to_idx)
        while pending:
            timeout = self._checkpoint_seconds_until_due()
            done, pending = concurrent.futures.wait(
                pending,
                timeout=timeout,
                return_when=concurrent.futures.FIRST_COMPLETED,
            )
            if not done:
                self.persist_runtime_checkpoint(reason="checkpoint_heartbeat")
                continue
            for future in done:
                idx = future_to_idx[future]
                sample = future_to_sample.pop(future)

                try:
                    result = float(future.result())
                    if np.isfinite(result):
                        values[idx] = result
                except Exception as exc:
                    self.logger.error(
                        format_two_column_log(
                            "[WorkerFactory] future exception consumed",
                            [("uuid", sample.uuid), ("error", exc)],
                        )
                    )
                    raise RuntimeError(f"Diver sample evaluation failed at index {idx}") from exc
                finally:
                    self._on_sample_completed(sample.info)
                    sample.close()

        return values

    def run_nested(self):
        if self._de_cfg is None:
            self.logger.error("Diver sampler is not configured.")
            sys.exit(2)
        if not hasattr(self, "factory") or self.factory is None:
            self.logger.error("Diver sampler has no WorkerFactory.")
            sys.exit(2)

        if self._de_result is not None and self._best_params is not None:
            self.logger.warning("Diver sampler resume detected -> using restored DE result")
            return

        optimizer = DifferentialEvolution(self._de_cfg)
        if self._de_state is not None:
            optimizer.import_state(self._de_state)

        def _checkpoint_diver_state(engine: DifferentialEvolution) -> None:
            self._de_state = engine.export_state()
            if self._runtime_checkpoint_enabled:
                self.persist_runtime_checkpoint(force=True, reason="Diver civilization checkpoint")

        self._de_result = optimizer.run(
            self._evaluate_population_loglike,
            logger=self.logger,
            checkpoint_callback=_checkpoint_diver_state,
        )

        self._best_params = self.map_point_into_distribution(self._de_result.best_vector)
        self._de_state = optimizer.export_state()
        self.logger.warning(
            "Diver finished: best LogL={:.6e}, civ={}, gen={}, evals={}, accept_rate={:.3f}".format(
                self._de_result.best_loglike,
                self._de_result.best_civilization,
                self._de_result.best_generation,
                self._de_result.evaluations,
                self._de_result.acceptance_rate,
            )
        )
        self.logger.warning(f"Diver best-fit parameters -> {self._best_params}")

    def finalize(self):
        pass

    def combine_data(self, df_full) -> None:
        # Diver results are already persisted through the global HDF5 writer pipeline.
        return
