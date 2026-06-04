#!/usr/bin/env python3 
import os, sys
from re import S 
from abc import ABCMeta, abstractmethod
from mpmath.functions.functions import re
import numpy as np
import pandas as pd 
from numpy import meshgrid
from pandas.core.series import Series
from sympy import sympify
from sympy.geometry import parabola
import time 
from random import randint
# from jarvishep.Sampling.sampler import SamplingVirtial
from jarvishep.log_kv import format_two_column_log
from jarvishep.Sampling.sampler import SamplingVirtial
import json
from scipy.special import gammainc
from jarvishep.sample import Sample
import concurrent.futures
import itertools
from copy import deepcopy


class Grid(SamplingVirtial):
    def __init__(self) -> None:
        super().__init__()
        self.load_schema_file()
        self.method = "Grid"
        self._P     = None
        self._index = 0 
        self.tasks  = set()
        self.info   = {}
        self.future_to_sample = {}
        self._runtime_pending_samples = []

    def load_schema_file(self):
        self.schema = self.path['GridSchema']

    def set_config(self, config_info) -> None:
        self.config = config_info
        self.set_bucket_alloc()
        self.init_generator()

    def supports_runtime_checkpointing(self) -> bool:
        return True

    def __iter__(self):
        if self._P is None:
            self.initialize()  # ensure the _P is generated before iteration 
        self._index = 0  # Ensure the index starting from 0
        return self

    def __next__(self):
        if self._P is None or not isinstance(self._P, np.ndarray):
            raise StopIteration  # Stop iteration, if _P is not defined or _P is not np.array
        while self._index < len(self._P):
            result = self._P[self._index]
            result = self.map_point_into_distribution(result)
            self._index += 1
            if self._selectionexp and not self.evaluate_selection(self._selectionexp, result):
                continue
            return result
        raise StopIteration

    def next_sample(self):
        return self.__next__()

    def map_point_into_distribution(self, row) -> np.ndarray:
        result = {}
        for ii in range(len(row)):
            result[self.vars[ii].name] = self.vars[ii].map_standard_random_to_distribution(row[ii])
        return result

    def set_logger(self, logger) -> None:
        super().set_logger(logger)
        self.logger.warning("Sampling method initializaing ...")

    def init_generator(self):
        self.load_variable()
        self._dimensions = len(self.vars)
        if "selection" in self.config["Sampling"]:
            self._selectionexp = self.config["Sampling"]["selection"]

    def initialize(self):
        self.logger.warning("Initializing the Grid Sampling")
        try:
            t0 = time.time()
            dims = np.array([var.parameters["num"] for var in self.vars])
            self._P = grid_sampling(dimensions=dims)
            self.info["NSamples"] = self._P.shape[0]
            self.info["t0"]       = time.time() - t0 
            self._index           = 0
            self.logger.info("Grid Sampler obtains {} samples in {:.2f} sec".format(self.info['NSamples'], self.info['t0']))
        except:
            self.logger.error("Grid Sampler meets error when trying scan the parameter space.")
            sys.exit(2)

    def _export_sampler_state(self):
        return {
            "grid_points": self._P,
            "index": int(self._index),
            "selectionexp": self._selectionexp,
            "info": deepcopy(self.info),
            "bucket_allocator": None if self.bucket_alloc is None else self.bucket_alloc.get_state(),
            "pending_samples": self._collect_pending_sample_infos(),
        }

    def _import_sampler_state(self, payload):
        self._P = payload.get("grid_points", self._P)
        self._index = int(payload.get("index", self._index or 0))
        self._selectionexp = payload.get("selectionexp", self._selectionexp)
        self.info = deepcopy(payload.get("info", self.info))
        bucket_state = payload.get("bucket_allocator")
        if bucket_state and self.bucket_alloc is not None:
            self.bucket_alloc.set_state(dict(bucket_state))
        self._runtime_pending_samples = list(payload.get("pending_samples", []))

    def run_nested(self):
        total_cores = os.cpu_count() or 1
        self.tasks = set()
        self.future_to_sample = {}
        pending_samples = list(getattr(self, "_runtime_pending_samples", []) or [])
        self._runtime_pending_samples = []
        exhausted = False
        base_sample_cfg = self.info['sample']

        for sample_info in pending_samples:
            sample = self._rebuild_sample_from_info(sample_info)
            try:
                future = self.factory.submit_task(sample.info)
            except Exception:
                self._on_sample_completed(sample.info)
                sample.close()
                raise
            self.tasks.add(future)
            self.future_to_sample[future] = sample

        while (not exhausted) or self.tasks:
            while not exhausted and len(self.tasks) < total_cores:
                try: 
                    param = self.next_sample()
                except StopIteration:
                    exhausted = True
                    break

                sample = Sample(param)
                sconfig = self.build_sample_config(
                    base_sample_cfg,
                    save_dir=self._next_bucket_dir_for_sample(),
                )
                sample.set_config(sconfig)
                try:
                    future = self.factory.submit_task(sample.info)
                except Exception:
                    self._on_sample_completed(sample.info)
                    sample.close()
                    raise
                self.tasks.add(future)
                self.future_to_sample[future] = sample

            if not self.tasks:
                continue

            done, _ = concurrent.futures.wait(
                self.tasks,
                return_when=concurrent.futures.FIRST_COMPLETED,
            )
            self.tasks.difference_update(done)
            for future in done:
                sample = self.future_to_sample.pop(future, None)
                try:
                    future.result()
                except Exception as exc:
                    suuid = sample.uuid if sample else "UNKNOWN"
                    self.logger.error(
                        format_two_column_log(
                            "[WorkerFactory] future exception consumed",
                            [("uuid", suuid), ("error", exc)],
                        )
                    )
                    raise
                finally:
                    if sample is not None:
                        self._on_sample_completed(sample.info)
                        sample.close()

    def finalize(self):
        pass


    def set_factory(self, factory) -> None:
        self.factory = factory
        self.logger.warning("WorkerFactory is ready for Grid sampler")

def grid_sampling(dimensions):
    """
    Generate grid samples based on dimensions and dimension-specific number of steps.

    Args:
        dimensions (list of tuple): A list of tuples where each tuple specifies the number of step for each dimension.

    Returns:
        numpy.ndarray: A 2D array where each row represents a grid sample.
    """

    # Generate grid points for each dimension based on the number of steps.
    # Clamp endpoints to open interval (eps, 1-eps) to avoid non-finite
    # values in ppf/logit style transforms at exact 0/1.
    eps = np.finfo(np.float64).eps
    grid_ranges = []
    for steps in dimensions:
        nsteps = int(steps)
        if nsteps <= 1:
            grid_ranges.append(np.array([0.5], dtype=np.float64))
            continue
        vals = np.linspace(0.0, 1.0, nsteps, dtype=np.float64)
        grid_ranges.append(np.clip(vals, eps, 1.0 - eps))
    # Generate the Cartesian product of all grid points
    grid_samples = np.array(list(itertools.product(*grid_ranges)))
    return grid_samples
