#!/usr/bin/env python3 
import os, sys
from re import S 
from abc import ABCMeta, abstractmethod
# from mpmath.functions.functions import re
import numpy as np
import time 
from jarvishep import benchmark
from jarvishep.log_kv import format_two_column_log
from jarvishep.Sampling.sampler import SamplingVirtial, BoolConversionError
import json
from jarvishep.sample import Sample
import concurrent.futures
import sympy as sp 
from copy import deepcopy

# from _pytest.mark import param


class RandomS(SamplingVirtial):
    def __init__(self) -> None:
        super().__init__()
        self.load_schema_file()
        self.method = "Random"
        self._P     = None
        self._index = None
        self.tasks  = set()
        self.info   = {}
        self._selectionexp = None
        self.future_to_sample = {}
        self._runtime_pending_samples = []

    def load_schema_file(self):
        self.schema = self.path['RandomSchema']

    def set_config(self, config_info) -> None:
        self.config = config_info
        self.set_bucket_alloc()
        self.init_generator()

    def supports_runtime_checkpointing(self) -> bool:
        return True

    def __iter__(self):
        if self._index is None:
            self.initialize()  # ensure the _P is generated before iteration 
        self._index = 0  # Ensure the index starting from 0
        return self

    def __next__(self):# Stop iteration, if _P is not defined or _P is not np.array
        if self._index < self._maxp:
            if self._selectionexp:
                is_selection = False
                while not is_selection:
                    temp    = np.random.rand(self._dimensions)
                    param   = self.map_point_into_distribution(temp)
                    if self._selectionexp: 
                        is_selection = self.evaluate_selection(self._selectionexp, param)
                self._index += 1
                return param
            else: 
                temp = np.random.random(self._dimensions)
                param = self.map_point_into_distribution(temp)
                self._index += 1 
                return param 
        else:
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
        self._maxp  = int(self.config['Sampling']['Point number'])
        if "selection" in self.config["Sampling"]:
            self._selectionexp = self.config["Sampling"]["selection"]

    def initialize(self):
        self.logger.warning("Initializing the Random Sampling")
        self._index = 0 
        if self._selectionexp:
            try:
                self.info["t0"]       = time.time() 
                temp    = np.random.rand(self._dimensions)
                param   = self.map_point_into_distribution(temp)
                self.evaluate_selection(self._selectionexp, param)
            except BoolConversionError:
                self.logger.error("Wrong selection condition in input YAML -> \n\t{}".format(self._selectionexp))
                sys.exit(2)
            except:
                self.logger.error("Random Sampler meets error when trying scan the parameter space.")
                sys.exit(2)
        # else: 
        #     try: 
        #         self.info['to']     = time.time() 

    def _export_sampler_state(self):
        return {
            "index": int(self._index if self._index is not None else 0),
            "maxp": int(self._maxp),
            "selectionexp": self._selectionexp,
            "numpy_random_state": np.random.get_state(),
            "info": deepcopy(self.info),
            "bucket_allocator": None if self.bucket_alloc is None else self.bucket_alloc.get_state(),
            "pending_samples": self._collect_pending_sample_infos(),
        }

    def _import_sampler_state(self, payload):
        self._index = int(payload.get("index", self._index or 0))
        self._maxp = int(payload.get("maxp", self._maxp or 0))
        self._selectionexp = payload.get("selectionexp", self._selectionexp)
        np_state = payload.get("numpy_random_state")
        if np_state is not None:
            np.random.set_state(np_state)
        self.info = deepcopy(payload.get("info", self.info))
        bucket_state = payload.get("bucket_allocator")
        if bucket_state and self.bucket_alloc is not None:
            self.bucket_alloc.set_state(dict(bucket_state))
        self._runtime_pending_samples = list(payload.get("pending_samples", []))
                

    def run_nested(self):
        timing_enabled = benchmark.TIMING_ENABLED
        timing_start = benchmark.monotonic_seconds() if timing_enabled else None
        try:
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
        finally:
            if timing_enabled and timing_start is not None:
                benchmark.record_stage(
                    "sampler_loop",
                    benchmark.monotonic_seconds() - timing_start,
                )

    def finalize(self):
        pass

    def set_factory(self, factory) -> None:
        self.factory = factory
        self.logger.warning("WorkerFactory is ready for Random sampler")
