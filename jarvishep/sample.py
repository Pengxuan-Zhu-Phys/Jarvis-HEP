#!/usr/bin/env python3

import os, sys 
import json
import numpy as np 
import sympy as sp 
from copy import deepcopy 
from time import sleep
from jarvishep.base import Base
from jarvishep import benchmark
from jarvishep.sample_logger import SampleLogger
from jarvishep.runtime_config import should_eager_materialize, should_materialize_on_failure
from uuid import uuid4
from jarvishep.inner_func import update_const, update_funcs


class _NullSampleLogger:
    """Discard sample logs until artifacts are materialized."""

    _options = (None, None, None, None, None, None, None, None, {})

    def bind(self, **_extra):
        return self

    def log(self, *_args, **_kwargs):
        return None

    def debug(self, *_args, **_kwargs):
        return None

    def info(self, *_args, **_kwargs):
        return None

    def warning(self, *_args, **_kwargs):
        return None

    def error(self, *_args, **_kwargs):
        return None

    def critical(self, *_args, **_kwargs):
        return None

    def close(self):
        return None


class Sample(Base):
    def __init__(self, params):
        self._params = dict(params)
        self._uuid = str(uuid4())
        self._likelihood = None  # Initialize likelihood with None
        self.processed = False
        self.observables = dict(params)
        self.observables['uuid'] = self.uuid
        self.logger = None 
        self.handlers = {}
        self._with_nuisance = False
        self._nuisance_status   = False
        self._u = None
        self._materialized = False

    @property
    def uuid(self):
        return self._uuid  

    @property
    def u(self):
        return self._u 

    @property
    def params(self):
        return self._params  

    @property
    def likelihood(self):
        return self._likelihood

    def update_uuid(self, uuid):
        self._uuid = uuid 
        self.observables['uuid'] = self.uuid

    @likelihood.setter
    def likelihood(self, value):
        self._likelihood = value  # 允许更新likelihood值

    @staticmethod
    def _resolve_sample_root(info: dict) -> str:
        bucket_parent = info.get("save_dir")
        if bucket_parent:
            return os.path.dirname(str(bucket_parent).rstrip(os.sep))

        sample_dirs = info.get("sample_dirs")
        if sample_dirs:
            return str(sample_dirs)

        task_root = info.get("task_result_dir", os.getcwd())
        return os.path.join(task_root, "SAMPLE")

    def set_config(self, config):
        timing_enabled = benchmark.TIMING_ENABLED
        timing_start = benchmark.monotonic_seconds() if timing_enabled else None
        try:
            self.info = config
            self.create_info()
            if self.info.get("nuisance", {}):
                self.combine_nuisance_card()
                self._with_nuisance = True
            if should_eager_materialize(self.info):
                self.materialize()
        finally:
            if timing_enabled and timing_start is not None:
                benchmark.record_stage(
                    "sample_setup",
                    benchmark.monotonic_seconds() - timing_start,
                )

    def create_info(self):
        bucket_parent = self.info.get("save_dir")
        if not bucket_parent:
            bucket_parent = self.info.get("sample_dirs")
        sample_root = self._resolve_sample_root(self.info)

        self.info.update({
            "uuid": self.uuid,
            "params": self.params,
            "observables": self.observables,
            "sample_dirs": sample_root,
            "save_dir": None,
            "run_log": None,
            "logger_name": f"Sample@{self.uuid}",
            "logger": _NullSampleLogger(),
            "handlers": self.handlers,
            "status": "Init",
            "_materialized": False,
            "_bucket_parent": bucket_parent if isinstance(bucket_parent, str) else None,
        })
        self.logger = self.info["logger"]

    def combine_nuisance_card(self):
        card = self.info['nuisance']
        uuid = "{}@{}".format(self.uuid, card['NAttempt'])
        params = card['active']['param']
        params.update({"uuid": uuid})
        self.info['params'].update(params)
        self.info['observables'] = self.info['params']
        self.info['NAttempt'] = card['NAttempt']
    
    def gather_nuisance(self):
        self.info['observables'].update({"uuid": self.uuid})
    
    def materialize(self, *, bucket_parent: str | None = None, minimal_log_message: str | None = None) -> str:
        if self._materialized:
            return str(self.info.get("save_dir"))

        if bucket_parent is None:
            bucket_parent = self.info.get("_bucket_parent")
        if bucket_parent is None:
            bucket_parent = self._resolve_sample_root(self.info)
            os.makedirs(bucket_parent, exist_ok=True)

        save_dir = os.path.join(str(bucket_parent), self.uuid)
        run_log = os.path.join(save_dir, "Sample_running.log")

        os.makedirs(save_dir, exist_ok=True)
        self.info["save_dir"] = save_dir
        self.info["run_log"] = run_log
        self.info["_materialized"] = True
        self._materialized = True

        self._open_sample_logger()
        if minimal_log_message:
            self.logger.info(minimal_log_message)
        return save_dir

    def _open_sample_logger(self):
        logger_name = f"Sample@{self.info['uuid']}"
        self.info['logger_name'] = logger_name

        self.logger = SampleLogger.open(
            self.info['run_log'],
            module=logger_name,
            extra={
                "to_console": True,
                "Jarvis": True,
                "_log_domain": "jarvis_hep",
            },
        )
        self.logger.info("Sample created into the Disk")
        self.info['logger'] = self.logger

    def set_logger(self):
        """Backward-compatible entry point used by checkpoint rebuild paths."""
        if self._materialized:
            self._open_sample_logger()
        else:
            self.info["logger"] = _NullSampleLogger()
            self.logger = self.info["logger"]
    
    def custom_format(self, record):
        module = record["extra"].get("module", "No module")
        if "raw" in record["extra"]:
            return "{message}"
        else:
            return f"\n·•· <cyan>{module}</cyan> \n\t-> <green>{record['time']:MM-DD HH:mm:ss.SSS}</green> - [<level>{record['level']}</level>] >>> \n<level>{{message}}</level>"

    
    def manage_directories(self, base_path):
        return super().manage_directories(base_path)
    
    def evaluate_output(self, outputs):
        custom_functions = update_funcs({})
        cunsom_constants = update_const({})
        result = {}
        for op in outputs: 
            if op in self.info['observables']:
                result[op] = self.info['observables'][op]
            else: 
                try:
                    symbols = {var: sp.symbols(var) for var in self.info['observables']}
                    locals_context = {**custom_functions, **cunsom_constants, **symbols}
                    expr = sp.sympify(op, locals=locals_context)
                    res = expr.subs(self.info['observables'])
                    result[op] = res
                except: 
                    raise ValueError
        return result

    def start(self):
        self.logger.info("Sample -> {} is ready for submittion".format(self.uuid))
        self.info['status'] = "Running"
        if self._with_nuisance: 
            self.info['nuisance']['status'] = "Running"
            self.logger.info("{}\nSample start {}-th nuisance attempt".format(">"*60, self.info['nuisance']['NAttempt']))

    def close(self):
        if getattr(self, "logger", None) is not None:
            self.logger.info(self._build_close_message())
        self.close_logger() 

    def _build_close_message(self) -> str:
        observables = self.info.get("observables", {}) if isinstance(getattr(self, "info", None), dict) else {}
        if not isinstance(observables, dict) or not observables:
            return "Sample closed"

        from jarvishep.Module.likelihood import LogLikelihood

        return (
            "Sample SUMMARY\n"
            "============================================================================\n"
            f"{LogLikelihood.format_summary(observables)}\n"
            "============================================================================\n"
            "Sample closed"
        )
        
    def close_logger(self):
        """Close per-sample logger handler safely.

        This removes the per-sample file sink (if registered) and clears self.logger.
        It is safe to call multiple times.
        """
        if getattr(self, 'logger', None) is not None:
            try:
                self.logger.close()
            except Exception:
                pass

        self.handlers = {}
        self.logger = None
        
        # keep info serializable; do not store logger objects there
        if hasattr(self, 'info') and isinstance(self.info, dict):
            if self.info.get('logger') is not None:
                self.info['logger'] = None

    @staticmethod
    def format_summary(values: dict,
                       key_width: int | None = None,
                       value_width: int = 60,
                       float_precision: int = 6) -> str:
        """Fast, pandas-free summary formatter.

        Produces an aligned two-column listing similar to `pandas.Series(values).to_string()`.
        - Keys are left-aligned.
        - Values are right-aligned.
        - Long strings are truncated with '...'.
        - Floats are formatted with up to `float_precision` decimal places, trimming trailing zeros.
        """

        def _format_value(v):
            # Keep None explicit
            if v is None:
                s = "None"
            # Numpy scalars -> python scalars
            elif isinstance(v, (np.generic,)):
                s = str(v.item())
            elif isinstance(v, bool):
                s = "True" if v else "False"
            elif isinstance(v, (int,)):
                s = str(v)
            elif isinstance(v, (float,)):
                if np.isnan(v):
                    s = "nan"
                elif np.isposinf(v):
                    s = "inf"
                elif np.isneginf(v):
                    s = "-inf"
                else:
                    av = abs(v)
                    # Use scientific notation for very small/large magnitudes to avoid rounding to 0
                    if (av != 0.0 and av < 1e-4) or av >= 1e6:
                        s = f"{v:.{float_precision}e}"
                    else:
                        s = f"{v:.{float_precision}f}".rstrip("0").rstrip(".")
            else:
                s = str(v)

            # Truncate overly long strings similar to pandas display
            if len(s) > value_width:
                if value_width <= 3:
                    s = s[:value_width]
                else:
                    s = s[: value_width - 3] + "..."
            return s

        # Preserve insertion order of `values` (dict preserves insertion order in Python 3.7+)
        keys = list(values.keys())
        if key_width is None:
            key_width = max((len(str(k)) for k in keys), default=1)
            # Keep it reasonable; very long keys reduce readability
            key_width = min(key_width, 60)

        lines = []
        for k in keys:
            ks = str(k)
            if len(ks) > key_width:
                # Truncate long keys to keep alignment stable
                if key_width <= 3:
                    ks = ks[:key_width]
                else:
                    ks = ks[: key_width - 3] + "..."
            vs = _format_value(values[k])
            # Match the look: key left, value right with plenty of spacing
            lines.append(f"{ks:<{key_width}}  {vs:>{value_width}}")
        return "\n".join(lines)

                
    def record(self):
        if not self._with_nuisance:
            return

        msg = (
            "Nuisance SUMMARY  -> {}-th\t attempt  \n".format(self.info.get("NAttempt", {}))
            + "=================================================================\n"
            + self.format_summary(self.info.get("observables", {}))
            + "\n================================================================="
        )

        self.logger.info(msg)


def ensure_sample_materialized(sample_info: dict, *, minimal_log_message: str | None = None) -> str | None:
    """Materialize sample artifacts on demand (e.g. @Sdir resolution)."""
    if not isinstance(sample_info, dict):
        return None
    if str(sample_info.get("sample_artifacts", "auto")).strip().lower() == "never":
        return None
    if sample_info.get("_materialized") and sample_info.get("save_dir"):
        return str(sample_info["save_dir"])

    seed_params = sample_info.get("params")
    if not isinstance(seed_params, dict) or not seed_params:
        seed_params = sample_info.get("observables", {})
    sample = Sample(seed_params if isinstance(seed_params, dict) else {})
    sample.info = sample_info
    if sample_info.get("uuid"):
        sample.update_uuid(str(sample_info["uuid"]))
    sample._materialized = bool(sample_info.get("_materialized"))
    sample.logger = sample_info.get("logger") or _NullSampleLogger()
    return sample.materialize(minimal_log_message=minimal_log_message)


def materialize_failure_artifacts(sample_info: dict, *, error: Exception | str | None = None) -> str | None:
    """Create per-sample artifacts for a failed sample (WP-1.1 minimal replay)."""
    if not isinstance(sample_info, dict):
        return None
    if not should_materialize_on_failure(sample_info):
        return None
    if sample_info.get("_materialized") and sample_info.get("save_dir"):
        return str(sample_info["save_dir"])

    message = "Sample failed"
    if error is not None:
        message = f"Sample failed -> {error}"

    return ensure_sample_materialized(sample_info, minimal_log_message=message)