#!/usr/bin/env python3

import os
import threading
import inspect
import json
import pickle
from copy import deepcopy
from datetime import datetime, timezone
from uuid import uuid4

import numpy as np
import pandas as pd
from prettytable import PrettyTable

from jarvishep.log_kv import format_two_column_log
from jarvishep.Sampling.nested_checkpoint_bridge import NestedLikelihoodBridge
from jarvishep.Sampling.Source.Dynesty.py.dynesty.pool import JarvisFactoryAsyncPool
from jarvishep.Sampling.sampler import SamplingVirtial
from jarvishep.sample import Sample


def prior_transform(u):
    uid = str(uuid4())
    ret = np.append(u, [uid])
    return ret


MultiNestFactoryPool = JarvisFactoryAsyncPool


def _linear_interp_extrapolate(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    def interp(values):
        arr = np.asarray(values, dtype=float)
        result = np.interp(arr, x, y)
        left = arr < x[0]
        if np.any(left):
            slope = (y[1] - y[0]) / (x[1] - x[0])
            result = np.asarray(result, dtype=float)
            result[left] = y[0] + (arr[left] - x[0]) * slope
        right = arr > x[-1]
        if np.any(right):
            slope = (y[-1] - y[-2]) / (x[-1] - x[-2])
            result = np.asarray(result, dtype=float)
            result[right] = y[-1] + (arr[right] - x[-1]) * slope
        if np.isscalar(values):
            return float(np.asarray(result))
        return result

    return interp


class MultiNest(SamplingVirtial):
    """Static nested sampling wrapper (dynesty NestedSampler)."""

    def __init__(self) -> None:
        super().__init__()
        self.load_schema_file()
        self.method = "MultiNest"
        self.sampler = None
        self.max_workers = os.cpu_count() or 1
        self._multinest_pool = None
        self._multinest_workers = 1
        self._factory_submit_limit = 1
        self._submit_gate = threading.BoundedSemaphore(value=1)
        self._execution_profile = {}
        self.df = None
        self.lnX_from_LogLike = None
        self._sampler_results_snapshot = None
        self._native_sampler_loaded = False

    def supports_runtime_checkpointing(self) -> bool:
        return True

    def load_schema_file(self):
        self.schema = self.path.get("MultiNestSchema", self.path.get("DynestySchema"))

    def set_config(self, config_info) -> None:
        self.config = config_info
        self.set_bucket_alloc()
        self.init_generator()

    def set_factory(self, factory) -> None:
        self.factory = factory
        self.logger.warning("WorkerFactory is ready for MultiNest sampler")

    def set_execution_profile(self, multinest_workers=None, max_pending_factory=None):
        profile = {}
        if multinest_workers is not None:
            profile["multinest_workers"] = int(multinest_workers)
        if max_pending_factory is not None:
            profile["max_pending_factory"] = int(max_pending_factory)
        self._execution_profile = profile

    @staticmethod
    def _sanitize_nested(value):
        if isinstance(value, dict):
            return {
                key: MultiNest._sanitize_nested(item)
                for key, item in value.items()
                if key not in {"logger", "handlers"}
            }
        if isinstance(value, list):
            return [MultiNest._sanitize_nested(item) for item in value]
        if isinstance(value, tuple):
            return tuple(MultiNest._sanitize_nested(item) for item in value)
        return value

    @staticmethod
    def _pickleable_or_none(value):
        if value is None:
            return None
        try:
            pickle.dumps(value, protocol=pickle.HIGHEST_PROTOCOL)
        except Exception:
            return None
        return value

    def _export_sampler_state(self):
        native_sampler = self._pickleable_or_none(self.sampler if self._native_sampler_loaded or self.sampler is not None else None)
        return {
            "dimensions": int(getattr(self, "_dimensions", 0) or 0),
            "nlive": int(getattr(self, "_nlive", 0) or 0),
            "rstate": getattr(self._rstate, "bit_generator", None).state if getattr(self, "_rstate", None) is not None else None,
            "runnested": dict(getattr(self, "_runnested", {}) or {}),
            "selectionexp": self._selectionexp,
            "execution_profile": dict(self._execution_profile),
            "native_sampler": native_sampler,
            "sampler_results": None if self.sampler is None else getattr(self.sampler, "results", None),
            "df": None if self.df is None else self.df.copy(deep=True),
            "info": self._sanitize_nested(deepcopy(self.info)) if isinstance(self.info, dict) else {},
        }

    def _import_sampler_state(self, payload):
        self._dimensions = int(payload.get("dimensions", self._dimensions or 0))
        self._nlive = int(payload.get("nlive", self._nlive or 0))
        rstate = payload.get("rstate")
        if rstate is not None:
            self._rstate = np.random.default_rng()
            self._rstate.bit_generator.state = rstate
        self._runnested = dict(payload.get("runnested", self._runnested or {}))
        self._selectionexp = payload.get("selectionexp", self._selectionexp)
        self._execution_profile = dict(payload.get("execution_profile", self._execution_profile or {}))
        native_sampler = payload.get("native_sampler", None)
        if native_sampler is not None:
            self.sampler = native_sampler
            self._native_sampler_loaded = True
        self._sampler_results_snapshot = payload.get("sampler_results", self._sampler_results_snapshot)
        self.df = payload.get("df", self.df)
        self.info = deepcopy(payload.get("info", self.info))

    def _ensure_sampler_results_proxy(self):
        if self.sampler is not None:
            return
        if self._sampler_results_snapshot is None:
            return
        proxy = type("_MultiNestCheckpointResultsProxy", (), {})()
        proxy.results = self._sampler_results_snapshot
        self.sampler = proxy

    def _build_likelihood_bridge(self) -> NestedLikelihoodBridge:
        sample_cfg = deepcopy(self.info.get("sample", {})) if isinstance(self.info, dict) else {}
        bucket_state = self.bucket_alloc.get_state() if getattr(self, "bucket_alloc", None) is not None else None
        submit_limit = max(1, int(getattr(self, "_factory_submit_limit", self._multinest_workers or 1) or 1))
        limit = 200
        width = 6
        if getattr(self, "bucket_alloc", None) is not None:
            try:
                state = self.bucket_alloc.get_state()
                limit = int(state.get("limit", limit))
                width = int(state.get("width", width))
            except Exception:
                pass
        bridge = NestedLikelihoodBridge(
            sampler_name=self.method,
            variables=self.vars,
            base_sample_cfg=sample_cfg,
            sample_cls=Sample,
            bucket_state=bucket_state,
            bucket_limit=limit,
            bucket_width=width,
            submit_limit=submit_limit,
            selection_expression=self._selectionexp,
        )
        bridge.attach_runtime(factory=self.factory, logger=self.logger)
        return bridge

    def _extract_bridge(self):
        candidates = [
            getattr(self.sampler, "loglikelihood", None),
            getattr(getattr(self.sampler, "loglikelihood", None), "loglikelihood", None),
            getattr(getattr(getattr(self.sampler, "loglikelihood", None), "loglikelihood", None), "func", None),
        ]
        for candidate in candidates:
            if candidate is None:
                continue
            if hasattr(candidate, "attach_runtime"):
                return candidate
            bridge = getattr(candidate, "func", None)
            if hasattr(bridge, "attach_runtime"):
                return bridge
        return None

    def _reattach_native_sampler_runtime(self) -> None:
        if self.sampler is None:
            return
        bridge = self._extract_bridge()
        if bridge is not None:
            bridge.attach_runtime(factory=self.factory, logger=self.logger)
        pool = self._multinest_pool
        if pool is None:
            return
        samplers = [self.sampler]
        inner = getattr(self.sampler, "sampler", None)
        if inner is not None:
            samplers.append(inner)
        for cursamp in samplers:
            try:
                cursamp.M = pool.map
                cursamp.pool = pool
                if hasattr(cursamp, "loglikelihood"):
                    cursamp.loglikelihood.pool = pool
            except Exception:
                pass

    def __next__(self):
        while True:
            u = np.random.random(self._dimensions)
            param = self.map_point_into_distribution(u)
            if self._selectionexp and not self.evaluate_selection(self._selectionexp, param):
                continue
            return param

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
        bounds = self.config["Sampling"]["Bounds"]
        self._nlive = int(bounds["nlive"])
        self._rstate = np.random.default_rng(bounds["rseed"])
        self._runnested = bounds["run_nested"]
        self._dimensions = len(self.vars)
        self._selectionexp = self.config["Sampling"].get("selection")

    def initialize(self):
        self.logger.warning("Initializing the MultiNest Sampling")

    def init_sampler_db(self):
        self.info["db"] = {
            "nested result": os.path.join(
                self.info["sample"]["task_result_dir"],
                "DATABASE",
                "multinest_result.csv",
            ),
            "summary": os.path.join(
                self.info["sample"]["task_result_dir"],
                "DATABASE",
                "multinest_summary.json",
            ),
        }

    def _shutdown_multinest_pool(self):
        if self._multinest_pool is not None:
            self._multinest_pool.shutdown(wait_for_tasks=True, cancel_futures=True)
            self._multinest_pool = None

    def _ensure_bucket_allocator(self):
        if getattr(self, "bucket_alloc", None) is not None:
            return

        sample_cfg = self.info.get("sample", {}) if isinstance(self.info, dict) else {}
        base_path = sample_cfg.get("sample_dirs")
        if not base_path:
            task_root = sample_cfg.get("task_result_dir", os.getcwd())
            base_path = os.path.join(task_root, "SAMPLE")
            sample_cfg["sample_dirs"] = base_path

        limit = 200
        width = 6
        cfg = getattr(self, "config", None)
        if isinstance(cfg, dict):
            scan_cfg = cfg.get("Scan", {})
            directory_cfg = None
            if isinstance(scan_cfg, dict):
                directory_cfg = scan_cfg.get("sample_directory", None)
            if not isinstance(directory_cfg, dict):
                directory_cfg = cfg.get("Directory_Setting", None)
            if isinstance(directory_cfg, dict):
                limit = int(directory_cfg.get("limit", limit))
                width = int(directory_cfg.get("width", width))

        from jarvishep.Sampling.bucketallocator import BucketAllocator

        self.bucket_alloc = BucketAllocator(
            base_path=base_path,
            limit=limit,
            width=width,
            start_bucket=1,
            on_bucket_sealed=self._on_bucket_sealed if self._archive_enabled() else None,
        )
        self.bucket_alloc.check_and_update()
        if self._archive_enabled():
            self._ensure_archive_manager()

    def _resolve_execution_profile(self):
        factory_workers = int(getattr(self.factory, "_max_workers", self.max_workers) or self.max_workers or 1)
        nlive = int(getattr(self, "_nlive", 1) or 1)
        base_workers = int(self.max_workers or 1)

        requested_workers = int(self._execution_profile.get("multinest_workers", base_workers))
        multinest_workers = max(1, min(requested_workers, base_workers, factory_workers, nlive))

        requested_pending = int(self._execution_profile.get("max_pending_factory", multinest_workers))
        max_pending_factory = max(1, min(requested_pending, factory_workers))

        self._multinest_workers = multinest_workers
        self._factory_submit_limit = max_pending_factory
        self._submit_gate = threading.BoundedSemaphore(value=max_pending_factory)

        if self._multinest_pool is None or self._multinest_pool.size != multinest_workers:
            self._shutdown_multinest_pool()
            self._multinest_pool = MultiNestFactoryPool(multinest_workers)

        self._log_execution_profile(
            multinest_workers=multinest_workers,
            max_pending_factory=max_pending_factory,
            factory_workers=factory_workers,
            nlive=nlive,
        )

    @staticmethod
    def _format_worker_summary(rows):
        metric_header = "Metric"
        value_header = "Value"
        metric_width = max([len(metric_header)] + [len(str(k)) for k, _ in rows])
        value_width = max([len(value_header)] + [len(str(v)) for _, v in rows])

        lines = [
            f"{metric_header:<{metric_width}}\t{value_header:<{value_width}}",
            f"{'-' * metric_width}\t{'-' * value_width}",
        ]
        for metric, value in rows:
            lines.append(f"{str(metric):<{metric_width}}\t{str(value):<{value_width}}")
        return "\n\t" + "\n\t".join(lines)

    def _log_execution_profile(self, multinest_workers, max_pending_factory, factory_workers, nlive):
        table = self._format_worker_summary(
            [
                ("multinest_workers", multinest_workers),
                ("factory_submit_limit", max_pending_factory),
                ("factory_workers", factory_workers),
                ("nlive", nlive),
            ]
        )
        self.logger.warning(f"MultiNest Worker Summary ->{table}")

    def _submit_with_backpressure(self, sample):
        if not self._submit_gate.acquire(blocking=False):
            self.logger.info(
                format_two_column_log(
                    "MultiNest factory submit window full; waiting",
                    [
                        ("factory_submit_limit", self._factory_submit_limit),
                        ("uuid", sample.uuid),
                    ],
                )
            )
            self._submit_gate.acquire()
        try:
            future = self.factory.submit_task(sample.info)
            return future.result()
        finally:
            self._submit_gate.release()

    def run_nested(self):
        self._ensure_bucket_allocator()
        self._resolve_execution_profile()
        self.init_sampler_db()

        try:
            if self._native_sampler_loaded and self.sampler is not None:
                self._reattach_native_sampler_runtime()
                run_kwargs = dict(self._runnested)
                try:
                    sig = inspect.signature(self.sampler.run_nested)
                    if "resume" in sig.parameters:
                        run_kwargs.setdefault("resume", True)
                except (TypeError, ValueError):
                    run_kwargs.setdefault("resume", True)
                self.sampler.run_nested(**run_kwargs)
                return

            bridge = self._build_likelihood_bridge()

            from jarvishep.Sampling.Source.Dynesty.py.dynesty import NestedSampler
            sampler_kwargs = {
                "loglikelihood": bridge,
                "prior_transform": prior_transform,
                "ndim": self._dimensions,
                "nlive": self._nlive,
                "pool": self._multinest_pool,
                "rstate": self._rstate,
                "queue_size": self._multinest_workers,
            }
            try:
                sig = inspect.signature(NestedSampler)
                if "log_file_path" in sig.parameters:
                    sampler_kwargs["log_file_path"] = self.info["logfile"]
                if "inner_logger" in sig.parameters:
                    sampler_kwargs["inner_logger"] = self.logger
            except (TypeError, ValueError):
                # Fallback for objects without introspectable signature.
                pass

            self.sampler = NestedSampler(**sampler_kwargs)
            self.sampler.logger = self.logger
            self._native_sampler_loaded = True
            self._reattach_native_sampler_runtime()
            self.sampler.run_nested(**self._runnested)
        finally:
            self._shutdown_multinest_pool()

    def _result_field(self, name, default=None):
        results = self.sampler.results
        try:
            return results[name]
        except Exception:
            return default

    def _build_multinest_table(self):
        samples = np.asarray(self._result_field("samples", np.empty((0, self._dimensions))))
        n = int(samples.shape[0]) if samples.ndim == 2 else 0

        uu = self._result_field("samples_u", None)
        if uu is None:
            uu = np.zeros_like(samples)
        uu = np.asarray(uu)

        uid = self._result_field("samples_uid", None)
        if uid is None or len(uid) != n:
            uid = np.array([f"multinest-{i:08d}" for i in range(n)], dtype=object)

        logwt = np.asarray(self._result_field("logwt", np.zeros(n)))
        logl = np.asarray(self._result_field("logl", np.zeros(n)))
        logvol = np.asarray(self._result_field("logvol", np.zeros(n)))
        logz = np.asarray(self._result_field("logz", np.zeros(n)))
        logzerr = np.asarray(self._result_field("logzerr", np.zeros(n)))
        ncall = np.asarray(self._result_field("ncall", np.zeros(n)))
        samples_n = np.asarray(self._result_field("samples_n", np.full(n, self._nlive)))
        samples_it = np.asarray(self._result_field("samples_it", np.arange(n)))
        samples_id = np.asarray(self._result_field("samples_id", np.arange(n)))
        information = self._result_field("information", 0.0)
        if np.isscalar(information):
            information = np.full(n, float(information))
        else:
            information = np.asarray(information)

        data = {
            "uuid": uid,
            "log_weight": logwt,
            "log_Like": logl,
            "log_PriorVolume": logvol,
            "log_Evidence": logz,
            "log_Evidence_err": logzerr,
            "samples_nlive": samples_n,
            "ncall": ncall,
            "samples_it": samples_it,
            "samples_id": samples_id,
            "information": information,
        }
        for ii in range(samples.shape[1] if samples.ndim == 2 else 0):
            data[f"samples_v[{ii}]"] = samples[:, ii]
        for ii in range(uu.shape[1] if uu.ndim == 2 else 0):
            data[f"samples_u[{ii}]"] = uu[:, ii]
        return data

    def save_multinest_results_to_csv(self):
        data = self._build_multinest_table()
        self.df = pd.DataFrame(data)
        out_csv = self.info["db"]["nested result"]
        out_dir = os.path.dirname(out_csv)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        self.df.to_csv(out_csv, index=False)
        self.logger.warning(f"Results saved to {self.info['db']['nested result']}")

        self._build_logl_to_logx_interpolator()

    @staticmethod
    def _as_scalar_number(value, default=np.nan):
        if value is None:
            return default
        if isinstance(value, (list, tuple)):
            if len(value) == 0:
                return default
            return MultiNest._as_scalar_number(value[-1], default=default)
        if isinstance(value, np.ndarray):
            if value.size == 0:
                return default
            return MultiNest._as_scalar_number(value.reshape(-1)[-1], default=default)
        try:
            return float(value)
        except Exception:
            return default

    @staticmethod
    def _as_int(value, default=0):
        fv = MultiNest._as_scalar_number(value, default=np.nan)
        if not np.isfinite(fv):
            return int(default)
        return int(fv)

    def _build_multinest_summary_payload(self):
        results = getattr(self.sampler, "results", None)
        nlive = self._as_int(self._result_field("nlive", None), default=getattr(self, "_nlive", 0))
        niter = self._as_int(self._result_field("niter", None), default=0)
        if niter <= 0:
            niter = int(self.df.shape[0]) if isinstance(self.df, pd.DataFrame) else 0

        ncall_field = self._result_field("ncall", None)
        if isinstance(ncall_field, np.ndarray):
            ncall = int(np.nansum(ncall_field)) if ncall_field.size else 0
        elif isinstance(ncall_field, (list, tuple)):
            ncall = int(np.nansum(np.asarray(ncall_field, dtype=float))) if len(ncall_field) else 0
        else:
            ncall = self._as_int(ncall_field, default=0)

        eff = self._as_scalar_number(self._result_field("eff", None), default=np.nan)
        if not np.isfinite(eff) and ncall > 0:
            eff = 100.0 * float(niter) / float(max(ncall, 1))

        logz = self._as_scalar_number(self._result_field("logz", None), default=np.nan)
        logzerr = self._as_scalar_number(self._result_field("logzerr", None), default=np.nan)

        payload = {
            "method": "MultiNest",
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "nlive": int(nlive),
            "niter": int(niter),
            "ncall": int(ncall),
            "eff_percent": float(eff) if np.isfinite(eff) else None,
            "logz": float(logz) if np.isfinite(logz) else None,
            "logzerr": float(logzerr) if np.isfinite(logzerr) else None,
            "result_rows": int(self.df.shape[0]) if isinstance(self.df, pd.DataFrame) else 0,
            "results_available": bool(results is not None),
        }
        return payload

    def _write_multinest_summary_json(self, payload):
        out_json = self.info.get("db", {}).get("summary")
        if not out_json:
            task_dir = self.info.get("sample", {}).get("task_result_dir", os.getcwd())
            out_json = os.path.join(task_dir, "DATABASE", "multinest_summary.json")
            self.info.setdefault("db", {})["summary"] = out_json
        out_dir = os.path.dirname(out_json)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(out_json, "w", encoding="utf-8") as f1:
            json.dump(payload, f1, indent=2)
        self.logger.warning(f"MultiNest summary json saved -> {out_json}")

    def _log_multinest_summary(self, payload):
        eff_text = "nan" if payload.get("eff_percent") is None else "{:.8f}".format(payload["eff_percent"])
        logz_text = "nan" if payload.get("logz") is None else "{:.8f}".format(payload["logz"])
        err_text = "nan" if payload.get("logzerr") is None else "{:.8f}".format(payload["logzerr"])

        table = PrettyTable()
        table.field_names = ["Metric", "Value"]
        table.align = "l"
        table.add_row(["nlive", payload.get("nlive", 0)])
        table.add_row(["niter", payload.get("niter", 0)])
        table.add_row(["ncall", payload.get("ncall", 0)])
        table.add_row(["eff(%)", eff_text])
        table.add_row(["logz", logz_text])
        table.add_row(["logzerr", err_text])

        table_str = "\n".join(f"\t{line}" for line in table.get_string().splitlines())
        self.logger.warning("MultiNest Summary ->\n{}".format(table_str))

    def _build_logl_to_logx_interpolator(self):
        self.lnX_from_LogLike = None
        if self.df is None or self.df.shape[0] < 2:
            return

        x = np.asarray(self.df.get("log_Like", []), dtype=float)
        y = np.asarray(self.df.get("log_PriorVolume", []), dtype=float)
        finite = np.isfinite(x) & np.isfinite(y)
        x = x[finite]
        y = y[finite]
        if x.size < 2:
            self.logger.warning("MultiNest interpolation skipped -> insufficient finite points")
            return

        # Ensure monotonic x and collapse duplicate logL values to avoid
        # divide-by-zero in scipy.interpolate on repeated knots.
        order = np.argsort(x)
        x = x[order]
        y = y[order]
        uniq_x, inv = np.unique(x, return_inverse=True)
        if uniq_x.size < 2:
            self.logger.warning("MultiNest interpolation skipped -> insufficient unique log_Like points")
            return

        agg_y = np.zeros(uniq_x.size, dtype=float)
        cnt = np.zeros(uniq_x.size, dtype=int)
        np.add.at(agg_y, inv, y)
        np.add.at(cnt, inv, 1)
        agg_y /= np.maximum(cnt, 1)

        self.lnX_from_LogLike = _linear_interp_extrapolate(uniq_x, agg_y)

    def finalize(self):
        self._ensure_sampler_results_proxy()
        if self.sampler is None or not hasattr(self.sampler, "results"):
            self.logger.warning("MultiNest finalize skipped -> no sampler results available")
            return
        self.save_multinest_results_to_csv()
        summary = self._build_multinest_summary_payload()
        self._log_multinest_summary(summary)
        self._write_multinest_summary_json(summary)

    def _resolve_full_dataframe_path(self, fulldf):
        candidates = []
        if isinstance(fulldf, str) and fulldf:
            candidates.append(fulldf)
            root, ext = os.path.splitext(fulldf)
            if ext.lower() == ".hdf5":
                candidates.append(f"{root}.0.csv")
                candidates.append(f"{root}.csv")

        task_dir = self.info.get("sample", {}).get("task_result_dir")
        if isinstance(task_dir, str) and task_dir:
            db_dir = os.path.join(task_dir, "DATABASE")
            db_info = os.path.join(db_dir, "running.json")
            if os.path.exists(db_info):
                try:
                    with open(db_info, "r", encoding="utf-8") as f1:
                        meta = json.load(f1)
                    converted = meta.get("converted", [])
                    if isinstance(converted, str):
                        converted = [converted]
                    if isinstance(converted, list):
                        # Prefer latest converted csv first.
                        candidates.extend(list(reversed([p for p in converted if isinstance(p, str) and p])))
                    pathroot = meta.get("pathroot")
                    active_no = meta.get("activeNO")
                    if isinstance(pathroot, str) and isinstance(active_no, int):
                        candidates.append(f"{pathroot}.{active_no}.csv")
                except Exception:
                    pass

        seen = set()
        for path in candidates:
            if not path or path in seen:
                continue
            seen.add(path)
            if os.path.exists(path):
                return path
        return None

    def combine_data(self, fulldf):
        if self.df is None:
            return
        full_df_path = self._resolve_full_dataframe_path(fulldf)
        if not full_df_path:
            self.logger.warning(
                format_two_column_log(
                    "MultiNest combine_data skipped",
                    [
                        ("reason", "no full sample table found"),
                        ("input", fulldf),
                    ],
                )
            )
            return
        df_full = pd.read_csv(full_df_path)
        merged_df = pd.merge(df_full, self.df, on="uuid", how="inner")
        merged_path = os.path.join(self.info["sample"]["task_result_dir"], "DATABASE", "multinest_full.csv")
        os.makedirs(os.path.dirname(merged_path), exist_ok=True)
        merged_df.to_csv(merged_path, index=False)
        if self.lnX_from_LogLike is not None and "LogL" in df_full.columns:
            df_full["log_PriorVolume"] = self.lnX_from_LogLike(df_full["LogL"])
