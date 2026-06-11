#!/usr/bin/env python3 
from __future__ import annotations

from jarvishep.base import Base
from abc import ABCMeta, abstractmethod
from jarvishep.Sampling.variables import Variable
import sympy as sp 
from jarvishep.inner_func import update_const, update_funcs
import numpy as np 
import sys, os, time 
import threading
from typing import Any, Dict, Tuple
from copy import deepcopy

from jarvishep.runtime_config import should_eager_materialize
from jarvishep.Sampling.Source.MCMC.runtime_checkpoint import (
    StateSaver,
    build_sampler_signature,
    build_state_payload,
    utc_now_iso,
    validate_sampler_signature,
    validate_state_payload,
)

class BoolConversionError(Exception):
    """Nothing but raise bool error"""
    pass
class SamplingVirtial(Base):
    __metaclass__ = ABCMeta
    def __init__(self) -> None:
        super().__init__()
        self.schema                     = None
        self.info: Dict[str, Any]       = {}
        self.vars: Tuple[Variable]      = None
        self.method                     = None
        self.max_workers                = 4
        self._selectionexp              = None
        self._loglike                   = None 
        self.nuisance_sampler           = None 
        self._with_nuisance             = False 
        self.bucket_alloc               = None 
        self.sample_archive_manager     = None
        self.total_core                 = os.cpu_count()
        self._runtime_checkpoint_root: str | None = None
        self._runtime_checkpoint_file: str | None = None
        self._runtime_state_saver: StateSaver | None = None
        self._runtime_checkpoint_logger = None
        self._runtime_checkpoint_interval_seconds = 30.0
        self._runtime_checkpoint_enabled = False
        self._runtime_checkpoint_auto_resume = True
        self._runtime_checkpoint_last_save_ts = 0.0
        self._runtime_checkpoint_loaded = False
        self._runtime_checkpoint_payload: Dict[str, Any] | None = None
        self._runtime_checkpoint_run_spec: Dict[str, Any] | None = None
        self._runtime_checkpoint_factory_blueprint: Dict[str, Any] | None = None
        self._runtime_checkpoint_run_spec_getter = None
        self._runtime_checkpoint_factory_blueprint_getter = None
        self._runtime_checkpoint_reason = ""
        self._runtime_checkpoint_resume_hint = False
        self._runtime_checkpoint_save_lock = threading.Lock()
        self._runtime_checkpoint_stop_event = None
        self._runtime_checkpoint_thread = None

    def supports_runtime_checkpointing(self) -> bool:
        return False

    def configure_runtime_checkpointing(
        self,
        checkpoint_root: str,
        *,
        interval_seconds: float = 30.0,
        auto_resume: bool = True,
        logger=None,
    ) -> None:
        root = os.path.abspath(str(checkpoint_root))
        os.makedirs(root, exist_ok=True)
        self._runtime_checkpoint_root = root
        self._runtime_checkpoint_file = os.path.join(root, "state.pkl")
        self._runtime_checkpoint_logger = logger
        self._runtime_checkpoint_interval_seconds = max(1.0, float(interval_seconds))
        self._runtime_checkpoint_enabled = True
        self._runtime_checkpoint_auto_resume = bool(auto_resume)
        self._runtime_state_saver = StateSaver(
            self._runtime_checkpoint_file,
            logger=self._checkpoint_logger(),
        )
        self._start_runtime_checkpoint_heartbeat()

    def _runtime_checkpoint_heartbeat_loop(self) -> None:
        stop_event = self._runtime_checkpoint_stop_event
        if stop_event is None:
            return
        while not stop_event.wait(float(self._runtime_checkpoint_interval_seconds)):
            try:
                self.persist_runtime_checkpoint(force=True, reason="checkpoint_heartbeat")
            except Exception as exc:
                self._checkpoint_log("warning", f"StateSaver heartbeat checkpoint failed -> {exc}")

    def _start_runtime_checkpoint_heartbeat(self) -> None:
        self.shutdown_runtime_checkpointing()
        if not self._runtime_checkpoint_enabled or self._runtime_state_saver is None:
            return
        self._runtime_checkpoint_stop_event = threading.Event()
        self._runtime_checkpoint_thread = threading.Thread(
            target=self._runtime_checkpoint_heartbeat_loop,
            name=f"{self.__class__.__name__}-StateSaverHeartbeat",
            daemon=True,
        )
        self._runtime_checkpoint_thread.start()

    def shutdown_runtime_checkpointing(self) -> None:
        stop_event = self._runtime_checkpoint_stop_event
        thread = self._runtime_checkpoint_thread
        if stop_event is not None:
            stop_event.set()
        if thread is not None and thread.is_alive():
            thread.join(timeout=2.0)
        self._runtime_checkpoint_stop_event = None
        self._runtime_checkpoint_thread = None

    def set_runtime_checkpoint_context(
        self,
        *,
        run_spec: Dict[str, Any] | None = None,
        factory_blueprint: Dict[str, Any] | None = None,
        run_spec_getter=None,
        factory_blueprint_getter=None,
    ) -> None:
        if run_spec is not None:
            self._runtime_checkpoint_run_spec = deepcopy(run_spec)
        if factory_blueprint is not None:
            self._runtime_checkpoint_factory_blueprint = deepcopy(factory_blueprint)
        if run_spec_getter is not None:
            self._runtime_checkpoint_run_spec_getter = run_spec_getter
        if factory_blueprint_getter is not None:
            self._runtime_checkpoint_factory_blueprint_getter = factory_blueprint_getter

    def set_runtime_checkpoint_resume_hint(self, enabled: bool) -> None:
        self._runtime_checkpoint_resume_hint = bool(enabled)

    def get_runtime_checkpoint_payload(self) -> Dict[str, Any] | None:
        return deepcopy(self._runtime_checkpoint_payload)

    def _checkpoint_logger(self):
        return self._runtime_checkpoint_logger or getattr(self, "logger", None)

    def _checkpoint_log(self, level: str, message: str) -> None:
        clogger = self._checkpoint_logger()
        if clogger is None:
            return
        method = getattr(clogger, str(level).lower(), None)
        if callable(method):
            method(message)

    def _checkpoint_seconds_until_due(self) -> float | None:
        if not self._runtime_checkpoint_enabled or self._runtime_state_saver is None:
            return None
        last_save = float(self._runtime_checkpoint_last_save_ts or 0.0)
        if last_save <= 0.0:
            return 0.0
        elapsed = time.monotonic() - last_save
        remaining = float(self._runtime_checkpoint_interval_seconds) - float(elapsed)
        return max(0.0, remaining)

    def _checkpoint_should_attempt(self) -> bool:
        remaining = self._checkpoint_seconds_until_due()
        return remaining is not None and remaining <= 0.0

    def _resolve_runtime_checkpoint_context(self, attr_name: str, getter_name: str) -> Dict[str, Any]:
        getter = getattr(self, getter_name, None)
        if callable(getter):
            try:
                value = getter()
                if isinstance(value, dict):
                    return deepcopy(value)
            except Exception:
                pass
        stored = getattr(self, attr_name, None)
        return deepcopy(stored) if isinstance(stored, dict) else {}

    def _checkpoint_safe_value(self, value: Any) -> Any:
        if value is None or isinstance(value, (str, bool, int, float)):
            return value
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, dict):
            safe: Dict[str, Any] = {}
            for key, item in value.items():
                if key in {"logger", "handlers"}:
                    continue
                safe[str(key)] = self._checkpoint_safe_value(item)
            return safe
        if isinstance(value, (list, tuple, set)):
            return [self._checkpoint_safe_value(item) for item in value]
        try:
            return deepcopy(value)
        except Exception:
            return str(value)

    def _sanitize_sample_info(self, sample_info: Dict[str, Any]) -> Dict[str, Any]:
        if not isinstance(sample_info, dict):
            return {}

        info = self._checkpoint_safe_value(sample_info)
        if not isinstance(info, dict):
            return {}

        # The runtime logger/handler objects must never enter the checkpoint.
        info["logger"] = None
        info["handlers"] = {}

        for key in ("uuid", "save_dir", "run_log", "logger_name", "status"):
            if key in sample_info and sample_info[key] is not None:
                value = sample_info[key]
                info[key] = value.item() if isinstance(value, np.generic) else str(value) if key in {"uuid", "save_dir", "run_log", "logger_name", "status"} else value
        return info

    def _collect_pending_sample_infos(self) -> list[Dict[str, Any]]:
        pending = []
        mapping = getattr(self, "future_to_sample", None)
        if isinstance(mapping, dict):
            values = mapping.values()
        else:
            tasks = getattr(self, "tasks", None)
            values = tasks.values() if isinstance(tasks, dict) else []
        for sample in values:
            info = getattr(sample, "info", None)
            if isinstance(info, dict):
                pending.append(self._sanitize_sample_info(info))
        return pending

    def _rebuild_sample_from_info(self, sample_info: Dict[str, Any]):
        from jarvishep.sample import Sample

        params = deepcopy(sample_info.get("params", {})) if isinstance(sample_info, dict) else {}
        sample = Sample(params)
        if isinstance(sample_info, dict):
            if sample_info.get("uuid"):
                sample.update_uuid(sample_info["uuid"])
            sample.info = deepcopy(sample_info)
            sample.info["logger"] = None
            sample.info["handlers"] = {}
            sample._with_nuisance = bool(sample_info.get("nuisance"))
            sample.handlers = {}
            existing_save_dir = sample_info.get("save_dir")
            if existing_save_dir:
                sample.materialize(
                    bucket_parent=os.path.dirname(str(existing_save_dir).rstrip(os.sep))
                )
            else:
                sample.set_logger()
        return sample

    def _export_sampler_state(self) -> Dict[str, Any]:
        raise NotImplementedError

    def _import_sampler_state(self, payload: Dict[str, Any]) -> None:
        raise NotImplementedError

    def export_runtime_state(self) -> Dict[str, Any]:
        if not self.supports_runtime_checkpointing():
            raise RuntimeError("Runtime checkpointing is not supported for this sampler.")
        sampler_state = self._export_sampler_state()
        if isinstance(sampler_state, dict) and "sampler_signature" not in sampler_state:
            sampler_state["sampler_signature"] = build_sampler_signature(self)
        return build_state_payload(
            run_spec=self._resolve_runtime_checkpoint_context(
                "_runtime_checkpoint_run_spec",
                "_runtime_checkpoint_run_spec_getter",
            ),
            factory_blueprint=self._resolve_runtime_checkpoint_context(
                "_runtime_checkpoint_factory_blueprint",
                "_runtime_checkpoint_factory_blueprint_getter",
            ),
            sampler_state=sampler_state,
            reason=self._runtime_checkpoint_reason,
        )

    def import_runtime_state(self, payload: Dict[str, Any]) -> None:
        if "sampler_state" in payload and isinstance(payload.get("sampler_state"), dict):
            self._runtime_checkpoint_payload = deepcopy(payload)
            self._runtime_checkpoint_run_spec = deepcopy(payload.get("run_spec", {}))
            self._runtime_checkpoint_factory_blueprint = deepcopy(payload.get("factory_blueprint", {}))
            payload = dict(payload.get("sampler_state", {}))
        signature = dict(payload.get("sampler_signature", {}))
        if signature:
            ok, reason = validate_sampler_signature(build_sampler_signature(self), signature)
            if not ok:
                raise ValueError(reason)
        self._import_sampler_state(payload)
        self._runtime_checkpoint_loaded = True
        if self._runtime_checkpoint_payload is None:
            self._runtime_checkpoint_payload = {
                "format": "jarvis-hep.simple-runtime",
                "version": 1,
                "timestamp_utc": utc_now_iso(),
                "sampler_state": deepcopy(payload),
                "run_spec": deepcopy(self._runtime_checkpoint_run_spec or {}),
                "factory_blueprint": deepcopy(self._runtime_checkpoint_factory_blueprint or {}),
                "integrity": {},
            }

    def restore_runtime_checkpoint_if_available(self, payload: Dict[str, Any] | None = None) -> bool:
        if not self._runtime_checkpoint_enabled or not self._runtime_checkpoint_auto_resume:
            return False
        if self._runtime_state_saver is None:
            return False
        if payload is None:
            payload = self._runtime_state_saver.load(default=None)
        if not payload:
            return False
        ok, reason = validate_state_payload(payload)
        if not ok:
            self._checkpoint_log(
                "warning",
                f"StateSaver checkpoint rejected -> {reason}",
            )
            return False
        self.import_runtime_state(payload)
        self._checkpoint_log(
            "warning",
            f"StateSaver checkpoint restored -> {self._runtime_checkpoint_file}",
        )
        return True

    def persist_runtime_checkpoint(self, *, force: bool = False, reason: str = "") -> bool:
        if not self._runtime_checkpoint_enabled or self._runtime_state_saver is None:
            return False
        with self._runtime_checkpoint_save_lock:
            now = time.monotonic()
            if (not force) and (
                (now - self._runtime_checkpoint_last_save_ts)
                < float(self._runtime_checkpoint_interval_seconds)
            ):
                return False

            self._runtime_checkpoint_reason = str(reason or "")
            payload = self.export_runtime_state()
            self._runtime_checkpoint_payload = payload
            self._runtime_state_saver.save(payload)
            self._runtime_checkpoint_last_save_ts = now
            if reason:
                self._checkpoint_log(
                    "info",
                    f"StateSaver checkpoint saved ({reason}) -> {self._runtime_checkpoint_file}",
                )
            else:
                self._checkpoint_log(
                    "info",
                    f"StateSaver checkpoint saved -> {self._runtime_checkpoint_file}",
                )
            return True

    def load_nuisance_sampler(self): 
        nui_config = (self.config.get('Sampling', {}) or {}).get('Nuisance', None)
        if nui_config: 
            if nui_config['Method'] == "Profile1D": 
                from jarvishep.Sampling.Source.Nuisance.profile1d import Profile1D
                self.nuisance_sampler = Profile1D()
                self.nuisance_sampler.set_logger(self.logger)
                self.nuisance_sampler.set_config(nui_config)
                self.info['sample']['nuisance'] = self.nuisance_sampler.get_info_card()
                self._with_nuisance = True

    @property
    def loglike(self):
        return self._loglike

    def set_max_workers(self, nn):
        self.max_workers = nn 

    @abstractmethod
    def load_schema_file(self) -> None:
        pass 

    @abstractmethod
    def set_config(self, config_info) -> None:
        pass 

    @abstractmethod
    def set_bucket_alloc(self) -> None: 
        limit = 200 
        width = 6 
        ba_config = None
        scan_cfg = self.config.get("Scan", {}) if isinstance(self.config, dict) else {}
        if isinstance(scan_cfg, dict):
            ba_config = scan_cfg.get("sample_directory", None)
        if not isinstance(ba_config, dict):
            ba_config = self.config.get("Directory_Setting", None)
        if isinstance(ba_config, dict):
            limit = ba_config.get("limit", 200)
            width = ba_config.get("width", 6)

        from jarvishep.Sampling.bucketallocator import BucketAllocator
        self.bucket_alloc = BucketAllocator(
            base_path=self.info['sample']['sample_dirs'],
            limit=limit,
            width=width,
            start_bucket=1,
            on_bucket_sealed=self._on_bucket_sealed if self._archive_enabled() else None,
        )
        self.bucket_alloc.check_and_update()
        if self._archive_enabled():
            self._ensure_archive_manager()
        
        
    @abstractmethod
    def init_generator(self) -> None:
        pass 
         

    @abstractmethod
    def set_factory(self) -> None: 
        pass 

    @abstractmethod
    def load_variable(self) -> None:
        variables = []
        for var_config in self.config['Sampling'].get("Variables", []):  # 使用get以防'Variables'键不存在
            var = Variable(
                name=var_config['name'],
                description=var_config['description'],
                distribution=var_config['distribution']['type'],
                parameters=var_config['distribution']['parameters']
            )
            variables.append(var)
        self.vars = tuple(variables)

    @abstractmethod
    def set_logger(self, logger) -> None:
        self.logger = logger 

    @abstractmethod 
    def combine_data(self, df_full) -> None:
        pass 
        

    
    @abstractmethod
    def set_likelihood(self, loglike):
        self._loglike = loglike
        self._loglike.logger = self.logger 
    
    @abstractmethod
    def logo_at_pos(self, pos, ax):
        from PIL import Image
        image_path = os.path.abspath(
            os.path.join(self.path.get("package_root", ""), "icons", "JarvisHEP.png")
        )
        with Image.open(image_path) as image:
            image = np.array(image.convert("RGBA"))
            ax.axes("off")
            ax.imshow(image, extent=[pos[0]-0.25, pos[0]+0.25, pos[1]-0.25, pos[1]+0.25], zorder=100)
            ax.text(pos[0]+0.27, pos[1]+0.15, "Jarvis-HEP", ha="left", va='top', color="#0F66C3", fontfamily="sans-serif", fontsize="small", fontstyle="normal", fontweight="bold")

    def evaluate_selection(self, expression, variables):
        """Evaluate selection expression and return a strict bool."""
        if expression is None:
            return True
        if not isinstance(variables, dict):
            raise BoolConversionError("Selection variables must be a dict.")

        custom_functions = update_funcs({})
        custom_constants = update_const({})
        symbols = {name: sp.symbols(name) for name in variables.keys()}
        locals_context = {**custom_functions, **custom_constants, **symbols}

        try:
            expr = sp.sympify(expression, locals=locals_context)
            evaluated = expr.subs(variables)
            return bool(evaluated)
        except Exception as exc:
            raise BoolConversionError(
                f"Cannot evaluate selection expression '{expression}' as boolean."
            ) from exc

    def check_evaluation(self):
        """Smoke-check selection evaluation on a deterministic probe point."""
        if not getattr(self, "_selectionexp", None):
            return True

        probe_values = {}
        for var in self.vars or []:
            name = getattr(var, "name", getattr(var, "_name", None))
            if not name:
                continue
            try:
                probe_values[name] = var.map_standard_random_to_distribution(0.5)
            except Exception:
                probe_values[name] = 0.5

        self.evaluate_selection(self._selectionexp, probe_values)
        return True

    def _should_eager_materialize(self) -> bool:
        sample_cfg = self.info.get("sample", {}) if isinstance(self.info, dict) else {}
        return should_eager_materialize(sample_cfg)

    def build_sample_config(self, base_sample_cfg: Dict[str, Any], save_dir: str | None = None) -> Dict[str, Any]:
        """Build per-sample config with minimal safe copying.

        Only deep-copy fields that are mutated downstream (e.g. ``nuisance``),
        while keeping read-only large fields shared to avoid hot-path deepcopy cost.
        """
        cfg = dict(base_sample_cfg or {})
        if save_dir is not None and self._should_eager_materialize():
            cfg["save_dir"] = save_dir

        nuisance = cfg.get("nuisance")
        if isinstance(nuisance, dict):
            cfg["nuisance"] = deepcopy(nuisance)
        return cfg

    def _archive_enabled(self) -> bool:
        sample_cfg = self.info.get("sample", {}) if isinstance(self.info, dict) else {}
        if "archive_samples" not in sample_cfg:
            # Keep backward compatibility for legacy/manual test setups that
            # do not carry the normalized runtime flag from Core.
            return False
        return bool(sample_cfg.get("archive_samples"))

    def _ensure_archive_manager(self):
        if not self._archive_enabled():
            return None
        if self.sample_archive_manager is None:
            from jarvishep.Sampling.sample_archive import SampleArchiveManager

            self.sample_archive_manager = SampleArchiveManager(
                sample_root=self.info["sample"]["sample_dirs"],
                logger=self.logger,
                enabled=True,
            )
            self.sample_archive_manager.start()
        return self.sample_archive_manager

    def _on_bucket_sealed(self, bucket_dir: str):
        manager = self._ensure_archive_manager()
        if manager is None:
            return
        manager.enqueue_bucket_dir(bucket_dir, blocking=True, timeout=30.0)

    def _next_bucket_dir_for_sample(self) -> str | None:
        if getattr(self, "bucket_alloc", None) is None:
            return None
        if not self._should_eager_materialize():
            return None
        return self.bucket_alloc.next_bucket_dir()

    def _on_sample_completed(self, sample_info: Dict[str, Any] | None) -> None:
        if getattr(self, "bucket_alloc", None) is None:
            return
        if not isinstance(sample_info, dict):
            return
        if not sample_info.get("_materialized"):
            return
        save_dir = sample_info.get("save_dir")
        if not save_dir:
            return
        bucket_dir = os.path.dirname(str(save_dir).rstrip(os.sep))
        try:
            self.bucket_alloc.mark_sample_finished(bucket_dir)
        except Exception:
            # Never break sampling cleanup on archive-accounting errors.
            pass

    def finalize_sample_archive(self):
        if not self._archive_enabled():
            return
        manager = self._ensure_archive_manager()
        if manager is None:
            return
        if self.bucket_alloc is not None:
            try:
                self.bucket_alloc.seal_current_bucket()
            except Exception:
                pass
        if self.bucket_alloc is not None:
            try:
                manager.enqueue_bucket_dir(
                    self.bucket_alloc.current_bucket_dir(),
                    blocking=True,
                    timeout=30.0,
                    force=True,
                )
            except Exception:
                pass
        manager.enqueue_all_existing_buckets()
        manager.shutdown(wait=True, timeout=300.0)
        self.sample_archive_manager = None
