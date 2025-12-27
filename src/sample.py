#!/usr/bin/env python3

import os, sys 
import json
import numpy as np 
import sympy as sp 
from copy import deepcopy 
from time import sleep
from base import Base
from loguru import logger
from uuid import uuid4
from inner_func import update_const, update_funcs
class Sample(Base):
    def __init__(self, params):
        self._params = params
        self._uuid = str(uuid4())
        self._likelihood = None  # Initialize likelihood with None
        self.processed = False
        self.observables = params 
        self.observables['uuid'] = self.uuid
        self.logger = None 
        self.handlers = {}
        self._with_nuisance = False
        self._nuisance_status   = False
        self._u = None

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

    def set_config(self, config):
        self.info = config
        self.create_info()
        if self.info.get("nuisance", {}):
            self.combine_nuisance_card()
            self._with_nuisance = True

    def create_info(self):
        save_dir = (self.info['save_dir'])
        # print(f"{self.uuid} -> save_dir is \n{save_dir}")
        self.info.update({
            "uuid": self.uuid,
            "params": self.params,
            "observables": self.observables, 
            "save_dir": os.path.join(save_dir, self.uuid), 
            "run_log":  os.path.join(save_dir, self.uuid, "Sample_running.log"),
            "logger":   None,
            "handlers": self.handlers,
            "status":   "Init"
        })
        self.set_logger()
    
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
    
    def set_logger(self): 
        logger_name = f"Sample@{self.info['uuid']}"
        self.info['logger_name'] = logger_name
        def filte_func(record):
            return logger_name in record['extra']['module']
        
        if not os.path.exists(self.info['save_dir']):
            os.makedirs(self.info['save_dir'])
        
        self.logger = logger.bind(module=logger_name, to_console=True, Jarvis=True)
        sample_handler = self.logger.add(self.info['run_log'], format=self.custom_format, level="DEBUG", rotation=None, retention=None, filter=filte_func)
        self.handlers['sample'] = sample_handler
        self.logger.info("Sample created into the Disk")
        self.info['logger'] = self.logger
        
    def custom_format(self, record):
        module = record["extra"].get("module", "No module")
        if "raw" in record["extra"]:
            return "{message}"
        else:
            return f"\n·•· <cyan>{module}</cyan> \n\t-> <green>{record['time']:MM-DD HH:mm:ss.SSS}</green> - [<level>{record['level']}</level>] >>> \n<level>{record['message']}</level>"

    
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
        self.info['status'] == "Running"
        if self._with_nuisance: 
            self.info['nuisance']['status'] == "Running"
            self.logger.info("{}\nSample start {}-th nuisance attempt".format(">"*60, self.info['nuisance']['NAttempt']))

    def close(self):
        self.close_logger() 
        
    def close_logger(self):
        """Close per-sample logger handler safely.

        This removes the per-sample file sink (if registered) and clears self.logger.
        It is safe to call multiple times.
        """
        if getattr(self, 'logger', None) is None:
            pass 

        # Remove per-sample handler if we have it
        handlers = getattr(self, 'handlers', None)
        if isinstance(handlers, dict):
            for kk, hh in handlers.items(): 
                try: 
                    logger.remove(hh)
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
                    