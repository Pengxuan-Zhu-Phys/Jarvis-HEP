#!/usr/bin/env python3
from __future__ import annotations

from base import Base
from abc import ABCMeta, abstractmethod

from typing import Any, Dict, List, Optional
import numpy as np

class NuisanceBase(Base):
    __metaclass__ = ABCMeta
    def __init__(self):
        super().__init__()
        self.config: Dict[str, Any] = {}
        self.name   = ""
        self.logger = None 
        self.max_attempt: int = 1
        self.names: List[str] = []
        self.var_defs: List[Dict[str, Any]] = []
        self._loglikelihoods: List[Dict[str, Any]] = []
        self._passconditions: List[Dict[str, Any]] = []
        self.seed: Optional[float] = None

    @property
    def loglikelihoods(self):
        return self._loglikelihoods
        
    @property
    def passconditions(self):
        return self._passconditions


    @abstractmethod
    def set_config(self, config: Dict[str, Any]) -> None:
        pass 

    @abstractmethod
    def get_info_card(self) -> Dict[str, Any]: 
        pass 

    @abstractmethod
    def renew_sample_info(self, sinfo: Dict[str, Any]) -> None: 
        pass 

    @abstractmethod
    def set_logger(self, logger) -> None: 
        # In Jarvis-HEP we store the logger name in loguru `extra['module']`.
        parent = "Jarvis-HEP.Nuisance"
        try:
            parent = logger._options[-1]['module']
        except Exception:
            pass

        # Derive a child name while inheriting all sinks/filters/extra flags.
        self.logger = logger.bind(module=f"{parent}.{self.name}")
        self.logger.warning("Profile1D nuisance sampler logger loaded")

    @abstractmethod
    def attach(self, point: Dict[str, Any], sample_info: Dict[str, Any]) -> Dict[str, Any]:
        # self.validate()  # Removed redundant validate call
        # Implement attach logic here
        return sample_info
