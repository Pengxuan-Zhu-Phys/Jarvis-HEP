#!/usr/bin/env python3 

import os, sys 
from base import Base
from abc import ABCMeta, abstractmethod

from typing import Any, Dict, List, Optional
import numpy as np
from Sampling.nuisance_sampler import NuisanceBase
from math import sqrt

class Profile1D(NuisanceBase):
    def __init__(self):
        super().__init__()
        self.name               =     "Profile1D"
        self.var_name:  str     = ""
        self.var_desc:  str     = ""
        self.zmin:      float   = 0.0
        self.zmax:      float   = 1.0
        self._phi               = float(0.5 + sqrt(5.0) / 2)
        self.mode:      str     = "max"

    @property 
    def invphi(self):
        return 1.0 / self._phi 
    
    @property
    def invphi2(self):
        return self.invphi * self.invphi

    def set_logger(self, logger):
        super().set_logger(logger)
        self.logger.warning("Profile-1D nuisance sampler is selected, using 'Golden-section search' method for best estimation search.")

    def set_config(self, config):
        try: 
            self.logger.warning("Loading nuisance parameter setting for Profile 1D sampler")
            vars = config['Variables']
            logLs = config['LogLikelihood']
            pconds= config['PassCondition']
            self.max_attempt = config['MaxAttempt']
            self.mode       = config['TargetMode']
            self.logger.info("Loading nuisance setting")
        except Exception as e: 
            self.logger.error(e)
        
        self.load_vars(vars)
        self.load_loglike(logLs)
        self.load_Conditions(pconds)
        
    def get_info_card(self):
        from copy import deepcopy
        card = {
            "vname":    self.var_name,
            "active":   {
                "param":    {}, 
                "limits":   [deepcopy(self.zmin), deepcopy(self.zmax)],
                "LogLikelihoods":   {},
                "current":  "left",
                "pass":     False, 
                "PassConditions":   {}
            },
            "status":   "Init",
            "NAttempt": 1, 
            "LogLikelihoods":   list(self.loglikelihoods.keys()), 
            "PassConditions":   list(self.passconditions.keys()),
            "history":  {
                "param":    [],
                "LogLs":    [], 
                "PassConditions": {kk: [] for kk in self.passconditions.keys()},
                "LogLikelihoods": {kk: [] for kk in self.loglikelihoods.keys()},
            }
        }
        res = self.get_cd(card['active']['limits'])
        card['active']['probes'] = [res["c"], res["d"]]
        card['active']['probes_logL']  = [None, None]
        card['active']['param'].update({self.var_name: res["c"]}) 
        
        return card 

    def get_cd(self, ab):
        a =     ab[0]
        b =     ab[1]
        c =     b - (b - a) * self.invphi
        d =     a + (b - a) * self.invphi
        return {"c": c, "d": d}
        
    def renew_sample_info(self, sinfo):
        dnuisance = sinfo['nuisance']
        if dnuisance['active']['pass']: 
            sinfo['status'] = "Accept"
            dnuisance['status'] = "Accept" 
        elif dnuisance["NAttempt"] == self.max_attempt: 
            sinfo['status'] = "Accept"
            dnuisance['status'] = "Accept" 
        else: 
            self.next(dnuisance) 
        
    def next(self, dns):
        logl = self.update_active_into_history(dns)
        if dns['NAttempt'] == 1: 
            if dns['active']['current'] == "left": 
                dns['active']['probes_logL'][0] = logl 
                dns['active']['current'] = "right"
                dns['active']['param'].update({self.var_name: dns['active']['probes'][1] }) 
            else: 
                dns['active']['probes_logL'][1] = logl 
                dns['active']['current'] = "left"
                dns['active']['param'].update({self.var_name: dns['active']['probes'][0] }) 
            dns["NAttempt"] += 1
        else: 
            if dns['active']['current'] == "left":
                dns['active']['probes_logL'][0] = logl
            else:
                dns['active']['probes_logL'][1] = logl
            
            self.update_limits_by_probe(dns)

            if dns['active']['probes_logL'][0] is None:
                dns['active']['current'] = "left"
                dns['active']['param'].update({self.var_name: dns['active']['probes'][0]})
            else:
                dns['active']['current'] = "right"
                dns['active']['param'].update({self.var_name: dns['active']['probes'][1]})

            dns["NAttempt"] += 1 

    def update_limits_by_probe(self, dns): 
        fn      = dns['active']['probes_logL'].copy()
        limits  = dns['active']['limits'].copy()
        probe   = dns['active']['probes'].copy()
        if self.mode == "max":
            gn = [-xx for xx in fn]
        else: 
            gn = fn 
            
        # if yc < yd:
        if gn[0] < gn[1]: 
            # Keep [a, d]
            # b, d, yd = d, c, yc
            limits[1]   = probe[1]
            probe[1]    = probe[0]
            fn[1]       = fn[0]
            cd          = self.get_cd(limits)
            probe[0]    = cd["c"]
            fn[0]       = None 
        else:
            # c = b - (b - a) * invphi
            # yc = g(c)
            limits[0]   = probe[0]
            probe[0]    = probe[1]
            fn[0]       = fn[1]
            # d = a + (b - a) * invphi
            cd          = self.get_cd(limits)
            probe[1]    = cd["d"]
            fn[1]       = None

        dns['active']['probes_logL']  = fn
        dns['active']['limits']       = limits
        dns['active']['probes']       = probe

    def update_active_into_history(self, dns):
        dns['history']['param'].append(dns['active']['param'])
        logl    = 0. 
        for kk, vv in dns['active']['LogLikelihoods'].items():
            dns['history']['LogLikelihoods'][kk].append(vv)
            logl += vv 
        for kk, vv in dns['active']['PassConditions'].items():
            dns['history']['PassConditions'][kk].append(vv)
        dns['history']['LogLs'].append(logl)
        dns['active']['LogLikelihoods'] = {}
        dns['active']['PassConditions'] = {}
        dns['active']['param'] = {}
        return logl     
    
    @property
    def loglikelihoods(self):
        return self._loglikelihoods
        
    @property
    def passconditions(self):
        return self._passconditions
        
    def load_vars(self, vars: List[Dict[str, Any]]):
        """Load nuisance variable definition.

        Notes (per design):
        - Profile1D *prefers* exactly ONE nuisance variable.
        - If more than one is provided, we accept ONLY the first and emit a warning.
        - For malformed definitions, we warn and fall back to safe defaults.
        """

        if not vars:
            # This is fatal: no nuisance variable to profile.
            self.logger.error("Profile1D expects at least ONE nuisance variable, got 0")
            raise ValueError("Profile1D expects at least ONE nuisance variable, got 0")

        if len(vars) != 1:
            self.logger.warning(
                f"Profile1D expects exactly ONE nuisance variable; got {len(vars)}. "
                "Only the first entry will be used."
            )

        var = vars[0]

        # defaults (safe)
        self.var_name = str(var.get("name", self.var_name or "ratio"))
        self.var_desc = str(var.get("description", ""))

        try:
            dist = var.get("distribution", {})
            params = dist.get("parameters", {})

            # parse bounds with fallbacks
            if "min" in params:
                self.zmin = float(params["min"])
            else:
                self.logger.warning(
                    f"Profile1D nuisance var '{self.var_name}' missing distribution.parameters.min; "
                    f"using default zmin={self.zmin}."
                )

            if "max" in params:
                self.zmax = float(params["max"])
            else:
                self.logger.warning(
                    f"Profile1D nuisance var '{self.var_name}' missing distribution.parameters.max; "
                    f"using default zmax={self.zmax}."
                )

        except Exception as e:
            # Per comment: do not throw here; warn and keep defaults.
            self.logger.warning(
                f"Invalid nuisance variable definition for Profile1D (using defaults z in [{self.zmin}, {self.zmax}]): {var}. "
                f"Error: {e}"
            )

        # validate range (warn + repair)
        if self.zmin > self.zmax:
            self.logger.warning(
                f"Profile1D nuisance var '{self.var_name}' has inverted range: min ({self.zmin}) > max ({self.zmax}). "
                "Swapping bounds."
            )
            self.zmin, self.zmax = self.zmax, self.zmin

        if self.zmin == self.zmax:
            self.logger.warning(
                f"Profile1D nuisance var '{self.var_name}' has degenerate range: min == max == {self.zmin}. "
                "Expanding max by +1.0."
            )
            self.zmax = self.zmin + 1.0

        self.logger.info(
            f"Loaded Profile1D nuisance variable '{self.var_name}' with range [{self.zmin}, {self.zmax}]"
        )
            
    def load_loglike(self, logls: List[Dict[str, Any]]):
        """
        Parse nuisance LogLikelihood expressions.

        This function ONLY:
        - parses LogLikelihood items from config
        - builds NuisanceExpressionRegistry objects (one per expression)
        - preserves order
        - does NOT evaluate anything
        """

        from Module.nuisance_LogLikelihood import NuisanceExpressionRegistry, build_default_locals

        self._loglikelihoods = {}

        if not logls:
            self.logger.warning(
                "No nuisance LogLikelihood specified for Profile1D sampler. "
            )
            return

        for item in logls:
            try:
                name = item.get("name", None)
                expr = item.get("expression", None)

                if name is None or expr is None:
                    self.logger.warning(
                        f"Invalid nuisance LogLikelihood item (missing name or expression): {item}"
                    )
                    continue

                reg = NuisanceExpressionRegistry()
                reg.set_config(
                    name=str(name),
                    expression=str(expr)
                )

                self._loglikelihoods[str(name)] = {
                    "expression": str(expr),
                    "deps": reg.deps,
                    "expr": reg.eval,
                    "checker": reg.can_eval
                }
                self.logger.info(
                    "Nuisance loglikelihood function loaded \n\t name -> {}\n\texpr -> {}\n\tdependences -> {}".format(name, expr, reg.deps)
                )
            except Exception as e:
                self.logger.warning(
                    f"Failed to parse nuisance LogLikelihood item {item}: {e}"
                )

        if not self.loglikelihoods:
            self.logger.warning(
                "All nuisance LogLikelihood entries are invalid; "
                "Profile1D will run with empty objective."
            )
        else:
            names = self.loglikelihoods.keys()
            self.logger.info(
                f"Loaded {len(self.loglikelihoods)} nuisance LogLikelihood expression(s): {names}"
            )
            
    def load_Conditions(self, conds: List[Dict[str, Any]]):
        """Parse nuisance PassCondition expressions.

        This function ONLY:
        - parses PassCondition items from config
        - builds NuisancePassCondition objects (one per expression)
        - preserves names as keys
        - does NOT evaluate anything here

        Stored entries provide callable access via dict input:
        - entry['expr'](values) -> bool
        - entry['checker'](keys) -> (ok, missing)
        """

        from Module.nuisance_passCondition import NuisancePassCondition

        self._passconditions = {}

        if not conds:
            self.logger.warning(
                "No nuisance PassCondition specified for Profile1D sampler."
            )
            return

        for item in conds:
            try:
                name = item.get("name", None)
                expr = item.get("expression", None)

                if name is None or expr is None:
                    self.logger.warning(
                        f"Invalid nuisance PassCondition item (missing name or expression): {item}"
                    )
                    continue

                reg = NuisancePassCondition()
                reg.set_config(
                    name=str(name),
                    expression=str(expr),
                )

                self._passconditions[str(name)] = {
                    "expression": str(expr),
                    "deps": reg.deps,
                    "expr": reg.eval,
                    "checker": reg.can_eval,
                }

                self.logger.info(
                    "Nuisance passcondition function loaded \n\t name -> {}\n\texpr -> {}\n\tdependences -> {}".format(name, expr, reg.deps)
                )

            except Exception as e:
                self.logger.warning(
                    f"Failed to parse nuisance PassCondition item {item}: {e}"
                )

        if not self.passconditions:
            self.logger.warning(
                "All nuisance PassCondition entries are invalid; "
                "Profile1D will run without pass conditions."
            )
        else:
            names = self.passconditions.keys()
            self.logger.info(
                f"Loaded {len(self.passconditions)} nuisance PassCondition expression(s): {names}"
            )
            