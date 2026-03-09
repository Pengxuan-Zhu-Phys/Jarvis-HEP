#!/usr/bin/env python3
from __future__ import annotations

from jarvishep.Sampling.Source.MCMC.config_contract import (
    bounds_get_float,
    bounds_get_int,
    bounds_get_list,
)
from jarvishep.Sampling.dream import DREAM


class DREAMLite(DREAM):
    """Backward-compatible lite profile on top of DREAM sampler runtime."""

    def __init__(self) -> None:
        super().__init__()
        self.method = "DREAMLite"
        self._dream_snooker_prob = 0.08
        self._dream_archive_size = 256
        self._dream_crossover_values = [0.5, 0.9]
        self._dream_crossover_adapt_interval = 128
        self._dream_scale_jitter = 0.05

    def load_schema_file(self):
        self.schema = self.path.get("DREAMLiteSchema", self.path.get("DREAMSchema", self.path.get("DEMCMCSchema", self.path["MCMCSchema"])))

    def init_generator(self):
        super().init_generator()
        smp = self.config["Sampling"]["Bounds"]
        self._dream_snooker_prob = bounds_get_float(
            smp,
            "dream_snooker_prob",
            aliases=("dream.snooker_prob",),
            default=self._dream_snooker_prob,
            minimum=0.0,
        )
        self._dream_archive_size = bounds_get_int(
            smp,
            "dream_archive_size",
            aliases=("dream.archive_size",),
            default=self._dream_archive_size,
            minimum=1,
        )

        values = bounds_get_list(
            smp,
            "dream_crossover_values",
            aliases=("dream.crossover_values",),
            default=self._dream_crossover_values,
        )
        if isinstance(values, (int, float)):
            values = [float(values)]
        self._dream_crossover_values = [float(x) for x in values] if values else [0.9]
        self._dream_crossover_adapt_interval = bounds_get_int(
            smp,
            "dream_crossover_adapt_interval",
            aliases=("dream.crossover_adapt_interval",),
            default=self._dream_crossover_adapt_interval,
            minimum=1,
        )
        self._dream_scale_jitter = bounds_get_float(
            smp,
            "dream_scale_jitter",
            aliases=("dream.scale_jitter",),
            default=self._dream_scale_jitter,
            minimum=0.0,
        )
