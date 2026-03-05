#!/usr/bin/env python3
from __future__ import annotations

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
        self._dream_snooker_prob = float(smp.get("dream_snooker_prob", self._dream_snooker_prob))
        self._dream_archive_size = int(smp.get("dream_archive_size", self._dream_archive_size))

        values = smp.get("dream_crossover_values", self._dream_crossover_values)
        if isinstance(values, (int, float)):
            values = [float(values)]
        self._dream_crossover_values = [float(x) for x in values] if values else [0.9]
        self._dream_crossover_adapt_interval = int(
            smp.get("dream_crossover_adapt_interval", self._dream_crossover_adapt_interval)
        )
        self._dream_scale_jitter = float(smp.get("dream_scale_jitter", self._dream_scale_jitter))
