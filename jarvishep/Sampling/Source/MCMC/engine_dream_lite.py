#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

from .engine_demcmc import DEMCMCChain


class DREAMLiteChain(DEMCMCChain):
    """DREAM-lite chain: DE proposals + snooker-style global moves."""

    def __init__(
        self,
        *args,
        dream_snooker_prob=0.1,
        dream_archive_size=256,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self._snooker_prob = min(1.0, max(0.0, float(dream_snooker_prob)))
        self._archive_size = max(16, int(dream_archive_size))
        self._archive = []
        self.last_move_type = "de"

    def _archive_population(self):
        if not self._archive:
            return []
        return [np.asarray(x, dtype=float) for x in self._archive]

    def _draw_snooker_move(self):
        others = self._other_population()
        archive = self._archive_population()
        pool = others + archive
        if len(pool) < 2:
            return self._draw_demove()

        ids = np.random.choice(len(pool), size=2, replace=False)
        xr1 = np.asarray(pool[int(ids[0])], dtype=float)
        xr2 = np.asarray(pool[int(ids[1])], dtype=float)

        direction = xr1 - xr2
        norm = float(np.linalg.norm(direction))
        if norm <= 1.0e-12:
            return self._draw_demove()
        direction = direction / norm

        jump_scale = float(np.random.normal(0.0, self._gamma))
        parallel = jump_scale * direction
        orthogonal = self._noise * np.random.normal(0.0, 1.0, size=self._dim)
        return parallel + orthogonal

    def _draw_demove(self):
        if np.random.rand() < self._snooker_prob:
            self.last_move_type = "snooker"
            return self._draw_snooker_move()
        self.last_move_type = "de"
        return super()._draw_demove()

    def update(self, new_loglikelihood, beta=1.0):
        accepted = super().update(new_loglikelihood, beta=beta)
        if accepted:
            self._archive.append(np.asarray(self.param, dtype=float))
            if len(self._archive) > self._archive_size:
                self._archive = self._archive[-self._archive_size :]
        return accepted
