#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

from .engine_demcmc import DEMCMCChain


class DREAMChain(DEMCMCChain):
    """DREAM chain: DE + snooker + adaptive crossover weights + archive."""

    def __init__(
        self,
        *args,
        dream_snooker_prob=0.1,
        dream_archive_size=512,
        dream_crossover_values=None,
        dream_crossover_adapt_interval=64,
        dream_scale_jitter=0.1,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self._snooker_prob = min(1.0, max(0.0, float(dream_snooker_prob)))
        self._archive_size = max(32, int(dream_archive_size))
        self._archive = []

        if dream_crossover_values is None:
            dream_crossover_values = [0.2, 0.5, 0.9]
        self._cr_values = [min(1.0, max(1.0e-3, float(x))) for x in dream_crossover_values]
        if not self._cr_values:
            self._cr_values = [0.9]
        self._cr_weights = np.ones(len(self._cr_values), dtype=float) / float(len(self._cr_values))
        self._cr_acc = np.zeros(len(self._cr_values), dtype=float)
        self._cr_trials = np.zeros(len(self._cr_values), dtype=float)
        self._cr_adapt_interval = max(8, int(dream_crossover_adapt_interval))
        self._scale_jitter = max(0.0, float(dream_scale_jitter))

        self._last_cr_idx = 0
        self.last_move_type = "de"
        self.last_crossover = self._cr_values[0]

    def _archive_population(self):
        if not self._archive:
            return []
        return [np.asarray(x, dtype=float) for x in self._archive]

    def _pick_cr_index(self):
        if len(self._cr_values) == 1:
            return 0
        probs = np.asarray(self._cr_weights, dtype=float)
        psum = float(np.sum(probs))
        if psum <= 0.0:
            probs = np.ones(len(self._cr_values), dtype=float) / float(len(self._cr_values))
        else:
            probs = probs / psum
        return int(np.random.choice(len(self._cr_values), p=probs))

    def _draw_de_step_with_cr(self, cr):
        others = self._other_population()
        archive = self._archive_population()
        pool = others + archive
        if len(pool) < 2:
            return np.random.normal(0.0, self.proposal_scale, size=self._dim)

        ids = np.random.choice(len(pool), size=2, replace=False)
        xr1 = np.asarray(pool[int(ids[0])], dtype=float)
        xr2 = np.asarray(pool[int(ids[1])], dtype=float)

        jitter = 1.0 + self._scale_jitter * float(np.random.normal(0.0, 1.0))
        gamma = max(1.0e-4, float(self._gamma) * max(0.25, jitter))
        step = gamma * (xr1 - xr2) + self._noise * np.random.normal(0.0, 1.0, size=self._dim)

        mask = np.random.rand(self._dim) < float(cr)
        if not np.any(mask):
            mask[np.random.randint(0, self._dim)] = True
        step = step * mask
        return step

    def _draw_snooker_move(self):
        others = self._other_population()
        archive = self._archive_population()
        pool = others + archive
        if len(pool) < 2:
            return self._draw_de_step_with_cr(0.9)

        ids = np.random.choice(len(pool), size=2, replace=False)
        xr1 = np.asarray(pool[int(ids[0])], dtype=float)
        xr2 = np.asarray(pool[int(ids[1])], dtype=float)

        direction = xr1 - xr2
        norm = float(np.linalg.norm(direction))
        if norm <= 1.0e-12:
            return self._draw_de_step_with_cr(0.9)
        direction = direction / norm

        jump_scale = float(np.random.normal(0.0, self._gamma))
        parallel = jump_scale * direction
        orthogonal = self._noise * np.random.normal(0.0, 1.0, size=self._dim)
        return parallel + orthogonal

    def _draw_demove(self):
        if np.random.rand() < self._snooker_prob:
            self.last_move_type = "snooker"
            self.last_crossover = 1.0
            return self._draw_snooker_move()

        self.last_move_type = "de"
        self._last_cr_idx = self._pick_cr_index()
        cr = float(self._cr_values[self._last_cr_idx])
        self.last_crossover = cr
        return self._draw_de_step_with_cr(cr)

    def _maybe_adapt_cr(self):
        if self.iterations <= 0:
            return
        if (self.iterations % self._cr_adapt_interval) != 0:
            return
        rates = (self._cr_acc + 1.0) / (self._cr_trials + 2.0)
        rates_sum = float(np.sum(rates))
        if rates_sum <= 0.0:
            return
        self._cr_weights = np.asarray(rates, dtype=float) / rates_sum
        self._cr_acc[:] = 0.0
        self._cr_trials[:] = 0.0

    def update(self, new_loglikelihood, beta=1.0):
        if self.last_move_type == "de":
            self._cr_trials[self._last_cr_idx] += 1.0

        accepted = super().update(new_loglikelihood, beta=beta)
        if accepted:
            self._archive.append(np.asarray(self.param, dtype=float))
            if len(self._archive) > self._archive_size:
                self._archive = self._archive[-self._archive_size :]
            if self.last_move_type == "de":
                self._cr_acc[self._last_cr_idx] += 1.0

        self._maybe_adapt_cr()
        return accepted
