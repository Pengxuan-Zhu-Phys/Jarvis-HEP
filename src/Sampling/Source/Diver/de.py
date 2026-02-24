#!/usr/bin/env python3
"""Pure-Python Diver differential-evolution engine.

This module provides a modern NumPy implementation inspired by the original
Fortran Diver code path (jDE/lambdajDE/current variants, boundary handling,
duplicate control, and smoothed fractional-improvement convergence).

The optimizer works in a unit hypercube and minimizes an objective value.
In Jarvis-HEP integration we use objective = -LogL.
"""

from __future__ import annotations

from collections import deque
from dataclasses import asdict, dataclass
from typing import Callable, Optional

import numpy as np


Array = np.ndarray


@dataclass(slots=True)
class DEConfig:
    """Configuration of the DE optimizer."""

    dim: int
    pop_size: int
    mode: str = "current"  # current | jDE | lambdajDE
    max_gen: int = 500
    max_civ: int = 1
    F: float = 0.7
    Cr: float = 0.9
    lambda_: float = 0.0
    current: bool = False
    expon: bool = False
    bndry: int = 1  # 0/1: none/brick-wall, 2: random reinit, 3: reflection
    convthresh: float = 1e-3
    convsteps: int = 10
    remove_duplicates: bool = False
    discard_unfit_points: bool = False
    max_acceptable_value: float = np.inf
    seed: Optional[int] = None


@dataclass(slots=True)
class DEResult:
    """Optimization result."""

    best_vector: Array
    best_cost: float
    best_loglike: float
    best_civilization: int
    best_generation: int
    evaluations: int
    acceptance_rate: float
    population: Array
    costs: Array


class DifferentialEvolution:
    """Differential evolution optimizer in a unit hypercube.

    - Objective is minimized.
    - Population vectors are always represented in [0, 1]^D.
    """

    # jDE/lambdajDE control parameters (Brest et al. style)
    TAU_F = 0.1
    F_L = 0.1
    F_U = 0.9  # Effective range [0.1, 1.0)
    TAU_CR = 0.1
    TAU_LAMBDA = 0.1

    def __init__(self, config: DEConfig) -> None:
        cfg = DEConfig(**asdict(config))
        cfg.mode = self._normalize_mode(cfg.mode)
        cfg.pop_size = max(int(cfg.pop_size), 5)
        cfg.dim = max(int(cfg.dim), 1)
        cfg.max_gen = max(int(cfg.max_gen), 1)
        cfg.max_civ = max(int(cfg.max_civ), 1)
        cfg.convsteps = max(int(cfg.convsteps), 1)
        cfg.Cr = float(np.clip(cfg.Cr, 0.0, 1.0))
        cfg.lambda_ = float(np.clip(cfg.lambda_, 0.0, 1.0))
        cfg.F = float(cfg.F)

        self.cfg = cfg
        self.rng = np.random.default_rng(cfg.seed)

        # Self-adaptive memories (active only in jDE/lambdajDE)
        self._fjde: Optional[Array] = None
        self._crjde: Optional[Array] = None
        self._lambdajde: Optional[Array] = None

        # Convergence state (SFIM-style)
        self._mean_cost: float = np.inf
        self._improvements: deque[float] = deque([1.0] * self.cfg.convsteps, maxlen=self.cfg.convsteps)

    @staticmethod
    def _normalize_mode(mode: str) -> str:
        m = (mode or "current").strip().lower()
        aliases = {
            "jde": "jDE",
            "j_de": "jDE",
            "jde+": "jDE",
            "lambdajde": "lambdajDE",
            "lambda_jde": "lambdajDE",
            "current": "current",
        }
        return aliases.get(m, "current")

    @staticmethod
    def _to_cost(loglike: Array) -> Array:
        cost = np.full(loglike.shape, np.inf, dtype=float)
        finite = np.isfinite(loglike)
        cost[finite] = -loglike[finite]
        return cost

    def _reset_convergence(self) -> None:
        self._mean_cost = np.inf
        self._improvements = deque([1.0] * self.cfg.convsteps, maxlen=self.cfg.convsteps)

    def _update_convergence(self, costs: Array) -> bool:
        finite = costs[np.isfinite(costs)]
        if finite.size == 0:
            fracdiff = 1.0
            cur = np.inf
        else:
            cur = float(np.mean(finite))
            if not np.isfinite(self._mean_cost) or self._mean_cost == 0.0:
                fracdiff = 1.0
            else:
                fracdiff = 1.0 - cur / self._mean_cost

        self._mean_cost = cur
        self._improvements.append(float(fracdiff))
        sfim = float(np.mean(self._improvements))
        return sfim < self.cfg.convthresh

    def _init_adaptive_memory(self) -> None:
        if self.cfg.mode in ("jDE", "lambdajDE"):
            self._fjde = self.F_L + self.rng.random(self.cfg.pop_size) * self.F_U
            self._crjde = self.rng.random(self.cfg.pop_size)
            if self.cfg.mode == "lambdajDE":
                self._lambdajde = self.rng.random(self.cfg.pop_size)
            else:
                self._lambdajde = None
        else:
            self._fjde = None
            self._crjde = None
            self._lambdajde = None

    def _sample_distinct(self, k: int, forbidden: set[int]) -> Array:
        pool = np.array([idx for idx in range(self.cfg.pop_size) if idx not in forbidden], dtype=int)
        if pool.size < k:
            raise ValueError(
                f"Population too small ({self.cfg.pop_size}) for mutation with excluded={sorted(forbidden)}"
            )
        return self.rng.choice(pool, size=k, replace=False)

    def _trial_parameters(self, idx: int) -> tuple[float, float, float]:
        if self.cfg.mode == "current":
            return self.cfg.F, self.cfg.Cr, self.cfg.lambda_

        assert self._fjde is not None and self._crjde is not None

        if self.rng.random() < self.TAU_F:
            trial_f = self.F_L + self.rng.random() * self.F_U
        else:
            trial_f = float(self._fjde[idx])

        if self.rng.random() < self.TAU_CR:
            trial_cr = float(self.rng.random())
        else:
            trial_cr = float(self._crjde[idx])

        if self.cfg.mode == "lambdajDE":
            assert self._lambdajde is not None
            if self.rng.random() < self.TAU_LAMBDA:
                trial_lambda = float(self.rng.random())
            else:
                trial_lambda = float(self._lambdajde[idx])
        else:
            trial_lambda = float(self.cfg.lambda_)

        return trial_f, trial_cr, trial_lambda

    def _mutate(self, population: Array, costs: Array, idx: int, trial_f: float, trial_lambda: float) -> Array:
        best_idx = int(np.argmin(costs)) if trial_lambda > 0.0 or self.cfg.mode == "lambdajDE" else idx

        if self.cfg.mode in ("jDE", "lambdajDE"):
            r1, r2, r3 = self._sample_distinct(3, {idx, best_idx})
            donor = (
                trial_lambda * population[best_idx]
                + (1.0 - trial_lambda) * population[r1]
                + trial_f * (population[r2] - population[r3])
            )
            return donor

        # Generic/current mode (including fixed-lambda variants)
        if self.cfg.current or np.isclose(trial_lambda, 1.0):
            ri = idx
        else:
            ri = int(self._sample_distinct(1, {idx, best_idx})[0])

        rj, rk = self._sample_distinct(2, {idx, best_idx, ri})
        donor = (
            trial_lambda * population[best_idx]
            + (1.0 - trial_lambda) * population[ri]
            + trial_f * (population[rj] - population[rk])
        )
        return donor

    def _crossover(self, target: Array, donor: Array, trial_cr: float) -> Array:
        trial = target.copy()

        if self.cfg.expon and self.cfg.mode == "current":
            start = int(self.rng.integers(0, self.cfg.dim))
            length = 1
            while length < self.cfg.dim and self.rng.random() <= trial_cr:
                length += 1
            idx = (start + np.arange(length)) % self.cfg.dim
            trial[idx] = donor[idx]
            return trial

        mask = self.rng.random(self.cfg.dim) <= trial_cr
        mask[int(self.rng.integers(0, self.cfg.dim))] = True
        trial[mask] = donor[mask]
        return trial

    def _apply_boundary(self, trial: Array, target: Array) -> tuple[Array, bool]:
        if np.all((trial >= 0.0) & (trial <= 1.0)):
            return trial, True

        if self.cfg.bndry == 1:
            # Brick wall: revert to target and mark invalid (trial gets +inf cost).
            return target.copy(), False

        if self.cfg.bndry == 2:
            # Random reinitialization.
            return self.rng.random(self.cfg.dim), True

        if self.cfg.bndry == 3:
            # Reflection back into hypercube.
            reflected = trial.copy()
            for _ in range(4):
                reflected = np.where(reflected > 1.0, 2.0 - reflected, reflected)
                reflected = np.where(reflected < 0.0, -reflected, reflected)
                if np.all((reflected >= 0.0) & (reflected <= 1.0)):
                    break
            return np.clip(reflected, 0.0, 1.0), True

        # Not enforced.
        return trial, True

    def _replace_duplicates(
        self,
        population: Array,
        costs: Array,
        evaluate_loglike: Callable[[Array], Array],
    ) -> tuple[int, int]:
        rounded = np.round(population, decimals=12)
        _, first_indices = np.unique(rounded, axis=0, return_index=True)
        first_set = set(int(i) for i in first_indices)
        duplicates = [idx for idx in range(self.cfg.pop_size) if idx not in first_set]

        if not duplicates:
            return 0, 0

        evals = 0
        for idx in duplicates:
            new_vec = self.rng.random(self.cfg.dim)
            logl = np.asarray(evaluate_loglike(new_vec[None, :]), dtype=float)
            if logl.size != 1:
                raise ValueError("evaluate_loglike must return one value per input vector")
            population[idx] = new_vec
            costs[idx] = self._to_cost(logl)[0]
            evals += 1

            if self.cfg.mode in ("jDE", "lambdajDE") and self._fjde is not None and self._crjde is not None:
                self._fjde[idx] = self.F_L + self.rng.random() * self.F_U
                self._crjde[idx] = self.rng.random()
                if self.cfg.mode == "lambdajDE" and self._lambdajde is not None:
                    self._lambdajde[idx] = self.rng.random()

        return evals, len(duplicates)

    def run(
        self,
        evaluate_loglike: Callable[[Array], Array],
        *,
        logger=None,
    ) -> DEResult:
        """Run differential evolution.

        Parameters
        ----------
        evaluate_loglike:
            Callable that receives vectors with shape ``(N, D)`` in unit-cube
            space and returns an array-like of log-likelihood values with shape
            ``(N,)``.
        logger:
            Optional logger with ``info`` method.
        """

        total_evals = 0
        total_accept = 0
        total_trials = 0

        best_vector = np.zeros(self.cfg.dim, dtype=float)
        best_cost = np.inf
        best_civ = 1
        best_gen = 0

        final_population = np.zeros((self.cfg.pop_size, self.cfg.dim), dtype=float)
        final_costs = np.full(self.cfg.pop_size, np.inf, dtype=float)

        for civ in range(1, self.cfg.max_civ + 1):
            self._init_adaptive_memory()
            self._reset_convergence()

            population = self.rng.random((self.cfg.pop_size, self.cfg.dim))
            logl0 = np.asarray(evaluate_loglike(population), dtype=float)
            if logl0.shape != (self.cfg.pop_size,):
                raise ValueError(
                    f"evaluate_loglike returned shape {logl0.shape}; expected {(self.cfg.pop_size,)}"
                )
            costs = self._to_cost(logl0)
            total_evals += self.cfg.pop_size

            civ_best_idx = int(np.argmin(costs))
            if costs[civ_best_idx] < best_cost:
                best_cost = float(costs[civ_best_idx])
                best_vector = population[civ_best_idx].copy()
                best_civ = civ
                best_gen = 0

            if logger is not None:
                finite = costs[np.isfinite(costs)]
                mean_logl = -float(np.mean(finite)) if finite.size else float("-inf")
                logger.info(
                    f"[Diver] civ={civ}/{self.cfg.max_civ} gen=0/{self.cfg.max_gen} "
                    f"best_logl={-costs[civ_best_idx]:.6e} mean_logl={mean_logl:.6e}"
                )

            for gen in range(1, self.cfg.max_gen + 1):
                trial_population = np.empty_like(population)
                trial_f = np.full(self.cfg.pop_size, self.cfg.F, dtype=float)
                trial_cr = np.full(self.cfg.pop_size, self.cfg.Cr, dtype=float)
                trial_lambda = np.full(self.cfg.pop_size, self.cfg.lambda_, dtype=float)
                valid_mask = np.ones(self.cfg.pop_size, dtype=bool)

                for i in range(self.cfg.pop_size):
                    f_i, cr_i, lambda_i = self._trial_parameters(i)
                    donor = self._mutate(population, costs, i, f_i, lambda_i)
                    trial = self._crossover(population[i], donor, cr_i)
                    bounded, valid = self._apply_boundary(trial, population[i])

                    trial_population[i] = bounded
                    trial_f[i] = f_i
                    trial_cr[i] = cr_i
                    trial_lambda[i] = lambda_i
                    valid_mask[i] = valid

                trial_logl = np.full(self.cfg.pop_size, -np.inf, dtype=float)
                if np.any(valid_mask):
                    eval_block = trial_population[valid_mask]
                    block_logl = np.asarray(evaluate_loglike(eval_block), dtype=float)
                    if block_logl.shape != (eval_block.shape[0],):
                        raise ValueError(
                            f"evaluate_loglike returned shape {block_logl.shape}; expected {(eval_block.shape[0],)}"
                        )
                    trial_logl[valid_mask] = block_logl
                    total_evals += int(eval_block.shape[0])

                trial_costs = self._to_cost(trial_logl)

                acceptable = trial_costs <= self.cfg.max_acceptable_value
                if self.cfg.discard_unfit_points:
                    selectable = acceptable
                else:
                    selectable = np.ones_like(acceptable, dtype=bool)

                improved = selectable & (trial_costs <= costs)
                accepted = int(np.count_nonzero(improved))

                next_population = np.where(improved[:, None], trial_population, population)
                next_costs = np.where(improved, trial_costs, costs)

                if self.cfg.mode in ("jDE", "lambdajDE") and self._fjde is not None and self._crjde is not None:
                    self._fjde = np.where(improved, trial_f, self._fjde)
                    self._crjde = np.where(improved, trial_cr, self._crjde)
                    if self.cfg.mode == "lambdajDE" and self._lambdajde is not None:
                        self._lambdajde = np.where(improved, trial_lambda, self._lambdajde)

                if self.cfg.remove_duplicates:
                    extra_evals, n_dups = self._replace_duplicates(next_population, next_costs, evaluate_loglike)
                    total_evals += extra_evals
                    if logger is not None and n_dups > 0:
                        logger.info(f"[Diver] civ={civ} gen={gen}: replaced {n_dups} duplicate vectors")

                population = next_population
                costs = next_costs

                total_accept += accepted
                total_trials += self.cfg.pop_size

                civ_best_idx = int(np.argmin(costs))
                civ_best_cost = float(costs[civ_best_idx])
                if civ_best_cost < best_cost:
                    best_cost = civ_best_cost
                    best_vector = population[civ_best_idx].copy()
                    best_civ = civ
                    best_gen = gen

                if logger is not None:
                    finite = costs[np.isfinite(costs)]
                    mean_logl = -float(np.mean(finite)) if finite.size else float("-inf")
                    acc_rate = accepted / float(self.cfg.pop_size)
                    logger.info(
                        f"[Diver] civ={civ}/{self.cfg.max_civ} gen={gen}/{self.cfg.max_gen} "
                        f"best_logl={-civ_best_cost:.6e} mean_logl={mean_logl:.6e} acc={acc_rate:.3f}"
                    )

                if self._update_convergence(costs):
                    if logger is not None:
                        logger.info(
                            f"[Diver] civ={civ} converged at gen={gen} (convthresh={self.cfg.convthresh})"
                        )
                    break

            final_population = population
            final_costs = costs

        acceptance_rate = (total_accept / total_trials) if total_trials > 0 else 0.0
        best_loglike = -best_cost if np.isfinite(best_cost) else float("-inf")

        return DEResult(
            best_vector=best_vector,
            best_cost=float(best_cost),
            best_loglike=float(best_loglike),
            best_civilization=int(best_civ),
            best_generation=int(best_gen),
            evaluations=int(total_evals),
            acceptance_rate=float(acceptance_rate),
            population=final_population,
            costs=final_costs,
        )
