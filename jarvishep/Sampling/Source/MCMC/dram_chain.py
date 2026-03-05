#!/usr/bin/env python3

import numpy as np

from .ammcmc_chain import AMMCMCChain


class DRAMChain(AMMCMCChain):
    """Delayed-Rejection Adaptive Metropolis (DRAM) chain runtime.

    Stage-0 uses AM proposal, then falls back to delayed-rejection stage(s)
    within the same iteration before marking final reject.
    """

    def __init__(
        self,
        initial_param,
        proposal_scale,
        n_iterations,
        adapt_enabled=True,
        adapt_start_iter=100,
        adapt_window=25,
        adapt_eps=1e-6,
        adapt_scale=2.38,
        dr_steps=2,
        dr_scale_factors=None,
    ):
        super().__init__(
            initial_param=initial_param,
            proposal_scale=proposal_scale,
            n_iterations=n_iterations,
            adapt_enabled=adapt_enabled,
            adapt_start_iter=adapt_start_iter,
            adapt_window=adapt_window,
            adapt_eps=adapt_eps,
            adapt_scale=adapt_scale,
        )
        self._dr_steps = max(1, int(dr_steps))
        if dr_scale_factors is None:
            dr_scale_factors = [1.0, 0.5]
        self._dr_scale_factors = self._normalize_scale_factors(dr_scale_factors)

        self._stage_proposals = {}
        self._stage_cov = {}
        self._stage_logl = {}
        self._stage_alpha = {}

    def _normalize_scale_factors(self, values):
        factors = [max(1e-6, float(x)) for x in values]
        if not factors:
            factors = [1.0, 0.5]
        if len(factors) < self._dr_steps:
            last = factors[-1]
            while len(factors) < self._dr_steps:
                last = max(1e-6, last * 0.5)
                factors.append(last)
        return factors[: self._dr_steps]

    def _clear_stage_cache(self):
        self._stage_proposals.clear()
        self._stage_cov.clear()
        self._stage_logl.clear()
        self._stage_alpha.clear()

    def _base_cov(self):
        if not self._adapt_enabled or self.iterations < self._adapt_start_iter:
            return (self.proposal_scale**2) * np.eye(self._dim, dtype=float)
        scale = (self._adapt_scale**2) / max(1.0, float(self._dim))
        cov = np.asarray(self._cov, dtype=float) * scale
        if cov.shape != (self._dim, self._dim):
            cov = (self.proposal_scale**2) * np.eye(self._dim, dtype=float)
        cov += self._adapt_eps * np.eye(self._dim, dtype=float)
        return cov

    def _stage_covariance(self, stage_index):
        factor = self._dr_scale_factors[min(int(stage_index), len(self._dr_scale_factors) - 1)]
        return np.asarray(self._base_cov(), dtype=float) * (factor**2)

    def _draw_stage_step(self, stage_index):
        cov = self._stage_covariance(stage_index)
        try:
            return np.random.multivariate_normal(np.zeros(self._dim, dtype=float), cov)
        except Exception:
            scale = self.proposal_scale * self._dr_scale_factors[
                min(int(stage_index), len(self._dr_scale_factors) - 1)
            ]
            return np.random.normal(0.0, scale, size=self._dim)

    def _accept_prob(self, old_logl, new_logl, beta=1.0):
        if old_logl is None:
            return 1.0
        delta = (float(new_logl) - float(old_logl)) * float(beta)
        if delta >= 0.0:
            return 1.0
        return float(np.exp(np.clip(delta, -700.0, 0.0)))

    def _log_mvn_density(self, point, mean, cov):
        x = np.asarray(point, dtype=float)
        mu = np.asarray(mean, dtype=float)
        c = np.asarray(cov, dtype=float)
        if c.shape != (self._dim, self._dim):
            c = (self.proposal_scale**2) * np.eye(self._dim, dtype=float)
        c = c + self._adapt_eps * np.eye(self._dim, dtype=float)
        sign, logdet = np.linalg.slogdet(c)
        if sign <= 0:
            c = c + (10.0 * self._adapt_eps) * np.eye(self._dim, dtype=float)
            sign, logdet = np.linalg.slogdet(c)
        inv = np.linalg.pinv(c)
        diff = x - mu
        quad = float(diff.T @ inv @ diff)
        return -0.5 * (self._dim * np.log(2.0 * np.pi) + logdet + quad)

    def _accept_prob_stage2(self, logl_stage2, beta):
        if self.last_loglikelihood is None:
            return 1.0
        if 0 not in self._stage_proposals or 0 not in self._stage_logl or 0 not in self._stage_alpha:
            return self._accept_prob(self.last_loglikelihood, logl_stage2, beta=beta)
        if 1 not in self._stage_proposals:
            return self._accept_prob(self.last_loglikelihood, logl_stage2, beta=beta)

        x = np.asarray(self.param, dtype=float)
        y1 = np.asarray(self._stage_proposals[0], dtype=float)
        y2 = np.asarray(self._stage_proposals[1], dtype=float)

        logl_x = float(self.last_loglikelihood)
        logl_y1 = float(self._stage_logl[0])
        logl_y2 = float(logl_stage2)

        alpha1_xy1 = float(self._stage_alpha[0])
        alpha1_y2y1 = float(self._accept_prob(logl_y2, logl_y1, beta=beta))

        cov1 = self._stage_cov.get(0, self._stage_covariance(0))
        log_q_y2_y1 = self._log_mvn_density(y1, y2, cov1)
        log_q_x_y1 = self._log_mvn_density(y1, x, cov1)

        numer = float(beta) * logl_y2 + log_q_y2_y1 + np.log(max(1e-12, 1.0 - alpha1_y2y1))
        denom = float(beta) * logl_x + log_q_x_y1 + np.log(max(1e-12, 1.0 - alpha1_xy1))
        ratio = float(np.exp(np.clip(numer - denom, -700.0, 0.0)))
        return min(1.0, ratio)

    def _finalize_accept(self, stage_index, accepted_logl):
        self.param = np.asarray(self._stage_proposals[int(stage_index)], dtype=float)
        self.proposed_param = np.asarray(self.param, dtype=float)
        self.last_loglikelihood = float(accepted_logl)
        self._history.append(np.asarray(self.param, dtype=float))
        self._accepted_since_adapt += 1
        self._maybe_adapt_cov()
        self.iterations += 1
        self._clear_stage_cache()

    def _finalize_reject(self, stage_index):
        self.proposed_param = np.asarray(
            self._stage_proposals.get(int(stage_index), self.param),
            dtype=float,
        )
        self.iterations += 1
        self._clear_stage_cache()

    def propose_stage(self, stage_index):
        if self.iterations >= self.n_iterations:
            raise StopIteration
        sid = int(stage_index)
        if sid == 0:
            self._clear_stage_cache()

        if self.iterations == 0 and sid == 0:
            proposal = np.random.rand(self._dim)
            self._stage_cov[sid] = (self.proposal_scale**2) * np.eye(self._dim, dtype=float)
        else:
            proposal = None
            for _ in range(2048):
                step = self._draw_stage_step(sid)
                cand = self.param + step
                if np.all((cand >= 0.0) & (cand <= 1.0)):
                    proposal = cand
                    break
            if proposal is None:
                scale = self.proposal_scale * self._dr_scale_factors[
                    min(sid, len(self._dr_scale_factors) - 1)
                ]
                proposal = np.clip(
                    self.param + np.random.normal(0.0, scale, size=self._dim),
                    0.0,
                    1.0,
                )
            self._stage_cov[sid] = self._stage_covariance(sid)

        proposal = np.asarray(proposal, dtype=float)
        self.proposed_param = proposal
        self._stage_proposals[sid] = proposal
        return proposal

    def consume_stage_result(self, stage_index, new_loglikelihood, beta=1.0):
        sid = int(stage_index)
        if sid not in self._stage_proposals:
            raise RuntimeError(f"DRAM stage proposal not found for stage={sid}")

        new_logl = float(new_loglikelihood)
        beta = float(beta)

        if sid == 0:
            alpha = float(self._accept_prob(self.last_loglikelihood, new_logl, beta=beta))
            self._stage_logl[sid] = new_logl
            self._stage_alpha[sid] = alpha
            accepted = bool(np.random.rand() < alpha)
            if accepted:
                self._finalize_accept(sid, new_logl)
                return {
                    "iteration_done": True,
                    "accepted": True,
                    "next_stage": None,
                    "stage_attempts": 1,
                    "logl": self.last_loglikelihood,
                }
            if self._dr_steps > 1 and self.iterations > 0:
                return {
                    "iteration_done": False,
                    "accepted": False,
                    "next_stage": 1,
                    "stage_attempts": 1,
                    "logl": self.last_loglikelihood,
                }
            self._finalize_reject(sid)
            return {
                "iteration_done": True,
                "accepted": False,
                "next_stage": None,
                "stage_attempts": 1,
                "logl": self.last_loglikelihood,
            }

        if sid == 1:
            alpha = float(self._accept_prob_stage2(new_logl, beta=beta))
        else:
            alpha = float(self._accept_prob(self.last_loglikelihood, new_logl, beta=beta))

        self._stage_logl[sid] = new_logl
        self._stage_alpha[sid] = alpha
        accepted = bool(np.random.rand() < alpha)

        if accepted:
            self._finalize_accept(sid, new_logl)
            return {
                "iteration_done": True,
                "accepted": True,
                "next_stage": None,
                "stage_attempts": sid + 1,
                "logl": self.last_loglikelihood,
            }

        if sid + 1 < self._dr_steps:
            return {
                "iteration_done": False,
                "accepted": False,
                "next_stage": sid + 1,
                "stage_attempts": sid + 1,
                "logl": self.last_loglikelihood,
            }

        self._finalize_reject(sid)
        return {
            "iteration_done": True,
            "accepted": False,
            "next_stage": None,
            "stage_attempts": sid + 1,
            "logl": self.last_loglikelihood,
        }

    def __iter__(self):
        return self

    def __next__(self):
        return self.propose_stage(0)

    def update(self, new_loglikelihood, beta=1.0):
        outcome = self.consume_stage_result(0, new_loglikelihood, beta=beta)
        return bool(outcome.get("accepted", False))
