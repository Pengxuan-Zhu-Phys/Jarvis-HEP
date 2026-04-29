#!/usr/bin/env python3
from __future__ import annotations

import numpy as np


class ESSChain:
    """Elliptical proposal chain for ESS-style MCMC."""

    def __init__(
        self,
        initial_param,
        proposal_scale,
        n_iterations,
        prior_cov=None,
    ):
        self.param = np.asarray(initial_param, dtype=float)
        self.proposal_scale = float(proposal_scale)
        self.n_iterations = int(n_iterations)
        self.iterations = 0
        self.last_loglikelihood = None
        self.proposed_param = None

        self._dim = int(self.param.shape[0])
        self._prior_cov = self._normalize_cov(prior_cov)
        self._chol = np.linalg.cholesky(self._prior_cov)

    def _normalize_cov(self, prior_cov):
        if prior_cov is None:
            cov = np.eye(self._dim, dtype=float)
        else:
            cov = np.asarray(prior_cov, dtype=float)
            if cov.ndim == 0:
                cov = np.eye(self._dim, dtype=float) * float(cov)
            if cov.ndim == 1:
                if cov.shape[0] != self._dim:
                    raise ValueError(
                        f"ESS prior_cov vector size mismatch: expect {self._dim}, got {cov.shape[0]}"
                    )
                cov = np.diag(np.asarray(cov, dtype=float))
            if cov.shape != (self._dim, self._dim):
                raise ValueError(
                    f"ESS prior_cov matrix size mismatch: expect {(self._dim, self._dim)}, got {cov.shape}"
                )
        cov = np.asarray(cov, dtype=float)
        cov += 1.0e-8 * np.eye(self._dim, dtype=float)
        eigmin = float(np.min(np.linalg.eigvalsh(cov)))
        if eigmin <= 0.0:
            cov += (abs(eigmin) + 1.0e-6) * np.eye(self._dim, dtype=float)
        return cov

    def __iter__(self):
        return self

    def __next__(self):
        if self.iterations >= self.n_iterations:
            raise StopIteration
        if self.iterations == 0:
            proposal = np.random.rand(self._dim)
            self.proposed_param = proposal
            return proposal

        nu = self._chol @ np.random.normal(0.0, 1.0, size=self._dim)
        theta = np.random.uniform(0.0, 2.0 * np.pi)
        centered = self.param - 0.5
        cand = centered * np.cos(theta) + nu * np.sin(theta)
        proposal = np.clip(cand + 0.5, 0.0, 1.0)
        self.proposed_param = proposal
        return proposal

    def update(self, new_loglikelihood, beta=1.0):
        new_loglikelihood = float(new_loglikelihood)
        beta = float(beta)
        accepted = False

        if self.iterations == 0 or self.last_loglikelihood is None:
            accepted = True
        else:
            delta = (new_loglikelihood - float(self.last_loglikelihood)) * beta
            if delta >= 0.0:
                accepted = True
            else:
                accepted = bool(np.random.rand() < np.exp(np.clip(delta, -700.0, 0.0)))

        if accepted:
            self.param = np.asarray(self.proposed_param, dtype=float)
            self.last_loglikelihood = new_loglikelihood

        self.iterations += 1
        return accepted
