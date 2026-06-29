#!/usr/bin/env python3
"""Sampler factory dispatch (V1-compatible method names)."""

from __future__ import annotations

from jarvishep2.Sampling.sampler import SamplingVirtial

STATELESS_METHODS = frozenset({"Bridson", "Random", "Grid", "CSV"})


class Distributor:
    RESUME_SUPPORT_STATUS = {
        "Bridson": "implemented",
        "Random": "implemented",
        "Grid": "implemented",
    }

    @classmethod
    def get_resume_status(cls, method: str) -> str:
        return cls.RESUME_SUPPORT_STATUS.get(method, "intentionally unsupported")

    @staticmethod
    def set_method(method: str) -> SamplingVirtial:
        match str(method).strip():
            case "Bridson":
                from jarvishep2.Sampling.bridson import Bridson

                return Bridson()
            case _:
                raise NotImplementedError(
                    f"Sampling.Method '{method}' is not implemented in Jarvis-HEP V2 yet"
                )


__all__ = ["Distributor", "STATELESS_METHODS"]