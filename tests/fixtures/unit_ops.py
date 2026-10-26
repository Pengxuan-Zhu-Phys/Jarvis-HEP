"""Portable registered operators for the HEP YAML/Worker unit contract tests."""

import asyncio
import time
from jarvis_operas import OperaFunction

CALLS = []


def limit(mass_GeV=200.0, logger=None, source_units=None, observables=None, delay=0.0):
    CALLS.append(mass_GeV)
    if delay:
        time.sleep(delay)
    if logger is not None:
        logger.info("unit fixture canonical mass {}", mass_GeV)
    return dict(
        value=0.038332,
        unit="pb",
        mass_received=mass_GeV,
        record_id="canonical-record",
        provenance={"dependencies": ["mass-record"]},
        status="ok",
        special={"a/b~c.d": [1.0, 2.0]},
    )


async def alimit(mass_GeV=200.0, logger=None):
    await asyncio.sleep(0)
    return limit(mass_GeV, logger=logger)


def consume(rate_pb):
    return {"value": rate_pb}


def logical(observables=None, **kwargs):
    return {"value": kwargs.get("mass_GeV", 0.0)}


def nullable(mass_GeV):
    CALLS.append(mass_GeV)
    return {"value": None}


def missing(mass_GeV):
    CALLS.append(mass_GeV)
    return {"status": "unavailable"}


def array(mass_GeV):
    return {"value": mass_GeV}


def decl(name, func, inputs, outputs):
    return OperaFunction(
        namespace="hep_unit_fixture",
        name=name,
        arity=None,
        return_dtype=None,
        numpy_impl=func,
        metadata={
            "params": {k: k for k in inputs},
            "unit_contract": {"version": "1.0", "inputs": inputs, "outputs": outputs},
        },
    )


MASS = {"mass_GeV": {"unit": "GeV", "shape": "scalar", "nullable": True}}
RATE = {"/value": {"unit": "pb", "shape": "scalar", "nullable": True}}
__JARVIS_OPERAS__ = {
    "limit": decl("limit", limit, MASS, RATE),
    "alimit": decl("alimit", alimit, MASS, RATE),
    "consume": decl(
        "consume", consume, {"rate_pb": {"unit": "pb", "shape": "scalar"}}, RATE
    ),
    "logical": decl("logical", logical, MASS, RATE),
    "nullable": decl("nullable", nullable, MASS, RATE),
    "missing": decl("missing", missing, MASS, RATE),
    "array": decl(
        "array",
        array,
        {"mass_GeV": {"unit": "GeV", "shape": "array"}},
        {"/value": {"unit": "GeV", "shape": "array"}},
    ),
}
