"""Opt-in GMFit acceptance declarations; project/native locations are external.

Only re-exports the real consumer. No physical implementation or unit arithmetic
is substituted, and loading this fixture never writes to the external project.
"""

import json
import os
from pathlib import Path
import sys

from jarvis_operas import OperaFunction

root = Path(os.environ["JARVIS_GMFIT_ROOT"])
sys.path.insert(0, os.environ["JARVIS_GMFIT_NATIVE_PATH"])
sys.path.insert(0, str(root / "deps/gmfit_hts_external/src"))

from gm_higgs.spheno_h1 import fit_h1_signals, opera_metadata  # noqa: E402


def source_report():
    report = json.loads(Path(os.environ["JARVIS_GMFIT_TRANSPORT_FILE"]).read_text())[
        "spheno_h1_inputs"
    ]
    return {
        "spheno_h1_inputs": report,
        "spheno_h1_mass": report["mass_GeV"],
        "spheno_h1_width": report["width_GeV"],
    }


__JARVIS_OPERAS__ = {
    "source_report": OperaFunction(
        namespace="gmfit_unit_ops",
        name="source_report",
        arity=None,
        return_dtype=None,
        numpy_impl=source_report,
    ),
    "fit_h1_signals": OperaFunction(
        namespace="gmfit_unit_ops",
        name="fit_h1_signals",
        arity=None,
        return_dtype=None,
        numpy_impl=fit_h1_signals,
        metadata=opera_metadata(),
    ),
}
