"""Opt-in acceptance with actual SPheno, native HiggsTools and the GMFit helper.

Set JARVIS_GMFIT_ROOT and JARVIS_GMFIT_NATIVE_PATH to existing external artifacts.
All generated files and JO registrations stay inside pytest's temporary project.
"""

import copy
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml

from jarvishep2.io.database import SimpleHDF5Writer
from test_operas_units_runtime import _setup

GENERATOR = """
import csv,json,subprocess,sys
from pathlib import Path
import sympy,yaml
root,native,scratch=sys.argv[1:]
root=Path(root)
sys.path.insert(0,native)
sys.path.insert(0,str(root/'deps/gmfit_hts_external/src'))
from gm_higgs.spheno_h1 import fit_h1_signals,transport_spheno_h1
scratch=Path(scratch)/'point'
scratch.mkdir(parents=True)
card=yaml.safe_load((root/'bin/GMFit_SPheno_h1_HS.yaml').read_text())
with (root/'data/gmfit_spheno_hs_points.csv').open() as stream:
    point=next(csv.DictReader(stream))
values={sympy.Symbol(k):float(v) for k,v in point.items() if k!='uuid'}
text=(root/'deps/LesHouches.in.GM_h1_HS').read_text()
variables=card['Calculators']['Modules'][0]['execution']['input'][0]['actions'][0]['variables']
for variable in variables:
    expression=sympy.sympify(str(variable['expression']),locals={str(k):k for k in values})
    text=text.replace(variable['placeholder'],format(float(expression.subs(values)),'.16e'))
(scratch/'input.slha').write_text(text)
# SPheno's fixed-size Fortran argv buffer cannot hold pytest's long filenames.
# Use relative names in a private cwd; the physical input and executable stay identical.
completed=subprocess.run([str(root/'deps/SPheno-4.0.5/bin/SPhenoGM'),'input.slha','spectrum.slha'],
                         cwd=scratch,text=True,capture_output=True,timeout=300)
(scratch/'spheno.log').write_text(completed.stdout+completed.stderr)
report=transport_spheno_h1(input_text=text,spectrum_text=(scratch/'spectrum.slha').read_text())
(scratch/'transport.json').write_text(json.dumps({'spheno_h1_inputs':report},allow_nan=False))
assert report['ready'],report['diagnostics']
result=fit_h1_signals(observables={'spheno_h1_inputs':report},mass_GeV=report['mass_GeV'],
                      width_GeV=report['width_GeV'],log_detail='none')
result.update(transport_status=report['status'],spheno_returncode=completed.returncode)
(scratch/'result.json').write_text(json.dumps(result,allow_nan=False))
assert result['ready'],result['diagnostics']
"""

DRIVER = """
import json,sys
from pathlib import Path
from jarvis_operas import persist_user_ops
from jarvishep2.runtime.core import Jarvis2Core
card,source,output=sys.argv[1:]
persist_user_ops(source,namespace='gmfit_unit_ops')
core=Jarvis2Core()
core.load_task_yaml(card)
core.config.update(skip_draw_flowchart=True,skip_render_plots=True)
outcome=core.run(write_run_summary=False)
Path(output).write_text(json.dumps(outcome.to_dict()))
raise SystemExit(outcome.exit_code)
"""


@pytest.mark.slow
@pytest.mark.skipif(
    not os.environ.get("JARVIS_GMFIT_ROOT"),
    reason="external GMFit native acceptance is opt-in",
)
def test_gmfit_three_units_through_real_worker_and_database(tmp_path):
    root = Path(os.environ["JARVIS_GMFIT_ROOT"]).resolve()
    native = Path(os.environ["JARVIS_GMFIT_NATIVE_PATH"]).resolve()
    assert (root / "deps/SPheno-4.0.5/bin/SPhenoGM").is_file()
    assert shutil.which("redis-server"), "real Redis is required"
    path, cfg, env = _setup(tmp_path)
    env.update(JARVIS_GMFIT_ROOT=str(root), JARVIS_GMFIT_NATIVE_PATH=str(native))
    producer = tmp_path / "producer"
    completed = subprocess.run(
        [sys.executable, "-c", GENERATOR, str(root), str(native), str(producer)],
        env=env,
        cwd=tmp_path,
        text=True,
        capture_output=True,
        timeout=300,
    )
    (tmp_path / "producer.log").write_text(completed.stdout + completed.stderr)
    assert completed.returncode == 0, completed.stderr[-4000:]
    transport = next(producer.glob("*/transport.json"))
    original_transport = transport.read_bytes()
    baseline = json.loads(transport.with_name("result.json").read_text())
    assert baseline["ready"] == 1 and baseline["transport_status"] == "ok"
    env["JARVIS_GMFIT_TRANSPORT_FILE"] = str(transport)
    original_card = yaml.safe_load((root / "bin/GMFit_SPheno_h1_HS.yaml").read_text())
    original_module = original_card["Operas"]["Modules"][0]
    cfg["Sampling"]["Bounds"]["point_number"] = 1
    cfg["EnvReqs"]["V2"]["workers"] = 1
    cfg["Operas"]["Modules"] = [
        {
            "name": "LoadReport",
            "operator": "gmfit_unit_ops.source_report",
            "output": [
                {"name": field, "entry": field}
                for field in ("spheno_h1_inputs", "spheno_h1_mass", "spheno_h1_width")
            ],
        }
    ]
    for tag, unit, factor in (
        ("gev", "GeV", 1.0),
        ("tev", "TeV", 0.001),
        ("mev", "MeV", 1000.0),
        ("mismatch", "GeV", 1.0),
    ):
        module = copy.deepcopy(original_module)
        module.update(
            name=tag,
            operator="gmfit_unit_ops.fit_h1_signals",
            required_modules=["LoadReport"],
            unit_policy="strict",
            call_mode="acall" if tag == "tev" else "call",
            kwargs={"log_detail": "none"},
        )
        for binding in module["input"]:
            source = binding.pop("entry")
            binding.update(expression=f"{source} * {factor}", unit=unit)
            if tag == "mismatch" and binding["name"] == "mass_GeV":
                binding["expression"] = f"{source} + 1"
        for binding in module["output"]:
            binding["name"] = tag + "_" + binding["name"]
        module["output"].extend(
            [
                {"name": tag + "_canonical", "entry": ""},
                {
                    "name": tag + "_xs_fb",
                    "entry": "quantities.xs_ggH.value",
                    "unit": "fb",
                },
            ]
        )
        cfg["Operas"]["Modules"].append(module)
    path.write_text(yaml.safe_dump(cfg))
    output = tmp_path / "outcome.json"
    source = Path(__file__).parent / "fixtures/gmfit_unit_ops.py"
    completed = subprocess.run(
        [sys.executable, "-c", DRIVER, str(path), str(source.resolve()), str(output)],
        cwd=tmp_path,
        env=env,
        text=True,
        capture_output=True,
        timeout=180,
    )
    (tmp_path / "host.log").write_text(completed.stdout + completed.stderr)
    assert completed.returncode == 0, (completed.stdout + completed.stderr)[-6000:]
    assert json.loads(output.read_text())["failed"] == 0
    database = tmp_path / "outputs/unit_scan/DATABASE/samples.hdf5"
    rows = SimpleHDF5Writer(str(database)).read_records()
    assert len(rows) == 1
    row = rows[0]
    canonical = row["gev_canonical"]
    for tag in ("gev", "tev", "mev"):
        actual = row[tag + "_canonical"]
        assert actual["ready"] == 1 and actual["nobs"] == baseline["nobs"]
        assert actual["chi2"] == pytest.approx(baseline["chi2"], rel=1e-10, abs=1e-10)
        assert actual["loglike"] == pytest.approx(
            baseline["loglike"], rel=1e-10, abs=1e-10
        )
        assert actual["prediction_id"] == canonical["prediction_id"]
        assert actual["quantities"] == canonical["quantities"]
        record = actual["quantities"]["xs_ggH"]
        assert record["unit"] == "pb"
        assert row[tag + "_xs_fb"] == pytest.approx(1000 * record["value"], rel=1e-13)
        assert row["_jarvis_units"]["observable_units"]["/" + tag + "_xs_fb"] == "fb"
    assert row["mismatch_canonical"]["ready"] == 0
    assert row["mismatch_canonical"]["loglike"] is None
    assert "disagree" in row["mismatch_canonical"]["diagnostics"][0]["message"]
    assert transport.read_bytes() == original_transport
    # Canonical snapshots are valid JSON cells larger than csv's default 128 KiB.
    previous_limit = csv.field_size_limit()
    try:
        csv.field_size_limit(16 * 1024 * 1024)
        with database.with_suffix(".csv").open() as stream:
            exported = list(csv.DictReader(stream))
    finally:
        csv.field_size_limit(previous_limit)
    assert json.loads(exported[0]["_jarvis_units"]) == row["_jarvis_units"]
