"""S2b: project archives carry exact user sampler sources without executing them."""

from __future__ import annotations

import json
import tarfile
from pathlib import Path

import pytest
import yaml
from jarvis_hep_sampler import SamplerRegistry

from jarvishep2.project_packager import (
    create_project_package,
    create_project_package_from_manifest,
    create_project_pack_manifest,
)
from jarvishep2.project_samplers import INDEX, inspect_packaged_samplers
from jarvishep2.official_project_library import (
    OfficialProjectFetchError,
    fetch_official_project,
)
from jarvishep2.sampling.plugin_sampler import preflight_plugin

from test_sampler_plugins import HALF_CUBE, _config


@pytest.mark.parametrize("profile", ["share", "repro", "full", "manifest"])
def test_pack_export_fetch_conflict_and_integrity(tmp_path, monkeypatch, profile):
    monkeypatch.setenv("JARVIS_SAMPLERS_HOME", str(tmp_path / "registry"))
    project = tmp_path / "project"
    (project / "bin").mkdir(parents=True)
    (project / ".jarvis-project.json").write_text("{}")
    source = tmp_path / "HalfCube.py"
    source.write_text(HALF_CUBE)
    first = SamplerRegistry().load(source).entry
    cfg = _config(project, "HalfCube", {"point_number": 8, "upper": 0.8})
    (project / "bin/task.yaml").write_text(
        yaml.safe_dump({"Scan": cfg["Scan"], "Sampling": cfg["Sampling"]})
    )
    preflight = preflight_plugin(
        cfg, task_root=str(project), scan_name=cfg["scan_name"]
    )
    source.write_text(HALF_CUBE.replace('"version": "1.0"', '"version": "2.0"'))
    second = SamplerRegistry().load(source).entry
    if profile == "manifest":
        manifest = create_project_pack_manifest(
            str(project), manifest_dir=str(tmp_path)
        )
        doc = yaml.safe_load(Path(manifest.manifest_path).read_text())
        doc["output"] = str(tmp_path / "manifest.tar.gz")
        Path(manifest.manifest_path).write_text(yaml.safe_dump(doc))
        report = create_project_package_from_manifest(manifest.manifest_path)
    else:
        report = create_project_package(str(project), profile)
    target = tmp_path / "unpacked"
    target.mkdir()
    with tarfile.open(report.archive_path) as archive:
        archive.extractall(target, filter="data")
    root = target if profile == "manifest" else target / "project"
    index = json.loads((root / INDEX).read_text())
    digest = second.source_sha256 if profile == "share" else first.source_sha256
    assert index["entries"][0]["source_sha256"] == digest
    payload = root / index["entries"][0]["file"]
    assert payload.is_file()
    # Reading/fetching assets must not execute the user's top-level statements.
    marker = tmp_path / "imports.txt"
    monkeypatch.setenv("JHS_TEST_MARKER", str(marker))
    notices = inspect_packaged_samplers(str(root))
    assert not marker.exists()
    if profile == "share":
        assert "already loaded" in notices[0]
    else:
        assert "different source" in notices[0] and "Jarvis sampler load" in notices[0]
    SamplerRegistry().remove("HalfCube")
    assert "not loaded" in inspect_packaged_samplers(str(root))[0]
    payload.write_text(payload.read_text() + "\n# corrupt\n")
    with pytest.raises(ValueError, match="hash mismatch"):
        inspect_packaged_samplers(str(root))
    assert preflight.source_copy.exists()
    # Repacking a fetched project replaces stale assets and does not duplicate names.
    stale = root / ".jarvis-pack/samplers/stale.py"
    stale.write_text("old source")
    SamplerRegistry().load(source)
    repacked = create_project_package(str(root), "share")
    with tarfile.open(repacked.archive_path) as archive:
        members = [m.name for m in archive.getmembers()]
        assert len(members) == len(set(members))
        assert not any(name.endswith("/stale.py") for name in members)


def test_same_method_multiple_scan_revisions_are_preserved(tmp_path, monkeypatch):
    monkeypatch.setenv("JARVIS_SAMPLERS_HOME", str(tmp_path / "registry"))
    project = tmp_path / "project"
    project.mkdir()
    (project / ".jarvis-project.json").write_text("{}")
    source = tmp_path / "HalfCube.py"
    digests = []
    for i in (1, 2):
        source.write_text(HALF_CUBE.replace('"version": "1.0"', f'"version": "{i}.0"'))
        SamplerRegistry().load(source)
        cfg = _config(project, "HalfCube", {"point_number": 8})
        cfg["scan_name"] = f"scan{i}"
        cfg["Scan"]["name"] = f"scan{i}"
        preflight = preflight_plugin(cfg, task_root=str(project), scan_name=f"scan{i}")
        digests.append(preflight.plugin.source_sha256)
    report = create_project_package(str(project), "repro")
    with tarfile.open(report.archive_path) as archive:
        index = json.load(archive.extractfile(f"project/{INDEX}"))
    assert {entry["source_sha256"] for entry in index["entries"]} == set(digests)


@pytest.mark.parametrize("corrupt", [False, True])
def test_actual_fetch_verifies_before_installing_project(
    tmp_path, monkeypatch, capsys, corrupt
):
    monkeypatch.setenv("JARVIS_SAMPLERS_HOME", str(tmp_path / "registry"))
    project = tmp_path / "original"
    (project / "bin").mkdir(parents=True)
    (project / ".jarvis-project.json").write_text("{}")
    source = tmp_path / "HalfCube.py"
    source.write_text(HALF_CUBE)
    SamplerRegistry().load(source)
    cfg = _config(project, "HalfCube", {"point_number": 8})
    (project / "bin/run.yaml").write_text(
        yaml.safe_dump({"Scan": cfg["Scan"], "Sampling": cfg["Sampling"]})
    )
    report = create_project_package(str(project), "share")
    archive_path = Path(report.archive_path)
    if corrupt:
        extracted = tmp_path / "edit"
        with tarfile.open(archive_path) as archive:
            archive.extractall(extracted, filter="data")
        payload = next((extracted / "original/.jarvis-pack/samplers").glob("*.py"))
        payload.write_text(payload.read_text() + "\n# damaged\n")
        with tarfile.open(archive_path, "w:gz") as archive:
            archive.add(extracted / "original", arcname="original")
    SamplerRegistry().remove("HalfCube")
    marker = tmp_path / "imported"
    monkeypatch.setenv("JHS_TEST_MARKER", str(marker))
    monkeypatch.setattr(
        "jarvishep2.official_project_library._load_catalog_payload",
        lambda **kwargs: {
            "schema_version": 1,
            "projects": [
                {
                    "name": "SamplerDemo",
                    "access": "public",
                    "archive_url": archive_path.as_uri(),
                    "archive_root": "original",
                    "entrypoint": "bin/run.yaml",
                }
            ],
        },
    )
    target = tmp_path / "fetched"
    if corrupt:
        with pytest.raises(OfficialProjectFetchError, match="hash mismatch"):
            fetch_official_project(
                "SamplerDemo", target_dir=str(target), progress=False
            )
        assert not target.exists()
    else:
        fetch_official_project("SamplerDemo", target_dir=str(target), progress=False)
        assert (target / INDEX).is_file()
        output = capsys.readouterr().out
        assert "Jarvis sampler load" in output and str(target) in output
    assert not marker.exists()
    assert SamplerRegistry().get("HalfCube") is None
