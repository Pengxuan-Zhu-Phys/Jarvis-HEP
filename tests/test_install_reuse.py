#!/usr/bin/env python3
"""clone_shadow install stamp: skip reinstall when calculator settings match."""

from __future__ import annotations

import os
import tempfile
import unittest
from unittest import mock

from jarvishep2.Module.calculator_spec import CalculatorSpec
from jarvishep2.Module.runtime_preparer import (
    INSTALL_STAMP_BASENAME,
    RuntimePreparer,
    build_install_fingerprint,
    force_calc_install_requested,
    read_install_stamp,
    write_install_stamp,
)


def _shadow_spec(tmpdir: str, *, install_cmds: list[str] | None = None) -> CalculatorSpec:
    source = os.path.join(tmpdir, "source")
    os.makedirs(source, exist_ok=True)
    with open(os.path.join(source, "marker.txt"), "w", encoding="utf-8") as handle:
        handle.write("src\n")
    runtime = os.path.join(tmpdir, "runtime", "@PackID")
    cmds = install_cmds or [f"cp -r ${{source}}/* ${{path}}"]
    return CalculatorSpec.from_config(
        "DemoCalc",
        {
            "name": "DemoCalc",
            "clone_shadow": True,
            "path": runtime,
            "source": source,
            "installation": cmds,
            "initialization": [],
            "execution": {"commands": ["true"]},
        },
    )


class InstallFingerprintTests(unittest.TestCase):
    def test_fingerprint_stable_for_same_settings(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            spec = _shadow_spec(tmp)
            a = build_install_fingerprint(spec)
            b = build_install_fingerprint(spec)
            self.assertEqual(a, b)
            self.assertEqual(len(a), 64)

    def test_fingerprint_changes_when_install_cmd_changes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            s1 = _shadow_spec(tmp, install_cmds=["cp -r ${source}/* ${path}"])
            s2 = _shadow_spec(tmp, install_cmds=["cp -R ${source}/. ${path}"])
            self.assertNotEqual(
                build_install_fingerprint(s1), build_install_fingerprint(s2)
            )


class InstallReuseTests(unittest.TestCase):
    def test_second_ensure_skips_run_stage(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            spec = _shadow_spec(tmp)
            prep = RuntimePreparer(spec)
            prep.acquire_pack_id("001")
            calls: list[str] = []

            def run_stage(commands, stage: str) -> None:
                calls.append(stage)
                # Mimic install: materialize marker into pack.
                runtime = prep.shadow_runtime_path()
                os.makedirs(runtime, exist_ok=True)
                with open(
                    os.path.join(runtime, "marker.txt"), "w", encoding="utf-8"
                ) as handle:
                    handle.write("installed\n")

            prep.ensure_shadow_installed(run_stage=run_stage)
            self.assertEqual(calls, ["install"])
            self.assertEqual(prep.last_install_action, "installed")
            stamp = read_install_stamp(prep.shadow_runtime_path())
            self.assertIsNotNone(stamp)
            self.assertEqual(stamp["pack_id"], "001")
            self.assertTrue(
                os.path.isfile(
                    os.path.join(prep.shadow_runtime_path(), INSTALL_STAMP_BASENAME)
                )
            )

            # Simulate new Worker process: fresh preparer, same disk pack.
            prep2 = RuntimePreparer(spec)
            prep2.acquire_pack_id("001")
            calls2: list[str] = []

            def run_stage2(commands, stage: str) -> None:
                calls2.append(stage)

            prep2.ensure_shadow_installed(run_stage=run_stage2)
            self.assertEqual(calls2, [])
            self.assertEqual(prep2.last_install_action, "reused")

    def test_force_reinstall_runs_again(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            spec = _shadow_spec(tmp)
            prep = RuntimePreparer(spec)
            prep.acquire_pack_id("002")
            n = {"install": 0}

            def run_stage(commands, stage: str) -> None:
                n[stage] = n.get(stage, 0) + 1
                runtime = prep.shadow_runtime_path()
                os.makedirs(runtime, exist_ok=True)

            prep.ensure_shadow_installed(run_stage=run_stage)
            prep2 = RuntimePreparer(spec, force_reinstall=True)
            prep2.acquire_pack_id("002")
            prep2.ensure_shadow_installed(run_stage=run_stage)
            self.assertEqual(n["install"], 2)
            self.assertEqual(prep2.last_install_action, "installed")

    def test_env_force_reinstall(self) -> None:
        self.assertFalse(force_calc_install_requested({}))
        self.assertTrue(force_calc_install_requested({"JARVIS_FORCE_CALC_INSTALL": "1"}))

    def test_changed_settings_trigger_reinstall(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            s1 = _shadow_spec(tmp, install_cmds=["echo one"])
            prep = RuntimePreparer(s1)
            prep.acquire_pack_id("001")
            calls: list[str] = []

            def run_stage(commands, stage: str) -> None:
                calls.append(str(commands[0].get("cmd")))
                os.makedirs(prep.shadow_runtime_path(), exist_ok=True)

            prep.ensure_shadow_installed(run_stage=run_stage)
            s2 = _shadow_spec(tmp, install_cmds=["echo two"])
            prep2 = RuntimePreparer(s2)
            prep2.acquire_pack_id("001")
            prep2.ensure_shadow_installed(run_stage=run_stage)
            self.assertEqual(calls, ["echo one", "echo two"])


if __name__ == "__main__":
    unittest.main()
