#!/usr/bin/env python3
"""SAMPLE bucket allocator + tar packing (V1 sample_directory parity)."""

from __future__ import annotations

import os
import tarfile
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import fakeredis

from jarvishep2.redis_queue import RedisQueue
from jarvishep2.runtime._runtime_supervisor import _RuntimeSupervisor
from jarvishep2.runtime_config import (
    get_archiver_config,
    get_sample_directory_config,
    normalize_archiver_block,
    pack_buckets_enabled,
)
from jarvishep2.sample_bucket import normalize_sample_directory, pack_bucket_dir
from jarvishep2.task_config import load_task_yaml


class SampleBucketDefaultsTests(unittest.TestCase):
    def test_defaults_write_directly_and_pack_buckets(self) -> None:
        self.assertNotIn("handoff", get_archiver_config({}))
        sample_dir = get_sample_directory_config({})
        self.assertTrue(sample_dir["enabled"])
        self.assertEqual(sample_dir["limit"], 200)
        self.assertTrue(pack_buckets_enabled({}))

    def test_envreqs_default_yaml_merges_sample_directory(self) -> None:
        with tempfile.TemporaryDirectory() as root:
            os.makedirs(os.path.join(root, "bin"))
            os.makedirs(os.path.join(root, "deps"))
            with open(os.path.join(root, "jarvis.project.yaml"), "w", encoding="utf-8") as handle:
                handle.write("project: bucket-defaults\n")
            defaults = os.path.join(root, "deps", "environment_default.yaml")
            with open(defaults, "w", encoding="utf-8") as handle:
                handle.write(
                    "EnvReqs:\n"
                    "  V2:\n"
                    "    workers: 2\n"
                    "    sample_directory:\n"
                    "      limit: 3\n"
                    "      width: 4\n"
                    "      pack: true\n"
                    "    archiver:\n"
                    "      pack_buckets: true\n"
                )
            task_path = os.path.join(root, "bin", "task.yaml")
            with open(task_path, "w", encoding="utf-8") as handle:
                handle.write(
                    "Scan:\n  name: bucket-test\n"
                    "EnvReqs:\n"
                    "  Check_default_dependencies:\n"
                    "    required: true\n"
                    "    default_yaml_path: '&J/deps/environment_default.yaml'\n"
                    "Sampling:\n  Method: Random\n"
                )
            config = load_task_yaml(task_path)
        self.assertEqual(config["EnvReqs"]["V2"]["sample_directory"]["limit"], 3)
        self.assertEqual(config["EnvReqs"]["V2"]["sample_directory"]["width"], 4)
        self.assertNotIn("sample_directory", config["Scan"])
        self.assertNotIn("Archiver", config.get("Calculators", {}))
        self.assertNotIn("handoff", config["EnvReqs"]["V2"]["archiver"])
        self.assertTrue(config["EnvReqs"]["V2"]["archiver"]["pack_buckets"])
        self.assertTrue(get_archiver_config(config)["pack_buckets"])

    def test_envreqs_v2_is_the_only_sample_directory_source(self) -> None:
        config = {
            "Scan": {
                "name": "ignored-scan-setting",
                "sample_directory": {"limit": 1},
            },
            "EnvReqs": {
                "V2": {"sample_directory": {"limit": 9, "width": 4}},
            },
        }
        sample_dir = get_sample_directory_config(config)
        self.assertEqual(sample_dir["limit"], 9)
        self.assertEqual(sample_dir["width"], 4)

    def test_envreqs_v2_sample_pack_gate_is_effective(self) -> None:
        config = {
            "EnvReqs": {"V2": {"sample_directory": {"pack": False}}},
        }
        self.assertFalse(pack_buckets_enabled(config))

    def test_string_booleans_are_rejected_before_pack_or_delete_policy(self) -> None:
        with self.assertRaisesRegex(ValueError, "sample_directory.enabled"):
            normalize_sample_directory({"enabled": "false"})
        with self.assertRaisesRegex(ValueError, "sample_directory.pack"):
            normalize_sample_directory({"pack": "false"})
        with self.assertRaisesRegex(ValueError, "archiver.delete_after_archive"):
            normalize_archiver_block({"delete_after_archive": "false"})
        with self.assertRaisesRegex(ValueError, "archiver.pack_buckets"):
            normalize_archiver_block({"pack_buckets": "false"})


class SampleBucketRedisTests(unittest.TestCase):
    def test_resume_initializes_above_existing_bucket_numbers(self) -> None:
        cases = (
            ("empty", None, None, 1, 1),
            ("directory", "000004/", None, 1, 5),
            ("archive", "000007.tar.gz", None, 1, 8),
            ("suffixed archive", "000009__0001.tar.gz", None, 1, 10),
            ("temporary archive", "000012__0001.tar.gz.tmp", None, 1, 13),
            ("manifest ahead", "000007.tar.gz", '{"bucket":"000016"}\n{"bucket":', 1, 17),
            ("configured floor", "000007.tar.gz", '{"bucket":"000016"}\n', 20, 20),
            ("unrelated file", "999999.txt", None, 1, 1),
        )
        for case, entry, manifest, start, expected in cases:
            with self.subTest(case=case), tempfile.TemporaryDirectory() as tmp:
                sample_root = Path(tmp) / "SAMPLE"
                sample_root.mkdir()
                if entry is not None:
                    path = sample_root / entry
                    if entry.endswith("/"):
                        path.mkdir()
                    else:
                        path.touch()
                if manifest is not None:
                    path = Path(tmp) / "DATABASE" / "archive_manifest.jsonl"
                    path.parent.mkdir()
                    path.write_text(manifest, encoding="utf-8")
                redis = RedisQueue({"host": "x", "port": 1, "db": 0})
                redis.r = fakeredis.FakeRedis(decode_responses=True)
                core = SimpleNamespace(
                    config={"EnvReqs": {"V2": {"store_samples": True, "sample_directory": {"start_bucket": start}}}},
                    info={}, redis=redis, _resume_policy="resume",
                    _resolve_sample_root=lambda: str(sample_root), _logger=Mock(),
                )
                _RuntimeSupervisor(core)._init_sample_buckets()
                allocation = redis.allocate_sample_bucket()
                assert allocation is not None
                self.assertEqual(allocation["bucket_id"], expected)

    def test_allocate_finish_seal_and_pack(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sample_root = os.path.join(tmp, "SAMPLE")
            q = RedisQueue({"host": "x", "port": 1, "db": 0})
            q.r = fakeredis.FakeRedis(decode_responses=True)
            q.init_sample_buckets(
                sample_root=sample_root,
                limit=2,
                width=6,
                pack=True,
                enabled=True,
            )
            a1 = q.allocate_sample_bucket()
            a2 = q.allocate_sample_bucket()
            assert a1 is not None and a2 is not None
            self.assertEqual(a1["bucket_id"], 1)
            self.assertEqual(a2["bucket_id"], 1)
            a3 = q.allocate_sample_bucket()
            assert a3 is not None
            self.assertEqual(a3["bucket_id"], 2)
            self.assertEqual(a3["sealed_id"], 1)

            # Worker finish alone must NOT enqueue pack (Archiver may still lag).
            self.assertFalse(q.finish_sample_bucket(1))
            self.assertFalse(q.finish_sample_bucket(1))
            self.assertIsNone(q.pull_ready_bucket())
            # Pack only after every assigned sample is archived.
            self.assertFalse(q.note_sample_archived(1))
            self.assertTrue(q.note_sample_archived(1))
            ready = q.pull_ready_bucket()
            assert ready is not None
            self.assertEqual(ready["bucket_id"], 1)

            for idx, alloc in enumerate((a1, a2), start=1):
                sample_dir = os.path.join(alloc["bucket_dir"], f"uuid-{idx}")
                os.makedirs(sample_dir)
                with open(os.path.join(sample_dir, "marker.txt"), "w", encoding="utf-8") as handle:
                    handle.write("ok")

            with self.assertRaises(InterruptedError):
                pack_bucket_dir(ready["bucket_dir"], cancelled=lambda: True)
            self.assertTrue(os.path.isdir(ready["bucket_dir"]))
            self.assertFalse(os.path.exists(ready["bucket_dir"] + ".tar.gz.tmp"))

            tar_path = pack_bucket_dir(
                ready["bucket_dir"],
                sample_root=sample_root,
                prune=True,
            )
            self.assertTrue(tar_path and tar_path.endswith("000001.tar.gz"))
            self.assertTrue(os.path.isfile(str(tar_path)))
            self.assertFalse(os.path.isdir(ready["bucket_dir"]))
            with tarfile.open(str(tar_path), "r:gz") as tar:
                names = tar.getnames()
            self.assertTrue(any(name.endswith("marker.txt") for name in names))

            # A pre-existing tar is not evidence that these files were packed.
            old_archive = Path(str(tar_path)).read_bytes()
            new_sample = Path(ready["bucket_dir"]) / "new-uuid"
            new_sample.mkdir(parents=True)
            (new_sample / "result.txt").write_text("new output", encoding="utf-8")
            new_tar = pack_bucket_dir(ready["bucket_dir"])
            self.assertEqual(new_tar, ready["bucket_dir"] + "__0001.tar.gz")
            self.assertEqual(Path(str(tar_path)).read_bytes(), old_archive)
            with tarfile.open(str(new_tar), "r:gz") as tar:
                self.assertIn("000001/new-uuid/result.txt", tar.getnames())
            self.assertFalse(new_sample.exists())


if __name__ == "__main__":
    unittest.main()
