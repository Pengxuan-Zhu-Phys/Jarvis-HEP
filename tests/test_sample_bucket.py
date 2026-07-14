#!/usr/bin/env python3
"""SAMPLE bucket allocator + tar packing (V1 sample_directory parity)."""

from __future__ import annotations

import os
import tarfile
import tempfile
import unittest

import fakeredis

from jarvishep2.redis_queue import RedisQueue
from jarvishep2.runtime_config import (
    get_cleanup_config,
    get_sample_directory_config,
    handoff_to_staging_enabled,
    pack_buckets_enabled,
)
from jarvishep2.sample_bucket import pack_bucket_dir
from jarvishep2.task_config import load_task_yaml


class SampleBucketDefaultsTests(unittest.TestCase):
    def test_defaults_are_direct_handoff_and_bucket_pack(self) -> None:
        self.assertEqual(get_cleanup_config({})["strategy"], "direct")
        self.assertFalse(handoff_to_staging_enabled({}))
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
                    "    cleanup:\n"
                    "      strategy: direct\n"
                    "    archiver:\n"
                    "      handoff: direct\n"
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
        self.assertEqual(config["Scan"]["sample_directory"]["limit"], 3)
        self.assertEqual(config["Scan"]["sample_directory"]["width"], 4)
        self.assertEqual(config["Calculators"]["Cleanup"]["strategy"], "direct")
        self.assertEqual(config["Calculators"]["Archiver"]["handoff"], "direct")
        self.assertTrue(config["Calculators"]["Archiver"]["pack_buckets"])


class SampleBucketRedisTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
