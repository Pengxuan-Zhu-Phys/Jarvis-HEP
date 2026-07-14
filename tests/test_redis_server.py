#!/usr/bin/env python3
"""Managed redis-server process title (Jarvis-Redis:<scan>)."""

from __future__ import annotations

import socket
import subprocess
import tempfile
import time
import unittest

from jarvishep2.proc_title import redis_title
from jarvishep2.redis_server import ManagedRedisServer, find_redis_server_binary, redis_port_open


def _free_port() -> int:
    sock = socket.socket()
    sock.bind(("127.0.0.1", 0))
    port = int(sock.getsockname()[1])
    sock.close()
    return port


@unittest.skipUnless(find_redis_server_binary(), "redis-server not installed")
class ManagedRedisServerTests(unittest.TestCase):
    def test_process_title_is_jarvis_redis_scan(self) -> None:
        port = _free_port()
        self.assertFalse(redis_port_open("127.0.0.1", port))
        with tempfile.TemporaryDirectory() as tmp:
            managed = ManagedRedisServer(
                host="127.0.0.1",
                port=port,
                scan_name="eggbox",
                work_dir=tmp,
            )
            started = managed.ensure(ready_timeout=10.0)
            self.assertTrue(started)
            self.assertTrue(managed.started_by_us)
            self.assertTrue(managed.is_reachable())
            assert managed.process is not None
            # Give redis a moment to rewrite its title (host:port suffix).
            time.sleep(0.2)
            cmd = subprocess.check_output(
                ["ps", "-p", str(managed.process.pid), "-o", "command="],
                text=True,
            ).strip()
            expected = redis_title(scan_name="eggbox")
            self.assertTrue(
                cmd.startswith(expected) or expected in cmd,
                msg=f"unexpected process title: {cmd!r}",
            )
            managed.stop()
            time.sleep(0.15)
            self.assertFalse(redis_port_open("127.0.0.1", port))

    def test_existing_listener_is_left_alone(self) -> None:
        port = _free_port()
        with tempfile.TemporaryDirectory() as tmp:
            first = ManagedRedisServer(
                host="127.0.0.1",
                port=port,
                scan_name="a",
                work_dir=tmp,
            )
            self.assertTrue(first.ensure())
            second = ManagedRedisServer(
                host="127.0.0.1",
                port=port,
                scan_name="b",
                work_dir=tmp,
            )
            self.assertFalse(second.ensure())
            self.assertFalse(second.started_by_us)
            first.stop()


if __name__ == "__main__":
    unittest.main()
