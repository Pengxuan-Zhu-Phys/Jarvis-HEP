from __future__ import annotations

import unittest

from jarvishep2.worker import Worker


class _FailingRedis:
    def force_release_calc(self, _name: str, _pack: str) -> bool:
        raise ConnectionError("redis unavailable")


class _SuccessfulRedis:
    def force_release_calc(self, _name: str, _pack: str) -> bool:
        return True


class WorkerReleaseTests(unittest.TestCase):
    def test_release_error_keeps_pack_visible_to_watchdog(self) -> None:
        worker = Worker(0, {"host": "127.0.0.1", "port": 6379, "db": 0}, {})
        worker._redis = _FailingRedis()
        worker._held_calc_packs["DemoCalc"] = "001"
        worker._force_release_pack("DemoCalc", "001")
        self.assertEqual(worker._held_calc_packs, {"DemoCalc": "001"})

    def test_release_success_removes_pack_after_confirmation(self) -> None:
        worker = Worker(0, {"host": "127.0.0.1", "port": 6379, "db": 0}, {})
        worker._redis = _SuccessfulRedis()
        worker._held_calc_packs["DemoCalc"] = "001"
        worker._force_release_pack("DemoCalc", "001")
        self.assertEqual(worker._held_calc_packs, {})


if __name__ == "__main__":
    unittest.main()
