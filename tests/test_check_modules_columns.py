#!/usr/bin/env python3
"""D9.5: check-modules CSV column inference."""

from __future__ import annotations

import unittest

from jarvishep2.testing.check_modules import (
    build_check_module_u_coords,
    coordinate_columns_for_check_modules,
)


class CheckModulesColumnTests(unittest.TestCase):
    def test_columns_from_sampling_variables(self) -> None:
        config = {
            "Sampling": {
                "Variables": [
                    {"name": "m0", "distribution": {"type": "Flat", "parameters": {"min": 0, "max": 1}}},
                    {"name": "m12", "distribution": {"type": "Flat", "parameters": {"min": 0, "max": 1}}},
                ]
            }
        }
        cols = coordinate_columns_for_check_modules(config, ["uuid", "m0", "m12", "extra"])
        self.assertEqual(cols, ["m0", "m12"])

    def test_columns_from_csv_headers_when_no_variables(self) -> None:
        cols = coordinate_columns_for_check_modules({}, ["uuid", "a", "b"])
        self.assertEqual(cols, ["a", "b"])

    def test_build_u_coords_order(self) -> None:
        row = {"m0": "0.1", "m12": "0.2", "uuid": "u1"}
        coords = build_check_module_u_coords(row, ["m0", "m12"])
        self.assertEqual(coords, [0.1, 0.2])


if __name__ == "__main__":
    unittest.main()
