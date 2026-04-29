#!/usr/bin/env python3
from __future__ import annotations

import unittest

from jarvishep.IOs.IOs import Parameter as CanonicalParameter
from jarvishep.IOs.parameter import Parameter as CompatParameter


class IOParameterCompatTests(unittest.TestCase):
    def test_legacy_parameter_module_reexports_canonical_class(self):
        self.assertIs(CompatParameter, CanonicalParameter)


if __name__ == "__main__":
    unittest.main()
