"""The core tier's list is true: every name in SLOW_MODULES.txt is a real test module, and core +
slow together are every test module (nothing is lost between the two tiers)."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_core  # noqa: E402


class Tiers(unittest.TestCase):
    def test_every_slow_name_is_a_test_module(self):
        self.assertFalse(run_core.slow_modules() - set(run_core.all_modules()))

    def test_core_and_slow_are_every_module(self):
        self.assertEqual(set(run_core.core_modules()) | run_core.slow_modules(),
                         set(run_core.all_modules()))
        self.assertFalse(set(run_core.core_modules()) & run_core.slow_modules())


if __name__ == "__main__":
    unittest.main()
