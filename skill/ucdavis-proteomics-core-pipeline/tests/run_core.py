#!/usr/bin/env python3
"""Run the CORE test tier: every test module except those listed in tests/SLOW_MODULES.txt.

    python3 tests/run_core.py            # from the skill folder or anywhere
    python3 tests/run_core.py -q         # quieter
    python3 tests/run_core.py --list     # print the modules it would run, run nothing

Plain unittest underneath (no extra packages). Prints the wall time at the end. Run this on every
change; run the FULL suite (python3 -m unittest discover -s tests, or the HIVE runner) before a
release. Exits non-zero if any test fails or a module cannot be imported.
"""
import os
import sys
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SLOW_FILE = os.path.join(HERE, "SLOW_MODULES.txt")


def slow_modules():
    with open(SLOW_FILE, encoding="utf-8") as fh:
        lines = [ln.split("#", 1)[0].strip() for ln in fh]
    return {ln for ln in lines if ln}


def all_modules():
    return sorted(f[:-3] for f in os.listdir(HERE) if f.startswith("test_") and f.endswith(".py"))


def core_modules():
    slow = slow_modules()
    unknown = slow - set(all_modules())
    if unknown:
        raise SystemExit("SLOW_MODULES.txt names no such test module: " + ", ".join(sorted(unknown)))
    return [m for m in all_modules() if m not in slow]


def main(argv):
    mods = core_modules()
    if "--list" in argv:
        print("\n".join(mods))
        return 0
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    suite = unittest.TestLoader().loadTestsFromNames(mods)
    t0 = time.time()
    result = unittest.TextTestRunner(verbosity=0 if "-q" in argv else 1).run(suite)
    print("core tier: %d modules (%d slow modules skipped), %d tests, %.0f s"
          % (len(mods), len(slow_modules()), result.testsRun, time.time() - t0))
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
