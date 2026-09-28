#!/usr/bin/env python3
"""
Every script must parse on the oldest Python 3 a user is likely to run it with: macOS's
/usr/bin/python3 is 3.9. report_style.py once put a quote of an f-string's own kind inside
its braces (legal only from Python 3.12), so make_analysis_html.py -- which imports it --
died with a SyntaxError there while the suite, run on 3.13, stayed green.

ast.parse(feature_version=(3, 9)) does not catch that, so this compiles with a real old
interpreter when one is installed, and skips when none is.
"""
import glob
import os
import shutil
import subprocess
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
CANDIDATES = ("/usr/bin/python3", "python3.9", "python3.10", "python3.11")


def old_python():
    """(path, version) of an installed Python 3 older than 3.12, or (None, None)."""
    for c in CANDIDATES:
        exe = c if os.path.isabs(c) else shutil.which(c)
        if not exe or not os.access(exe, os.X_OK):
            continue
        try:
            out = subprocess.run([exe, "-c", "import sys; print(sys.version_info[:2])"],
                                 capture_output=True, text=True, timeout=30)
        except (OSError, subprocess.SubprocessError):
            continue
        v = out.stdout.strip()
        if out.returncode == 0 and v.startswith("(3, ") and int(v[4:].rstrip(")")) < 12:
            return exe, v
    return None, None


class OldPythonSyntax(unittest.TestCase):
    def test_every_script_parses_before_3_12(self):
        exe, v = old_python()
        if not exe:
            self.skipTest("no Python 3.9-3.11 installed to check with")
        code = ("import ast, sys\nbad = []\nfor p in sys.argv[1:]:\n"
                "    try:\n        ast.parse(open(p, encoding='utf-8').read(), p)\n"
                "    except SyntaxError as e:\n        bad.append(f'{p}:{e.lineno}: {e.msg}')\n"
                "print('\\n'.join(bad))\n")
        files = sorted(glob.glob(os.path.join(SCRIPTS, "*.py")))
        out = subprocess.run([exe, "-c", code, *files], capture_output=True, text=True,
                             timeout=120)
        self.assertEqual(out.returncode, 0, out.stderr)
        self.assertEqual(out.stdout.strip(), "", f"does not parse on Python {v} ({exe})")


if __name__ == "__main__":
    unittest.main()
