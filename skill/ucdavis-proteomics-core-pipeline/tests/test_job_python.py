#!/usr/bin/env python3
"""
A job's Python is the pipeline env's, so the Sage job's LFQ check can run (T22 of the
2.14.2 plan).

_job_python() baked in the GENERATOR's sys.executable. A search generated under /usr/bin/python3
(no pyarrow) gave the job that interpreter, and the Sage job's sage_lfq_check.py recorded
"unchecked: pyarrow is not installed in this Python" (a Core submission, 2026-10-29). The job now takes
the first executable of:
  1. ${PROTEOMICS_PIPELINE_HOME:-$HOME/.proteomics-pipeline}/micromamba/envs/proteomics-pipeline/
     bin/python -- the pipeline env, as acquire_tools.sh finds it, expanded IN THE JOB (HIVE's
     $HOME, not the laptop's);
  2. the generator's sys.executable;
  3. python3.
The generated line is run in bash under the job's own `set -euo pipefail`. Hermetic.
"""
import os
import shlex
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "scripts"))

import run_search  # noqa: E402

ENV_PY = ("micromamba", "envs", "proteomics-pipeline", "bin", "python")


@unittest.skipUnless(shutil.which("bash"), "the line is bash")
class JobPython(unittest.TestCase):
    def setUp(self):
        self.d = os.path.realpath(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.d, True)

    def pick(self, env, line=None):
        """$_py after the generated line, in bash under set -euo pipefail."""
        script = f"set -euo pipefail\n{line or run_search._job_python()}\nprintf '%s' \"$_py\"\n"
        r = subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=env,
                           timeout=30)
        self.assertEqual(r.returncode, 0, r.stderr)
        return r.stdout

    def make_env_python(self, home):
        py = os.path.join(home, *ENV_PY)
        os.makedirs(os.path.dirname(py))
        with open(py, "w") as fh:
            fh.write("#!/bin/sh\n")
        os.chmod(py, os.stat(py).st_mode | stat.S_IXUSR)
        return py

    def test_job_prefers_the_pipeline_env_python(self):
        home = os.path.join(self.d, "pp home")              # a space: still one word
        py = self.make_env_python(home)
        self.assertEqual(self.pick({"PATH": "/usr/bin:/bin", "HOME": self.d,
                                    "PROTEOMICS_PIPELINE_HOME": home}), py)

    def test_without_the_variable_the_jobs_home_is_used(self):
        py = self.make_env_python(os.path.join(self.d, ".proteomics-pipeline"))
        self.assertEqual(self.pick({"PATH": "/usr/bin:/bin", "HOME": self.d}), py)

    def test_with_no_env_it_is_the_generators_python(self):
        """(pin) The generator's interpreter, as before, when there is no pipeline env."""
        self.assertEqual(self.pick({"PATH": "/usr/bin:/bin", "HOME": self.d,
                                    "PROTEOMICS_PIPELINE_HOME": os.path.join(self.d, "none")}),
                         sys.executable)

    def test_then_python3(self):
        """(pin) python3 when neither exists."""
        line = run_search._job_python().replace(shlex.quote(sys.executable),
                                                shlex.quote(os.path.join(self.d, "gone")))
        self.assertEqual(self.pick({"PATH": "/usr/bin:/bin", "HOME": self.d}, line), "python3")

    def test_an_unset_home_does_not_end_the_job(self):
        """(pin) set -u: a job without $HOME falls through instead of failing on the
        expansion the env path adds."""
        self.assertEqual(self.pick({"PATH": "/usr/bin:/bin"}), sys.executable)

    def test_the_path_is_expanded_in_the_job_not_here(self):
        line = run_search._job_python()
        self.assertIn("${PROTEOMICS_PIPELINE_HOME:-", line)
        # the fallback is this interpreter by its own path, which on HIVE IS under
        # ~/.proteomics-pipeline; only the env path must be left for the job to expand
        env_part = line.replace(shlex.quote(sys.executable), "<this python>")
        self.assertNotIn(os.path.expanduser("~") + "/.proteomics-pipeline", env_part)


if __name__ == "__main__":
    unittest.main()
