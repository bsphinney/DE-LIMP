#!/usr/bin/env python3
"""
The 5-step chain's recovery notes (RECOVERY.md, .recovery.json) go in the analysis session, and
never into a folder that is not one.

submit.sh records the submission with checkpoint.py in the folder two levels above --out --
<session>/output/search -> <session>. A validation chain whose --out was <home>/<work>/search
(2026-10-15) therefore wrote RECOVERY.md and .recovery.json into the HIVE home directory, where
they could replace another run's notes. Now a folder that is not a session
(session.is_session_dir: session.json, or the input/ + output/ layout init makes) gets nothing:
the notes go beside --out, and the generator and submit.sh say where.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)
from job_env import job_env  # noqa: E402

import session  # noqa: E402

FAKE_SBATCH = r"""#!/bin/bash
for last; do :; done
[ -f "$last" ] || { echo "sbatch: error: Unable to open file $last" >&2; exit 1; }
n=$(( $(cat "$SB_COUNTER" 2>/dev/null || echo 800) + 1 )); echo "$n" > "$SB_COUNTER"
echo "$n"
"""


def _exe(path, body):
    with open(path, "w") as fh:
        fh.write(body)
    os.chmod(path, 0o755)
    return path


class RecoveryNotes(unittest.TestCase):
    def chain(self, d, out):
        raws = []
        for i in range(6):
            raws.append(os.path.join(d, "raw", f"s{i}.d"))
            os.makedirs(raws[-1])
        fasta = os.path.join(d, "db.fasta")
        with open(fasta, "w") as fh:
            fh.write(">sp|P1|X\nPEPTIDEKR\n")
        cfg = os.path.join(d, "params.cfg")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                            "--engine", "diann", "--acquisition", "DIA", "--instrument",
                            "timsTOF HT", "--precursor-mz-range", "299.5", "1200.5", "--out", cfg],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(p.returncode, 0, p.stderr)
        g = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                            "--diann", "/bin/true", "--raw", *raws, "--fasta", fasta,
                            "--out", out, "--cfg", cfg, "--no-notify", "--no-fran"],
                           capture_output=True, text=True, timeout=120)
        self.assertEqual(g.returncode, 0, g.stderr)
        fake = os.path.join(d, "fakebin")
        os.makedirs(fake, exist_ok=True)
        _exe(os.path.join(fake, "sbatch"), FAKE_SBATCH)
        env = job_env(d, PATH=fake + os.pathsep + os.environ.get("PATH", ""),
                      SB_COUNTER=os.path.join(d, "counter"))
        s = subprocess.run(["bash", os.path.join(out, "submit.sh")], capture_output=True,
                           text=True, env=env, timeout=120)
        self.assertEqual(s.returncode, 0, s.stdout + s.stderr)
        return g, s

    def test_an_out_outside_a_session_keeps_its_notes_beside_it(self):
        with tempfile.TemporaryDirectory() as d:
            home = os.path.join(d, "home")
            out = os.path.join(home, "work", "search")
            g, s = self.chain(d, out)
            self.assertFalse(session.is_session_dir(home))
            for f in ("RECOVERY.md", ".recovery.json"):
                self.assertTrue(os.path.exists(os.path.join(out, f)), f)
                self.assertFalse(os.path.exists(os.path.join(home, f)),
                                 f"{f} must not be written into a folder that is not a session")
            self.assertEqual(json.loads(g.stdout)["recovery_notes"],
                             os.path.join(out, "RECOVERY.md"))
            self.assertIn("is not an analysis session folder", g.stderr)
            self.assertIn(f"recovery notes written to {out}/RECOVERY.md (beside --out", s.stdout)
            self.assertEqual(json.load(open(os.path.join(out, ".recovery.json")))["session"], out)

    def test_a_session_gets_them_where_it_always_did(self):
        with tempfile.TemporaryDirectory() as d:
            sess = os.path.join(d, "2026-10-15_x")
            for sub in session.SUBDIRS:
                os.makedirs(os.path.join(sess, sub))
            out = os.path.join(sess, "output", "search")
            g, s = self.chain(d, out)
            self.assertTrue(os.path.exists(os.path.join(sess, "RECOVERY.md")))
            self.assertFalse(os.path.exists(os.path.join(out, "RECOVERY.md")))
            self.assertNotIn("is not an analysis session folder", g.stderr)
            self.assertIn(f"recovery notes written to {sess}/RECOVERY.md —", s.stdout)

    def test_what_a_session_is(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertFalse(session.is_session_dir(d))
            os.makedirs(os.path.join(d, "output"))
            self.assertFalse(session.is_session_dir(d), "output/ alone is not a session")
            os.makedirs(os.path.join(d, "input"))
            self.assertTrue(session.is_session_dir(d))
            e = os.path.join(d, "e")
            os.makedirs(e)
            open(os.path.join(e, "session.json"), "w").write("{}")
            self.assertTrue(session.is_session_dir(e))


if __name__ == "__main__":
    unittest.main()
