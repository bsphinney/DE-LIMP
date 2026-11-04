#!/usr/bin/env python3
"""finalize warns when logs/commands.log is missing or behind the analysis (2.14.2 C2, E2E
F11; LEDGER Ruling 84 LOW).

SKILL.md golden rule 5: every command is appended to logs/commands.log as it is run, and the
bundle copies it. The E2E agent wrote none until after its first finalize, and neither finalize
nor `session.py status` said so; the bundle then had no commands.log. Now finalize says it --
a warning, never a refusal -- on stderr, in its JSON's finish_warnings (stdout), in MANIFEST.txt
as [WARN] and in the stamp, when the log is missing or empty, or older than the newest record
the analysis's steps wrote (session.step_records: the workflow manifest, the FASTA sidecar, the
search records, the DE record, the report's record). What these tests pin:
  * missing: the warning, and finalize still finishes (exit 0, stamped);
  * older than the search record: the warning names that record;
  * as new as every step record: no warning.
Stdlib only (FinalizeBase: a no-DE DIA session, finished, fake skill home).
"""
import json
import os
import sys
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import test_finish_gates as tfg            # noqa: E402

PART = "Command log (logs/commands.log)"


class CommandsLogWarning(tfg.FinalizeBase):
    def run_finalize(self):
        r = self.finalize()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)        # a warning, never a gap
        res = json.loads(r.stdout)
        warn = [ln for ln in self.manifest().splitlines() if ln.startswith("[WARN]")]
        stamp = json.loads(tfg.read(self.stamp_path()))
        return r, res, warn, stamp

    def log(self, age_s=0):
        tfg.write(self.p["commands_log"], "[hive --job] python3 scripts/run_search.py ...\n")
        if age_s:
            t = time.time() - age_s
            os.utime(self.p["commands_log"], (t, t))

    def test_a_missing_log_is_warned_of_and_finalize_finishes(self):
        self.finish()
        self.assertFalse(os.path.exists(self.p["commands_log"]))
        r, res, warn, stamp = self.run_finalize()
        self.assertEqual(len(warn), 1, self.manifest())
        self.assertTrue(warn[0].startswith(f"[WARN]    {PART}"), warn)
        self.assertIn("logs/commands.log is missing", warn[0])
        self.assertIn("golden rule 5", warn[0])
        self.assertTrue([w for w in res["finish_warnings"] if "commands.log is missing" in w],
                        res["finish_warnings"])
        self.assertIn("[session] WARNING: logs/commands.log is missing", r.stderr)
        self.assertTrue([w for w in stamp["warnings"] if "commands.log" in w], stamp)
        self.assertTrue(res["finalized"])

    def test_a_log_older_than_the_search_record_is_warned_of(self):
        self.finish()
        self.log(age_s=3600)
        later = time.time() + 120          # the search record is the newest step record
        os.utime(self.p["search_prov"], (later, later))
        _r, res, warn, _stamp = self.run_finalize()
        self.assertEqual(len(warn), 1, self.manifest())
        self.assertIn("is older than output/search/search_provenance.json", warn[0])
        self.assertIn("the newest step record", warn[0])
        self.assertTrue([w for w in res["finish_warnings"] if "is older than" in w])

    def test_a_current_log_is_not_warned_of(self):
        """(pin)"""
        self.finish()
        self.log()
        future = time.time() + 60          # newer than every step record the fixture wrote
        os.utime(self.p["commands_log"], (future, future))
        _r, res, warn, _stamp = self.run_finalize()
        self.assertEqual(warn, [], self.manifest())
        self.assertFalse([w for w in res["finish_warnings"] if "commands.log" in w])


if __name__ == "__main__":
    unittest.main()
