#!/usr/bin/env python3
"""`locate --choose <key>=<file>` accepts a run whose name carries the label followed by
letters (a digest suffix: label S1, run `...S1CT...`), 2.15 lane X6 (after-run F1). Such a run
is never a candidate on its own (the label must be a delimited token, so KG1 does not find
KG13), but staff naming it outright is a decision, not a guess. A digit after the label (KG13
for KG1) stays refused. Stdlib only, no network."""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "scripts"))
sys.path.insert(0, HERE)

import test_core_submission as tcs             # noqa: E402

tims, write_summary, raw, run, gates = tcs.tims, tcs.write_summary, tcs.raw, tcs.run, tcs.gates


class ChooseSuffixedLabel(tcs.TestLocate):
    # inherits setUp/helpers only; its own tests are not re-run (see load_tests)
    def suffixed(self, uid="KG1CT", n=9):
        return tims(self.tmp, "03122025", uid, n)

    def test_choose_takes_a_suffixed_run_beside_a_plain_one(self):
        files = self.clean_tree()
        sfx = self.suffixed()
        summary, _ = write_summary(self.tmp)
        rc, out, p = self.locate(summary, "--choose", f"KG1={sfx}")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertIn(sfx, self.files_txt())
        self.assertNotIn(files["KG1"], self.files_txt())

    def test_choose_takes_a_suffixed_run_when_no_run_carries_the_plain_label(self):
        files = self.clean_tree()
        os.remove(files["KG1"]) if os.path.isfile(files["KG1"]) else __import__("shutil").rmtree(files["KG1"])
        sfx = self.suffixed()
        summary, _ = write_summary(self.tmp)
        rc, out, p = self.locate(summary, "--choose", f"KG1={sfx}")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertIn(sfx, self.files_txt())

    def test_a_digit_after_the_label_is_still_refused(self):
        files = self.clean_tree()
        rc, out, _ = self.locate(write_summary(self.tmp)[0], "--choose", f"KG1={files['KG13']}")
        self.assertEqual(rc, 2)
        self.assertIn("not one of the", out["error"])


def load_tests(loader, tests, pattern):
    suite = unittest.TestSuite()
    suite.addTests(loader.loadTestsFromNames(
        [f"{__name__}.ChooseSuffixedLabel.{n}" for n in
         ("test_choose_takes_a_suffixed_run_beside_a_plain_one",
          "test_choose_takes_a_suffixed_run_when_no_run_carries_the_plain_label",
          "test_a_digit_after_the_label_is_still_refused")]))
    return suite


if __name__ == "__main__":
    unittest.main()
