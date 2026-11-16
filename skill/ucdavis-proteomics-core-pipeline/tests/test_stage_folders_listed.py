#!/usr/bin/env python3
"""`stage` lists the folders it does not choose (2.15 lane X6, after-run F15 and F16).
F15: with --service-dir, and for a folder named for the submission without its underscore,
the dry run used to say nothing about the group's existing folder for this submission. It still
stages into a new folder (the documented way to ask for one), but lists the others and how to
stage into one (--project-dir).
F16: two folders match the PI's surname: stage still refuses to choose, and now lists both
with the exact --service-dir to pass. Stdlib only, no network."""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "scripts"))
sys.path.insert(0, HERE)

import test_core_submission as tcs             # noqa: E402

write_summary, record = tcs.write_summary, tcs.record


class StageListsFolders(tcs.TestStage):
    def test_a_folder_named_for_the_submission_without_the_underscore_is_listed(self):
        group = self.folder("on_campus", "Placeholder lab")
        old = self.folder("on_campus", "Placeholder lab", "Earlier-analysis-PROT0807")
        summary, _ = write_summary(self.tmp)
        for extra in ((), ("--service-dir", "Placeholder lab")):
            rc, out, p = self.stage(summary, *extra)
            self.assertEqual(rc, 0, p.stderr)
            self.assertEqual(out["service_dir"], os.path.join(group, "PROT_0807"))
            listed = out["group_resolution"]["existing_folders_not_used"]
            self.assertEqual([x["dir"] for x in listed], [old])
            self.assertTrue(any("Earlier-analysis-PROT0807/" in w and "--project-dir" in w
                                for w in out["warnings"]), out["warnings"])

    def test_a_longer_number_is_not_this_submissions_folder(self):
        self.folder("on_campus", "Placeholder lab")
        self.folder("on_campus", "Placeholder lab", "Other-PROT08070")
        summary, _ = write_summary(self.tmp)
        rc, out, _ = self.stage(summary)
        self.assertEqual(rc, 0)
        self.assertNotIn("existing_folders_not_used", out["group_resolution"])

    def test_service_dir_still_lists_a_folder_that_holds_the_raw_files(self):
        group = self.folder("on_campus", "Placeholder lab")
        mine = self.folder("on_campus", "Placeholder lab", "Sampled03122025")
        for f in self.files:
            os.symlink(f, os.path.join(mine, os.path.basename(f)))
        summary, _ = write_summary(self.tmp)
        rc, out, _ = self.stage(summary, "--service-dir", "Placeholder lab")
        self.assertEqual((rc, out["service_dir"]), (0, os.path.join(group, "PROT_0807")))
        self.assertEqual([x["dir"] for x in out["group_resolution"]["existing_folders_not_used"]],
                         [mine])

    def test_two_surname_folders_are_listed_with_the_service_dir_to_pass(self):
        self.folder("on_campus", "Doe Jane")
        self.folder("on_campus", "Doe (John)")
        summary, _ = write_summary(self.tmp, record(pi_last="Doe", pi_first="Jane"))
        rc, out, _ = self.stage(summary)
        self.assertEqual(rc, 2)
        opts = out["service_dir_options"]
        self.assertEqual(sorted(o["folder"] for o in opts), ["Doe (John)", "Doe Jane"])
        for o in opts:
            self.assertIn(f"--service-dir '{o['folder']}'", o["rerun_with"])
        self.assertIn("Doe Jane", out["hint"])
        self.assertIn("Doe (John)", out["hint"])


def load_tests(loader, tests, pattern):
    names = [n for n in dir(StageListsFolders) if n.startswith("test_") and
             n in StageListsFolders.__dict__]
    return loader.loadTestsFromNames([f"{__name__}.StageListsFolders.{n}" for n in names])


if __name__ == "__main__":
    unittest.main()
