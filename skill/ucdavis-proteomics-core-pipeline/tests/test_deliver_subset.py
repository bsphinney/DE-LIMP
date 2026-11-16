#!/usr/bin/env python3
"""`core_submission.py deliver` for a declared plate SUBSET (2.15, lane X6).

A Core HT plate stages every file of the submission; staff analyse a subset (the study samples,
the NIST / QC runs left out) in a session that is not under the staged work dir. Before 2.15
deliver refused that session ("not under the staged work dir ... does not match the N staged
files"), and --force did not bypass it, so no delivery plan or email draft could be made.
Pinned here:
  * a session whose raw files are a SUBSET of the staged files belongs to the submission; the
    plan, delivery.json and the email draft say which staged files were left out;
  * a session with a file the stage record does not have (foreign, or a superset) still refuses,
    and the message names the sanctioned path;
  * a session that searched exactly the staged files is delivered as before, with no left-out line;
  * --force bypasses only a MISSING stage record, and its help says so.
Stdlib only, no network.
"""
import json
import os
import shutil
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import core_submission as cs                   # noqa: E402
import test_core_submission as tcs             # noqa: E402

write, read, run = tcs.write, tcs.read, tcs.run


class DeliverSubset(tcs.DeliverBase):
    def elsewhere(self, raws):
        """A copy of the finished session OUTSIDE the staged work dir, which searched `raws`."""
        other = os.path.join(self.tmp, "somewhere", "session")
        shutil.copytree(self.session, other, symlinks=True)
        write(other, "input/raw_files.txt", "# Raw MS files\n" + "\n".join(raws) + "\n")
        write(other, "logs/qc_bracket.json", json.dumps(tcs.qc_record("good", session=other)))
        return other

    def deliver_from(self, session, *extra):
        return run(["deliver", "--summary", self.summary, "--session", session, "--label", "t1",
                    "--include-raw", "no"] + list(extra), self.env)

    def test_a_subset_session_gets_a_plan_that_names_the_left_out_files(self):
        other = self.elsewhere(self.raw_files[:1])
        rc, out, p = self.deliver_from(other)
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        left = out["session_left_out"]
        self.assertEqual(left["n_staged"], 2)
        self.assertEqual(left["names"], [os.path.basename(self.raw_files[1])])
        self.assertTrue(any("1 staged file not in this analysis" in w
                            and os.path.basename(self.raw_files[1]) in w for w in out["warnings"]),
                        out["warnings"])
        self.assertIn("subset", out["session_ownership"])

    def test_the_subset_is_recorded_in_delivery_json_and_the_email_does_not_claim_all(self):
        other = self.elsewhere(self.raw_files[:1])
        rc, out, p = self.deliver_from(other, "--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        dj = json.loads(read(os.path.join(other, "delivery.json")))
        self.assertEqual(dj["session_left_out"]["names"], [os.path.basename(self.raw_files[1])])
        md = os.path.join(self.tmp, "EMAIL_DRAFT.md")
        run(["email-draft", "--summary", self.summary, "--delivery",
             os.path.join(other, "delivery.json"), "--out", md], self.env)
        text = read(md)
        gone = os.path.basename(self.raw_files[1])
        self.assertIn(gone, text.split("-->")[0])                       # the staff note
        # ...and the CLIENT is told, by name, in the body (C2 review bio H1)
        body = text.split("-->", 1)[1]
        want = (f"1 of the 2 submitted files were analysed. Not analysed, and not in this share: "
                f"{gone}. Ask us if you want them processed.")
        self.assertIn(want, body)
        self.assertEqual(body.count(want), 1)
        # the client README (.md and .html) says the same, in the same words
        for name in ("README.md", "README.html"):
            readme = read(os.path.join(dj["delivery_dir"], name))
            self.assertIn(gone, readme, name)
            self.assertIn("1 of the 2 submitted files were analysed", readme, name)

    def test_left_out_words(self):
        self.assertIsNone(cs.left_out_text(None))
        self.assertIsNone(cs.left_out_text({"names": [], "n_staged": 3}))
        lo = cs.left_out_text({"names": ["a.d", "b.d"], "n_staged": 96})
        self.assertTrue(lo["sentence"].startswith("94 of the 96 submitted files were analysed."))
        self.assertIn("not in this share: a.d, b.d.", lo["sentence"])
        self.assertIn("94 analysed files only", lo["raw"])
        many = cs.left_out_text({"names": [f"r{i}.d" for i in range(25)], "n_staged": 96})
        self.assertIn("r19.d, and 5 more.", many["sentence"])
        self.assertNotIn("r20.d", many["sentence"])
        allraw = cs.left_out_text({"names": ["a.d"], "n_staged": 2, "raw_all_linked": True})
        self.assertIn("their raw files are in the raw folder", allraw["sentence"])
        self.assertNotIn("not in this share", allraw["sentence"])
        self.assertIn("not analysed included", allraw["raw"])

    def test_the_raw_folder_line_does_not_claim_every_file(self):
        other = self.elsewhere(self.raw_files[:1])
        rc, out, p = self.deliver_from(other, "--include-raw", "yes", "--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        dj = json.loads(read(os.path.join(other, "delivery.json")))
        self.assertFalse(dj["session_left_out"]["raw_all_linked"])
        readme = read(os.path.join(dj["delivery_dir"], "README.md"))
        self.assertIn("the raw instrument files of the 1 analysed file only", readme)
        self.assertNotIn("your raw instrument files", readme)
        md = os.path.join(self.tmp, "EMAIL_DRAFT2.md")
        run(["email-draft", "--summary", self.summary, "--delivery",
             os.path.join(other, "delivery.json"), "--out", md], self.env)
        body = read(md).split("-->", 1)[1]
        self.assertNotIn("Your raw instrument files are in the raw folder", body)
        self.assertIn("The raw folder of the same share holds the raw instrument files of the 1 "
                      "analysed file only.", body)

    def test_a_session_that_searched_exactly_the_staged_files_has_no_left_out_line(self):
        other = self.elsewhere(self.raw_files)
        rc, out, p = self.deliver_from(other)
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertIsNone(out["session_left_out"])
        self.assertFalse(any("not in this analysis" in w for w in out["warnings"]))

    def test_a_session_with_a_file_the_stage_record_lacks_is_refused_with_the_sanctioned_path(self):
        for raws in (self.raw_files[:1] + ["/nfs/other/lab_run.d"],             # foreign + subset
                     self.raw_files + ["/nfs/other/lab_run.d"]):                # superset
            other = self.elsewhere(raws)
            rc, out, _ = self.deliver_from(other)
            self.assertEqual(rc, 2)
            self.assertIn("does not belong", out["error"])
            self.assertIn("lab_run.d", out["error"])
            self.assertIn("stage --apply", out["hint"])
            shutil.rmtree(os.path.dirname(other))

    def test_force_does_not_bypass_membership_when_a_stage_record_exists(self):
        other = self.elsewhere(["/nfs/other/lab_run.d"])
        rc, out, _ = self.deliver_from(other, "--force")
        self.assertEqual(rc, 2)
        self.assertIn("--force", out["hint"])
        self.assertIn("only", out["hint"])

    def test_include_raw_yes_links_the_analysed_files_only_and_all_links_every_staged_file(self):
        other = self.elsewhere(self.raw_files[:1])
        raw_dir = os.path.join(self.share, "raw")
        links = lambda: sorted(os.listdir(raw_dir)) if os.path.isdir(raw_dir) else []   # noqa: E731
        rc, out, p = run(["deliver", "--summary", self.summary, "--session", other, "--label", "t1",
                          "--include-raw", "yes"], self.env)                      # dry run
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertEqual([x["name"] for x in out["raw"]["links"]],
                         [os.path.basename(self.raw_files[0])])
        self.assertTrue(any("--include-raw all" in w for w in out["warnings"]))
        rc, out, p = self.deliver_from(other, "--include-raw", "yes", "--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertEqual(links(), [os.path.basename(self.raw_files[0])])
        rc, out, p = run(["deliver", "--summary", self.summary, "--session", other, "--label", "t2",
                          "--include-raw", "all"], self.env)
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertEqual(sorted(x["name"] for x in out["raw"]["links"]),
                         sorted(os.path.basename(f) for f in self.raw_files))
        self.assertFalse(any("--include-raw all links" in w for w in out["warnings"]))

    def test_help_says_what_force_and_allow_partial_cover(self):
        src = read(os.path.join(SCRIPTS, "core_submission.py"))
        force = src.split('de.add_argument("--force"', 1)[1].split("de.add_argument", 1)[0]
        self.assertIn("does NOT bypass", force)
        self.assertIn("stage record", force)
        partial = src.split('lo.add_argument("--allow-partial"', 1)[1].split("lo.add_argument", 1)[0]
        self.assertIn("not part of this analysis", partial)


if __name__ == "__main__":
    unittest.main()
