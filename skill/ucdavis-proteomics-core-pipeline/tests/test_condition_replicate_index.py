#!/usr/bin/env python3
"""Names with a number glued to them: a possible difference for staff, never a finding in the
report (Core review of the client IP-MS guide, 2026-10-20).

Before: submission_report's sheet check read ATG5 as ATG replicate 5 and HSP90 as HSP replicate
90, so a sheet with ATG5 x3, ATG7 x3 and WT x3 became "ATG" and "WT", and the report told the
client the analysed groups split the sheet's "ATG" ("conditions differ"), which was false. Three
reviews showed that no rule reading names alone can tell ATG5 / ATG7 (genes) from Parental12 /
Parental51 (sample IDs), so the sheet's groups stay exactly 2.13.2's (test_replicate_parity).
What changed is the note: a part of "Conditions analysed vs the sample sheet" that rests on a
group made by taking a number glued to the names off them is withheld from the report and said
to staff instead, as a possibility (submission_report.sheet_cautions), on both readings: what
differs if the names are replicates, and what differs if they are separate conditions. Staff
see it in the commands' JSON and the decisions log -- never in a file a collaborator gets: not
the report, not the brief, not the session zip. Every other note is 2.13.2's.

Stdlib only.
"""
import csv
import json
import os
import stat
import subprocess
import sys
import tempfile
import unittest
import zipfile
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import submission_report as sr  # noqa: E402
from test_submission_report import fixture, notes_by_id  # noqa: E402
from job_env import job_env  # noqa: E402


def twice(names):
    return [x for x in names for _ in (0, 1)]


def record(names):
    r = fixture()
    r["submission_data"]["samples"] = [
        {"unique_id": f"SMP{100 + i}", "sample_name": f"sample {i}", "condition_name": c}
        for i, c in enumerate(names, 1)]
    return sr.sanitize(r)


def session(tmp, names, analysed, s=None, runs=None):
    """A session whose conditions.csv holds `analysed` (by hand), each run assigned to its
    sample, with the sheet `names` attached. `s`, `runs`: an existing session and its runs."""
    s = s or os.path.join(tmp, "2026-10-20_conditions")
    os.makedirs(os.path.join(s, "input"), exist_ok=True)
    os.makedirs(os.path.join(s, "output"), exist_ok=True)
    rec = record(names)
    runs = dict(zip([x["unique_id"] for x in rec["samples"]], runs or [
        f"10202026_DIA_{x['unique_id']}_1_{i}" for i, x in enumerate(rec["samples"], 1)]))
    if not os.path.isfile(os.path.join(s, "input", "raw_files.txt")):
        with open(os.path.join(s, "input", "raw_files.txt"), "w") as fh:
            fh.write("# raw\n" + "".join(f"/data/OCT26/{r}.d\n" for r in runs.values()))
    with open(os.path.join(s, "input", "conditions.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["File.Name", "Group"])
        for x, g in zip(rec["samples"], analysed):
            w.writerow([runs[x["unique_id"]], g])
    with open(os.path.join(s, "input", "locate_decisions.json"), "w") as fh:
        json.dump({"assignments": [{"sample": u, "file": f"{r}.d"} for u, r in runs.items()]}, fh)
    sr.attach(s, rec)
    return s


def run(script, *args, env=None):
    return subprocess.run([sys.executable, os.path.join(SCRIPTS, script), *args],
                          capture_output=True, text=True, timeout=300, env=env)


def caution(names, analysed):
    """The staff caution's text after its lead, for a session with this sheet and analysis."""
    with tempfile.TemporaryDirectory() as tmp:
        s = session(tmp, names, analysed)
        cautions = sr.sheet_cautions(sr.load(s), s)
    return cautions[0]["text"][len(sr.CAUTION_LEAD):] if cautions else None


# (the sheet, the analysis's groups by hand: the names are separate conditions)
GLUED = (
    (["ATG5"] * 3 + ["ATG7"] * 3 + ["WT"] * 3, ["ATG5"] * 3 + ["ATG7"] * 3 + ["WT"] * 3, "ATG"),
    (twice(["WT3", "WT4", "WT6"]) + ["KO"] * 6, twice(["WT3", "WT4", "WT6"]) + ["KO"] * 6, "WT"),
    (twice(["Parental12", "Parental51"]) + ["Dox(+)"] * 4,
     twice(["Parental12", "Parental51"]) + ["Dox(+)"] * 4, "Parental"),
)


class GluedNumbersAreForStaff(unittest.TestCase):
    def test_the_groups_are_2_13_2s(self):
        for names, _analysed, base in GLUED:
            self.assertIn(base, set(sr.sheet_groups(record(names)).values()), names)

    def test_the_client_note_is_withheld_and_staff_are_told(self):
        for names, analysed, base in GLUED:
            with tempfile.TemporaryDirectory() as tmp:
                s = session(tmp, names, analysed)
                rec = sr.load(s)
                self.assertNotIn("conditions_differ", notes_by_id(rec, s), names)
                cautions = sr.sheet_cautions(rec, s)
                self.assertEqual([c["id"] for c in cautions], ["conditions_glued_names"], names)
                text = cautions[0]["text"]
                self.assertIn(f"as “{base}”), the sheet's “{base}” is split across groups",
                              text)
                self.assertIn("read as separate conditions, nothing differs for these samples",
                              text)
                self.assertIn("staff only, a possibility, never in the report", text)
                # not in the report's Submission section either
                self.assertNotIn("split across", sr.report_section(session=s)["md"])

    def test_a_note_not_resting_on_glued_numbers_is_2_13_2s(self):
        """Marked replicate numbers (WT_1): the client note stays, and no caution."""
        names = [f"{g}_{n}" for g in ("WT", "KO") for n in (1, 2, 3)]
        with tempfile.TemporaryDirectory() as tmp:
            s = session(tmp, names, ["WT"] * 6)
            rec = sr.load(s)
            self.assertIn("group “WT” pools the sheet's KO, WT",
                          notes_by_id(rec, s)["conditions_differ"])
            self.assertEqual(sr.sheet_cautions(rec, s), [])

    def test_glued_and_other_parts_are_split_between_client_and_staff(self):
        names = ["ATG5"] * 3 + ["ATG7"] * 3 + ["WT"] * 3 + ["KO"] * 3
        analysed = ["ATG5"] * 3 + ["ATG7"] * 3 + ["WT"] * 6        # WT pools the sheet's KO
        with tempfile.TemporaryDirectory() as tmp:
            s = session(tmp, names, analysed)
            rec = sr.load(s)
            note = notes_by_id(rec, s)["conditions_differ"]
            self.assertIn("group “WT” pools the sheet's KO, WT", note)
            self.assertNotIn("ATG", note)
            self.assertIn("“ATG” is split across", sr.sheet_cautions(rec, s)[0]["text"])

    def test_the_caution_says_what_differs_on_the_other_reading(self):
        """From the data, never from the kind of part: "nothing differs" only when nothing does,
        and a pool is a pool on either reading."""
        # IP and input rows: WT3 is split across WT_IP and WT_Input on either reading
        said = caution(twice(["WT3", "WT4", "WT6"]) + ["KO"] * 6,
                       ["WT_IP", "WT_Input"] * 3 + ["KO_IP", "KO_Input"] * 3)
        self.assertNotIn("nothing differs", said)
        self.assertIn("the sheet's “WT3” is split across groups WT_IP, WT_Input", said)
        # one group for all: it pools conditions either way
        said = caution(["ATG5"] * 3 + ["ATG7"] * 3 + ["WT"] * 3, ["All"] * 9)
        self.assertIn("group “All” pools the sheet's ATG, WT; read as separate conditions, "
                      "group “All” pools the sheet's ATG5, ATG7, WT: this pooling holds on "
                      "either reading.", said)
        # a split as replicates is a pool as separate conditions
        said = caution(["ATG5"] * 2 + ["ATG7"] * 2 + ["ATG9"] * 2 + ["WT"] * 2,
                       ["A"] * 4 + ["B"] * 2 + ["WT"] * 2)
        self.assertIn("the sheet's “ATG” is split across groups A, B; read as separate "
                      "conditions, group “A” pools the sheet's ATG5, ATG7.", said)
        # a pool of two glued groups: one sentence naming both
        said = caution(["ATG5"] * 2 + ["ATG7"] * 2 + ["HSP70"] * 2 + ["HSP90"] * 2, ["All"] * 8)
        self.assertEqual(said.count("Read as replicates"), 1, said)
        self.assertIn("(ATG5, ATG7 as “ATG”; HSP70, HSP90 as “HSP”)", said)

    def test_the_caution_is_logged_once_and_never_in_the_brief(self):
        names, analysed, _base = GLUED[0]
        with tempfile.TemporaryDirectory() as tmp:
            s = session(tmp, names, analysed)
            for logged in ("yes", "already"):            # `notes` run twice: logged once
                p = run("submission_report.py", "notes", "--session", s)
                self.assertEqual(p.returncode, 0, p.stderr)
                out = json.loads(p.stdout)
                self.assertEqual(out["staff_cautions_logged"], logged)
            self.assertNotIn("conditions_differ", [n["id"] for n in out["notes"]])
            self.assertEqual([c["id"] for c in out["staff_cautions"]], ["conditions_glued_names"])
            brief = []
            import analysis_prompt
            staff = analysis_prompt.submission_brief(brief.append, s)
            brief = "\n".join(brief)
            self.assertNotIn("“ATG” is split across", brief)
            self.assertNotIn("glued", brief)
            self.assertEqual(staff["staff_cautions"], out["staff_cautions"])
            self.assertEqual(staff["staff_cautions_logged"], "already")
            with open(os.path.join(s, "logs", "decisions.md"), encoding="utf-8") as fh:
                self.assertEqual(fh.read().count("“ATG” is split across"), 1)

    def test_a_log_that_cannot_be_written_is_said_never_a_crash(self):
        """decisions.md is 0640: a colleague running `notes` on the session cannot write it."""
        names, analysed, _base = GLUED[0]
        with tempfile.TemporaryDirectory() as tmp:
            s = session(tmp, names, analysed)
            log = os.path.join(s, "logs", "decisions.md")
            os.makedirs(os.path.dirname(log), exist_ok=True)
            with open(log, "w", encoding="utf-8") as fh:
                fh.write("# Decisions log\n")
            os.chmod(log, stat.S_IRUSR | stat.S_IRGRP)
            try:
                if os.access(log, os.W_OK):
                    self.skipTest("this user can write a read-only file (root)")
                for cmd in (["notes", "--session", s], ["show", "--session", s]):
                    p = run("submission_report.py", *cmd)
                    self.assertEqual(p.returncode, 0, p.stderr)
                    out = json.loads(p.stdout)
                    self.assertEqual([c["id"] for c in out["staff_cautions"]],
                                     ["conditions_glued_names"])
                    self.assertTrue(out["staff_cautions_logged"].startswith("no: "), out)
                    self.assertIn("the staff caution was not logged in", p.stderr)
            finally:
                os.chmod(log, stat.S_IRUSR | stat.S_IWUSR)

    def test_no_file_in_the_session_zip_carries_the_caution(self):
        """notes and the brief, then finalize --zip: the caution is in logs/decisions.md and the
        commands' JSON, and in no file of the zip a collaborator gets."""
        import test_deposit_package as tdp
        names, analysed, _base = GLUED[0]
        runs = [f"HeLa_atg_{i:02d}" for i in range(1, len(names) + 1)]
        with tempfile.TemporaryDirectory() as tmp, mock.patch.object(tdp, "RUNS", runs):
            p = tdp.dia_session(tmp)
            s = session(tmp, names, analysed, s=p["session_dir"], runs=runs)
            env = job_env(tmp, PDF_BROWSER="/usr/bin/false")   # never a browser, never Slack
            r = run("submission_report.py", "notes", "--session", s, env=env)
            self.assertEqual(r.returncode, 0, r.stderr)
            text = json.loads(r.stdout)["staff_cautions"][0]["text"]
            # what this session's caution says (the skill's own scripts, in the zip, hold the
            # wording but not these)
            marks = ("ATG5, ATG7 as “ATG”", "the sheet's “ATG” is split across groups")
            self.assertTrue(all(m in text for m in marks), text)
            brief = os.path.join(p["output_dir"], "ANALYSIS_PROMPT.md")
            r = run("analysis_prompt.py", "--de-dir", p["de_dir"], "--out", brief,
                    "--submission", s, env=env)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual([c["id"] for c in json.loads(r.stdout)["staff_cautions"]],
                             ["conditions_glued_names"])
            r = run("session.py", "finalize", "--dir", s, "--zip", env=env)
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(os.path.join(s, "logs", "decisions.md"), encoding="utf-8") as fh:
                logged = fh.read()                              # the text exists, on disk
            self.assertTrue(all(m in logged for m in marks), logged)
            with zipfile.ZipFile(json.loads(r.stdout)["zip"]) as z:
                listed = z.namelist()
                self.assertTrue(any(n.endswith("output/ANALYSIS_PROMPT.md") for n in listed))
                self.assertFalse([n for n in listed if n.endswith("decisions.md")], listed)
                for n in listed:
                    data = z.read(n).decode("utf-8", "replace")
                    for m in marks:
                        self.assertNotIn(m, data, n)


if __name__ == "__main__":
    unittest.main()
