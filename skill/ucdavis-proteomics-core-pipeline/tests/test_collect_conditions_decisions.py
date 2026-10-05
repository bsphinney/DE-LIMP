#!/usr/bin/env python3
"""Which run IS which sample is the user's answer: collect_conditions --map never decides it.

A staff search (2026-09-30): the runs were numbered, two runs carried the same number and one
number had no run. The agent mapped by plate well "so the duplicate doesn't matter" and went
straight to the compute confirmation; the user was never asked, and a later re-analysis showed
sample identity really was in doubt. SKILL.md called a label that hits several files "usually
fine", and needs_confirmation did not even count it.

Guards:
  * a label that names several runs, or that two samples share, puts those runs under
    decisions_required with the question, and writes them with NO group -- so --validate
    refuses the CSV until the user answers;
  * --confirm-multi '<label>' is the answer "it names every one of those runs"; it is recorded
    in <csv>.decisions.json, which provenance.py copies with the CSV.
"""
import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
SCRIPT = os.path.join(SCRIPTS, "collect_conditions.py")

# Ten numbered samples; two runs carry 61 and none carries 62.
NUMBERED = [f"exp_{n}_plate" for n in (55, 56, 57, 58, 59, 60, 63, 64)] + \
    ["exp_61_plate_a", "exp_61_plate_b"]
PREFIXED = [f"{g}_{i}" for g in ("ctrl", "trt") for i in (1, 2, 3)]


def run_map(tmp, runs, *extra, mapping=None, sheet=None):
    out = os.path.join(tmp, "conditions.csv")
    args = [sys.executable, SCRIPT, "--map", out, "--runs", ",".join(runs), *extra]
    if sheet is not None:
        path = os.path.join(tmp, "sheet.csv")
        with open(path, "w", newline="") as fh:
            csv.writer(fh).writerows(sheet)
        args += ["--from-file", path]
    else:
        args += ["--mapping-json", json.dumps(mapping)]
    p = subprocess.run(args, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    with open(out, newline="") as fh:
        return json.loads(p.stdout), {r["File.Name"]: r["Group"] for r in csv.DictReader(fh)}, out


def validate(csv_path):
    return subprocess.run([sys.executable, SCRIPT, "--validate", csv_path],
                          capture_output=True, text=True)


class Decisions(unittest.TestCase):
    def test_a_number_naming_two_runs_is_asked_and_left_without_a_group(self):
        mapping = {"groups": {"A": ["55", "56", "57", "58", "59"],
                              "B": ["60", "61", "62", "63", "64"]}}
        with tempfile.TemporaryDirectory() as tmp:
            rep, groups, out = run_map(tmp, NUMBERED, mapping=mapping)
            self.assertTrue(rep["needs_confirmation"])
            amb = rep["ambiguities"]
            self.assertEqual(amb["multi_match_identifiers"], {"61": ["exp_61_plate_a", "exp_61_plate_b"]})
            self.assertEqual(amb["unmatched_identifiers"], ["62"])
            self.assertEqual(amb["awaiting_decision_runs"], ["exp_61_plate_a", "exp_61_plate_b"])
            self.assertEqual(amb["unassigned_runs"], [])
            (d,) = rep["decisions_required"]
            self.assertEqual((d["identifier"], d["kind"]), ("61", "label_names_several_runs"))
            self.assertIn("exp_61_plate_a, exp_61_plate_b", d["question"])
            self.assertIn("--confirm-multi '61'", d["question"])
            self.assertEqual((groups["exp_61_plate_a"], groups["exp_61_plate_b"]), ("", ""))
            self.assertEqual(groups["exp_60_plate"], "B")
            self.assertIn("never yours", rep["guidance"])
            self.assertNotEqual(validate(out).returncode, 0, "a CSV with an open question validates")

    def test_a_group_prefix_is_assigned_only_once_the_user_says_so(self):
        mapping = {"groups": {"control": ["ctrl"], "treated": ["trt"]}}
        with tempfile.TemporaryDirectory() as tmp:
            rep, groups, out = run_map(tmp, PREFIXED, mapping=mapping)
            self.assertEqual(set(groups.values()), {""})
            self.assertEqual({d["identifier"] for d in rep["decisions_required"]}, {"ctrl", "trt"})
            rep, groups, out = run_map(tmp, PREFIXED, "--confirm-multi", "ctrl",
                                       "--confirm-multi", "TRT", mapping=mapping)
            self.assertFalse(rep["needs_confirmation"], rep["ambiguities"])
            self.assertEqual(rep["decisions_required"], [])
            self.assertEqual(groups, {**{f"ctrl_{i}": "control" for i in (1, 2, 3)},
                                      **{f"trt_{i}": "treated" for i in (1, 2, 3)}})
            with open(out + ".decisions.json") as fh:
                rec = json.load(fh)
            self.assertEqual(rec["confirmed_multi_match"]["ctrl"], ["ctrl_1", "ctrl_2", "ctrl_3"])
            self.assertEqual(rep["decisions_file"], os.path.abspath(out + ".decisions.json"))
            # the same mapping without the answer: the old answer does not stay behind
            run_map(tmp, PREFIXED, mapping=mapping)
            self.assertFalse(os.path.exists(out + ".decisions.json"))

    def test_a_label_two_samples_share_is_asked_not_merged(self):
        sheet = [["sample", "group"], ["S1", "A"], ["S1", "A"], ["S2", "A"], ["S3", "B"], ["S4", "B"]]
        with tempfile.TemporaryDirectory() as tmp:
            rep, groups, _ = run_map(tmp, ["S1", "S2", "S3", "S4"], sheet=sheet)
            self.assertEqual(rep["ambiguities"]["duplicate_identifiers"], {"S1": ["A", "A"]})
            (d,) = rep["decisions_required"]
            self.assertEqual(d["kind"], "duplicate_identifier")
            self.assertEqual(groups["S1"], "")
            self.assertEqual(rep["ambiguities"]["conflicting_runs"], {})
            # a confirmation does not settle two samples sharing a label
            rep, groups, _ = run_map(tmp, ["S1", "S2", "S3", "S4"], "--confirm-multi", "S1", sheet=sheet)
            self.assertEqual(groups["S1"], "")

    def test_an_answer_for_a_label_that_names_one_run_is_reported(self):
        mapping = {"groups": {"control": ["ctrl_1", "ctrl_2"], "treated": ["trt_1", "trt_2"]}}
        with tempfile.TemporaryDirectory() as tmp:
            rep, _, _ = run_map(tmp, PREFIXED, "--confirm-multi", "ctrl_1", mapping=mapping)
            self.assertEqual(rep["confirm_multi_unused"], ["ctrl_1"])

    def test_provenance_keeps_the_answers_with_the_csv(self):
        mapping = {"groups": {"control": ["ctrl"], "treated": ["trt"]}}
        with tempfile.TemporaryDirectory() as tmp:
            _, _, out = run_map(tmp, PREFIXED, "--confirm-multi", "ctrl", "--confirm-multi", "trt",
                                mapping=mapping)
            wm = os.path.join(tmp, "workflow.manifest.json")
            with open(wm, "w") as fh:
                fh.write("{}")
            bundle = os.path.join(tmp, "repro")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"), "--outdir", bundle,
                                "--workflow-manifest", wm, "--conditions", out, "--engine", "diann",
                                "--acquisition", "DIA", "--instrument", "timsTOF HT"],
                               capture_output=True, text=True, cwd=tmp)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertTrue(os.path.exists(os.path.join(bundle, "inputs", "conditions.csv.decisions.json")))


def _well(i):
    return "ABCDEFGH"[i // 12] + str(i % 12 + 1)


class DelimitersAreTokenBoundaries(unittest.TestCase):
    """2.9.1 (a staff search, 2026-10-02): --mapping-json identifiers 'DIA-XY1_' and 'DIA-XY2_',
    meant for XY1 and XY2 only, matched every XY10..XY19 and XY20..XY24 run -- the trailing '_' was
    normalised away and what was left matched as a prefix: 13 conflicting runs and a broken
    conditions.csv. A delimiter in an identifier is a token boundary, by the one definition of a
    delimited token (core_submission.token_pattern, which `locate` matches sample ids with), and a
    run carrying the identifier as a whole token is preferred over runs that only extend it.
    Anything still ambiguous is asked, never guessed."""

    RUNS = [f"10012026__60SPD_DIA-XY{i + 1}_S3-{_well(i)}_1_{24001 + i}" for i in range(24)]

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, SCRIPTS)
        import collect_conditions
        cls.cc = collect_conditions

    def test_a_trailing_delimiter_ends_the_identifier(self):
        mapping = {"groups": {"ctrl": ["DIA-XY1_", "DIA-XY3_"], "trt": ["DIA-XY2_", "DIA-XY4_"]}}
        with tempfile.TemporaryDirectory() as tmp:
            res, groups, _ = run_map(tmp, self.RUNS, mapping=mapping)
        r = self.RUNS
        self.assertEqual(res["assigned"], {r[0]: "ctrl", r[2]: "ctrl", r[1]: "trt", r[3]: "trt"})
        self.assertEqual(res["ambiguities"]["conflicting_runs"], {})
        self.assertEqual(res["ambiguities"]["multi_match_identifiers"], {})
        self.assertEqual([groups[x] for x in r[9:]], [""] * 15)      # XY10..XY24: no group

    def test_the_one_definition_of_a_delimited_token(self):
        import core_submission
        self.assertIs(self.cc.token_pattern, core_submission.token_pattern)

    def test_a_separator_between_letters_and_digits_is_optional_as_in_locate(self):
        """token_pattern: LRS96 finds LRS-96 (measured on the Core share) -- and never LRS-960."""
        m = self.cc.match_to_runs
        self.assertEqual(m("XY-1", self.RUNS), [self.RUNS[0]])
        self.assertEqual(m("dia xy 2_", self.RUNS), [self.RUNS[1]])

    def test_a_whole_token_is_preferred_over_a_prefix(self):
        m = self.cc.match_to_runs
        runs = ["exp_Ctrl_1", "exp_Ctrl1", "exp_Ctrl2"]
        self.assertEqual(m("Ctrl", runs), ["exp_Ctrl_1"])
        self.assertEqual(m("Ctrl", runs[1:]), runs[1:])      # only a prefix: every one of them

    def test_a_prefix_naming_several_runs_is_still_asked(self):
        mapping = {"groups": {"ctrl": ["DIA-XY"], "trt": ["DIA-XY2_"]}}
        with tempfile.TemporaryDirectory() as tmp:
            res, groups, _ = run_map(tmp, self.RUNS, mapping=mapping)
        self.assertTrue(res["needs_confirmation"])
        self.assertEqual(len(res["ambiguities"]["multi_match_identifiers"]["DIA-XY"]), 24)
        self.assertIn("label_names_several_runs", [q["kind"] for q in res["decisions_required"]])
        self.assertEqual(groups[self.RUNS[0]], "")


if __name__ == "__main__":
    unittest.main()
