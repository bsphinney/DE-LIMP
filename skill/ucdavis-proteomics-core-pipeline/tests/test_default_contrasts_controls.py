#!/usr/bin/env python3
"""With the control groups recorded, the default contrasts are group minus control (2.14.2, I18).

Before: run_de.R's default contrast is every group against the FIRST level, alphabetically
(`paste0(l[-1], "-", l[1])`), so a TurboID study with groups Bait and IgG got IgG-Bait -- control
minus bait, every fold change upside down -- even with `experiment_type.py set --controls IgG`
recorded. normalization_check.py passed --contrasts only when one was given, and the pull-down
note then said "no --controls were recorded" when they had been. Now:

  * experiment_type.default_contrasts(groups, rec): `<group>-<control>` for every recorded
    control and every other group; None when no controls were recorded (group names alone stay
    advisory);
  * normalization_check.step8_commands() (and so required()) adds those contrasts to all three
    printed DEs when the asking design has no --contrasts of its own; given contrasts win;
  * pulldown_design()'s note tells "controls recorded, but on the left" from "none recorded".

Stdlib only; hermetic.
"""
import csv
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import experiment_type as et  # noqa: E402
import normalization_check as nc  # noqa: E402


def session_with(d, groups, controls=None):
    """A session recording an experiment type (proximity labelling) and `controls`, with a
    conditions.csv of two runs per group -> (session, conditions, report)."""
    s = os.path.join(d, "session")
    os.makedirs(os.path.join(s, "input"))
    et.record(s, "proximity", "user", controls=",".join(controls) if controls else None)
    cond = os.path.join(s, "input", "conditions.csv")
    with open(cond, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["File.Name", "Group"])
        for g in groups:
            w.writerows([(f"{g}_1", g), (f"{g}_2", g)])
    return s, cond, os.path.join(s, "output", "search", "report.parquet")


class DefaultContrasts(unittest.TestCase):
    def test_recorded_control_puts_bait_on_the_left(self):
        self.assertEqual(et.default_contrasts(["Bait", "IgG"], {"controls": ["IgG"]}), ["Bait-IgG"])

    def test_every_group_against_every_recorded_control(self):
        self.assertEqual(et.default_contrasts(["Bait", "Beads", "Fusion", "IgG"],
                                              {"controls": ["IgG", "Beads"]}),
                         ["Bait-Beads", "Fusion-Beads", "Bait-IgG", "Fusion-IgG"])

    def test_no_recorded_controls_gives_none(self):
        self.assertIsNone(et.default_contrasts(["A", "B"], {"controls": []}))
        self.assertIsNone(et.default_contrasts(["A", "B"], {"controls": None}))
        self.assertIsNone(et.default_contrasts(["A", "B"], None))
        # a group NAMED like a control is still advisory only: nothing recorded, nothing chosen
        self.assertIsNone(et.default_contrasts(["Bait", "IgG"], {"controls": None}))
        # recorded controls that are not among the groups choose nothing either
        self.assertIsNone(et.default_contrasts(["A", "B"], {"controls": ["IgG"]}))


class Step8Commands(unittest.TestCase):
    def test_step8_commands_carry_the_default_contrasts(self):
        with tempfile.TemporaryDirectory() as d:
            s, cond, report = session_with(d, ["Bait", "IgG"], controls=["IgG"])
            cmds = nc.step8_commands(report, cond, s)
        self.assertTrue(all("--contrasts Bait-IgG" in c for c in cmds["check"][:2]), cmds["check"])
        self.assertIn("--contrasts Bait-IgG", cmds["final"])
        self.assertNotIn("--contrasts", cmds["check"][2])      # the check itself takes none

    def test_required_prints_them(self):
        with tempfile.TemporaryDirectory() as d:
            s, cond, report = session_with(d, ["Bait", "IgG"], controls=["IgG"])
            r = nc.required(report, cond)
        self.assertTrue(r["required"])
        self.assertEqual(r["message"].count("--contrasts Bait-IgG"), 3)

    def test_given_contrasts_win(self):
        """(pin) given --contrasts were always carried as given."""
        with tempfile.TemporaryDirectory() as d:
            s, cond, report = session_with(d, ["Bait", "IgG"], controls=["IgG"])
            cmds = nc.step8_commands(report, cond, s, design={
                "contrasts": ["IgG-Bait"], "flags": ["--contrasts", "IgG-Bait"]})
        for c in cmds["check"][:2] + [cmds["final"]]:
            self.assertIn("--contrasts IgG-Bait", c)
            self.assertNotIn("Bait-IgG", c)

    def test_run_des_own_default_is_not_a_given_contrast(self):
        """A design_request records run_de.R's alphabetical default in `contrasts` with no
        --contrasts among its flags: that is the default being replaced, not a choice."""
        with tempfile.TemporaryDirectory() as d:
            s, cond, report = session_with(d, ["Bait", "IgG"], controls=["IgG"])
            cmds = nc.step8_commands(report, cond, s, design={"contrasts": ["IgG-Bait"],
                                                              "flags": []})
        self.assertIn("--contrasts Bait-IgG", cmds["final"])

    def test_no_recorded_controls_adds_nothing(self):
        """(pin) without recorded controls the commands carry no --contrasts, as before."""
        with tempfile.TemporaryDirectory() as d:
            s, cond, report = session_with(d, ["Bait", "IgG"])
            cmds = nc.step8_commands(report, cond, s)
        self.assertNotIn("--contrasts", cmds["next"])


class TheNote(unittest.TestCase):
    def test_reversed_contrast_note_differs_from_none_recorded(self):
        d = et.pulldown_design({"type": "ip", "label": "IP", "controls": ["IgG"]}, ["Bait", "IgG"],
                               ["IgG-Bait"])
        self.assertIn("control on the left", d["note"])
        self.assertNotIn("no --controls were recorded", d["note"])
        self.assertEqual(d["note"], (
            "the control groups were recorded (IgG), but every contrast has a control on the left "
            "(IgG-Bait), so its fold changes are control minus bait: use --contrasts Bait-IgG"))

    def test_none_recorded_keeps_its_note(self):
        """(pin) no control group at all: the 2.14.1 note."""
        d = et.pulldown_design({"type": "ip", "label": "IP", "source": "the user"}, ["A", "B"],
                               ["B-A"])
        self.assertIn("no --controls were recorded", d["note"])

    def test_bait_minus_control_has_no_note(self):
        """(pin)"""
        d = et.pulldown_design({"type": "ip", "label": "IP", "source": "the user",
                                "controls": ["IgG"]}, ["Bait", "IgG"], ["Bait-IgG"])
        self.assertIsNone(d["note"])
        self.assertEqual(d["vs_control"], ["Bait-IgG"])


class NoContrastsIsRefused(unittest.TestCase):
    """C2 statistician M1: "fold changes point the right way" held only for step 8's printed
    commands. A DE typed without --contrasts (from a template, both check DEs too) in a session
    with recorded controls fell back to run_de.R's alphabetical default -- IgG-Bait, control
    minus bait -- and nothing warned. Now it is refused, with the exact --contrasts
    (normalization_check.contrasts_refusal, asked by run_de.R)."""

    def test_the_rule_names_the_exact_contrasts(self):
        with tempfile.TemporaryDirectory() as d:
            s, cond, report = session_with(d, ["Bait", "IgG"], controls=["IgG"])
            r = nc.contrasts_refusal(report, cond)
        self.assertTrue(r["refused"])
        self.assertEqual(r["contrasts"], "Bait-IgG")
        self.assertIn("REFUSED: this session records the control group(s) IgG", r["message"])
        self.assertIn("--contrasts Bait-IgG", r["message"])

    def test_several_groups_and_controls(self):
        with tempfile.TemporaryDirectory() as d:
            s, cond, report = session_with(d, ["Bait", "Beads", "IgG"],
                                           controls=["IgG", "Beads"])
            r = nc.contrasts_refusal(report, cond)
        self.assertEqual(r["contrasts"], ",".join(et.default_contrasts(
            ["Bait", "Beads", "IgG"], {"controls": ["IgG", "Beads"]})))

    def test_no_recorded_controls_or_no_session_keeps_the_default(self):
        """(pin) group names alone stay advisory; a DE outside a session runs as before"""
        with tempfile.TemporaryDirectory() as d:
            s, cond, report = session_with(d, ["Bait", "IgG"])
            self.assertFalse(nc.contrasts_refusal(report, cond)["refused"])
        with tempfile.TemporaryDirectory() as d:
            cond = os.path.join(d, "conditions.csv")
            with open(cond, "w") as fh:
                fh.write("File.Name,Group\nA_1,A\nA_2,A\nB_1,B\nB_2,B\n")
            self.assertFalse(nc.contrasts_refusal(os.path.join(d, "report.parquet"),
                                                  cond)["refused"])

    def test_the_cli_answers_run_de(self):
        import json
        import subprocess
        with tempfile.TemporaryDirectory() as d:
            s, cond, report = session_with(d, ["Bait", "IgG"], controls=["IgG"])
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "normalization_check.py"),
                                "contrasts", "--input", report, "--metadata", cond],
                               capture_output=True, text=True, timeout=60)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(json.loads(p.stdout)["contrasts"], "Bait-IgG")

    def test_run_de_refuses_without_contrasts(self):
        import shutil
        import subprocess
        from test_run_de_contaminants import r_has
        if not (shutil.which("Rscript") and r_has("limma", "jsonlite")):
            self.skipTest("needs R with limma and jsonlite")
        with tempfile.TemporaryDirectory() as d:
            s, cond, report = session_with(d, ["Bait", "IgG"], controls=["IgG"])
            p = subprocess.run(["Rscript", os.path.join(SCRIPTS, "run_de.R"), "--input", report,
                                "--metadata", cond, "--quantities", "normalised",
                                "--check-input", "--outdir", os.path.join(d, "de")],
                               capture_output=True, text=True, timeout=600)
        self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIn("REFUSED: this session records the control group(s) IgG", p.stderr)
        self.assertIn("--contrasts Bait-IgG", p.stderr)


if __name__ == "__main__":
    unittest.main()
