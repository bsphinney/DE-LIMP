#!/usr/bin/env python3
"""A pull-down's control groups are recognised by more of their usual names (Core review of the
client IP-MS guide, 2026-10-20).

Before: experiment_type.IP_CONTROL_NAME knew only IgG and beads, so groups named IgG1, Isotype,
EmptyBeads, Parental, Untagged, Mock, Ctrl or KO were tagged by hand (set --controls). Now the
names are read in tiers, the first that names any group giving the controls:

  * IgG with a subclass or species (IgG1, IgG2a, mIgG, RbIgG), beads (EmptyBeads, BeadsOnly),
    isotype: a control in any study -- and with no experiment type recorded, the names alone
    still decide that the design is a pull-down;
  * Parental, Untagged, Mock, Empty / EV, KO / Knockout, then Ctrl / Control: also ordinary
    groups in a whole-proteome study (KO vs WT), so read only when the RECORDED type is an
    enrichment, and as a guess (controls_from "enrichment_name"): the brief and the session
    README tag it as not confirmed, and the normalisation check stops (exit 3) until the
    controls are recorded -- `decide` cannot stand in for them.

A name is a whole token of the group name: gene symbols that merely contain the letters never
match. Stdlib only.
"""
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import experiment_type as et  # noqa: E402
import normalization_check as nc  # noqa: E402
import session_docs  # noqa: E402
import test_deposit_package as tdp  # noqa: E402
from test_normalization_check import session  # noqa: E402
from test_normalization_check_hold import BALANCED, fake_de  # noqa: E402

ANY_STUDY = ("IgG", "Old_IgG", "IgG1", "IgG2a", "IgG2b", "IgG3", "IgG4", "mIgG", "RbIgG",
             "Beads", "Bead-only", "BeadsOnly", "EmptyBeads", "Empty_beads", "Empty-beads",
             "Isotype", "Isotype_IgG")
ENRICHMENT_ONLY = ("Parental", "Untagged", "Mock", "Empty", "EmptyVector", "EV", "KO", "Bait_KO",
                   "Knockout", "Knock-out", "Knock_out", "Ctrl", "Control", "Neg_Ctrl", "Controls")
NOT_A_CONTROL = ("SKO1", "SKOR1", "MOCS1", "CTRB1", "Controller", "Mockup", "Emptor", "IgGFc",
                 "Every", "WT", "Bait", "Treated")


def ip(**kw):
    return dict({"type": "ip", "label": et.DEFAULTS["ip"][0], "source": "the user"}, **kw)


def whole_proteome():
    return {"type": "whole_proteome", "label": et.DEFAULTS["whole_proteome"][0],
            "source": "the user"}


class TheNames(unittest.TestCase):
    def test_any_study_names_make_a_pulldown_with_no_type_recorded(self):
        for g in ANY_STUDY:
            pd = et.pulldown_design(None, ["Bait", g], [f"Bait-{g}"])
            self.assertEqual((pd["pulldown"], pd["controls"], pd["controls_from"]),
                             (True, [g], "name"), g)

    def test_every_control_name_is_recognised_in_a_recorded_enrichment(self):
        for rec in (ip(), dict(ip(), type="proximity")):
            for g in ANY_STUDY + ENRICHMENT_ONLY:
                pd = et.pulldown_design(rec, ["Bait", g], [f"Bait-{g}"])
                self.assertEqual((pd["controls"], pd["vs_control"]), ([g], [f"Bait-{g}"]),
                                 (rec["type"], g))
                self.assertEqual(pd["controls_from"],
                                 "name" if g in ANY_STUDY else "enrichment_name", g)

    def test_letters_inside_a_word_are_not_a_control(self):
        for rec in (None, ip()):
            for g in NOT_A_CONTROL:
                self.assertEqual(et.controls_in(["Bait", g], rec), [], g)


class NotEveryKOIsAnIP(unittest.TestCase):
    DESIGNS = ((["WT", "KO"], ["KO-WT", "WT-KO"]), (["Infected", "Mock"], ["Infected-Mock"]),
               (["Engineered", "Parental"], ["Engineered-Parental"]),
               (["Treated", "Control"], ["Treated-Control"]), (["OE", "Empty"], ["OE-Empty"]),
               (["OE", "EV"], ["OE-EV"]), (["Drug", "Ctrl"], ["Drug-Ctrl"]))

    def test_a_whole_proteome_design_is_not_a_pulldown_by_its_names(self):
        for groups, contrasts in self.DESIGNS:
            pd = et.pulldown_design(None, groups, contrasts)
            self.assertEqual((pd["pulldown"], pd["controls"]), (False, []), groups)

    def test_a_recorded_whole_proteome_raises_nothing(self):
        for groups, contrasts in self.DESIGNS:
            pd = et.pulldown_design(whole_proteome(), groups, contrasts)
            self.assertEqual((pd["pulldown"], pd["controls"], pd["note"]), (False, [], None), groups)
            # session_docs' README reads the same rule
            self.assertEqual(et.controls_in(groups, whole_proteome()), [], groups)


class WhichControlsInAnIP(unittest.TestCase):
    def test_igg_beside_a_treatment_is_the_only_control(self):
        """Bait_Ctrl is the untreated bait IP, not a negative control: IgG is."""
        groups = ["Bait_Ctrl", "Bait_Drug", "IgG"]
        contrasts = ["Bait_Ctrl-IgG", "Bait_Drug-IgG", "Bait_Drug-Bait_Ctrl"]
        pd = et.pulldown_design(ip(), groups, contrasts)
        self.assertEqual(pd["controls"], ["IgG"])
        self.assertEqual(pd["vs_control"], ["Bait_Ctrl-IgG", "Bait_Drug-IgG"])
        self.assertIsNone(pd["note"])

    def test_a_specific_name_beats_a_generic_ctrl(self):
        groups = ["Bait_Ctrl", "Bait_Drug", "Parental"]
        pd = et.pulldown_design(ip(), groups, ["Bait_Ctrl-Parental", "Bait_Drug-Bait_Ctrl"])
        self.assertEqual((pd["controls"], pd["vs_control"]), (["Parental"], ["Bait_Ctrl-Parental"]))

    def test_an_enrichment_only_name_is_used_and_raised_to_confirm(self):
        pd = et.pulldown_design(ip(), ["Bait", "Parental"], ["Bait-Parental"])
        self.assertTrue(pd["pulldown"])
        self.assertEqual(pd["controls_from"], "enrichment_name")
        self.assertIn("confirm with the user", pd["note"])
        self.assertIn("Parental", pd["note"])

    def test_the_recorded_controls_still_win(self):
        with tempfile.TemporaryDirectory() as d:
            s = session(d, "ip", controls="Ctrl")
            pd = et.pulldown_design(et.load(s), ["Bait", "Ctrl", "IgG"], ["Bait-Ctrl", "Bait-IgG"])
            self.assertEqual((pd["controls"], pd["controls_from"], pd["vs_control"], pd["note"]),
                             (["Ctrl"], "recorded", ["Bait-Ctrl"], None))

    def test_no_control_named_says_what_names_count(self):
        note = et.pulldown_design(ip(), ["A", "B"], ["A-B"])["note"]
        self.assertIn("record the control groups", note)
        self.assertIn("parental", note)

    def test_recorded_controls_that_name_no_group_are_said_so(self):
        note = et.pulldown_design(ip(controls=["IgG"]), ["Bait", "Mock"], ["Bait-Mock"])["note"]
        self.assertIn("the recorded control groups (IgG) are none of this design's groups", note)
        self.assertNotIn("no --controls were recorded", note)

    def test_a_control_that_is_no_contrasts_reference_is_said_so(self):
        """KO named, but the contrast is KO-Bait: not "none is named like one"."""
        note = et.pulldown_design(ip(), ["Bait", "KO"], ["KO-Bait"])["note"]
        self.assertIn("KO read as the control group(s), but no contrast has one as its reference",
                      note)
        self.assertNotIn("none is named", note)


class TheGuessIsHeldForTheUser(unittest.TestCase):
    """Which contrasts are enrichments decides which volcano anomalies the normalisation check
    expects: controls only guessed from their names are confirmed before anything is decided."""

    def inputs(self, control, controls=None):
        """A session recorded as an IP, a header-only report (B1-B3 not assessed: no pyarrow
        needed), conditions.csv and the two check DEs -> the `run` arguments."""
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        s = session(d, "ip", **({"controls": controls} if controls else {}))
        report = os.path.join(d, "report.tsv")
        with open(report, "w") as fh:
            fh.write("Run\tPrecursor.Id\tProtein.Group\tPrecursor.Quantity\t"
                     "Precursor.Normalised\tQ.Value\n")
        cond = os.path.join(d, "conditions.csv")
        with open(cond, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["File.Name", "Group"])
            w.writerows([f"{g}_{i}", g] for g in ("Bait", control) for i in (1, 2, 3))
        c = f"Bait-{control}"
        return (report, cond, fake_de(os.path.join(d, "n"), c, BALANCED),
                fake_de(os.path.join(d, "r"), c, BALANCED), s)

    def check(self, control, controls=None):
        report, cond, dn, dr, s = self.inputs(control, controls)
        return nc.run(report, cond, dn, dr, session=s)

    def test_a_guessed_control_needs_a_decision(self):
        rec, out = self.check("Parental")
        self.assertFalse(rec["tripped"])
        self.assertTrue(rec["needs_decision"])
        self.assertIsNone(rec["decision"])
        self.assertIn("read from their names only", rec["needs_decision_why"])
        with open(os.path.join(out, nc.SUMMARY)) as fh:
            summary = fh.read()
        self.assertIn("**Ask the user:** the control group(s) Parental", summary)
        self.assertNotIn("**Design:**", summary)               # said once, not twice

    def test_decide_cannot_stand_in_for_the_controls(self):
        """A decision would be reported as the type's default with the data check agreeing."""
        _rec, out = self.check("Parental")
        with self.assertRaisesRegex(ValueError, "nothing to decide yet.*set --controls"):
            nc.decide(os.path.join(out, nc.RECORD), et.NORMALISED, "someone", "a guess")
        self.assertFalse(nc.gate(os.path.join(out, nc.RECORD), et.RAW)[0])

    def test_the_cli_stops_and_says_why(self):
        report, cond, dn, dr, s = self.inputs("Parental")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "normalization_check.py"), "run",
                            "--report", report, "--conditions", cond, "--de-normalised", dn,
                            "--de-raw", dr, "--session", s],
                           capture_output=True, text=True, timeout=120)
        self.assertEqual(p.returncode, nc.TRIPPED_EXIT, p.stdout + p.stderr)
        self.assertIn("CONTROLS NOT CONFIRMED", p.stderr)
        self.assertTrue(json.loads(p.stdout)["needs_decision"])

    def test_recording_the_controls_clears_it(self):
        """The way out: record the confirmed controls, run the check again -- the type's default."""
        report, cond, dn, dr, s = self.inputs("Parental")
        self.assertTrue(nc.run(report, cond, dn, dr, session=s)[0]["needs_decision"])
        et.record(s, "ip", "user", controls="Parental")
        rec = nc.run(report, cond, dn, dr, session=s)[0]
        self.assertEqual((rec["needs_decision"], rec["controls_unconfirmed"]), (False, False))
        self.assertEqual(rec["decision"]["source"], "default")
        self.assertEqual(rec["design"]["controls_from"], "recorded")

    def test_igg_or_recorded_controls_are_decided_as_before(self):
        for control, controls in (("IgG", None), ("Parental", "Parental")):
            rec = self.check(control, controls)[0]
            self.assertFalse(rec["needs_decision"], control)
            self.assertEqual(rec["decision"]["source"], "default", control)


class TheGuessIsSaidWhereItIsRead(unittest.TestCase):
    """The brief and the session README name where the controls came from (rule 2)."""

    def session(self, controls=None, other="Parental"):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        p = tdp.dia_session(d)
        et.record(p["session_dir"], "ip", "user", controls=controls)
        c = f"Bait-{other}"
        prov = dict(tdp.DE_PROV, groups={"Bait": 3, other: 3}, contrasts=[c],
                    significant_per_contrast={c: 7}, n_samples=6)
        with open(os.path.join(p["de_dir"], "de_provenance.json"), "w") as fh:
            json.dump(prov, fh)
        return d, p

    def brief(self, d, p):
        out = os.path.join(d, "brief.md")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "analysis_prompt.py"), "--out",
                            out, "--de-dir", p["de_dir"], "--session", p["session_dir"]],
                           capture_output=True, text=True, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(out) as fh:
            return fh.read()

    def test_guessed_controls_are_tagged(self):
        d, p = self.session()
        said = "read from their names — NOT confirmed"
        self.assertIn(f"Control group(s), {said}", self.brief(d, p))
        readme = session_docs.agents_md(session_docs.gather(p["session_dir"]))
        self.assertIn(f"Pull-down controls, {said}", readme)

    def test_no_control_group_prints_no_empty_control_line(self):
        d, p = self.session(other="Lysate")
        brief = self.brief(d, p)
        self.assertIn("## Pull-down design", brief)
        self.assertNotIn("Control group(s)", brief)
        self.assertNotIn("Contrasts against them ()", brief)
        self.assertIn("no contrast is against a control group", brief)    # the design's note

    def test_recorded_controls_are_not_called_a_name_guess(self):
        d, p = self.session(controls="Parental")
        self.assertIn("Control group(s), as recorded", self.brief(d, p))
        readme = session_docs.agents_md(session_docs.gather(p["session_dir"]))
        self.assertIn("Pull-down controls, as recorded", readme)
        self.assertNotIn("the name is the only evidence", readme)


if __name__ == "__main__":
    unittest.main()
