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
    enrichment, and as a guess (controls_from "enrichment_name"): the brief tags it as not
    confirmed, and the normalisation check holds the decision for the user.

A name is a whole token of the group name: gene symbols that merely contain the letters never
match. Stdlib only.
"""
import csv
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import experiment_type as et  # noqa: E402
import normalization_check as nc  # noqa: E402
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

    def test_a_control_that_is_no_contrasts_reference_is_said_so(self):
        """KO named, but the contrast is KO-Bait: not "none is named like one"."""
        note = et.pulldown_design(ip(), ["Bait", "KO"], ["KO-Bait"])["note"]
        self.assertIn("KO read as the control group(s), but no contrast has one as its reference",
                      note)
        self.assertNotIn("none is named", note)


class TheGuessIsHeldForTheUser(unittest.TestCase):
    """Which contrasts are enrichments decides which volcano anomalies the normalisation check
    expects: controls only guessed from their names are confirmed before anything is decided."""

    def check(self, control, controls=None):
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
        return nc.run(report, cond, fake_de(os.path.join(d, "n"), c, BALANCED),
                      fake_de(os.path.join(d, "r"), c, BALANCED), session=s)[0]

    def test_a_guessed_control_needs_a_decision(self):
        rec = self.check("Parental")
        self.assertFalse(rec["tripped"])
        self.assertTrue(rec["needs_decision"])
        self.assertIsNone(rec["decision"])
        self.assertIn("read from their names only", rec["needs_decision_why"])

    def test_igg_or_recorded_controls_are_decided_as_before(self):
        for control, controls in (("IgG", None), ("Parental", "Parental")):
            rec = self.check(control, controls)
            self.assertFalse(rec["needs_decision"], control)
            self.assertEqual(rec["decision"]["source"], "default", control)


if __name__ == "__main__":
    unittest.main()
