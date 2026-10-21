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
import contextlib
import csv
import io
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
import audit_results  # noqa: E402
import make_methods  # noqa: E402
from test_normalization_check import no_staff_list, session  # noqa: E402
from test_normalization_check_hold import BALANCED, fake_de  # noqa: E402

ANY_STUDY = ("IgG", "Old_IgG", "IgG1", "IgG2a", "IgG2b", "IgG3", "IgG4", "mIgG", "RbIgG",
             "Beads", "Bead-only", "BeadsOnly", "EmptyBeads", "Empty_beads", "Empty-beads",
             "Isotype", "Isotype_IgG", "NoAb", "No_antibody")
ENRICHMENT_ONLY = ("Parental", "Untagged", "Mock", "Empty", "EmptyVector", "EV", "KO", "Bait_KO",
                   "Knockout", "Knock-out", "Knock_out", "Ctrl", "Control", "Neg_Ctrl", "Controls",
                   "NegCtrl")
NOT_A_CONTROL = ("SKO1", "SKOR1", "MOCS1", "CTRB1", "Controller", "Mockup", "Emptor", "IgGFc",
                 "Every", "Nobody", "WT", "Bait", "Treated")


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
            self.assertEqual((pd["pulldown"], pd["controls"], pd["note"]), (False, [], None),
                             groups)
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
        # recorded as none: no control, no guess, nothing to raise
        pd = et.pulldown_design(ip(controls=[]), ["Bait_Ctrl", "Bait_Drug"],
                                ["Bait_Drug-Bait_Ctrl"])
        self.assertEqual((pd["controls"], pd["controls_from"], pd["note"]), ([], "recorded", None))

    def test_a_control_that_is_no_contrasts_reference_is_said_so(self):
        """KO named, but the contrast is KO-Bait: not "none is named like one"."""
        note = et.pulldown_design(ip(), ["Bait", "KO"], ["KO-Bait"])["note"]
        self.assertIn("KO read as the control group(s), but no contrast has one as its reference",
                      note)
        self.assertNotIn("none is named", note)


class TheGuessIsHeldForTheUser(unittest.TestCase):
    """Which contrasts are enrichments decides which volcano anomalies the normalisation check
    expects: controls only guessed from their names are confirmed before anything is decided."""

    def inputs(self, control, controls=None, bait="Bait", contrast=None):
        """A session recorded as an IP, a header-only report (B1-B3 not assessed: no pyarrow
        needed), conditions.csv and the two check DEs -> the `run` arguments."""
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        s = session(d, "ip", **({"controls": controls} if controls else {}))
        report = os.path.join(d, "report.tsv")
        with open(report, "w") as fh:
            fh.write("Run\tPrecursor.Id\tProtein.Group\tPrecursor.Quantity\t"
                     "Precursor.Normalised\tQ.Value\n")
        cond = os.path.join(s, "input", "conditions.csv")      # set-controls checks names here
        with open(cond, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["File.Name", "Group"])
            w.writerows([f"{g}_{i}", g] for g in (bait, control) for i in (1, 2, 3))
        c = contrast or f"{bait}-{control}"
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
        with self.assertRaisesRegex(ValueError, "nothing to decide yet.*set-controls"):
            nc.decide(os.path.join(out, nc.RECORD), et.NORMALISED, "someone", "a guess")
        # and nothing the user reads sends them to `decide` for it
        ok, msg, _rec = nc.gate(os.path.join(out, nc.RECORD), et.RAW)
        self.assertFalse(ok)
        self.assertIn("waiting for the control groups", msg)
        self.assertIn("set-controls", msg)
        self.assertNotIn("decide", msg)
        with open(os.path.join(out, nc.SUMMARY)) as fh:
            summary = fh.read()
        self.assertIn("**Decision: not possible yet.**", summary)
        self.assertNotIn("normalization_check.py decide", summary)

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
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "experiment_type.py"),
                            "set-controls", "--session", s, "--controls", "Parental"],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(p.returncode, 0, p.stderr)
        rec = nc.run(report, cond, dn, dr, session=s)[0]
        self.assertEqual((rec["needs_decision"], rec["controls_unconfirmed"]), (False, False))
        self.assertEqual(rec["decision"]["source"], "default")
        self.assertEqual(rec["design"]["controls_from"], "recorded")

    def test_a_design_with_no_control_group_says_none(self):
        """Bait_Drug vs Bait_Ctrl with no IgG: Ctrl is the untreated bait, not a control."""
        report, cond, dn, dr, s = self.inputs("Bait_Ctrl", bait="Bait_Drug")
        rec = nc.run(report, cond, dn, dr, session=s)[0]
        self.assertTrue(rec["controls_unconfirmed"])
        et.record_controls(s, et.NO_CONTROLS)
        self.assertEqual(et.load(s)["controls"], [])
        rec = nc.run(report, cond, dn, dr, session=s)[0]
        self.assertEqual((rec["needs_decision"], rec["decision"]["source"]), (False, "default"))
        pd = rec["design"]
        self.assertEqual((pd["controls"], pd["controls_from"], pd["vs_control"], pd["note"],
                          pd["no_control"]), ([], "recorded", [], None, True))
        # the default stands, for a reason that names no control -- and the Methods say so
        self.assertEqual(rec["default"]["why"], et.NO_CONTROL["ip"][1])
        self.assertEqual(et.load(s)["label"], et.NO_CONTROL["ip"][0])
        methods = make_methods.de_normalisation_sentence(
            {"normalisation": "non-normalised",
             "normalization_check": dict(rec, status="decided")})
        # one consistent sentence: no "with IgG or bead controls" beside "no negative control"
        self.assertIn("pull-down, with no negative control group (per the user)", methods)
        for never in ("IgG or bead controls", "the controls carry"):
            self.assertNotIn(never, methods)

    def test_the_gate_refuses_controls_changed_since_the_check(self):
        report, cond, dn, dr, s = self.inputs("Parental", controls="Parental")
        rec, out = nc.run(report, cond, dn, dr, session=s)
        path = os.path.join(out, nc.RECORD)
        self.assertTrue(nc.gate(path, rec["decision"]["quantities"])[0])
        # a later `set` that names no --controls keeps them: the check still stands
        et.record(s, "ip", "user")
        self.assertEqual(et.load(s)["controls"], ["Parental"])
        self.assertTrue(nc.gate(path, rec["decision"]["quantities"])[0])
        # set-controls none: the check's volcano read Parental as a control -- run it again
        et.record_controls(s, et.NO_CONTROLS)
        ok, msg, _rec = nc.gate(path, rec["decision"]["quantities"])
        self.assertFalse(ok)
        self.assertIn("the control groups changed since this check ran", msg)

    def test_a_tripped_check_with_guessed_controls_is_not_decided(self):
        """The trips were computed with the guessed controls: confirm them first."""
        report, cond, dn, dr, s = self.inputs("Parental")
        no_staff_list(self, os.path.dirname(s))
        _rec, out = nc.run(report, cond, dn, dr, session=s)
        path = os.path.join(out, nc.RECORD)
        with open(path) as fh:
            rec = json.load(fh)
        rec.update(tripped=True, trips=["B1 a stand-in trip"])
        with open(path, "w") as fh:
            json.dump(rec, fh)
        with self.assertRaisesRegex(ValueError, "its trips with them"):
            nc.decide(path, et.NORMALISED, "canalyst", "the staff member's reason")
        self.assertIn("waiting for the control groups", nc.gate(path, et.NORMALISED)[1])

    def test_a_check_that_trips_still_asks_for_the_controls(self):
        rec, _out = self.check("Parental")
        rec.update(tripped=True, trips=["B1 a stand-in trip"])
        summary = nc.summary_md(rec)
        self.assertIn("**Decision: not possible yet.**", summary)
        self.assertNotIn("normalization_check.py decide", summary)

    def test_a_guess_holds_even_when_no_contrast_uses_it(self):
        report, cond, dn, dr, s = self.inputs("KO", contrast="KO-Bait")
        rec = nc.run(report, cond, dn, dr, session=s)[0]
        self.assertEqual((rec["design"]["vs_control"], rec["controls_unconfirmed"],
                          rec["needs_decision"]), ([], True, True))
        self.assertIn("no contrast has one as its reference", rec["needs_decision_why"])

    def test_the_final_des_contrasts_against_a_control_must_be_the_checks(self):
        report, cond, dn, dr, s = self.inputs("Parental", controls="Parental")
        rec, out = nc.run(report, cond, dn, dr, session=s)
        path, q = os.path.join(out, nc.RECORD), rec["decision"]["quantities"]
        self.assertTrue(nc.gate(path, q, design={"contrasts": ["Bait-Parental"]})[0])
        ok, msg, _rec = nc.gate(path, q, design={"contrasts": ["Parental-Bait"]})
        self.assertFalse(ok)
        self.assertIn("contrasts against a control (none) are not the check's (Bait-Parental)",
                      msg)

    def test_the_gate_refuses_a_changed_or_missing_type(self):
        report, cond, dn, dr, s = self.inputs("IgG")
        rec, out = nc.run(report, cond, dn, dr, session=s)
        path, q = os.path.join(out, nc.RECORD), rec["decision"]["quantities"]
        et.record(s, "proximity", "user")
        self.assertIn("the experiment type changed since this check ran (ip -> proximity)",
                      nc.gate(path, q)[1])
        os.remove(et.path_for(s))
        self.assertIn("there is no experiment-type record at", nc.gate(path, q)[1])

    def test_a_moved_session_is_the_same_session(self):
        """The gate reads the session the check sits in: moving it changes nothing. A check
        whose session is gone says so, not that the type changed."""
        report, cond, dn, dr, s = self.inputs("IgG")
        out = os.path.join(s, "output", "norm_check")
        os.makedirs(out)
        rec, _out = nc.run(report, cond, dn, dr, session=s, out=out)
        q = rec["decision"]["quantities"]
        moved = s + "_moved"
        shutil.move(s, moved)
        self.assertTrue(nc.gate(os.path.join(moved, "output", "norm_check", nc.RECORD), q)[0])
        elsewhere = os.path.join(os.path.dirname(s), "check_copy")
        shutil.copytree(os.path.join(moved, "output", "norm_check"), elsewhere)
        ok, msg, _rec = nc.gate(os.path.join(elsewhere, nc.RECORD), q)
        self.assertFalse(ok)
        self.assertIn(f"the session this check ran in ({s}) is not there any more", msg)

    def test_a_check_run_without_the_session_is_not_compared_with_one(self):
        """Saved in an IP session's output/norm_check, but run without --session: the gate
        compares no session (as before 2.14.1). Where a session is compared (the delivery), it
        says the check never read its type -- not that the type changed."""
        report, cond, dn, dr, s = self.inputs("IgG")
        out = os.path.join(s, "output", "norm_check")
        os.makedirs(out)
        rec = nc.run(report, cond, dn, dr, out=out)[0]
        self.assertIsNone(rec.get("session"))
        msg = nc.gate(os.path.join(out, nc.RECORD), et.RAW)[1] or ""
        self.assertIn("needs a decision", msg)              # the usual path, no comparison
        for never in ("ran without --session", "changed since this check ran"):
            self.assertNotIn(never, msg)
        said = nc.controls_changed(rec, session=s, legacy=False)
        self.assertIn("this check ran without --session, so it did not read the experiment type "
                      f"recorded in {s} (ip)", said)
        self.assertNotIn("changed", said)

    def test_a_check_from_before_2_14_1_is_compared_too(self):
        """It recorded controls and vs_control but no controls_from: a guess since is refused."""
        report, cond, dn, dr, s = self.inputs("Parental", controls="Parental")
        rec, out = nc.run(report, cond, dn, dr, session=s)
        path, q = os.path.join(out, nc.RECORD), rec["decision"]["quantities"]
        with open(path) as fh:
            old = json.load(fh)
        del old["design"]["controls_from"], old["design"]["no_control"]
        with open(path, "w") as fh:
            json.dump(old, fh)
        self.assertTrue(nc.gate(path, q)[0])                   # the same controls: it stands
        with open(et.path_for(s)) as fh:
            r = json.load(fh)
        r.update(controls=None)                                # a 2.13.2 session: none recorded
        with open(et.path_for(s), "w") as fh:
            json.dump(r, fh)
        self.assertIn("a guess this check (from before 2.14.1) never saw", nc.gate(path, q)[1])

    def test_the_delivery_holds_an_analysis_whose_controls_changed(self):
        report, cond, dn, dr, s = self.inputs("Parental", controls="Parental")
        rec, _out = nc.run(report, cond, dn, dr, session=s)
        tables = os.path.join(s, "output", "tables")
        os.makedirs(tables)
        with open(os.path.join(tables, "de_provenance.json"), "w") as fh:
            json.dump({"normalization_check": dict(rec, status="decided",
                                                   quantities_applied=et.RAW),
                       "design_request": {"contrasts": ["Bait-Parental"]}}, fh)
        self.assertTrue(nc.delivery_gate(s)["proceed"])
        et.record_controls(s, et.NO_CONTROLS)
        g = nc.delivery_gate(s)
        self.assertEqual((g["proceed"], g["status"]), (False, "stale"))
        self.assertIn("the control groups changed since this check ran", g["reason"])

    def test_an_unreadable_record_is_said_at_the_gate(self):
        report, cond, dn, dr, s = self.inputs("IgG")
        rec, out = nc.run(report, cond, dn, dr, session=s)
        with open(et.path_for(s), "w") as fh:
            fh.write("{not json")
        self.assertIn("the experiment-type record could not be read",
                      nc.gate(os.path.join(out, nc.RECORD), rec["decision"]["quantities"])[1])

    def test_a_staff_choice_equal_to_the_default_is_said_too(self):
        report, cond, dn, dr, s = self.inputs("IgG")
        no_staff_list(self, os.path.dirname(s))
        _rec, out = nc.run(report, cond, dn, dr, session=s)
        rec = nc.decide(os.path.join(out, nc.RECORD), et.RAW, "canalyst", "as the default says")
        findings = []
        block = dict(rec, status="decided", quantities_applied=et.RAW)
        audit_results.audit_normalization_check(findings, {"normalization_check": block})
        f = [x for x in findings if x["check"] == "normalization_check"]
        self.assertEqual(f[0]["status"], "PASS")
        self.assertIn("chose non-normalised (raw) quantities because as the default says",
                      f[0]["message"])

    def test_a_staff_choice_on_an_agreeing_check_is_said_as_theirs(self):
        """Before: the Methods and AUDIT.md said "its default; the data check agreed"."""
        report, cond, dn, dr, s = self.inputs("IgG")
        no_staff_list(self, os.path.dirname(s))
        rec, out = nc.run(report, cond, dn, dr, session=s)
        self.assertEqual((rec["tripped"], rec["decision"]["source"]), (False, "default"))
        rec = nc.decide(os.path.join(out, nc.RECORD), et.NORMALISED, "canalyst",
                        "the eluates were loaded to equal protein")
        block = dict(rec, status="decided", quantities_applied=et.NORMALISED)
        methods = make_methods.de_normalisation_sentence(
            {"normalisation": "normalised", "normalization_check": block})
        self.assertIn("chose normalised quantities because the eluates were loaded", methods)
        findings = []
        audit_results.audit_normalization_check(findings, {"normalization_check": block})
        f = [x for x in findings if x["check"] == "normalization_check"]
        self.assertEqual([x["status"] for x in f], ["WARN"])
        self.assertNotIn("normalised quantities, its default", f[0]["message"])
        self.assertIn("chose normalised quantities", f[0]["message"])

    def test_set_controls_names_the_sessions_groups(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        s = session(d, "ip")
        with open(os.path.join(s, "input", "conditions.csv"), "w", newline="") as fh:
            csv.writer(fh).writerows([["File.Name", "Group"], ["r1", "Bait"], ["r2", "Parental"]])
        with self.assertRaisesRegex(ValueError, "parental: not a group.*Bait, Parental"):
            et.record_controls(s, "parental")
        err = io.StringIO()
        os.makedirs(os.path.join(s, "output", "norm_check"))
        open(os.path.join(s, "output", "norm_check", "normalization_check.json"), "w").close()
        with contextlib.redirect_stderr(err):
            et.record_controls(s, "Parental")
        self.assertIn("run it again", err.getvalue())                # a check that read others
        # and a re-recorded type keeps the answer, as it keeps the bait's
        et.record(s, "ip", "user", stated="an IP")
        self.assertEqual(et.load(s)["controls"], ["Parental"])

    def test_set_controls_is_for_an_enrichment_and_a_readable_conditions_csv(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        s = session(d, "whole_proteome")
        with self.assertRaisesRegex(ValueError, "not an enrichment"):
            et.record_controls(s, "KO")
        et.record(s, "ip", "user")
        with open(os.path.join(s, "input", "conditions.csv"), "w", newline="") as fh:
            csv.writer(fh).writerows([["File.Name", "Condition"], ["r1", "Bait"]])
        with self.assertRaisesRegex(ValueError, "has no Group column"):
            et.record_controls(s, "Bait")
        with open(os.path.join(s, "input", "conditions.csv"), "w", newline="") as fh:
            csv.writer(fh).writerows([["File.Name", "Group"], ["r1", "Bait"], ["r2", "IgG"]])
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            et.record_controls(s, et.NO_CONTROLS)
        self.assertIn("IgG is named like one", err.getvalue())      # 'none' beside an IgG group

    def test_a_record_from_before_2_14_1_reads_with_todays_label(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        s = session(d, "ip")
        with open(et.path_for(s)) as fh:
            r = json.load(fh)
        r["label"] = ("immunoprecipitation / affinity purification / pull-down with IgG or bead "
                      "controls")
        with open(et.path_for(s), "w") as fh:
            json.dump(r, fh)
        self.assertEqual(et.load(s)["label"], et.DEFAULTS["ip"][0])

    def test_set_controls_changes_nothing_else(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        s = session(d, "ip", bait="BAITX", stated="an IP of BAITX")
        with self.assertRaisesRegex(ValueError, "no input/conditions.csv"):
            et.record_controls(s, "Parental")
        with open(os.path.join(s, "input", "conditions.csv"), "w", newline="") as fh:
            csv.writer(fh).writerows([["File.Name", "Group"], ["r1", "Bait"], ["r2", "Parental"],
                                      ["r3", "Mock"]])
        before = et.load(s)
        et.record_controls(s, "Parental, Mock")
        after = et.load(s)
        self.assertEqual(after["controls"], ["Parental", "Mock"])
        for k in ("type", "source", "stated", "bait", "recorded_at"):
            self.assertEqual(after[k], before[k], k)
        with self.assertRaisesRegex(ValueError, "or 'none'"):
            et.record_controls(s, " ")

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
        agents = session_docs.agents_md(session_docs.gather(p["session_dir"]))
        self.assertIn(f"Pull-down controls, {said}", agents)
        self.assertIn("the name is the only evidence these are controls", agents)
        # and the type's label names no particular control
        self.assertNotIn("IgG or bead", et.load(p["session_dir"])["label"])

    def test_no_control_group_prints_no_empty_control_line(self):
        d, p = self.session(other="Lysate")
        brief = self.brief(d, p)
        self.assertIn("## Pull-down design", brief)
        self.assertNotIn("Control group(s)", brief)
        self.assertNotIn("Contrasts against them ()", brief)
        self.assertIn("no contrast is against a control group", brief)    # the design's note

    def test_a_design_recorded_with_no_control_group_says_so(self):
        d, p = self.session(controls="none", other="Ctrl")
        brief = self.brief(d, p)
        self.assertIn("## Pull-down design — no negative control group", brief)
        self.assertIn("Read each contrast as a comparison between pull-downs", brief)
        for never in ("enrichment over a control", "enrichment over the control",
                      "enriched over the control", "### Specificity"):
            self.assertNotIn(never, brief)
        self.assertIn("### Cross-Comparison Biomarkers", brief)
        agents = session_docs.agents_md(session_docs.gather(p["session_dir"]))
        self.assertIn("**No negative control group** (recorded by the user", agents)
        readme = session_docs.readme_md(session_docs.gather(p["session_dir"]))
        self.assertIn("Negative control: none -- the user recorded", readme)

    def test_recorded_controls_are_not_called_a_name_guess(self):
        d, p = self.session(controls="Parental")
        self.assertIn("Control group(s), as recorded", self.brief(d, p))
        agents = session_docs.agents_md(session_docs.gather(p["session_dir"]))
        self.assertIn("Pull-down controls, as recorded", agents)
        self.assertNotIn("the name is the only evidence", agents)


if __name__ == "__main__":
    unittest.main()
