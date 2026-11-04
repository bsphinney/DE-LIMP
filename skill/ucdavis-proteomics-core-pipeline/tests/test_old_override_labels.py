#!/usr/bin/env python3
"""A session from before 2.14.2, bundled again with 2.14.2: its overrides' "set by Core staff"
labels are checked (C1 review, F6).

Before 2.14.2 every override was labelled "set by Core staff" (staff.ROLE), whoever set it (S4,
fixed by T18: staff.role_for). Bundling such a session again with 2.14.2 put a commands.log saying
--override-by "<redacted>" beside a params sidecar (and a workflow manifest) still saying Core
staff -- the label T18 removed, in a client deliverable.

2.14.2 writes the role it gave beside the setter in the staff-only record
(estimate_params.record_override_setter `role`). provenance.py now, copying a params sidecar or a
workflow manifest whose staff-only entry records no role, re-derives it from that entry's setter
with staff.role_for, rebuilds the labels with estimate_params.override_source's wording, and says
so in a MANIFEST.txt NOTE (no name in it). A record whose entry has a role is copied as it is; so
is one whose setter is on the Core staff list now. With no staff-only entry to check, the label is
"the user", said. The session's own files are never changed.

Hermetic: estimate_params.py / resolve_defaults.py / provenance.py as subprocesses on a temp dir,
with no Core staff list unless a test writes one (as tests/test_overrides_are_the_users.py).
"""
import json
import os
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import staff  # noqa: E402
from test_overrides_are_the_users import NoStaffList  # noqa: E402

REASON = "instrument ~20 ppm off calibration"
OVERRIDES = '{"--mass-acc": 25, "--mass-acc-ms1": 25}'


class Base(NoStaffList):
    def run_ok(self, *argv):
        p = subprocess.run([sys.executable, *argv], capture_output=True, text=True, timeout=300)
        self.assertEqual(p.returncode, 0, p.stdout[-1500:] + p.stderr[-1500:])
        return p

    def no_list(self):
        os.environ[staff.STAFF_FILE_ENV] = os.path.join(self.tmp, "no_staff_list.txt")

    def as_2_14_1(self, record_path):
        """The staff-only record as 2.14.1 wrote it: no `role` in its entries."""
        with open(record_path) as fh:
            d = json.load(fh)
        for e in d["entries"]:
            e.pop("role", None)
        with open(record_path, "w") as fh:
            json.dump(d, fh)

    def sidecar(self, by="jdoe"):
        """A params.cfg whose sidecar says "set by Core staff" (`by` was on a list then)."""
        cfg = os.path.join(self.tmp, "wf", "params.cfg")
        os.makedirs(os.path.dirname(cfg), exist_ok=True)
        self.staff_list(by)
        self.run_ok(os.path.join(SCRIPTS, "estimate_params.py"), "--engine", "diann",
                    "--acquisition", "DIA", "--instrument", "timsTOF HT", "--precursor-mz-range",
                    "299.5", "1200.5", "--out", cfg, "--overrides", OVERRIDES, "--override-by",
                    by, "--override-reason", REASON)
        with open(cfg + ".rationale.json") as fh:
            self.assertEqual(json.load(fh)["overrides"]["--mass-acc"]["set_by"], staff.ROLE)
        return cfg

    def bundle(self, *args):
        out = os.path.join(self.tmp, "repro")
        self.run_ok(os.path.join(SCRIPTS, "provenance.py"), "--outdir", out, *args)
        with open(os.path.join(out, "MANIFEST.txt")) as fh:
            notes = [ln for ln in fh.read().splitlines()
                     if ln.startswith("[NOTE]") and "re-labelled" in ln]
        return os.path.join(out, "inputs"), notes

    def read(self, path):
        with open(path) as fh:
            return json.load(fh)


class ParamsSidecar(Base):
    def test_a_2_14_1_sidecar_whose_setter_is_not_verified_says_the_user(self):
        cfg = self.sidecar()
        self.as_2_14_1(cfg + staff.STAFF_SUFFIX)
        self.no_list()
        with open(cfg + ".rationale.json") as fh:
            before = fh.read()
        inputs, notes = self.bundle("--params", cfg)
        side = self.read(os.path.join(inputs, "params.cfg.rationale.json"))
        for flag in ("--mass-acc", "--mass-acc-ms1"):
            o = side["overrides"][flag]
            self.assertEqual(o["set_by"], staff.USER_ROLE)
            self.assertEqual(o["source"], f"user override: {flag} = 25, set by the user; reason: "
                                          f"{REASON}")
            self.assertEqual(side["rationale"][flag]["source"], o["source"])
        self.assertNotIn("Core staff", json.dumps(side))
        self.assertEqual(len(notes), 1, notes)
        self.assertIn("params.cfg.rationale.json", notes[0])
        self.assertIn("the user", notes[0])
        self.assertNotIn("jdoe", notes[0])
        with open(cfg + ".rationale.json") as fh:
            self.assertEqual(fh.read(), before, "the session's own sidecar is not changed")

    def test_a_2_14_1_sidecar_whose_setter_is_on_the_list_now_is_kept(self):
        cfg = self.sidecar()
        self.as_2_14_1(cfg + staff.STAFF_SUFFIX)
        inputs, notes = self.bundle("--params", cfg)          # jdoe's list still in place
        side = self.read(os.path.join(inputs, "params.cfg.rationale.json"))
        self.assertEqual(side["overrides"]["--mass-acc"]["set_by"], staff.ROLE)
        self.assertEqual(notes, [])

    def test_a_2_14_2_sidecar_keeps_its_recorded_role(self):
        cfg = self.sidecar()                                    # its entry records the role
        self.no_list()
        inputs, notes = self.bundle("--params", cfg)
        side = self.read(os.path.join(inputs, "params.cfg.rationale.json"))
        self.assertEqual(side["overrides"]["--mass-acc"]["set_by"], staff.ROLE)
        self.assertEqual(notes, [])

    def test_no_staff_record_to_check_says_the_user(self):
        cfg = self.sidecar()
        os.unlink(cfg + staff.STAFF_SUFFIX)
        self.no_list()
        inputs, notes = self.bundle("--params", cfg)
        side = self.read(os.path.join(inputs, "params.cfg.rationale.json"))
        self.assertEqual(side["overrides"]["--mass-acc"]["set_by"], staff.USER_ROLE)
        self.assertEqual(len(notes), 1, notes)
        self.assertIn("params.cfg.staff.json", notes[0])


    def test_a_staff_record_with_no_override_entry_is_named_as_such(self):
        cfg = self.sidecar()
        rec = cfg + staff.STAFF_SUFFIX
        d = self.read(rec)
        d["entries"] = [{"event": "something_else", "by": "jdoe"}]   # exists, no override entry
        with open(rec, "w") as fh:
            json.dump(d, fh)
        self.no_list()
        inputs, notes = self.bundle("--params", cfg)
        side = self.read(os.path.join(inputs, "params.cfg.rationale.json"))
        self.assertEqual(side["overrides"]["--mass-acc"]["set_by"], staff.USER_ROLE)
        self.assertEqual(len(notes), 1, notes)
        self.assertIn("params.cfg.staff.json) has no override entry", notes[0])
        self.assertNotIn("is not here", notes[0])
        self.assertNotIn("jdoe", notes[0])

    def test_a_missing_staff_record_says_it_is_not_here(self):
        cfg = self.sidecar()
        os.unlink(cfg + staff.STAFF_SUFFIX)
        self.no_list()
        _, notes = self.bundle("--params", cfg)
        self.assertIn("params.cfg.staff.json) is not here", notes[0])


class WorkflowManifest(Base):
    def test_a_2_14_1_manifests_ppm_overrides_say_the_user(self):
        dest = os.path.join(self.tmp, "wf")
        self.staff_list("jdoe")
        self.run_ok(os.path.join(SCRIPTS, "resolve_defaults.py"), "--engine", "diann",
                    "--acquisition", "DIA", "--instrument", "timsTOF HT", "--dest", dest,
                    "--ms1-ppm", "25", "--ms2-ppm", "20", "--override-by", "jdoe",
                    "--override-reason", REASON)
        man = os.path.join(dest, "workflow.manifest.json")
        self.assertIn("set by Core staff", self.read(man)["search"]["ppm_source"])
        self.as_2_14_1(man + staff.STAFF_SUFFIX)
        self.no_list()
        inputs, notes = self.bundle("--workflow-manifest", man)
        s = self.read(os.path.join(inputs, "workflow.manifest.json"))["search"]
        self.assertEqual(s["ppm_source"],
                         f"MS1: user override: --ms1-ppm = 25, set by the user; reason: {REASON}; "
                         f"MS2: user override: --ms2-ppm = 20, set by the user; reason: {REASON}")
        for flag in ("--ms1-ppm", "--ms2-ppm"):
            self.assertEqual(s["overrides"][flag]["set_by"], staff.USER_ROLE)
            self.assertIn("set by the user;", s["overrides"][flag]["source"])
        self.assertEqual(len(notes), 1, notes)
        self.assertIn("workflow.manifest.json", notes[0])


if __name__ == "__main__":
    unittest.main()
