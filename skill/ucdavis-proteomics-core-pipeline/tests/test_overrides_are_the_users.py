#!/usr/bin/env python3
"""
An override is the user's choice, never a "validated SOP" (staff report, 2026-09-25).

A timsTOF cohort ~20 ppm off calibration was searched at 25 ppm, a one-off workaround the user
chose. resolve_defaults.py --ms1-ppm/--ms2-ppm recorded it as "site SOP override" and
estimate_params.py --overrides as "user-override (validated SOP)", so the provenance -- and the
Methods and run record that quote it -- would have called it a validated SOP.

Every override is now tagged by one wording (estimate_params.override_source): a user override,
its value, who set it by role, and why (--override-reason, else "not given -- confirm"). The
sidecar and the manifest carry the same as a record (`overrides`). They reach the client, so who
set it by name (--override-by, else the login that ran the command, said so) is only in the
staff-only <file>.staff.json beside them (2.10 review, 2026-10-01). Sage overrides, merged into
the config with no rationale before, get one.

The role is "Core staff" only when the setter is a login on a usable Core staff list
(staff.role_for); anyone else -- a login not on the list, a name, or anyone at all when there is
no list to check against -- is "the user" (S4: every override used to say "set by Core staff").
Each test chooses its staff list: none (CORE_STAFF_FILE at a path that does not exist) unless it
writes one, so the tests do not depend on the machine they run on.

Hermetic: estimate_params.py / resolve_defaults.py as subprocesses on a temp dir; no network.
"""
import getpass
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import estimate_params as ep  # noqa: E402
import staff  # noqa: E402

REASON = "instrument ~20 ppm off calibration; user chose 25 ppm for this cohort"
ME = getpass.getuser()


class NoStaffList(unittest.TestCase):
    """CORE_STAFF_FILE at a path that does not exist, unless a test writes a list (staff_list)."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        env = {staff.STAFF_FILE_ENV: os.path.join(self.tmp, "no_staff_list.txt")}
        patch = mock.patch.dict(os.environ, env)
        patch.start()
        self.addCleanup(patch.stop)
        for k in ("SKILL_SLACK_TEST_LOOPBACK", "SKILL_CORE_ADMINS"):
            os.environ.pop(k, None)

    def staff_list(self, *logins):
        """A trusted staff list, built as tests/test_staff_decisions.py builds it: this account
        owns it and is the Core admin, under the test switch."""
        d = os.path.join(self.tmp, "config")
        os.makedirs(d, mode=0o755, exist_ok=True)
        os.chmod(d, 0o755)
        path = os.path.join(d, "core_staff.txt")
        with open(path, "w") as fh:
            fh.write("# the Core's staff\n" + "".join(f"{x}  # someone\n" for x in logins))
        os.chmod(path, 0o644)
        os.environ.update({staff.STAFF_FILE_ENV: path, "SKILL_SLACK_TEST_LOOPBACK": "1",
                           "SKILL_CORE_ADMINS": ME})
        return path


def estimate(d, engine, *extra, ok=True):
    out = os.path.join(d, "params.cfg" if engine == "diann" else "sage.json")
    p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                        "--engine", engine, "--acquisition", "DIA", "--instrument", "timsTOF HT",
                        "--precursor-mz-range", "299.5", "1200.5", "--out", out, *extra],
                       capture_output=True, text=True, timeout=60)
    if not ok:
        return p
    if p.returncode != 0:
        raise AssertionError(p.stderr)
    with open(out + ".rationale.json") as fh:
        return out, json.load(fh)


class EstimateParams(NoStaffList):
    def test_an_override_names_its_value_who_and_why_and_never_an_sop(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "diann", "--overrides",
                                 '{"--mass-acc": 25, "--mass-acc-ms1": 25}',
                                 "--override-by", "core analyst", "--override-reason", REASON)
            with open(cfg) as fh:
                lines = fh.read().splitlines()
            who = staff.read(cfg + staff.STAFF_SUFFIX)
        self.assertIn("--mass-acc 25", lines)
        for flag in ("--mass-acc", "--mass-acc-ms1"):
            src = side["rationale"][flag]["source"]
            self.assertEqual(src, f"user override: {flag} = 25, set by the user; reason: "
                                  f"{REASON}")
            self.assertNotIn("SOP", src)
            self.assertEqual(side["overrides"][flag],
                             {"value": 25, "set_by": "the user", "reason": REASON, "source": src,
                              "staff_record": "params.cfg.staff.json"})
        self.assertEqual([(e["by"], e["by_how"], e["reason"]) for e in who],
                         [("core analyst", "given", REASON)])
        self.assertNotIn("core analyst", json.dumps(side))

    def test_without_by_or_reason_it_says_so(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "diann", "--overrides",
                                 '{"--mass-acc": 25, "--mass-acc-ms1": 25}')
            who = staff.read(cfg + staff.STAFF_SUFFIX)
        src = side["rationale"]["--mass-acc"]["source"]
        self.assertIn("set by the user", src)
        self.assertTrue(src.endswith(ep.OVERRIDE_REASON_MISSING), src)
        self.assertIsNone(side["overrides"]["--mass-acc"]["reason"])
        # the login that ran it -- who TYPED it -- said to be that, in the staff-only record
        self.assertEqual([(e["by"], e["by_how"]) for e in who], [(staff.login(), "login")])

    def test_nothing_overridden_records_nothing(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "diann")
            self.assertFalse(os.path.exists(cfg + staff.STAFF_SUFFIX))
        self.assertEqual(side["overrides"], {})
        self.assertFalse(any("user override" in str(v.get("source", ""))
                             for v in side["rationale"].values() if isinstance(v, dict)))

    def test_by_or_reason_without_an_override_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p = estimate(d, "diann", "--override-reason", REASON, ok=False)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("describe --overrides", p.stderr)

    def test_a_sage_override_gets_a_rationale_and_replaces_the_derived_one(self):
        tol = {"ppm": [-25.0, 25.0]}
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "sage", "--overrides", json.dumps({"precursor_tol": tol}),
                                 "--override-by", "core analyst", "--override-reason", REASON)
            with open(cfg) as fh:
                self.assertEqual(json.load(fh)["precursor_tol"], tol)
        r = side["rationale"]
        self.assertNotIn("precursor_tol_ppm", r, "the derived value is no longer what ran")
        self.assertEqual(r["precursor_tol"]["value"], tol)
        self.assertIn('user override: precursor_tol = {"ppm": [-25.0, 25.0]}, set by the user',
                      r["precursor_tol"]["source"])
        self.assertIn("fragment_tol_ppm", r, "an entry nothing overrode stays")


class WhoSetIt(NoStaffList):
    """S4: "set by Core staff" only for a login on a usable Core staff list."""

    def said(self, side, flag):
        o = side["overrides"][flag]
        return o["set_by"], side["rationale"][flag]["source"]

    def test_listed_staff_override_says_core_staff(self):
        self.staff_list("jdoe", "asmith")
        with tempfile.TemporaryDirectory() as d:
            _, side = estimate(d, "diann", "--overrides", '{"--mass-acc": 25, "--mass-acc-ms1": 25}',
                               "--override-by", "jdoe", "--override-reason", REASON)
        role, src = self.said(side, "--mass-acc")
        self.assertEqual(role, "Core staff")
        self.assertIn("set by Core staff; reason:", src)

    def test_unlisted_login_says_the_user(self):
        self.staff_list("jdoe")
        with tempfile.TemporaryDirectory() as d:
            _, side = estimate(d, "sage", "--overrides",
                               json.dumps({"precursor_tol": {"ppm": [-25.0, 25.0]}}),
                               "--override-by", "asmith", "--override-reason", REASON)
        role, src = self.said(side, "precursor_tol")
        self.assertEqual(role, "the user")
        self.assertIn("set by the user; reason:", src)
        self.assertNotIn("Core staff", src)

    def test_no_list_says_the_user(self):
        with tempfile.TemporaryDirectory() as d:
            _, side = estimate(d, "diann", "--overrides", '{"--mass-acc": 25, "--mass-acc-ms1": 25}',
                               "--override-by", "jdoe")
        role, src = self.said(side, "--mass-acc-ms1")
        self.assertEqual(role, "the user", "with no list, nothing can be verified")
        self.assertIn("set by the user;", src)

    def test_the_login_that_ran_it_is_checked_when_no_by_is_given(self):
        self.staff_list(ME)
        with tempfile.TemporaryDirectory() as d:
            _, side = estimate(d, "diann", "--overrides", '{"--mass-acc": 25, "--mass-acc-ms1": 25}')
        self.assertEqual(self.said(side, "--mass-acc")[0], "Core staff")

    def test_a_caller_that_forgets_the_role_gets_the_less_privileged_label(self):
        self.assertIn("set by the user;", ep.override_source("--mass-acc", 25))
        self.assertEqual(ep.override_record({"--mass-acc": 25})["--mass-acc"]["set_by"],
                         staff.USER_ROLE)
        # the builders, handed override_meta without a role
        _, rat = ep.build_sage("DDA", "timstof", None, {"precursor_tol": {"ppm": [-25, 25]}},
                               override_meta={"reason": "x"})
        self.assertIn("set by the user;", rat["precursor_tol"]["source"]
                      if isinstance(rat["precursor_tol"], dict) else str(rat["precursor_tol"]))
        # an explicit role still wins
        self.assertIn("set by Core staff;",
                      ep.override_source("--mass-acc", 25, role=staff.ROLE))

    def test_role_for(self):
        self.assertEqual(staff.role_for(None), staff.USER_ROLE)
        self.assertEqual(staff.role_for("jdoe"), staff.USER_ROLE)       # no list
        self.staff_list("jdoe")
        self.assertEqual(staff.role_for("jdoe"), staff.ROLE)
        self.assertEqual(staff.role_for("JDoe"), staff.ROLE)            # a login, lowercased
        self.assertEqual(staff.role_for("asmith"), staff.USER_ROLE)
        self.assertEqual(staff.role_for("Jane Doe"), staff.USER_ROLE)   # a name, not a login
        self.assertEqual(staff.USER_ROLE, "the user")

    def test_resolve_defaults_listed_staff_says_core_staff(self):
        self.staff_list("jdoe")
        dest = os.path.join(self.tmp, "wf")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "resolve_defaults.py"),
                            "--engine", "diann", "--acquisition", "DIA",
                            "--instrument", "timsTOF HT", "--dest", dest, "--ms1-ppm", "25",
                            "--ms2-ppm", "25", "--override-by", "jdoe"],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(os.path.join(dest, "workflow.manifest.json")) as fh:
            s = json.load(fh)["search"]
        self.assertEqual(s["overrides"]["--ms1-ppm"]["set_by"], "Core staff")
        self.assertIn("--ms1-ppm = 25, set by Core staff;", s["ppm_source"])


class ResolveDefaults(NoStaffList):
    def resolve(self, d, *extra):
        dest = os.path.join(d, "wf")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "resolve_defaults.py"),
                            "--engine", "diann", "--acquisition", "DIA",
                            "--instrument", "timsTOF HT", "--dest", dest, *extra],
                           capture_output=True, text=True, timeout=60)
        return p, dest

    def test_ppm_overrides_are_the_users_with_who_and_why(self):
        with tempfile.TemporaryDirectory() as d:
            p, dest = self.resolve(d, "--ms1-ppm", "25", "--ms2-ppm", "25",
                                   "--override-by", "core analyst", "--override-reason", REASON)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(os.path.join(dest, "workflow.manifest.json")) as fh:
                s = json.load(fh)["search"]
            who = staff.read(os.path.join(dest, "workflow.manifest.json" + staff.STAFF_SUFFIX))
        self.assertEqual((s["ms1_ppm"], s["ms2_ppm"]), (25.0, 25.0))
        self.assertEqual(s["ppm_source"],
                         f"MS1: user override: --ms1-ppm = 25, set by the user; reason: "
                         f"{REASON}; MS2: user override: --ms2-ppm = 25, set by the user; "
                         f"reason: {REASON}")
        self.assertNotIn("SOP", s["ppm_source"])
        self.assertEqual(s["overrides"]["--ms2-ppm"]["set_by"], "the user")
        self.assertEqual(s["overrides"]["--ms2-ppm"]["staff_record"],
                         "workflow.manifest.json.staff.json")
        self.assertEqual(s["overrides"]["--ms2-ppm"]["reason"], REASON)
        self.assertNotIn("core analyst", json.dumps(s))
        self.assertEqual([(e["by"], e["by_how"]) for e in who], [("core analyst", "given")])

    def test_no_override_no_record_and_by_alone_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p, dest = self.resolve(d)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(os.path.join(dest, "workflow.manifest.json")) as fh:
                s = json.load(fh)["search"]
            self.assertEqual(s["overrides"], {})
            p, _ = self.resolve(d, "--override-by", "core analyst")
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("neither was given", p.stderr)


if __name__ == "__main__":
    unittest.main()
