#!/usr/bin/env python3
"""
A Core decision is a person's, and their name stays with the staff (2.10 review, 2026-10-01).

  * `--by` (normalization_check.py decide / ack-legacy; qc_bracket.py ack shares staff.py) must
    be a HIVE login in the Core's staff list (CORE_STAFF_FILE, trusted only when a Core admin owns
    it and nobody else can write it); with no usable list an agent's name, or anything not shaped
    like a login (a full name), is still refused, and a login passes with a warning saying it was
    not checked;
  * client documents -- the check record copied into de_provenance.json, methods.txt, the
    Methods, AUDIT.md, the cfg sidecar and manifest in reproducibility/inputs/, commands.log, the
    session zip -- say "Core staff"; the name and login are only in *.staff.json records, which
    deliver never ships and the zip leaves out;
  * a DE from before 2.10 (no check record) passes the delivery gate once staff record
    `ack-legacy` for that DE record, and only for it;
  * MaxLFQ's non-normalised DE has a next step: dpc, or a --no-norm report given as --report-raw,
    to which a raw decision then holds the final DE.

Hermetic: temp dirs, subprocesses of the skill's own scripts, synthetic logins (jdoe, asmith).
"""
import getpass
import json
import os
import shutil
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

import experiment_type as et  # noqa: E402
import normalization_check as nc  # noqa: E402
import staff  # noqa: E402
from test_normalization_check import record, session  # noqa: E402

CHECK = os.path.join(SCRIPTS, "normalization_check.py")
ME = getpass.getuser()
REASON = "the IgG controls carry little protein"


class StaffEnv(unittest.TestCase):
    """Each test chooses its staff list: none (a path that does not exist) unless it writes one."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.env = mock.patch.dict(os.environ, {
            staff.STAFF_FILE_ENV: os.path.join(self.tmp, "no_staff_list.txt")})
        self.env.start()
        self.addCleanup(self.env.stop)
        for k in ("SKILL_SLACK_TEST_LOOPBACK", "SKILL_CORE_ADMINS"):
            os.environ.pop(k, None)

    def staff_list(self, *logins, mode=0o644):
        """A staff list this test's account owns, made trusted the way slack_collab's tests are:
        the test switch and this account as the Core admin."""
        d = os.path.join(self.tmp, "config")
        os.makedirs(d, mode=0o755, exist_ok=True)
        os.chmod(d, 0o755)
        path = os.path.join(d, "core_staff.txt")
        with open(path, "w") as fh:
            fh.write("# the Core's staff\n" + "".join(f"{x}  # someone\n" for x in logins))
        os.chmod(path, mode)
        os.environ.update({staff.STAFF_FILE_ENV: path, "SKILL_SLACK_TEST_LOOPBACK": "1",
                           "SKILL_CORE_ADMINS": ME})
        return path


class WhoMaySign(StaffEnv):
    AGENTS = ("claude", "Claude Code", "claude-opus", "assistant", "the skill", "AI agent", "bot",
              "my-bot", "x.ai", "ai", "ChatGPT", "Copilot")

    def test_without_a_list_an_agents_name_is_refused_as_a_whole_part_only(self):
        for name in self.AGENTS:
            c = staff.check_by(name)
            self.assertFalse(c["ok"], name)
            self.assertIn("names an agent", c["error"])
        for person in ("talbot", "kaiser", "claudette", "abbott", "haig", "jdoe"):
            c = staff.check_by(person)
            self.assertTrue(c["ok"], person)
            self.assertTrue(c["checked"].startswith("NOT VERIFIED:"), c)

    def test_a_trusted_list_is_authoritative_both_ways(self):
        self.staff_list("jdoe", "claude.b")
        c = staff.check_by("claude.b")                # a person called Claude, on the list
        self.assertTrue(c["ok"], c)
        self.assertTrue(c["checked"].startswith("on the Core staff list"))
        for name in ("claude", "talbot", "my-bot"):   # not on it: refused, whatever it looks like
            c = staff.check_by(name)
            self.assertFalse(c["ok"], name)
            self.assertIn("is not on the Core staff list", c["error"])

    def test_a_name_is_not_a_login_list_or_no_list(self):
        for name in ("Staff Example", "Jane Doe", "jdoe (Jane)", "1jdoe", ""):
            c = staff.check_by(name)
            self.assertFalse(c["ok"], name)
            self.assertIsNone(c["checked"], name)
        self.assertIn("is not a HIVE login", staff.check_by("Staff Example")["error"])
        self.staff_list("jdoe")
        self.assertIn("is not a HIVE login", staff.check_by("Jane Doe")["error"])

    def test_without_a_list_a_login_passes_said_to_be_unchecked(self):
        c = staff.check_by("JDoe")
        self.assertEqual((c["ok"], c["by"]), (True, "jdoe"))
        self.assertTrue(c["checked"].startswith("NOT VERIFIED:"), c)
        self.assertIn("does not exist", c["checked"])
        self.assertIn("the Core staff list is not set up", c["warning"])

    def test_a_trusted_list_decides(self):
        path = self.staff_list("jdoe", "asmith")
        c = staff.check_by("JDoe")
        self.assertEqual((c["ok"], c["by"], c["checked"]),
                         (True, "jdoe", f"on the Core staff list ({path})"))
        c = staff.check_by("bjones")
        self.assertFalse(c["ok"])
        self.assertIn("is not on the Core staff list", c["error"])
        self.assertIsNone(c["warning"])
        with self.assertRaises(ValueError):
            staff.require("bjones")

    def test_a_list_others_can_write_is_not_used(self):
        self.staff_list("jdoe", mode=0o664)
        c = staff.check_by("bjones")
        self.assertTrue(c["ok"])
        self.assertIn("is not trusted", c["checked"])
        self.assertIn("other than its owner can write it", c["checked"])

    def test_a_command_log_keeps_the_command_and_loses_the_name(self):
        log = ('python3 normalization_check.py decide --check x.json --quantities raw --by jdoe '
               '--reason "low IgG"\n'
               'python3 qc_bracket.py ack --session S --by "Jane Doe" --note "ok"\n'
               "python3 estimate_params.py --override-by='Jane Doe' --override-reason x\n"
               'python3 x.py --by="a \\"quoted\\" name" --next\n'
               "Rscript run_de.R --block Mouse --input r.parquet\n")
        out = staff.redact(log)
        self.assertNotIn("jdoe", out)
        self.assertNotIn("Jane Doe", out)
        self.assertNotIn("quoted", out)
        self.assertEqual(out.count('"<Core staff>"'), 4)
        self.assertIn('--reason "low IgG"', out)
        self.assertIn("--block Mouse --input r.parquet", out)


def check_record(path, sess=None, **extra):
    """A tripped check record (normalised default) at `path`, undecided."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    rec = dict(record(et.NORMALISED), schema=nc.SCHEMA, schema_version=1, tripped=True,
               trips=["B1 normalisation factors 4-fold apart"], decision=None,
               session=sess, **extra)
    with open(path, "w") as fh:
        json.dump(rec, fh)
    return path


class TheCoresStaffFile(StaffEnv):
    # The Core's list on HIVE, byte for byte as written 2026-10-02
    # (/quobyte/proteomics-grp/.config/core_staff.txt: brettsp, 644, in a 2750 folder of his).
    REAL = (b"# Core staff HIVE logins allowed to sign off (QC ack, normalization decide, "
            b"ack-legacy). Admin-owned; one login per line.\nbrettsp\nmsalemi\ngabrig\n")

    def test_the_hive_format_is_read_exactly(self):
        d = os.path.join(self.tmp, "config")
        os.makedirs(d)
        os.chmod(d, 0o2750)
        path = os.path.join(d, "core_staff.txt")
        with open(path, "wb") as fh:
            fh.write(self.REAL)
        os.chmod(path, 0o644)
        os.environ.update({staff.STAFF_FILE_ENV: path, "SKILL_SLACK_TEST_LOOPBACK": "1",
                           "SKILL_CORE_ADMINS": ME})
        self.assertEqual(staff.load_staff(), (["brettsp", "msalemi", "gabrig"], None))
        for login in ("brettsp", "msalemi", "gabrig", "MSALEMI"):
            c = staff.check_by(login)
            self.assertTrue(c["ok"], c)
            self.assertEqual(c["checked"], f"on the Core staff list ({path})")
        for login in ("talbot", "claude", "core", "staff"):       # the comment's words included
            self.assertFalse(staff.check_by(login)["ok"], login)


class DecideIsThePersons(StaffEnv):
    def in_session(self):
        s = os.path.join(self.tmp, "S")
        return s, check_record(os.path.join(s, "output", "norm_check", nc.RECORD), s)

    def test_claude_cannot_decide(self):
        _, path = self.in_session()
        for name in ("claude", "Claude (on behalf of the user)"):
            with self.assertRaises(ValueError) as cm:
                nc.decide(path, et.RAW, name, REASON)
            self.assertIn("never Claude's", str(cm.exception))
        self.assertIsNone(nc.load(path)["decision"])

    def test_with_a_list_only_a_listed_login_decides(self):
        self.staff_list("jdoe")
        _, path = self.in_session()
        with self.assertRaises(ValueError):
            nc.decide(path, et.RAW, "bjones", REASON)
        self.assertEqual(nc.decide(path, et.RAW, "jdoe", REASON)["decision"]["quantities"], et.RAW)

    def test_the_record_names_the_role_and_the_staff_record_the_person(self):
        self.staff_list("jdoe")
        s, path = self.in_session()
        d = nc.decide(path, et.RAW, "jdoe", "Because " + REASON + ".")["decision"]
        self.assertEqual(d["by"], staff.ROLE)
        self.assertEqual(d["statement"], "Core staff chose non-normalised (raw) quantities because "
                                         "the IgG controls carry little protein.")
        self.assertEqual(d["staff_record"], os.path.join("logs", nc.STAFF_RECORD))
        for f in (path, os.path.join(os.path.dirname(path), nc.SUMMARY)):
            with open(f) as fh:
                text = fh.read()
            self.assertNotIn("jdoe", text, f)
            self.assertIn("Core staff chose non-normalised (raw) quantities", text)
        entries = staff.read(os.path.join(s, "logs", nc.STAFF_RECORD))
        self.assertEqual([(e["event"], e["by"], e["account"], e["quantities"]) for e in entries],
                         [("decide", "jdoe", ME, et.RAW)])
        self.assertTrue(entries[0]["staff_check"].startswith("on the Core staff list"))
        ok, msg, _ = nc.gate(path, et.RAW)
        self.assertTrue(ok)
        self.assertNotIn("jdoe", msg)

    def test_outside_a_session_the_staff_record_sits_beside_the_check(self):
        path = check_record(os.path.join(self.tmp, "loose", nc.RECORD))
        d = nc.decide(path, et.RAW, "jdoe", REASON)["decision"]
        self.assertEqual(d["staff_record"], nc.STAFF_RECORD)
        self.assertEqual(staff.read(os.path.join(self.tmp, "loose", nc.STAFF_RECORD))[0]["by"],
                         "jdoe")

    def test_the_methods_and_audit_say_core_staff(self):
        import make_methods
        _, path = self.in_session()
        rec = nc.decide(path, et.RAW, "jdoe", REASON)
        rec["experiment_type"] = {"type": "ip", "label": "IP / pull-down", "source": "user"}
        prov = {"normalisation": "none", "normalization_check": dict(rec, status="decided",
                                                                     quantities_applied=et.RAW)}
        sentence = make_methods.de_normalisation_sentence(prov)
        self.assertIn("and Core staff chose non-normalised (raw) quantities because the IgG "
                      "controls carry little protein.", sentence)
        self.assertNotIn("jdoe", sentence)
        tables = os.path.join(self.tmp, "tables")
        os.makedirs(tables)
        with open(os.path.join(tables, "de_provenance.json"), "w") as fh:
            json.dump(prov, fh)
        a = subprocess.run([sys.executable, os.path.join(SCRIPTS, "audit_results.py"), "--out",
                            "AUDIT.md", "--de-dir", tables], capture_output=True, text=True,
                           timeout=120, cwd=self.tmp)
        self.assertEqual(a.returncode, 0, a.stderr)
        for f in ("AUDIT.md", "AUDIT.json"):
            with open(os.path.join(self.tmp, f)) as fh:
                text = fh.read()
            self.assertNotIn("jdoe", text, f)
            self.assertIn("Core staff chose non-normalised (raw) quantities", text, f)

    def test_the_cli_refuses_an_agent(self):
        _, path = self.in_session()
        p = subprocess.run([sys.executable, CHECK, "decide", "--check", path, "--quantities", "raw",
                            "--by", "claude", "--reason", REASON], capture_output=True, text=True,
                           timeout=60)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("names an agent", p.stderr)


class AckLegacy(StaffEnv):
    def legacy(self, prov=None):
        s = os.path.join(self.tmp, "S")
        os.makedirs(os.path.join(s, "output", "tables"), exist_ok=True)
        with open(nc.de_provenance_path(s), "w") as fh:
            json.dump(prov if prov is not None else {"method": "dpc", "skill_version": "2.9.0"}, fh)
        return s

    def test_a_pre_2_10_de_passes_once_staff_say_it_stands(self):
        s = self.legacy()
        g = nc.delivery_gate(s)
        self.assertFalse(g["proceed"])
        self.assertIn("predates", g["reason"])
        self.assertIn("ack-legacy --session <S> --by <their HIVE login>", g["hint"])
        g = nc.ack_legacy(s, "jdoe", "searched and analysed under 2.9; quantities reviewed")
        self.assertEqual((g["proceed"], g["status"]), (True, "legacy_acknowledged"))
        self.assertIn("Core staff acknowledged", g["reason"])
        self.assertNotIn("jdoe", g["reason"])
        e = staff.read(os.path.join(s, "logs", nc.STAFF_RECORD))
        self.assertEqual([(x["event"], x["by"], x["account"]) for x in e],
                         [("ack-legacy", "jdoe", ME)])
        self.assertEqual(e[0]["de_provenance_sha256"], nc.sha256(nc.de_provenance_path(s)))

    def test_the_acknowledgement_covers_that_de_record_only(self):
        s = self.legacy()
        nc.ack_legacy(s, "jdoe", "reviewed")
        self.legacy({"method": "maxlfq"})                 # the DE ran again (still no check)
        g = nc.delivery_gate(s)
        self.assertFalse(g["proceed"])
        self.assertIn("ack-legacy", g["hint"])

    def test_a_de_under_the_check_cannot_be_acknowledged_away(self):
        for st in ("not_run", "check_input"):
            s = self.legacy({"normalization_check": {"status": st}})
            with self.assertRaises(ValueError) as cm:
                nc.ack_legacy(s, "jdoe", "reviewed")
            self.assertIn("only for a DE from before skill 2.10", str(cm.exception))
            self.assertFalse(nc.delivery_gate(s)["proceed"])

    def test_no_de_no_agent_no_unlisted_login_one_line_note(self):
        with self.assertRaises(ValueError) as cm:
            nc.ack_legacy(os.path.join(self.tmp, "empty"), "jdoe", "reviewed")
        self.assertIn("no DE record", str(cm.exception))
        s = self.legacy()
        with self.assertRaises(ValueError):
            nc.ack_legacy(s, "claude", "reviewed")
        with self.assertRaises(ValueError):
            nc.ack_legacy(s, "jdoe", "two\nlines")
        self.staff_list("jdoe")
        with self.assertRaises(ValueError):
            nc.ack_legacy(s, "bjones", "reviewed")
        self.assertFalse(nc.delivery_gate(s)["proceed"])
        p = subprocess.run([sys.executable, CHECK, "ack-legacy", "--session", s, "--by", "jdoe",
                            "--note", "reviewed"], capture_output=True, text=True, timeout=60)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertTrue(json.loads(p.stdout)["gate"]["proceed"])

    def test_core_submission_deliver_lets_the_acknowledged_de_through(self):
        from test_core_submission import DeliverBase

        class Legacy(DeliverBase):
            def runTest(self):
                pass
        t = Legacy()
        t.setUp()
        try:
            with open(os.path.join(t.output, "tables", "de_provenance.json"), "w") as fh:
                json.dump({"method": "dpc"}, fh)
            # a stray copy of a staff record among the tables is never shipped
            with open(os.path.join(t.output, "tables", nc.STAFF_RECORD), "w") as fh:
                fh.write("{}")
            rc, out, p = t.deliver()
            self.assertNotEqual(rc, 0, p.stdout)
            self.assertIn("ack-legacy", out["hint"])
            nc.ack_legacy(t.session, "jdoe", "a 2.9 analysis; its quantities were reviewed")
            rc, out, p = t.deliver("--apply")
            self.assertEqual(rc, 0, p.stdout + p.stderr)
            self.assertEqual(out["normalization_gate"]["status"], "legacy_acknowledged")
            shipped = [os.path.join(dp, f) for dp, _, fs in os.walk(t.delivery) for f in fs]
            self.assertTrue(shipped)
            self.assertFalse([f for f in shipped if staff.is_staff_only(f)], shipped)
        finally:
            t.tearDown()


class NamesStayWithStaff(StaffEnv):
    def test_deliver_never_plans_a_staff_record(self):
        import core_submission as cs
        out = os.path.join(self.tmp, "output")
        os.makedirs(os.path.join(out, "reproducibility", "inputs"))
        os.makedirs(os.path.join(out, "tables"))
        for f in ("tables/normalization_check.staff.json",
                  "reproducibility/inputs/params.cfg.staff.json", "tables/DE.csv"):
            with open(os.path.join(out, f), "w") as fh:
                fh.write("x")
        plan = {i["name"]: i for i in cs.plan_delivery(out)}
        for f in ("tables/normalization_check.staff.json",
                  "reproducibility/inputs/params.cfg.staff.json"):
            self.assertEqual((plan[f]["status"], plan[f]["reason"]),
                             ("SKIPPED", cs.STAFF_ONLY_REASON), f)
        self.assertEqual(plan["tables/DE.csv"]["status"], "OK")

    def test_the_session_zip_leaves_staff_records_out_and_redacts_commands_log(self):
        env = dict(os.environ, RECORD_RUN="off", SKILL_CONFIG_DIR=os.path.join(self.tmp, "cfg"),
                   CLAUDE_CONFIG_DIR=os.path.join(self.tmp, "cc"))
        env.pop("CLAUDE_CODE_SESSION_ID", None)
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "init", "--name",
                            "Staff zip", "--base", os.path.join(self.tmp, "base")],
                           capture_output=True, text=True, env=env)
        self.assertEqual(r.returncode, 0, r.stderr)
        sdir = json.loads(r.stdout)["paths"]["session_dir"]
        os.makedirs(os.path.join(sdir, "input", "wf"), exist_ok=True)
        for rel in ("logs/" + nc.STAFF_RECORD, "input/wf/params.cfg.staff.json"):
            with open(os.path.join(sdir, rel), "w") as fh:
                json.dump({"by": "jdoe"}, fh)
        with open(os.path.join(sdir, "logs", "commands.log"), "a") as fh:
            fh.write("python3 normalization_check.py decide --check c.json --quantities raw "
                     "--by jdoe --reason \"low IgG\"\n")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "finalize",
                            "--dir", sdir, "--zip", "--no-deposit", "--no-notify"],
                           capture_output=True, text=True, env=env)
        self.assertEqual(r.returncode, 0, r.stderr)
        z = zipfile.ZipFile(sdir + ".zip")
        self.assertFalse([n for n in z.namelist() if staff.is_staff_only(n)], z.namelist())
        log = [n for n in z.namelist() if n.endswith("logs/commands.log")]
        self.assertEqual(len(log), 1, z.namelist())
        text = z.read(log[0]).decode()
        self.assertIn('--by "<Core staff>" --reason "low IgG"', text)
        self.assertNotIn("jdoe", text)
        excluded = json.loads(r.stdout)["zip_excluded"]
        self.assertEqual(
            excluded["staff-only records (*.staff.json: who decided, by name -- kept on disk)"], 2)
        with open(os.path.join(sdir, "logs", "commands.log")) as fh:
            self.assertIn("--by jdoe", fh.read(), "the session's own log keeps the command as run")

    def test_an_override_names_the_role_and_the_staff_record_the_person(self):
        cfg = os.path.join(self.tmp, "wf", "params.cfg")
        os.makedirs(os.path.dirname(cfg))
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                            "--engine", "diann", "--acquisition", "DIA", "--instrument",
                            "timsTOF HT", "--precursor-mz-range", "299.5", "1200.5", "--out", cfg,
                            "--overrides", '{"--mass-acc": 25}', "--override-by", "Jane Doe",
                            "--override-reason", "off calibration"],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertNotIn("Jane Doe", p.stdout)
        with open(cfg + ".rationale.json") as fh:
            side = fh.read()
        self.assertNotIn("Jane Doe", side)
        said = json.loads(side)
        self.assertNotIn(ME, json.dumps([said["overrides"], said["rationale"]]))
        o = said["overrides"]["--mass-acc"]
        self.assertEqual((o["set_by"], o["staff_record"]), ("Core staff", "params.cfg.staff.json"))
        self.assertIn("set by Core staff; reason: off calibration", o["source"])
        e = staff.read(cfg + staff.STAFF_SUFFIX)
        self.assertEqual([(x["by"], x["by_how"], x["account"]) for x in e],
                         [("Jane Doe", "given", ME)])
        # the bundle the client gets: the sidecar by role, the commands log redacted, no staff file
        log = os.path.join(self.tmp, "commands.log")
        with open(log, "w") as fh:
            fh.write("python3 estimate_params.py --override-by \"Jane Doe\" --out wf/params.cfg\n")
        out = os.path.join(self.tmp, "repro")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"), "--outdir",
                            out, "--params", cfg, "--commands", log],
                           capture_output=True, text=True, timeout=300)
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        inputs = os.path.join(out, "inputs")
        self.assertEqual(sorted(os.listdir(inputs)),
                         ["commands.log", "params.cfg", "params.cfg.rationale.json"])
        for f in os.listdir(inputs):
            with open(os.path.join(inputs, f)) as fh:
                self.assertNotIn("Jane Doe", fh.read(), f)


class MaxLFQRawHasANextStep(StaffEnv):
    def test_the_refusal_names_dpc_and_the_no_norm_route(self):
        s = session(self.tmp, "ip")
        report, cond = (os.path.join(s, "output", "search", "report.parquet"),
                        os.path.join(s, "input", "conditions.csv"))
        msg = nc.required(report, cond, "maxlfq")["message"]
        for words in ("--method dpc", "--no-norm", "no_norm_report.parquet", "--report-raw",
                      "SKILL.md step 8"):
            self.assertIn(words, msg)
        self.assertNotIn("--report-raw", nc.required(report, cond, "dpc")["message"])
        cmds = nc.step8_commands(report, cond, s, "maxlfq", report_raw="/x/no_norm_report.parquet")
        raw_de = [c for c in cmds["check"] if "--quantities raw" in c][0]
        self.assertIn("--input /x/no_norm_report.parquet", raw_de)
        self.assertIn(f"--input {report}", cmds["check"][0])
        self.assertIn("--report-raw /x/no_norm_report.parquet", cmds["check"][-1])

    def test_a_raw_decision_holds_the_final_de_to_the_no_norm_report(self):
        a, b = os.path.join(self.tmp, "report.parquet"), os.path.join(self.tmp, "no_norm.parquet")
        for f, x in ((a, "normal"), (b, "no-norm")):
            with open(f, "w") as fh:
                fh.write(x)
        path = check_record(os.path.join(self.tmp, "nc", nc.RECORD), report=a,
                            report_sha256=nc.sha256(a), report_raw=b,
                            report_raw_sha256=nc.sha256(b))
        nc.decide(path, et.RAW, "jdoe", REASON)
        self.assertTrue(nc.gate(path, et.RAW, report=b)[0])
        ok, msg, _ = nc.gate(path, et.RAW, report=a)
        self.assertFalse(ok)
        self.assertIn("--report-raw", msg)
        nc.decide(path, et.NORMALISED, "jdoe", "whole lysate after all")
        self.assertTrue(nc.gate(path, et.NORMALISED, report=a)[0])
        self.assertFalse(nc.gate(path, et.NORMALISED, report=b)[0])


if __name__ == "__main__":
    unittest.main()
