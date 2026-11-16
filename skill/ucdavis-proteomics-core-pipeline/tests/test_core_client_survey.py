"""
The Core's feedback survey on Core client work that has no CoreOmics submission yet (2.15).

Before: submission_report.core_run -- the one test for the survey -- said yes only for a session
that answers a CoreOmics submission. Core staff doing client work before a submission existed (a
pilot, a favour, a rush job) could not say "this is a Core client delivery", and nothing told them
the survey was off: an external pilot went out without it, and staff noticed after delivery.

Now:
  * `session.py core-client --session <S> --why-file -` marks such a session as Core client work.
    Only Core staff may: staff.core_staff -- the one test of Core staff, the admin-owned staff
    list (from a laptop, HIVE asked) -- must say so, else it is refused and nothing is written.
    session.json records the mark, by role and time; the login goes to the staff-only record and
    the decisions log has the reason. A mark counts only with its staff record.
  * core_run then answers yes, with no submission number: the report, README and AGENTS.md carry
    the survey, its link without `prot`.
  * a Core staff member whose session's survey is off is told so in one line at the report step
    (make_analysis_html.py) and at finalize / status, with the two ways to add it. A user outside
    the Core never sees that line, and a session whose records say nothing of Core work costs no
    staff check at all.
  * attaching the submission later and rebuilding brings the survey (with its number).

The staff check is injected by patching staff.core_staff -- in the process, so the CLIs whose
output depends on it run in-process here (2.15 review M1: no environment switch for it). Every
command run in a subprocess gets an environment that cannot reach HIVE. Fake records only
(PROT_9756 is the skill's fixture). Stdlib only, no network.
"""
import contextlib
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, HERE)
sys.path.insert(0, SCRIPTS)

import make_analysis_html as mah  # noqa: E402
import session as session_mod  # noqa: E402
import session_docs  # noqa: E402
import staff  # noqa: E402
import submission_report as sr  # noqa: E402
from test_submission_report import FIXTURE, attach_fixture, make_session  # noqa: E402

SURVEY = "https://feedback-ucd-proteomics.azurewebsites.net/"
STAFF = ("staffexample", "on the Core staff list (a test's)", None)
NOT_STAFF = (None, None, "staffexample is not on the Core staff list (a test's)")
UNASKABLE = (None, None, "the Core staff list is on HIVE and could not be asked from here (a "
                         "test's)")
WHY = "a pilot for an external lab before its submission exists"


def no_hive_env():
    """os.environ for a subprocess that must never reach HIVE: no HIVE login, no hive.env."""
    env = {k: v for k, v in os.environ.items() if not k.startswith("HIVE_")}
    env["HIVE_ENV_FILE"] = "/nonexistent/hive.env"
    return env


def run(script, *args, stdin=None):
    return subprocess.run([sys.executable, os.path.join(SCRIPTS, script)] + list(args),
                          capture_output=True, text=True, timeout=180, env=no_hive_env(),
                          input=stdin)


def cli(module, argv, answer, stdin=None):
    """`module`.main() with `argv`, in this process, staff.core_staff answering `answer`.
    -> (exit code, stdout, stderr)."""
    out, err, code = io.StringIO(), io.StringIO(), 0
    calls = []

    def asked():
        calls.append(1)
        return answer
    with contextlib.ExitStack() as st:
        st.enter_context(mock.patch.object(staff, "core_staff", asked))
        st.enter_context(mock.patch.object(sys, "argv", [module.__file__] + list(argv)))
        if stdin is not None:
            st.enter_context(mock.patch.object(sys, "stdin", io.TextIOWrapper(
                io.BytesIO(stdin.encode("utf-8")), encoding="utf-8")))
        st.enter_context(contextlib.redirect_stdout(out))
        st.enter_context(contextlib.redirect_stderr(err))
        try:
            module.main()
        except SystemExit as e:
            if isinstance(e.code, str):
                err.write(e.code)
            code = e.code if isinstance(e.code, int) else (0 if e.code is None else 1)
    cli.calls = len(calls)
    return code, out.getvalue(), err.getvalue()


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def with_analysis(tmp, internal="Example lab pilot"):
    """A session with an analysis to report. `internal`: registered as internal Core work (init
    --internal), as the issue's external pilot was -- the records of a Core staff session."""
    s = make_session(tmp)
    with open(os.path.join(s, "output", "AI_Analysis_Report.md"), "w") as fh:
        fh.write("# Analysis\n\nA standfirst.\n\n## Overview\n\nText.\n")
    if internal:
        with open(session_mod.paths_for(s)["session_json"], "w") as fh:
            json.dump({"analysis_for": {"internal": internal}}, fh)
    # these tests are about the survey line, so the session keeps its commands log as golden rule 5
    # says (written last, so it is not older than any step record) and no log warning joins the notes
    log = session_mod.paths_for(s)["commands_log"]
    os.makedirs(os.path.dirname(log), exist_ok=True)
    with open(log, "w") as fh:
        fh.write("# a test's command\n")
    return s


def log_command(s, line="# a test's command"):
    """Append to the session's commands log, as golden rule 5 says each command is, so a
    step record written by the command just run is not newer than the log."""
    with open(session_mod.paths_for(s)["commands_log"], "a") as fh:
        fh.write(line + "\n")


def report(s, answer=NOT_STAFF):
    out = os.path.join(s, "output", "Analysis_Report.html")
    r = cli(mah, ["--session", s, "--out", out, "--no-pdf"], answer)
    log_command(s, "python3 make_analysis_html.py --session ... --no-pdf")   # golden rule 5
    return r, read(out), read(os.path.splitext(out)[0] + ".md")


def mark(s, answer=STAFF, why=WHY):
    return cli(session_mod, ["core-client", "--session", s, "--why-file", "-"], answer,
               stdin=why + "\n")


class TestMark(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.addCleanup(self._td.cleanup)
        self.tmp = self._td.name

    def unmarked_line(self):
        with tempfile.TemporaryDirectory() as t2, \
                mock.patch.object(staff, "core_staff", lambda: STAFF):
            return sr.survey_off_note(with_analysis(t2))

    def test_core_staff_mark_the_session_and_it_is_recorded(self):
        s = with_analysis(self.tmp)
        p = session_mod.paths_for(s)
        self.assertIsNone(sr.core_run(s))
        rc, out, err = mark(s)
        self.assertEqual(rc, 0, out + err)
        self.assertEqual(cli.calls, 1)                          # one staff check (review M6)
        j = json.loads(out)
        self.assertEqual(j["survey"], "on")
        block = json.loads(read(p["session_json"]))[sr.CORE_CLIENT]
        self.assertEqual(block, j["core_client"])
        # who and when: the role in session.json (it is zipped); the login in the staff record
        self.assertEqual(block["by"], "Core staff")
        self.assertRegex(block["marked"], r"^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ$")
        self.assertEqual(block["staff_check"], STAFF[1])
        self.assertEqual(block["recorded_in"], ["logs/decisions.md", "logs/core_client.staff.json"])
        entries = staff.read(os.path.join(s, "logs", "core_client.staff.json"))
        self.assertEqual([(e["event"], e["by"], e["at"]) for e in entries],
                         [("core-client", "staffexample", block["marked"])])
        log = read(os.path.join(s, "logs", "decisions.md"))
        self.assertIn("Marked this session as Core client work with no CoreOmics submission", log)
        self.assertIn(WHY, log)
        self.assertNotIn("staffexample", json.dumps(block))
        self.assertEqual(sr.core_run(s), {"prot": None, "core_client": True})

    def test_anyone_not_on_the_core_staff_list_is_refused(self):
        s = with_analysis(self.tmp)
        p = session_mod.paths_for(s)
        before = read(p["session_json"])
        for answer in (NOT_STAFF, UNASKABLE):
            with self.subTest(answer=answer[2]):
                rc, out, err = mark(s, answer)
                self.assertEqual(rc, 2, out + err)
                self.assertEqual(json.loads(out)["error"],
                                 f"session.py core-client is for UC Davis Proteomics Core staff: "
                                 f"{answer[2]}. Nothing was written.")
                self.assertEqual(read(p["session_json"]), before)
                self.assertFalse(os.path.exists(os.path.join(s, "logs", "decisions.md")))
                self.assertFalse(os.path.exists(os.path.join(s, "logs", "core_client.staff.json")))
                self.assertIsNone(sr.core_run(s))

    def test_the_reason_is_required_and_comes_from_a_file_or_standard_input(self):
        s = with_analysis(self.tmp)
        rc, out, err = cli(session_mod, ["core-client", "--session", s], STAFF)
        self.assertNotEqual(rc, 0)
        self.assertNotIn(sr.CORE_CLIENT, read(session_mod.paths_for(s)["session_json"]))
        rc, out, err = mark(s, why="a rush job; the lab's form is $(not run) `coming`")
        self.assertEqual(rc, 0, out + err)
        self.assertIn("the lab's form is $(not run) `coming`", read(os.path.join(s, "logs",
                                                                              "decisions.md")))

    def test_a_session_that_answers_a_submission_needs_no_mark(self):
        s = with_analysis(self.tmp, internal=None)
        attach_fixture(s)
        rc, out, err = mark(s)
        self.assertEqual(rc, 0, out + err)
        j = json.loads(out)
        self.assertIsNone(j["core_client"])
        self.assertIn("PROT_9756", j["note"])
        self.assertNotIn(sr.CORE_CLIENT, json.loads(read(session_mod.paths_for(s)["session_json"])))

    def test_internal_work_keeps_its_board_title_when_marked(self):
        # the external pilot of the issue was registered as internal work (init --internal), the
        # only way staff had to put it on the board
        s = with_analysis(self.tmp)
        rc, out, err = mark(s)
        self.assertEqual(rc, 0, out + err)
        sj = json.loads(read(session_mod.paths_for(s)["session_json"]))
        self.assertEqual(sj["analysis_for"], {"internal": "Example lab pilot"})
        self.assertIn(sr.CORE_CLIENT, sj)
        self.assertIn("internal", json.loads(out)["note"])

    def test_clear_takes_the_mark_away_and_says_so(self):
        s = with_analysis(self.tmp)
        self.assertEqual(mark(s)[0], 0)
        rc, out, err = cli(session_mod, ["core-client", "--session", s, "--clear", "--why-file",
                                         "-"], STAFF, stdin="marked by mistake: internal QC\n")
        self.assertEqual(rc, 0, out + err)
        self.assertEqual(json.loads(out)["survey"], "off")
        self.assertNotIn(sr.CORE_CLIENT, json.loads(read(session_mod.paths_for(s)["session_json"])))
        self.assertIn("Removed the Core client work mark", read(os.path.join(s, "logs",
                                                                           "decisions.md")))
        self.assertIsNone(sr.core_run(s))

    def test_the_mark_fails_closed(self):
        s = with_analysis(self.tmp)
        p = session_mod.paths_for(s)
        self.assertEqual(mark(s)[0], 0)
        good = json.loads(read(p["session_json"]))[sr.CORE_CLIENT]
        self.assertEqual(sr.core_run(s), {"prot": None, "core_client": True})
        for block in ({}, "yes", True, [good], dict(good, marked=""),
                      dict(good, marked="2026-01-01T00:00:00Z"),       # not the staff record's
                      {k: v for k, v in good.items() if k != "marked"}):
            with self.subTest(block=block):
                with open(p["session_json"], "w") as fh:
                    json.dump({sr.CORE_CLIENT: block}, fh)
                self.assertIsNone(sr.core_run(s))
        for damaged in ("{", "[1]", ""):
            with self.subTest(damaged=damaged):
                with open(p["session_json"], "w") as fh:
                    fh.write(damaged)
                self.assertIsNone(sr.core_run(s))

    def test_a_mark_written_by_hand_is_ignored_and_said(self):
        # review M1: anyone could type {"core_client": {"marked": "..."}} into session.json
        s = with_analysis(self.tmp)
        p = session_mod.paths_for(s)
        with open(p["session_json"], "w") as fh:
            json.dump({"analysis_for": {"internal": "x"},
                       sr.CORE_CLIENT: {"marked": "2026-11-04T18:00:00Z", "by": "Core staff",
                                        "core_member": True}}, fh)
        self.assertIsNone(sr.core_run(s))
        with mock.patch.object(staff, "core_staff", lambda: STAFF):
            line = sr.survey_off_note(s)
        self.assertIn("The Core's feedback survey will NOT be included: session.json's Core client "
                      "work mark is ignored (it has no matching entry in logs/core_client.staff.json,"
                      " and session.py core-client writes both), so its report, README and "
                      "AGENTS.md go out without it. If it was marked on HIVE, pull "
                      "logs/core_client.staff.json with session.json and rebuild the report. "
                      "Otherwise, for work for a Core client, attach the submission", line)
        self.assertNotIn("If it was marked on HIVE", self.unmarked_line())
        # a staff record that ends with --clear: the mark it cleared no longer counts either
        self.assertEqual(mark(s)[0], 0)
        good = json.loads(read(p["session_json"]))[sr.CORE_CLIENT]
        staff.append(os.path.join(s, "logs", "core_client.staff.json"),
                     {"event": "core-client --clear", "by": "staffexample", "at": "later"})
        with open(p["session_json"], "w") as fh:
            json.dump({sr.CORE_CLIENT: good}, fh)
        self.assertIsNone(sr.core_run(s))


class TestMarkGuards(unittest.TestCase):
    """review M5 (a named submission is attached, not marked), L1 (all three records or none)."""

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.addCleanup(self._td.cleanup)
        self.tmp = self._td.name

    def test_a_session_that_names_a_submission_is_refused(self):
        s = with_analysis(self.tmp, internal=None)
        p = session_mod.paths_for(s)
        with open(p["session_json"], "w") as fh:
            json.dump({"analysis_for": {"submissions": ["PROT_9756"]}}, fh)
        rc, out, err = mark(s)
        self.assertEqual(rc, 2, out + err)
        self.assertIn("session.json says this session is for PROT_9756 (session.py init "
                      "--submission): attach that submission instead", json.loads(out)["error"])
        self.assertNotIn(sr.CORE_CLIENT, read(p["session_json"]))
        self.assertFalse(os.path.exists(os.path.join(s, "logs", "core_client.staff.json")))

    def test_a_damaged_staff_record_writes_nothing(self):
        s = with_analysis(self.tmp)
        p = session_mod.paths_for(s)
        os.makedirs(os.path.join(s, "logs"), exist_ok=True)
        bad = os.path.join(s, "logs", "core_client.staff.json")
        with open(bad, "w") as fh:
            fh.write("{")
        before = read(p["session_json"])
        rc, out, err = mark(s)
        self.assertEqual(rc, 2, out + err)
        self.assertIn("core-client could not write logs/core_client.staff.json (JSONDecodeError",
                      json.loads(out)["error"])
        self.assertTrue(json.loads(out)["error"].endswith("; nothing was changed"))
        self.assertEqual((read(p["session_json"]), read(bad)), (before, "{"))
        self.assertFalse(os.path.exists(os.path.join(s, "logs", "decisions.md")))

    def test_a_failure_after_the_first_record_undoes_it(self):
        import log_decision
        s = with_analysis(self.tmp)
        p = session_mod.paths_for(s)
        before = read(p["session_json"])
        with mock.patch.object(log_decision, "append", side_effect=OSError(28, "No space left")):
            rc, out, err = mark(s)
        self.assertEqual(rc, 2, out + err)
        self.assertIn("core-client could not write logs/decisions.md (OSError", json.loads(out)["error"])
        self.assertTrue(json.loads(out)["error"].endswith("; nothing was changed"))
        self.assertEqual(read(p["session_json"]), before)
        self.assertFalse(os.path.exists(os.path.join(s, "logs", "core_client.staff.json")))
        self.assertIsNone(sr.core_run(s))


class TestMarkGuardsRound2(unittest.TestCase):
    """review2 L1: any failure undoes all three records; an unreadable decisions.md refuses with
    JSON (exit 2), nothing written."""

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.addCleanup(self._td.cleanup)
        self.tmp = self._td.name

    def test_any_exception_undoes_the_records(self):
        import log_decision
        s = with_analysis(self.tmp)
        before = read(session_mod.paths_for(s)["session_json"])
        with mock.patch.object(log_decision, "append", side_effect=RuntimeError("boom")):
            rc, out, err = mark(s)
        self.assertEqual(rc, 2, out + err)
        self.assertIn("core-client could not write logs/decisions.md (RuntimeError: boom); "
                      "nothing was changed", json.loads(out)["error"])
        self.assertEqual(read(session_mod.paths_for(s)["session_json"]), before)
        self.assertFalse(os.path.exists(os.path.join(s, "logs", "core_client.staff.json")))
        self.assertIsNone(sr.core_run(s))

    def test_an_unreadable_decisions_log_is_refused_with_json(self):
        s = with_analysis(self.tmp)
        os.makedirs(os.path.join(s, "logs", "decisions.md"))         # a folder: cannot be read
        before = read(session_mod.paths_for(s)["session_json"])
        rc, out, err = mark(s)
        self.assertEqual(rc, 2, out + err)
        self.assertIn("core-client cannot read logs/decisions.md", json.loads(out)["error"])
        self.assertEqual(read(session_mod.paths_for(s)["session_json"]), before)
        self.assertFalse(os.path.exists(os.path.join(s, "logs", "core_client.staff.json")))


class TestTheReasonNeverTravels(unittest.TestCase):
    """review H3: the reason (who the work is for) reached commands.log, and from there the
    delivered bundle, the zip and the run registry. staff.redact, the one redactor, takes it out."""

    LOG = ("python3 scripts/log_decision.py --session S --what x --why 'kept: not core-client'\n"
           "python3 scripts/session.py core-client --session S --why \"Dr Example's pilot\"\n"
           "bash scripts/hive_exec.sh \"python3 ~/x/session.py core-client --session S "
           "--why-file -\" <<'EOF'\n"
           "the Example lab's pilot, before its submission\n"
           "second line of it\n"
           "EOF\n"
           "bash scripts/hive_exec.sh \"python3 ~/x/session.py core-client --session S --why "
           "\\\"Example Lab, typed inline\\\"\"\n"
           "bash scripts/hive_exec.sh \"python3 ~/x/qc_bracket.py ack --by \\\"Jane Example\\\"\"\n"
           "python3 scripts/staff.py check --by staffexample\n")

    def test_staff_redact_takes_the_reason_out(self):
        got = staff.redact(self.LOG)
        self.assertNotIn("Example", got)
        self.assertNotIn("second line of it", got)
        self.assertIn(f'--why {staff.REDACTED_REASON_ARG}\n', got)
        # --why-file - (standard input) is kept: the heredoc's body is the reason
        self.assertIn('--why-file -" <<\'EOF\'', got)
        self.assertIn(f"\n{staff.REDACTED_REASON}\nEOF\n", got)
        # a value in \"...\" inside a double-quoted command (hive_exec.sh): all of it
        self.assertIn(f'--why {staff.REDACTED_REASON_ARG}"\n', got)
        self.assertIn('--by "<Core staff>""\n', got)
        self.assertIn("--why 'kept: not core-client'", got)      # another command's: as it was
        self.assertIn('--by "<Core staff>"\n', got)              # the names, as before

    def test_the_run_registry_copy_is_redacted(self):
        import record_run
        with tempfile.TemporaryDirectory() as tmp:
            src = os.path.join(tmp, "commands.log")
            with open(src, "w") as fh:
                fh.write(self.LOG)
            dst = os.path.join(tmp, "registry", "scripts", "commands.log")
            record_run.copy_file(src, dst, record_run.Deadline(60))
            self.assertNotIn("Example", read(dst))
            other = os.path.join(tmp, "registry", "scripts", "notes.txt")   # any other file: bytes
            record_run.copy_file(src, other, record_run.Deadline(60))
            self.assertEqual(read(other), self.LOG)

    def test_the_bundle_copy_is_redacted(self):
        # provenance.py copies commands.log into the delivered bundle's inputs/
        with tempfile.TemporaryDirectory() as tmp:
            log, out = os.path.join(tmp, "commands.log"), os.path.join(tmp, "repro")
            with open(log, "w") as fh:
                fh.write(self.LOG)
            r = run("provenance.py", "--outdir", out, "--commands", log)
            self.assertEqual(r.returncode, 0, r.stderr[-1500:])
            got = read(os.path.join(out, "inputs", "commands.log"))
            self.assertNotIn("Example", got)
            self.assertIn(f"\n{staff.REDACTED_REASON}\nEOF\n", got)

    def test_the_session_zip_copy_is_redacted(self):
        # session.py finalize zips logs/commands.log (the one log it zips)
        import finish_fixture
        import test_deposit_package as tdp
        import zipfile
        from job_env import job_env
        with tempfile.TemporaryDirectory() as tmp:
            p = tdp.dia_session(tmp)
            s = p["session_dir"]
            finish_fixture.finished(p["output_dir"])
            os.makedirs(os.path.dirname(p["commands_log"]), exist_ok=True)
            with open(p["commands_log"], "w") as fh:
                fh.write(self.LOG)
            with mock.patch.dict(os.environ, job_env(os.path.dirname(s)), clear=True):
                rc, out, err = cli(session_mod, ["finalize", "--dir", s, "--zip"], NOT_STAFF)
            self.assertEqual(rc, 0, err[-2000:])
            with zipfile.ZipFile(json.loads(out)["zip"]) as z:
                logs = [n for n in z.namelist() if n.endswith("logs/commands.log")]
                self.assertEqual(len(logs), 1, z.namelist())
                got = z.read(logs[0]).decode("utf-8")
            self.assertNotIn("Example", got)
            self.assertIn(f"\n{staff.REDACTED_REASON}\nEOF\n", got)


class TestPodcastConsent(unittest.TestCase):
    """Brett's ruling 92 (2026-11-04): marked Core client work is a Core run for the podcast's
    standing approval; never internal work."""

    def test_marked_work_is_a_core_run_never_internal(self):
        import make_podcast as mp
        with tempfile.TemporaryDirectory() as tmp:
            s = with_analysis(tmp)
            self.assertEqual(mark(s)[0], 0)
            script = os.path.join(s, "output", "podcast", "podcast_script.md")
            os.makedirs(os.path.dirname(script))
            with mock.patch.object(staff, "core_staff", lambda: STAFF):
                ok, rec, why = mp.core_run_consent(script, "submission")
                self.assertEqual((ok, why), (mp.CORE_CLOUD_OK, None))
                self.assertEqual(rec["kind"], "core_client")
                self.assertEqual(rec["by"], staff.ROLE)
                ok, rec, why = mp.core_run_consent(script, "internal")
                self.assertIsNone(ok)
                self.assertIn("is Core client work (session.py core-client), a client's data",
                              why)
            with mock.patch.object(staff, "core_staff", lambda: NOT_STAFF):
                ok, rec, why = mp.core_run_consent(script, "submission")
                self.assertIsNone(ok)
                self.assertIn("--core-run is for Core staff", why)


class TestSurveyOnMarkedWork(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.addCleanup(self._td.cleanup)
        self.tmp = self._td.name

    def test_the_report_readme_and_agents_carry_the_survey_without_a_number(self):
        s = with_analysis(self.tmp)
        self.assertEqual(mark(s)[0], 0)
        (rc, out, err), page, md = report(s)
        self.assertEqual(rc, 0, err)
        url = SURVEY + "?src=report"
        self.assertIn(f'<a href="{url}">a 5-minute survey</a></p></main>', page)
        self.assertTrue(md.endswith(f"[a 5-minute survey]({url}).\n"), md[-300:])
        self.assertNotIn("prot=", page)
        written = session_docs.write_docs(s)
        self.assertIn(f"[a 5-minute survey]({SURVEY}?src=readme)", read(written["readme"]))
        self.assertIn(f"**Feedback:** the UC Davis Proteomics Core's survey for this project is "
                      f"{SURVEY}?src=readme -- ", read(written["agents"]))

    def test_status_owes_no_coreomics_delivery_for_marked_work(self):
        s = with_analysis(self.tmp)
        self.assertEqual(mark(s)[0], 0)
        items = {i["item"]: i for i in session_mod.finish_status(s)}
        self.assertTrue(items["delivery"]["done"])
        self.assertEqual(items["delivery"]["state"], "not_needed")
        self.assertIn("Core client work with no CoreOmics submission", items["delivery"]["detail"])


class TestSurveyOffLine(unittest.TestCase):
    """One plain line for Core staff when the survey will be left out; never for anyone else."""

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.addCleanup(self._td.cleanup)
        self.tmp = self._td.name

    def line(self, s, step="report", answer=STAFF):
        with mock.patch.object(staff, "core_staff", lambda: answer):
            return sr.survey_off_note(s, step)

    def test_core_staff_are_told_at_the_report_and_status(self):
        s = with_analysis(self.tmp)
        want = self.line(s)
        q = s                                   # no characters to quote in the test path
        self.assertEqual(want, (
            "The Core's feedback survey will NOT be included: this session answers no CoreOmics "
            "submission, so its report, README and AGENTS.md go out without it. For work for a "
            "Core client, attach the submission (python3 scripts/submission_report.py attach "
            f"--session {q} --record <its fetch folder>) or mark the session as Core client work "
            f"(python3 scripts/session.py core-client --session {q} --why-file - <<'EOF', then "
            "the reason, then EOF), then rebuild the report."))
        (rc, out, err), page, md = report(s, STAFF)
        self.assertEqual(rc, 0, err)
        self.assertIn("[make_analysis_html] NOTE: " + want, err)
        self.assertEqual(cli.calls, 1)
        self.assertNotIn(SURVEY, page + md)
        then = self.line(s, "status")
        self.assertTrue(then.endswith(
            "then rebuild the report (make_analysis_html.py; a report built on a laptop needs "
            "session.json and logs/core_client.staff.json pulled again first) and run finalize "
            "again."), then)
        rc, out, err = cli(session_mod, ["status", "--session", s, "--json"], STAFF)
        self.assertEqual(json.loads(out)["notes"], [then])
        rc, out, err = cli(session_mod, ["status", "--session", s], STAFF)
        self.assertIn("  note     " + then, out)

    def test_the_line_says_why_from_the_records(self):
        s = with_analysis(self.tmp)
        p = session_mod.paths_for(s)
        with open(p["session_json"], "w") as fh:
            json.dump({"analysis_for": {"submissions": ["PROT_9756"]}}, fh)
        self.assertIn("session.json names PROT_9756 (session.py init --submission), but no "
                      "CoreOmics record is attached, so", self.line(s))
        attach_fixture(s)
        self.assertIsNone(self.line(s))                       # the survey is on: nothing to say
        with open(p["submission_record"], "w") as fh:
            fh.write("{")
        self.assertIn("its CoreOmics record (input/submission.json) cannot be read, so",
                      self.line(s))

    def test_finalize_says_it_and_to_run_finalize_again(self):
        import finish_fixture
        import test_deposit_package as tdp
        from job_env import job_env
        p = tdp.dia_session(self.tmp)
        s = p["session_dir"]
        with open(p["session_json"], "w") as fh:
            json.dump({"analysis_for": {"internal": "Example lab pilot"}}, fh)
        finish_fixture.finished(p["output_dir"])
        want = self.line(s, "finalize")
        for answer, said in ((STAFF, want), (NOT_STAFF, None)):
            with self.subTest(staff=bool(answer[0])), \
                    mock.patch.dict(os.environ, job_env(os.path.dirname(s)), clear=True):
                rc, out, err = cli(session_mod, ["finalize", "--dir", s], answer)
            self.assertEqual(rc, 0, err[-2000:])
            self.assertEqual(json.loads(out)["survey"], said)
            if said:
                self.assertIn("[session] NOTE: " + said, err)
                self.assertIn("and run finalize again.", said)
            else:
                self.assertNotIn("feedback survey", err)

    def test_a_user_outside_the_core_never_sees_it(self):
        s = with_analysis(self.tmp)
        for answer in (NOT_STAFF, UNASKABLE):
            with self.subTest(answer=answer[2]):
                self.assertIsNone(self.line(s, answer=answer))
                (rc, out, err), page, md = report(s, answer)
                self.assertEqual(rc, 0, err)
                for text in (out, err):
                    self.assertNotIn("feedback survey", text)
                    self.assertNotIn("core-client", text)
                rc, out, err = cli(session_mod, ["status", "--session", s, "--json"], answer)
                self.assertEqual(json.loads(out)["notes"], [])
                rc, out, err = cli(session_mod, ["status", "--session", s], answer)
                self.assertNotIn("feedback survey", out + err)

    def test_a_session_that_says_nothing_of_core_work_costs_no_check(self):
        # review M6: anyone else's report, finalize and status ask nothing of HIVE
        s = with_analysis(self.tmp, internal=None)
        (rc, out, err), page, md = report(s, STAFF)
        self.assertEqual((rc, cli.calls), (0, 0), err)
        rc, out, err = cli(session_mod, ["status", "--session", s, "--json"], STAFF)
        self.assertEqual((json.loads(out)["notes"], cli.calls), ([], 0))

    def test_a_marked_session_has_no_line(self):
        s = with_analysis(self.tmp)
        self.assertEqual(mark(s)[0], 0)
        self.assertIsNone(self.line(s))
        (rc, out, err), page, md = report(s, STAFF)
        self.assertNotIn("will NOT be included", err)
        self.assertEqual(cli.calls, 0)                    # core_run says yes: nothing to ask


class TestAttachLater(unittest.TestCase):
    """(c) The submission attached after a report was built: a rebuild carries the survey."""

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.addCleanup(self._td.cleanup)
        self.tmp = self._td.name

    def test_unmarked_then_attached(self):
        s = with_analysis(self.tmp, internal=None)
        (rc, out, err), page, md = report(s)
        self.assertNotIn(SURVEY, page + md)
        self.assertNotIn(SURVEY, read(session_docs.write_docs(s)["readme"]))
        r = run("submission_report.py", "attach", "--session", s, "--record", FIXTURE)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        (rc, out, err), page, md = report(s, STAFF)
        self.assertEqual(rc, 0, err)
        self.assertNotIn("will NOT be included", err)
        self.assertIn(f'href="{SURVEY}?prot=PROT_9756&amp;src=report"', page)
        self.assertIn(f"{SURVEY}?prot=PROT_9756&src=report", md)
        written = session_docs.write_docs(s)
        self.assertIn(f"{SURVEY}?prot=PROT_9756&src=readme", read(written["readme"]))
        self.assertIn(f"{SURVEY}?prot=PROT_9756&src=readme", read(written["agents"]))

    def test_marked_then_attached_gains_the_number(self):
        s = with_analysis(self.tmp)
        self.assertEqual(mark(s)[0], 0)
        (rc, out, err), page, md = report(s)
        self.assertIn(f"{SURVEY}?src=report", md)
        attach_fixture(s)
        self.assertEqual(sr.core_run(s), {"prot": "PROT_9756"})
        (rc, out, err), page, md = report(s)
        self.assertIn(f"{SURVEY}?prot=PROT_9756&src=report", md)
        self.assertNotIn(f"{SURVEY}?src=report", md)

    def test_a_report_built_before_the_mark_is_said_at_status_and_finalize(self):
        # review H4c: the report went out without the survey, silently; nothing re-adds it
        s = with_analysis(self.tmp)
        report(s)                                             # before the mark: no survey line
        self.assertEqual(mark(s)[0], 0)
        rc, out, err = cli(session_mod, ["status", "--session", s, "--json"], STAFF)
        self.assertEqual(json.loads(out)["notes"], [(
            "The report (output/Analysis_Report.html) does not carry the Core's feedback survey, "
            "although this session is Core client work (session.py core-client): it was built "
            "before that was recorded. Rebuild it (make_analysis_html.py; a report built on a "
            "laptop needs session.json, input/submission.json and logs/core_client.staff.json "
            "pulled again first), then run finalize again.")])
        report(s)                                             # rebuilt: it carries the survey
        rc, out, err = cli(session_mod, ["status", "--session", s, "--json"], STAFF)
        self.assertEqual(json.loads(out)["notes"], [])


class TestOneStaffTest(unittest.TestCase):
    """review M1 / H4a / M6: one definition of Core staff, asked once, never by an env switch."""

    def test_one_test_of_core_staff_everywhere(self):
        # review M1, H4a: core-client, the survey line, the podcast and the PRM follow-up all ask
        # staff.core_staff; from a laptop it asks HIVE once, through hive_exec.sh, stdin closed
        import make_podcast as mp
        import skyline_peaks as sp
        with mock.patch.object(staff, "core_staff", lambda: STAFF):
            self.assertEqual(mp.core_staff(), STAFF)
            self.assertEqual(sp.core_staff(), ("staffexample", None))
        with mock.patch.object(staff, "core_staff", lambda: NOT_STAFF):
            self.assertEqual(sp.core_staff(), (None, NOT_STAFF[2]))
        calls = []

        def fake_run(argv, **kw):
            calls.append((argv, kw))
            return subprocess.CompletedProcess(argv, 0, stderr="", stdout=(
                "a login banner\n" + json.dumps({"ok": True, "by": "staffexample",
                                                  "checked": "on the Core staff list (/x)"})))
        with mock.patch.object(staff, "_on_hive", lambda: False), \
                mock.patch.object(staff.subprocess, "run", fake_run):
            self.assertEqual(staff.core_staff(), ("staffexample",
                                                  "on the Core staff list (/x), asked on HIVE", None))
        (argv, kw), = calls
        self.assertEqual((os.path.basename(argv[1]), argv[2]), ("hive_exec.sh",
                                                                staff.STAFF_CHECK_ON_HIVE))
        self.assertIs(kw.get("stdin"), subprocess.DEVNULL)
        self.assertIn("$(id -un)", staff.STAFF_CHECK_ON_HIVE)     # HIVE's own name for the key

    def test_a_cleared_mark_leaves_no_survey_unflagged(self):
        # review2 L3: a report built while marked kept the survey after --clear, and finalize /
        # status said "will NOT be included" -- false for the report that goes out
        with tempfile.TemporaryDirectory() as tmp:
            s = with_analysis(tmp)
            self.assertEqual(mark(s)[0], 0)
            report(s, STAFF)                                     # carries the survey
            rc, out, err = cli(session_mod, ["core-client", "--session", s, "--clear",
                                             "--why-file", "-"], STAFF, stdin="by mistake\n")
            self.assertEqual(rc, 0, out + err)
            rc, out, err = cli(session_mod, ["status", "--session", s, "--json"], STAFF)
            self.assertEqual(json.loads(out)["notes"], [(
                "The report (output/Analysis_Report.html) carries the Core's feedback survey, but "
                "this session is no longer a Core run (it answers no CoreOmics submission and "
                "carries no Core client work mark -- the mark was cleared, or the record taken "
                "away, after the report was built). Rebuild it (make_analysis_html.py; a report "
                "built on a laptop needs session.json, input/submission.json and "
                "logs/core_client.staff.json pulled again first), then run finalize again.")])
            report(s, STAFF)                                     # rebuilt: no survey now
            rc, out, err = cli(session_mod, ["status", "--session", s, "--json"], STAFF)
            (only,) = json.loads(out)["notes"]
            self.assertTrue(only.startswith("The Core's feedback survey will NOT be included"))

    def test_no_environment_switch_turns_the_staff_check_on(self):
        # review M1: CORE_ACCESS_CHECK, a tests-only switch, was honoured in production
        self.assertFalse(hasattr(sr, "ACCESS_CHECK_ENV"))
        self.assertFalse(hasattr(sr, "core_member"))
        for name in ("submission_report.py", "session.py", "staff.py"):
            self.assertNotIn("CORE_ACCESS_CHECK", read(os.path.join(SCRIPTS, name)), name)


if __name__ == "__main__":
    unittest.main()
