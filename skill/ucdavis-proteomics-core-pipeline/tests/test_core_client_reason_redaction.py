#!/usr/bin/env python3
"""A core-client reason never travels, and the redactor never drops a command it cannot account
for (2.15 review2 NH1 + the H3 residual; ruling 97).

The reason staff give `session.py core-client` says who the work is for. Every command goes into
logs/commands.log, and copies of it travel: the reproducibility bundle's inputs/ (delivered), the
session zip, the run registry. staff.redact_log / redact_report:
  * remove the reason by VALUE (logs/core_client.staff.json keeps it as given), in any quoting;
  * on each core-client command (`\\`-continued lines joined) replace the --why / --why-file value
    and a heredoc body -- the body only when its end line is found;
  * withhold a command whole ("[core-client mark: reason withheld]") when its reason cannot be cut
    out exactly, and say so (a MANIFEST [NOTE], the registry's record) -- never silently, and never
    a line they cannot account for dropped.
Before, an opener whose heredoc end was not logged took every later line with it (NH1), and a
`\\` continuation, a pipe, a here-string, $'..', escaped spaces, a --why spanning lines or a
delimiter with a '-' let the reason through. Invented text only. Stdlib only, no network.
"""
import json
import os
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, HERE)
sys.path.insert(0, SCRIPTS)

import staff  # noqa: E402

REASON = "the Example lab pilot before its submission"
CC = "python3 scripts/session.py core-client --session S"
HX = 'bash scripts/hive_exec.sh "python3 ~/proteomics-pipeline/scripts/session.py core-client --session S'
NEXT = "python3 scripts/run_de.py --session S"

#: the leak forms review2 lists, each followed by an ordinary command that must survive as is
FORMS = {
    "a \\ continuation before --why-file": f"{CC} \\\n  --why-file - <<'EOF'\nthe Example lab pilot\nEOF",
    "a \\ continuation before the heredoc": f"{CC} --why-file - \\\n  <<'EOF'\nthe Example lab pilot\nEOF",
    "a \\ continuation before --why's value": f"{CC} --why \\\n  \"the Example lab pilot\"",
    "echo into it": f'echo "the Example lab pilot" | {CC} --why-file -',
    "printf into hive_exec.sh": f"printf '%s\\n' 'the Example lab pilot' | {HX} --why-file -\"",
    "a here-string": f'{CC} --why-file - <<< "the Example lab pilot"',
    "a $'..' value": f"{CC} --why $'the Example\\nlab pilot'",
    "escaped spaces": f"{CC} --why the\\ Example\\ lab",
    "quotes run together": f"{CC} --why 'the Example'\" lab\"",
    "a --why spanning lines": f'{CC} --why "the Example lab\npilot before its submission"',
    "a quoted delimiter with a -": f"{CC} --why-file - <<'REASON-END'\nthe Example lab pilot\nREASON-END",
    "a bare delimiter with a -": f"{CC} --why-file - <<REASON-END\nthe Example lab pilot\nREASON-END",
    "--why= inside hive_exec.sh": f'{HX} --why=\\"the Example lab\\""',
}


def redact(text, reasons=()):
    return staff.redact_report(text, reasons)


class EveryLeakFormIsClosed(unittest.TestCase):
    def test_nh1_an_opener_with_no_end_keeps_every_other_command(self):
        log = ["python3 scripts/session.py init --session S --internal 'x'",
               f"{CC} --why-file - <<'EOF'",
               "python3 scripts/run_search.py --session S",
               "sbatch job.sh",
               NEXT,
               "python3 scripts/make_analysis_html.py --session S",
               "python3 scripts/session.py finalize --dir S --zip"]
        for reasons in ((), (REASON,)):
            with self.subTest(reasons=bool(reasons)):
                got, notes = redact("\n".join(log), reasons)
                lines = got.split("\n")
                self.assertEqual(len(lines), 7, got)
                self.assertEqual(lines[1], staff.WITHHELD_COMMAND)
                self.assertEqual(lines[:1] + lines[2:], log[:1] + log[2:])
                self.assertEqual(len(notes), 1)
                self.assertIn("heredoc has no end line", notes[0])

    def test_each_form_review2_lists(self):
        for name, cmd in FORMS.items():
            for reasons in ((), (REASON,)):
                with self.subTest(form=name, reasons=bool(reasons)):
                    got, notes = redact(cmd + "\n" + NEXT, reasons)
                    self.assertNotIn("Example", got)
                    self.assertTrue(got.endswith("\n" + NEXT), got)        # the next one kept
                    self.assertNotIn("\n\n", got)

    def test_what_is_withheld_whole_is_said(self):
        for name in ("echo into it", "printf into hive_exec.sh", "a here-string"):
            with self.subTest(form=name):
                got, notes = redact(FORMS[name] + "\n" + NEXT)
                self.assertEqual(got, staff.WITHHELD_COMMAND + "\n" + NEXT)
                self.assertEqual(len(notes), 1)

    def test_the_documented_forms_stay_readable(self):
        got, notes = redact(f"{HX} --why-file -\" <<'EOF'\n{REASON}\nsecond line\nEOF\n{NEXT}")
        self.assertEqual(got, f"{HX} --why-file -\" <<'EOF'\n{staff.REDACTED_REASON}\nEOF\n{NEXT}")
        self.assertEqual(notes, [])
        got, notes = redact(f'{CC} --why "{REASON}"')
        self.assertEqual(got, f"{CC} --why {staff.REDACTED_REASON_ARG}")
        got, notes = redact(f'{HX} --why \\"{REASON}\\""')                    # hive_exec \"..\"
        self.assertEqual(got, f'{HX} --why {staff.REDACTED_REASON_ARG}"')
        self.assertEqual(notes, [])

    def test_a_heredoc_end_is_never_borrowed_from_a_later_command(self):
        log = (f"{CC} --why-file - <<'EOF'\n{NEXT}\n"
               "python3 scripts/other.py --x <<'EOF'\nits own text\nEOF\nlast")
        got, notes = redact(log)
        self.assertEqual(got.split("\n")[1:], log.split("\n")[1:])
        self.assertEqual(got.split("\n")[0], staff.WITHHELD_COMMAND)

    def test_the_single_quoted_by_form_keeps_no_login(self):
        # review2 L5: --by '\''login'\'' inside a single-quoted hive_exec.sh command
        got = staff.redact("bash scripts/hive_exec.sh 'python3 ~/p/qc_bracket.py ack --by "
                           "'\\''staffexample'\\'' --session S'\n" + NEXT)
        self.assertNotIn("staffexample", got)
        self.assertIn('--by "<Core staff>" --session S\'', got)


class TheReasonIsTakenOutByValue(unittest.TestCase):
    def test_anywhere_in_the_log_in_any_quoting(self):
        r = "Dr O'Neil's lab pilot, before the submission"
        log = "\n".join([
            f"printf '%s\\n' '{r}' > /tmp/w.txt",                     # as typed: plain
            "printf '%s\\n' 'Dr O'\\''Neil'\\''s lab pilot, before the submission' > /tmp/w",
            'echo "Dr O\'Neil\'s lab pilot, before the submission" >> /tmp/w',
            "echo Dr\\ O\\'Neil\\'s\\ lab\\ pilot,\\ before\\ the\\ submission",
            f"{CC} --why-file /tmp/w.txt", NEXT])
        got, notes = redact(log, (r,))
        self.assertNotIn("Neil", got)
        self.assertEqual(len(got.split("\n")), 6)
        self.assertTrue(got.endswith(NEXT))
        got, _ = redact(log)                                     # without the record: the form
        self.assertIn("Neil", got.split("\n")[0])                # alone cannot know a printf's text

    def test_each_line_of_a_multi_line_reason(self):
        r = "for the Example lab\nbefore their CoreOmics form is in"
        log = f"{CC} --why-file - <<'EOF'\nfor the Example lab\n"      # the end never logged
        got, notes = redact(log + NEXT, (r,))
        self.assertNotIn("Example", got)
        self.assertTrue(got.endswith(NEXT))

    def test_a_short_reason_is_not_replaced_everywhere(self):
        got, _ = redact(f"python3 scripts/run.py --name pilot\n{CC} --why pilot", ("pilot",))
        self.assertEqual(got.split("\n")[0], "python3 scripts/run.py --name pilot")
        self.assertEqual(got.split("\n")[1], f"{CC} --why {staff.REDACTED_REASON_ARG}")

    def test_redact_log_reads_the_staff_record_beside_the_log(self):
        with tempfile.TemporaryDirectory() as tmp:
            logs = os.path.join(tmp, "logs")
            os.makedirs(logs)
            path = os.path.join(logs, "commands.log")
            with open(path, "w") as fh:
                fh.write(f"printf '%s' '{REASON}' > /tmp/w\n{NEXT}\n")
            rec = os.path.join(logs, staff.CORE_CLIENT_RECORD)
            staff.append(rec, {"event": "core-client", "by": "staffexample", "at": "t",
                               "why": REASON, "why_text": REASON})
            text, notes = staff.redact_log(path)
            self.assertNotIn("Example", text)
            self.assertEqual(notes, ["a reason recorded for a core-client mark was found in 1 "
                                     "other line(s) and replaced there by [reason-redacted]"])
            with open(rec, "w") as fh:
                fh.write("{")
            text, notes = staff.redact_log(path)
            self.assertIn("could not be read", notes[0])         # said, not swallowed


class Round3(unittest.TestCase):
    """ruling 98: R1 (a quoted value is never joined across lines), R2 (by value outside a
    core-client command: bounded words, 2+ words and 12+ characters, and said), R3 (a short reason
    is warned about)."""

    def test_r1_review2s_four_lines_are_all_accounted_for(self):
        log = [f'{CC} --why "For Zorblat',
               'python3 scripts/run_search.py --session S --title "pilot',
               NEXT,
               "python3 scripts/session.py finalize --dir S"]
        for reasons in ((), ("For Zorblat",)):
            with self.subTest(reasons=bool(reasons)):
                got, notes = redact("\n".join(log), reasons)
                self.assertEqual(got.split("\n"), [staff.WITHHELD_COMMAND] + log[1:])
                self.assertEqual(notes, ["a core-client command's quoted value does not end on its "
                                         "line: the command was withheld, the lines after it kept "
                                         "as logged"])

    def test_r1_a_multi_line_value_keeps_its_other_lines_and_by_value_cleans_them(self):
        why = "for the Example lab\nbefore their CoreOmics form is in"
        log = f'{CC} --why "for the Example lab\nbefore their CoreOmics form is in"\n{NEXT}'
        got, notes = redact(log, (why,))
        self.assertEqual(got.split("\n"), [staff.WITHHELD_COMMAND, '[reason-redacted]"', NEXT])
        self.assertEqual(len(notes), 2)                  # withheld; one other line by value

    def test_r2_an_unrelated_path_or_title_is_untouched(self):
        log = ("python3 scripts/run_search.py --out /x/S/pilotproject2026/search\n"
               'python3 scripts/make_analysis_html.py --title "Results for the lab"\n'
               'echo "see Zorblat" \'lab notes\'\n'
               f"{CC} --why pilotproject2026")
        for reasons in (("pilotproject2026",), ("for the",), ("Zorblat lab",)):
            with self.subTest(reason=reasons[0]):
                got, notes = redact(log, reasons)
                self.assertEqual(got.split("\n")[:3], log.split("\n")[:3])
                self.assertEqual(notes, [])
        got, _ = redact(log, ("pilotproject2026",))
        self.assertEqual(got.split("\n")[3], f"{CC} --why {staff.REDACTED_REASON_ARG}")

    def test_r2_a_hit_outside_a_core_client_command_is_said(self):
        log = f"printf '%s' '{REASON}' > /tmp/w\n{NEXT}"
        got, notes = redact(log, (REASON,))
        self.assertEqual(got, f"printf '%s' '{staff.REDACTED_REASON_ARG}' > /tmp/w\n{NEXT}")
        self.assertEqual(notes, ["a reason recorded for a core-client mark was found in 1 other "
                                 "line(s) and replaced there by [reason-redacted]"])
        # inside a word, a path or a file name it is not the reason
        got, notes = redact(f"cp /x/{REASON.replace(' ', '_')}.txt /y\n{NEXT}", (REASON,))
        self.assertEqual(notes, [])

    def test_r3_a_short_reason_is_warned_about(self):
        from test_core_client_survey import mark, with_analysis
        with tempfile.TemporaryDirectory() as tmp:
            for k, (why, warned) in enumerate((("QC pilot", True), ("pilotproject2026", True),
                                               ("a pilot for the Example lab", False))):
                with self.subTest(why=why):
                    os.makedirs(os.path.join(tmp, str(k)))
                    s = with_analysis(os.path.join(tmp, str(k)))
                    rc, out, err = mark(s, why=why)
                    self.assertEqual(rc, 0, out + err)
                    w = json.loads(out)["warning"]
                    self.assertEqual(bool(w), warned, w)
                    if warned:
                        self.assertIn("[session] WARNING: the reason is under 12 characters", err)

    def test_r3_the_limit_is_documented(self):
        ref = os.path.join(os.path.dirname(HERE), "references", "core-submissions.md")
        with open(ref, encoding="utf-8") as fh:
            text = " ".join(fh.read().split())
        self.assertIn("**Known limit.** Lines after a withheld command are kept as logged", text)
        import session as session_mod
        from test_core_client_survey import cli, STAFF
        rc, out, err = cli(session_mod, ["core-client", "--help"], STAFF)
        self.assertIn("Known limit", " ".join(out.split()))


class TheWritersSayWhatWasWithheld(unittest.TestCase):
    LOG = f"python3 scripts/session.py init --session S\n{CC} --why-file - <<'EOF'\n{NEXT}\n"

    def test_the_bundle_manifest(self):
        from test_core_client_survey import run
        with tempfile.TemporaryDirectory() as tmp:
            log, out = os.path.join(tmp, "commands.log"), os.path.join(tmp, "repro")
            with open(log, "w") as fh:
                fh.write(self.LOG)
            r = run("provenance.py", "--outdir", out, "--commands", log)
            self.assertEqual(r.returncode, 0, r.stderr[-1500:])
            with open(os.path.join(out, "inputs", "commands.log")) as fh:
                self.assertEqual(len(fh.read().strip().split("\n")), 3)
            with open(os.path.join(out, "MANIFEST.txt")) as fh:
                man = fh.read()
            self.assertIn("[NOTE]    commands.log: a core-client command's heredoc has no end", man)

    def test_the_session_zip_manifest(self):
        import finish_fixture
        import session as session_mod
        import test_deposit_package as tdp
        from job_env import job_env
        from test_core_client_survey import NOT_STAFF, cli
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
            with open(p["manifest_txt"]) as fh:
                man = fh.read()
            self.assertRegex(man, r"\[NOTE\]\s+logs/commands.log in the zip\s+-- a core-client "
                                  r"command's heredoc has no end")

    def test_the_run_registry_record(self):
        import record_run
        with tempfile.TemporaryDirectory() as tmp:
            src = os.path.join(tmp, "commands.log")
            with open(src, "w") as fh:
                fh.write(self.LOG)
            notes = record_run.copy_file(src, os.path.join(tmp, "reg", "commands.log"),
                                         record_run.Deadline(60))
            self.assertEqual(len(notes), 1)
            self.assertEqual(record_run.copy_file(src, os.path.join(tmp, "reg", "x.txt"),
                                                  record_run.Deadline(60)), [])


class TheMarkKeepsTheReasonAsGiven(unittest.TestCase):
    def test_why_text_in_the_staff_record(self):
        from test_core_client_survey import mark, with_analysis
        with tempfile.TemporaryDirectory() as tmp:
            s = with_analysis(tmp)
            why = "for the Example lab\nbefore their form is in"
            self.assertEqual(mark(s, why=why)[0], 0)
            e = staff.read(os.path.join(s, "logs", staff.CORE_CLIENT_RECORD))[-1]
            self.assertEqual((e["why_text"], e["why"]), (why, " ".join(why.split())))
            self.assertEqual(staff.reasons_for(os.path.join(s, "logs", "commands.log"))[0],
                             [why, " ".join(why.split())])


if __name__ == "__main__":
    unittest.main()
