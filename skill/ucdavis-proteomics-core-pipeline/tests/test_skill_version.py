#!/usr/bin/env python3
"""
skill_version.py is the ONE reader of the skill's version (.claude-plugin/plugin.json beside
scripts/). On HIVE it read "unknown" -- and make_deposit wrote "0.0.0" into sdrf.tsv -- whenever
scripts/ went up without .claude-plugin/ (access.md's "put ./scripts once"). Every reader now asks
skill_version.py, and a missing plugin.json is a tagged UNKNOWN, never a made-up number.

notify_slack.py (runs from stdin on HIVE, no siblings), skill_version.sh (bash: report_issue.sh
sources it, and its --check-hive reads HIVE's copy with it) and skill_version.R (run_de.R's, in R)
carry mirrors; this file keeps each one equal to the original. It also keeps skill_version.py's
core_admins() equal to skill_version.sh's skill_core_admins, the two readers of core_admins.txt.
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL, "scripts")
sys.path.insert(0, SCRIPTS)

sys.path.insert(0, HERE)
import skill_version as sv      # noqa: E402
import notify_slack             # noqa: E402
from job_env import job_env     # noqa: E402

PLUGIN = os.path.join(SKILL, ".claude-plugin", "plugin.json")


def fixture(tmp, name, plugin_text):
    """<tmp>/<name>/scripts, with ../.claude-plugin/plugin.json holding plugin_text (None: absent)."""
    sd = os.path.join(tmp, name, "scripts")
    os.makedirs(sd)
    if plugin_text is not None:
        os.makedirs(os.path.join(tmp, name, ".claude-plugin"))
        with open(os.path.join(tmp, name, ".claude-plugin", "plugin.json"), "w",
                  encoding="utf-8") as fh:
            fh.write(plugin_text)
    return sd


CASES = {"present": '{"name": "x", "version": "9.8.7"}',
         "absent": None,
         "no_version": '{"name": "x"}',
         "blank_version": '{"version": "  "}',
         "not_json": "version: 1.2.3",
         "not_an_object": '["9.8.7"]'}


class TheOneReader(unittest.TestCase):
    def test_reads_the_installed_skill(self):
        with open(PLUGIN, encoding="utf-8") as fh:
            want = json.load(fh)["version"]
        self.assertEqual(sv.skill_version(), want)
        self.assertEqual(sv.label(want), f"v{want}")

    def test_anything_but_a_version_is_the_tagged_unknown(self):
        with tempfile.TemporaryDirectory() as tmp:
            for name, text in CASES.items():
                with self.subTest(name):
                    got = sv.skill_version(fixture(tmp, name, text))
                    self.assertEqual(got, "9.8.7" if name == "present" else sv.UNKNOWN)
        self.assertIn("unknown", sv.UNKNOWN)
        self.assertIn("plugin.json", sv.UNKNOWN)
        self.assertEqual(sv.label(sv.UNKNOWN), sv.UNKNOWN, "never 'v(unknown ...)'")

    def test_the_ascii_tag_is_the_same_tag(self):
        """sdrf.tsv gets UNKNOWN_ASCII: an SDRF validator may reject the em dash."""
        self.assertEqual(sv.UNKNOWN_ASCII, "unknown - plugin.json not found")
        self.assertTrue(sv.UNKNOWN_ASCII.isascii())
        self.assertEqual(sv.label(sv.UNKNOWN, ascii=True), sv.UNKNOWN_ASCII)
        self.assertEqual(sv.label("2.8.0", ascii=True), "v2.8.0")
        with open(os.path.join(SCRIPTS, "skill_version.py"), encoding="utf-8") as fh:
            src = fh.read()
        self.assertEqual(src.count("plugin.json not found"), 1, "one definition of the tag")

    def test_no_reader_spells_its_own_path_to_plugin_json(self):
        """record_run, provenance, make_deposit and core_submission each read plugin.json on
        their own, with their own fallbacks ("unknown", None, "0.0.0", "?")."""
        allowed = {"skill_version.py", "notify_slack.py",   # the reader and its mirror
                   "bump_version.py",                       # writes the version
                   "record_run.py"}                         # copies the file to HIVE, no read
        for name in sorted(os.listdir(SCRIPTS)):
            if not name.endswith(".py") or name in allowed:
                continue
            with open(os.path.join(SCRIPTS, name), encoding="utf-8") as fh:
                src = fh.read()
            self.assertNotIn('"plugin.json"', src, f"{name} reads plugin.json itself")
        with open(os.path.join(SCRIPTS, "record_run.py"), encoding="utf-8") as fh:
            src = fh.read()
        body = src[src.index("def skill_version("):]
        body = body[:body.index("\ndef ")]
        self.assertIn("import skill_version", body)
        self.assertNotIn("plugin.json\")", body)

    def test_sdrf_and_the_prep_script_never_carry_a_made_up_version(self):
        with open(os.path.join(SCRIPTS, "make_deposit.py"), encoding="utf-8") as fh:
            src = fh.read()
        self.assertNotIn('"0.0.0"', src)
        self.assertNotIn(" v{skill_ver}", src)


class NotifySlackMirror(unittest.TestCase):
    def test_same_constant(self):
        self.assertEqual(notify_slack.SKILL_VERSION_UNKNOWN, sv.UNKNOWN)

    def test_same_answer_for_every_fixture(self):
        with tempfile.TemporaryDirectory() as tmp:
            for name, text in CASES.items():
                with self.subTest(name):
                    sd = fixture(tmp, name, text)
                    self.assertEqual(notify_slack._skill_version(sd), sv.skill_version(sd))
        self.assertEqual(notify_slack._skill_version(), sv.skill_version())

    def test_from_stdin_there_is_no_here(self):
        self.assertEqual(notify_slack._skill_version(None), sv.UNKNOWN)

    def test_the_footer_never_reads_v_unknown(self):
        self.assertEqual(notify_slack._skill_label("2.8.0"), "ucdavis-proteomics-core-pipeline v2.8.0")
        self.assertEqual(notify_slack._skill_label(sv.UNKNOWN),
                         f"ucdavis-proteomics-core-pipeline {sv.UNKNOWN}")
        self.assertEqual(notify_slack._skill_label(None), "ucdavis-proteomics-core-pipeline")


class BashMirror(unittest.TestCase):
    """skill_version.sh: the one bash reading, sourced by report_issue.sh (Windows `python3` is
    often the Microsoft Store stub) and used by --check-hive on this computer and on HIVE."""
    SCRIPT = os.path.join(SCRIPTS, "skill_version.sh")

    def bash(self, snippet, sd=None):
        script = self.SCRIPT
        if sd is not None:                  # a copy of scripts/ in a fixture
            shutil.copy(self.SCRIPT, sd)
            script = os.path.join(sd, "skill_version.sh")
        with tempfile.TemporaryDirectory() as tmp:
            r = subprocess.run(["bash", "-c", f'. "$1"; {snippet}', "bash", script],
                               capture_output=True, timeout=60,
                               env=job_env(tmp, LC_ALL="C.UTF-8", LANG="C.UTF-8"))
        self.assertEqual(r.returncode, 0, r.stderr.decode("utf-8", "replace"))
        return r.stdout.decode("utf-8").rstrip("\n")

    def test_same_constant(self):
        self.assertEqual(self.bash('printf "%s" "$SKILL_VERSION_UNKNOWN"'), sv.UNKNOWN)

    def test_same_answer_for_every_fixture(self):
        with tempfile.TemporaryDirectory() as tmp:
            for name, text in CASES.items():
                with self.subTest(name):
                    sd = fixture(tmp, name, text)
                    self.assertEqual(self.bash("skill_version", sd), sv.skill_version(sd))
                    self.assertEqual(self.bash('skill_version "$(dirname "$1")"', sd),
                                     sv.skill_version(sd))
        self.assertEqual(self.bash("skill_version"), sv.skill_version())

    def test_run_as_a_command_it_prints_the_version(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = subprocess.run(["bash", self.SCRIPT], capture_output=True, text=True,
                               env=job_env(tmp), timeout=60)
        self.assertEqual((r.returncode, r.stdout.strip()), (0, sv.skill_version()))

    def test_same_label(self):
        self.assertEqual(self.bash('skill_label 2.8.0; skill_label "$SKILL_VERSION_UNKNOWN"'),
                         f"{sv.label('2.8.0')}\n{sv.label(sv.UNKNOWN)}")


#: core_admins.txt files that differ only in what both readers must strip or skip.
ADMIN_FILES = {
    "comments_and_blank_lines": "# admins\n\n#msalemi\n   # indented comment\nbrettsp\n",
    "blanks_at_both_ends": "  brettsp  \n\tgabrig\t\n \x0b\x0cmsalemi \x0c\n",
    "crlf_and_a_lone_cr": "brettsp\r\n\r\ngab\rrig\r\n \r \n",
    "comments_after_a_name": "brettsp # the director\ngabrig#staff\n#\n",
    "no_final_newline": "brettsp\ngabrig",
    "a_blank_or_comma_inside_and_twice": "Bad Name\na,b\nbrettsp\nbrettsp\n",
    "empty": "",
    "only_comments_and_blanks": "# nobody\n\n   \n\t\n",
}


class CoreAdminsMirror(unittest.TestCase):
    """scripts/core_admins.txt has two readers: skill_version.sh's skill_core_admins (bash:
    --check-hive, --publish-release) and skill_version.py's core_admins() (Python: notes.py,
    slack_collab.py). They must give the same admins for every file, or the two sides of the
    skill trust different people."""
    SCRIPT = os.path.join(SCRIPTS, "skill_version.sh")

    def bash_admins(self, path, locale):
        with tempfile.TemporaryDirectory() as tmp:
            r = subprocess.run(["bash", "-c", '. "$1"; skill_core_admins "$2"', "bash",
                                self.SCRIPT, path or ""],
                               capture_output=True, timeout=60,
                               env=job_env(tmp, LC_ALL=locale, LANG=locale))
        return r.stdout.decode("utf-8").splitlines()

    def test_same_admins_for_every_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            for name, text in ADMIN_FILES.items():
                p = os.path.join(tmp, name + ".txt")
                with open(p, "w", encoding="utf-8", newline="") as fh:
                    fh.write(text)
                for locale in ("C.UTF-8", "C"):
                    with self.subTest(name, locale=locale):
                        self.assertEqual(self.bash_admins(p, locale), sv.core_admins(p))
            missing = os.path.join(tmp, "missing.txt")
            self.assertEqual(self.bash_admins(missing, "C"), [])
            self.assertEqual(sv.core_admins(missing), [])

    def test_same_admins_for_the_shipped_list(self):
        self.assertEqual(self.bash_admins(None, "C.UTF-8"), sv.core_admins())
        self.assertEqual(sv.core_admins(), ["brettsp"])


class ReportIssueMirror(unittest.TestCase):
    SCRIPT = os.path.join(SCRIPTS, "report_issue.sh")

    def test_it_reads_the_version_through_skill_version_sh(self):
        """One bash reading, not a second copy of it (CLAUDE.md rule 3)."""
        with open(self.SCRIPT, encoding="utf-8") as fh:
            code = "\n".join(l for l in fh.read().splitlines() if not l.lstrip().startswith("#"))
        self.assertIn('. "$HERE/skill_version.sh"', code)
        self.assertIn('skill_version "$HERE"', code)
        self.assertNotIn("plugin.json", code)
        self.assertNotIn('"version"', code)

    def _entry(self, sd, siblings=("skill_version.sh",)):
        """The entry report_issue.sh writes when run from <sd> (a copy of scripts/)."""
        for f in ("report_issue.sh",) + tuple(siblings):
            shutil.copy(os.path.join(SCRIPTS, f), sd)
        out = os.path.join(os.path.dirname(sd), "issues")
        # no HIVE login, nothing of a real job: the entry goes to `out` directly
        env = job_env(os.path.dirname(sd), HOME=os.path.dirname(sd), SKILL_ISSUES_DIR=out,
                      SKILL_ISSUES_LOCAL_DIR=out + "_local", TMPDIR=os.path.dirname(sd),
                      LC_ALL="C.UTF-8", LANG="C.UTF-8")
        os.makedirs(out)
        r = subprocess.run(["bash", os.path.join(sd, "report_issue.sh"), "--title", "t",
                            "--what", "w"], capture_output=True, env=env, timeout=60)
        self.assertEqual(r.returncode, 0, r.stderr.decode("utf-8", "replace"))
        (name,) = [n for n in os.listdir(out) if n.endswith(".md")]
        with open(os.path.join(out, name), encoding="utf-8") as fh:
            return fh.read()

    def test_same_answer_as_the_reader(self):
        with tempfile.TemporaryDirectory() as tmp:
            for name in CASES:
                with self.subTest(name):
                    sd = fixture(tmp, name, CASES[name])
                    text = self._entry(sd)
                    want = sv.skill_version(sd)
                    self.assertIn(f"- **Skill version:** {want}\n", text)
                    self.assertIn(f"- **Skill:** {sv.label(want)}, mode ", text)

    def test_without_its_sibling_the_entry_is_still_written(self):
        """Recording a problem must never become a second problem."""
        with tempfile.TemporaryDirectory() as tmp:
            text = self._entry(fixture(tmp, "alone", CASES["present"]), siblings=())
        self.assertIn("- **Skill version:** (unknown -- skill_version.sh not found", text)


def _r_has(*pkgs):
    if not shutil.which("Rscript"):
        return False
    expr = ("quit(status = if (all(vapply(c(%s), requireNamespace, logical(1), quietly = TRUE))) 0 "
            "else 1)" % ", ".join(f'"{p}"' for p in pkgs))
    return subprocess.run(["Rscript", "-e", expr], capture_output=True).returncode == 0


@unittest.skipUnless(_r_has("jsonlite"), "needs Rscript + jsonlite")
class RunDeRMirror(unittest.TestCase):
    """run_de.R reads the version through skill_version.R, the R mirror: same constant, same
    answer for every fixture, with jsonlite and without it (its regex fallback)."""
    R_FILE = os.path.join(SCRIPTS, "skill_version.R")

    def r(self, expr):
        env = dict(os.environ, LC_ALL="en_US.UTF-8", LANG="en_US.UTF-8")
        p = subprocess.run(["Rscript", "-e", f'source("{self.R_FILE}", encoding = "UTF-8"); {expr}'],
                           capture_output=True, env=env)
        self.assertEqual(p.returncode, 0, p.stderr.decode("utf-8", "replace"))
        return json.loads(p.stdout.decode("utf-8"))

    def test_same_constant(self):
        self.assertEqual(self.r("cat(jsonlite::toJSON(SKILL_VERSION_UNKNOWN, auto_unbox = TRUE))"),
                         sv.UNKNOWN)

    def test_same_answer_for_every_fixture_both_parsers(self):
        with tempfile.TemporaryDirectory() as tmp:
            sds = {name: fixture(tmp, name, text) for name, text in CASES.items()}
            for use in ("TRUE", "FALSE"):
                calls = ", ".join(f'"{n}" = skill_version("{sd}", use_jsonlite = {use})'
                                  for n, sd in sds.items())
                got = self.r(f"cat(jsonlite::toJSON(list({calls}), auto_unbox = TRUE))")
                for name, sd in sds.items():
                    with self.subTest(name, jsonlite=use):
                        self.assertEqual(got[name], sv.skill_version(sd))
        got = self.r(f'cat(jsonlite::toJSON(skill_version("{SCRIPTS}"), auto_unbox = TRUE))')
        self.assertEqual(got, sv.skill_version())

    def test_same_label(self):
        got = self.r('cat(jsonlite::toJSON(list(skill_label("2.8.0"), '
                     'skill_label(SKILL_VERSION_UNKNOWN)), auto_unbox = TRUE))')
        self.assertEqual(got, [sv.label("2.8.0"), sv.label(sv.UNKNOWN)])

    def test_run_de_reads_it_through_the_mirror(self):
        with open(os.path.join(SCRIPTS, "run_de.R"), encoding="utf-8") as fh:
            code = [l for l in fh.read().splitlines() if not l.lstrip().startswith("#")]
        src = "\n".join(code)
        self.assertIn('.sibling("skill_version.R")', src)
        self.assertIn("skill_label(skill_version(.script_dir))", src)
        self.assertNotIn("plugin.json", src.replace("skill_version.R", ""))


if __name__ == "__main__":
    unittest.main(verbosity=2)
