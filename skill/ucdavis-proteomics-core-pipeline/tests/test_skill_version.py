#!/usr/bin/env python3
"""
skill_version.py is the ONE reader of the skill's version (.claude-plugin/plugin.json beside
scripts/). On HIVE it read "unknown" -- and make_deposit wrote "0.0.0" into sdrf.tsv -- whenever
scripts/ went up without .claude-plugin/ (access.md's "put ./scripts once"). Every reader now asks
skill_version.py, and a missing plugin.json is a tagged UNKNOWN, never a made-up number.

notify_slack.py (runs from stdin on HIVE, no siblings) and report_issue.sh (bash) carry mirrors;
this file keeps each one equal to the original.
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


class ReportIssueMirror(unittest.TestCase):
    SCRIPT = os.path.join(SCRIPTS, "report_issue.sh")

    def test_same_constant(self):
        with open(self.SCRIPT, encoding="utf-8") as fh:
            m = re.search(r'^VER_UNKNOWN="([^"]*)"$', fh.read(), re.M)
        self.assertIsNotNone(m, "report_issue.sh has no VER_UNKNOWN")
        self.assertEqual(m.group(1), sv.UNKNOWN)

    def _entry(self, sd):
        """The entry report_issue.sh writes when run from <sd> (a copy of scripts/)."""
        shutil.copy(self.SCRIPT, sd)
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


if __name__ == "__main__":
    unittest.main(verbosity=2)
