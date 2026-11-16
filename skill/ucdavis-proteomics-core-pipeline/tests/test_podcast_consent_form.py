#!/usr/bin/env python3
"""The consent for anyone else's data has one form (2.15 Ruling 96.5).

references/podcast.md ("Privacy and consent") records it as
`--cloud-ok "<who agreed>, <date>, <paid|free> key"`: who, when, and on which key tier, since the
tier decides what Google may do with the text. make_podcast.py's refusals and --help asked for
"who agreed, when" -- no tier. Now every message there says CLOUD_OK_FORM, and every form
podcast.md shows is the same one. Stdlib only; nothing is sent anywhere.
"""
import contextlib
import io
import os
import re
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL_DIR = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL_DIR, "scripts")
sys.path.insert(0, SCRIPTS)

import make_podcast as mp  # noqa: E402

FORM = '"<who agreed>, <date>, <paid|free> key"'


def read(*parts):
    with open(os.path.join(SKILL_DIR, *parts), encoding="utf-8") as fh:
        return fh.read()


class OneConsentForm(unittest.TestCase):
    def test_the_constant_is_podcast_md_s_form(self):
        self.assertEqual(mp.CLOUD_OK_FORM, FORM)
        self.assertIn(f"--cloud-ok {FORM}", read("references", "podcast.md"))

    def test_every_form_podcast_md_shows_is_the_one(self):
        doc = read("references", "podcast.md")
        shown = re.findall(r'--cloud-ok ("<[^"]*")', doc)
        self.assertTrue(shown)
        self.assertEqual(sorted(set(shown)), [FORM], shown)

    def test_make_podcast_asks_for_nothing_else(self):
        src = read("scripts", "make_podcast.py")
        self.assertNotIn("who agreed, when", src)
        self.assertNotIn("note of who agreed and when", src)
        # every quoted --cloud-ok example in the script is the constant
        for m in re.finditer(r'--cloud-ok \\?"([^"\\]*)', src):
            self.assertEqual(m.group(1), FORM.strip('"'), m.group(0))

    def test_the_help_of_render_verify_and_job_says_it(self):
        for cmd in ("render", "verify", "job"):
            out = io.StringIO()
            with contextlib.redirect_stdout(out), self.assertRaises(SystemExit):
                mp.main([cmd, "--help"])
            text = " ".join(out.getvalue().split())
            self.assertIn("<who agreed>, <date>, <paid|free> key", text, cmd)

    def test_each_no_consent_refusal_asks_for_the_form(self):
        """render, verify and job without consent, and the two --core-run refusals that send
        the reader to --cloud-ok."""
        src = read("scripts", "make_podcast.py")
        for where in ("not sending anything. --tts", "not sending anything. verify sends",
                      "so it needs consent: --core-run for a Core run, or",
                      "for anyone else's data, ask them \"\n"):
            i = src.index(where)
            self.assertIn("CLOUD_OK_FORM", src[i:i + 700], where)
        # the two --core-run refusals (internal work, no submission) that point at --cloud-ok
        self.assertEqual(src.count("pass --cloud-ok {CLOUD_OK_FORM}."), 2)


if __name__ == "__main__":
    unittest.main()
