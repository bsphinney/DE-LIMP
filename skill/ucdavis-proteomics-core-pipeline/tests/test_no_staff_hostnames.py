#!/usr/bin/env python3
"""No staff computer's name in the published skill (security review, 2.12).

A skill issue record is named <date>_<user>_<session>+<time>-<hostname>-<pid>-<random>.md, and a
staff laptop's hostname (the Core's Windows naming: a building code, a room and the person's
initials) went into a test's docstring with the file name. The repository is public. Cite an
issue record as `<date>_<user>_<session>+….md`. The patterns are built from parts so this file
does not match itself.
"""
import os
import re
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)

HOSTNAME = re.compile("CBS" + r"-GC\d{3,}|-(?:" + "|".join(("MI" + "CH", "G" + "AB", "MI" + "NI"))
                      + r")-\d+-\d+")
TEXT = (".py", ".md", ".R", ".sh", ".json", ".txt", ".yaml", ".tsv", ".html")


class NoStaffHostnames(unittest.TestCase):
    def test_the_pattern_catches_the_real_shape(self):
        self.assertTrue(HOSTNAME.search("2026-10-08_x_PROT_0000+123523-" + "CBS" + "-GC0000-"
                                        + "MI" + "CH-282-415932196.md"))
        self.assertFalse(HOSTNAME.search("2026-10-08_x_PROT_0000+….md"))

    def test_none_in_the_skill(self):
        hits = []
        for dp, dns, fns in os.walk(SKILL):
            dns[:] = [d for d in dns if d not in ("__pycache__", ".pytest_cache")]
            for fn in fns:
                if not fn.endswith(TEXT):
                    continue
                path = os.path.join(dp, fn)
                with open(path, encoding="utf-8", errors="replace") as fh:
                    for i, line in enumerate(fh, 1):
                        if HOSTNAME.search(line):
                            hits.append(f"{os.path.relpath(path, SKILL)}:{i}")
        self.assertEqual(hits, [], "a staff computer's name: cite the record as "
                                   "<date>_<user>_<session>+….md")


if __name__ == "__main__":
    unittest.main()
