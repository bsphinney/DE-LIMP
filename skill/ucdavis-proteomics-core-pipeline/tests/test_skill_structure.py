#!/usr/bin/env python3
"""SKILL.md is a router; the detail lives in references/ (the slim-down approved 2026-10-16).

Every analysis reads SKILL.md in full before anything else, so it stays short: each step keeps
its command, the rules that must never be missed, and the reference section to read before the
step. The checks here keep that shape from drifting back:

  * SKILL.md stays under 500 lines;
  * every reference SKILL.md names exists, and every quoted section after a reference name
    ("access.md \"Step 0: is this skill current?\"") is the start of one of its headings;
  * every reference over 100 lines opens with a contents list whose line numbers are right
    (regenerate them with `python3 tests/test_skill_structure.py --toc` after an edit);
  * a reference never sends the reader to read itself, and every reference another reference
    tells you to read is also named in SKILL.md, so no rule is two hops from the router.
"""
import os
import re
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
REFS = os.path.join(SKILL, "references")
TOC_START = "<!-- contents: tests/test_skill_structure.py --toc -->"
TOC_END = "<!-- /contents -->"
MAX_ROUTER_LINES = 500
TOC_OVER = 100


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def references():
    return sorted(f for f in os.listdir(REFS) if f.endswith(".md"))


def headings(text):
    """(line number, heading text) of every ## heading outside code fences."""
    out, fence = [], False
    for n, ln in enumerate(text.split("\n"), 1):
        if ln.lstrip().startswith("```"):
            fence = not fence
        elif not fence and ln.startswith("## "):
            out.append((n, ln[3:].strip()))
    return out


def strip_toc(text):
    if TOC_START not in text:
        return text
    a = text.index(TOC_START)
    b = text.index(TOC_END, a) + len(TOC_END)
    return text[:a].rstrip("\n") + "\n\n" + text[b:].lstrip("\n")


def with_toc(text):
    """The text with a fresh contents list after its title paragraph."""
    text = strip_toc(text)
    lines = text.split("\n")
    first = next(i for i, ln in enumerate(lines) if ln.startswith("## "))
    # after the title and its opening paragraph(s): just before the first ## heading
    at = first
    while at > 0 and lines[at - 1] == "":
        at -= 1
    heads = headings(text)
    added = len(heads) + 4                     # a blank line, start, label, entries..., end
    entries = [f"- {n + added}: {h}" for n, h in heads]
    block = [TOC_START, "**Contents** (section: line)"] + entries + [TOC_END]
    new = lines[:at] + [""] + block + lines[at:]
    out = "\n".join(new)
    # the entries must name the heading lines of the result
    for (n, h), entry in zip(headings(out), entries):
        assert entry == f"- {n}: {h}", (entry, n, h)
    return out


class RouterStaysShort(unittest.TestCase):
    def test_skill_md_is_under_500_lines(self):
        n = read(os.path.join(SKILL, "SKILL.md")).count("\n")
        self.assertLess(n, MAX_ROUTER_LINES, "SKILL.md is the router: move detail to references/")


class RouterPointsAtRealSections(unittest.TestCase):
    def setUp(self):
        self.skill = read(os.path.join(SKILL, "SKILL.md"))
        self.refs = {f: read(os.path.join(REFS, f)) for f in references()}

    def test_every_named_reference_exists(self):
        named = set(re.findall(r"(?<![\w/.-])([a-z][a-z_-]+\.md)\b", self.skill))
        named -= {"SKILL.md"}
        missing = sorted(f for f in named
                         if f not in self.refs
                         and not re.search(r"(?:output|logs|podcast|tables)/[^ ]*" + re.escape(f), self.skill)
                         and f.upper() != f)
        # files the pipeline writes (AUDIT.md, README.md ...) are upper case or output paths
        self.assertEqual([m for m in missing if m not in self.written_files()], [])

    @staticmethod
    def written_files():
        return {"meaning_check.md", "decisions.md", "review_proteomics.md"}

    def test_every_quoted_section_is_a_heading(self):
        bad = []
        flat = " ".join(self.skill.split())
        for ref, title in re.findall(r"([a-z][a-z_-]+\.md),? \"([^\"]+)\"", flat):
            if ref not in self.refs:
                continue
            heads = [h.replace("`", "") for _, h in headings(self.refs[ref])]
            want = title.replace("`", "").rstrip(".")
            # the title is the start of a heading, ending where a word ends: "4." names
            # "4. Resolving ...", never "4a. DIA ..."
            if not any(h == want or (h.startswith(want) and h[len(want)] in " .:(—,")
                       for h in heads):
                bad.append(f'{ref} "{title}"')
        self.assertEqual(bad, [], "SKILL.md names a section its reference does not have")


class StepLevelRemindersStayAtTheirStep(unittest.TestCase):
    """2.13.1/2.13.2 rules that bite at a particular step stay AT that step in the router, not only
    in a general sentence or a reference (independent review of the slim-down, H1-H3, M1-M3)."""

    def step(self, head, nxt):
        skill = read(os.path.join(SKILL, "SKILL.md"))
        i = skill.index(head)
        return " ".join(skill[i:skill.index(nxt, i)].split())

    def test_each_step_keeps_its_reminder(self):
        cases = {
            ("### 0a.", "### 0. "): ["These move themselves to a compute node", "setup.sh",
                                     "`session.py finalize`", "`core_submission.py locate`/`stage`/`deliver`",
                                     "provenance.py", "never the 2-minute default"],
            ("### 0. ", "### 0b."): ["**Check again**", "--session-start", "then `--session-start` again"],
            ("### 1c.", "### 2. "): ['save_transcript.py --hive "$S" --remember'],
            ("### 2. ", "### 3. "): ["--resolution-source detected", "--ms2-analyzer ITMS", "in the background",
                                     "dotnet-core-sdk/8.0.4"],
            ("### 7b.", "### 7c."): ["every 5 minutes or more"],
            ("### 8. ", "### 8a."): ["on_compute.py --header", "lint-jobscript.sh", "--fasta-meta"],
            ("### 10. ", "### 11a."): ["in the background on HIVE"],
            ("### 12. ", "### 12b."): ["in the background on HIVE", "[--call-id <token>]"],
            ("### 12b.", "## Keeping the record"): ["in the background", "70 minutes", "bioshare ensure"],
        }
        for (head, nxt), wants in cases.items():
            text = self.step(head, nxt)
            for w in wants:
                with self.subTest(step=head, want=w):
                    self.assertIn(w, text)


class LongReferencesHaveContents(unittest.TestCase):
    def test_contents_lists_are_present_and_right(self):
        stale = []
        for f in references():
            text = read(os.path.join(REFS, f))
            if strip_toc(text).count("\n") <= TOC_OVER:
                continue
            if TOC_START not in text or with_toc(text) != text:
                stale.append(f)
        self.assertEqual(stale, [], "run: python3 tests/test_skill_structure.py --toc")


class NoSecondHop(unittest.TestCase):
    def test_a_reference_never_sends_you_to_read_itself(self):
        bad = []
        for f in references():
            for m in re.finditer(r"[Rr]ead `references/([a-z_-]+\.md)`", read(os.path.join(REFS, f))):
                if m.group(1) == f:
                    bad.append(f)
        self.assertEqual(bad, [])

    def test_what_a_reference_says_to_read_is_named_in_skill_md(self):
        skill = read(os.path.join(SKILL, "SKILL.md"))
        bad = []
        for f in references():
            for m in re.finditer(r"[Rr]ead `references/([a-z_-]+\.md)`", read(os.path.join(REFS, f))):
                if m.group(1) not in skill:
                    bad.append(f"{f} -> {m.group(1)}")
        self.assertEqual(bad, [])


def regenerate():
    for f in references():
        p = os.path.join(REFS, f)
        text = read(p)
        if strip_toc(text).count("\n") <= TOC_OVER:
            new = strip_toc(text)
        else:
            new = with_toc(text)
        if new != text:
            with open(p, "w", encoding="utf-8") as fh:
                fh.write(new)
            print("contents updated:", f)


if __name__ == "__main__":
    if sys.argv[1:] == ["--toc"]:
        regenerate()
    else:
        unittest.main(verbosity=2)
