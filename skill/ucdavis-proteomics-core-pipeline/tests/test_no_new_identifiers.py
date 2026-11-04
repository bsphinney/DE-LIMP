#!/usr/bin/env python3
"""No NEW real identifier in the published skill (2.14.2, the C2 scrub).

The repository is public. 2.14.2 removed every client name from the skill's text (PIs, labs,
companies, people outside the Core) and every submission number, staff login and staff name that
2.14.x had added. The identifiers that were already public on main -- staff logins and names and
old submission numbers, almost all in test fixtures and comments -- stay for now: the 2.15
slim-down rewrites those same lines, and touching them here would make it conflict on every one.

So this guard fails on any identifier that is not one of those frozen, already-public lines:

* the rule is tests/identifier_rules.py (problems() is the slim-down's own rule, copied; plus
  client_names(), the client names that rule does not know, read from a PRIVATE list outside the
  repository -- $SKILL_PRIVATE_NAMES or ~/.config/ucdavis-proteomics/client_names.txt. Without
  it the client-name part is skipped, and test_the_client_name_check_ran says so; the
  submission-number and staff parts always run);
* tests/identifier_baseline.json freezes the hits that were left as KEYS, never their text: for
  each line, identifier_rules.baseline_key(path, line) -- the sha256 of "<path>\0<line>" -- and
  how often it occurs; for each file name, the key of (path, name). A list of ~830 staff and
  submission lines verbatim in one new public file would republish them (2.14.2, before the
  push). It holds staff and submission-number hits only; a client name is never excused by it
  (test_the_baseline_holds_no_identifier_text, test_a_new_identifier_is_still_found).

A new line, or an edited old one, that carries an identifier fails here. Write it the way
identifier_rules.py says instead ("a Core submission", PROT_9xxx, `<your HIVE user>`, staff1,
"a Core member", "a client lab"). Never add to the baseline: it only shrinks.

When the 2.15 slim-down lands, its full guard (tests/test_no_real_identifiers.py, no baseline,
no allowance for old lines) replaces this test, and identifier_baseline.json is deleted with it.
"""
import collections
import json
import os
import re
import shutil
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import identifier_rules as rules  # noqa: E402

BASELINE_REL = os.path.join("tests", "identifier_baseline.json")


# The baseline as frozen from the release head (2.14.2): entries, occurrences, named paths. It
# only shrinks -- lower these with it. A key gives no text to review, so a grown baseline (a key
# added to quiet the guard) is refused by size.
FROZEN_LINES, FROZEN_OCCURRENCES, FROZEN_NAMES = 826, 917, 3
KEY = re.compile(r"[0-9a-f]{64}")


def load_baseline(path=None):
    """(Counter of line key -> count, set of name keys) from the baseline file."""
    with open(path or os.path.join(SKILL, BASELINE_REL), encoding="utf-8") as fh:
        data = json.load(fh)
    lines = collections.Counter()
    for key, count in data["lines"]:
        lines[key] += count
    return lines, set(data["names"])


def new_identifiers(skill=SKILL, baseline=None, private=None):
    """Every identifier under `skill` the baseline does not excuse, one line each: a staff or
    submission-number hit whose key (identifier_rules.baseline_key) is not in it, or occurs more
    often than its count; a file name whose key is not in it; and every client name -- the
    baseline never excuses one. `private`: a client-name list ({} = none; default: the
    private list, when there is one)."""
    frozen, frozen_names = baseline if baseline is not None else load_baseline()
    if private is None:
        private = rules.client_list()[0] or {}      # read once; {} = absent (the client part skips)
    seen, new = collections.Counter(), []
    for rel, i, text in rules.walk(skill, skip={BASELINE_REL}):
        where = f"{rel}:{i}" if i else f"{rel}: in the file name"
        staff = rules.line_hits(rel, text, {})      # the slim-down's rule, the public text allowed
        if staff:
            key = rules.baseline_key(rel, text)
            if i == 0:
                excused = key in frozen_names
            else:
                seen[key] += 1
                excused = seen[key] <= frozen.get(key, 0)
            if not excused:
                new.append(f"{where}: {'; '.join(staff)}")
        client = rules.client_names(text, private) if private else []
        if client:
            new.append(f"{where}: " + "; ".join(f"a client name (private list entry {n})"
                                              for n in client))
    return new


class TheRules(unittest.TestCase):
    def test_the_slim_down_rule_catches_the_real_shapes(self):
        dataset = ("va08172026", "Sil")      # out of order: two parts side by side read as one word
        real = ["PROT_" + "07" + "56", "PROT" + "0793", "prot-" + "756", "PROT-" + "0812",
                "coreomics_prot" + "0756", "PROT " + "0801", "submission " + "756",
                "submission #" + "0793", "brett" + "sp", "AD3+" + "gab" + "rig", "MSAL" + "EMI",
                "Mich" + "elle's Claude", "info@" + "frag" + "matics.com", "the " + "Dick" + "son lab",
                dataset[1] + dataset[0], "0022066" + "cd85f"]
        for s in real:
            self.assertTrue(rules.line_hits("scripts/x.py", f"# {s} 2026-09-23"), s)

    def test_the_allowed_fakes_pass(self):
        for s in ("PROT_0807", "prot-807", "PROT_0000", "PROT_9001", "PROT_99999", "submission 807",
                  "Total_prot_10ug", "Uniprot_1234", "bstaff", "mstaff", "staff1", "tester",
                  "a Core member", "a client lab", "the bat-plate study", "Burkholderia",
                  "Silva Campitelli de Barros BC", "amschaal/coreomics_fs"):
            self.assertEqual(rules.line_hits("tests/x.py", f"x {s} y"), [], s)

    def private_list(self, text):
        """A private list holding `text`, as $SKILL_PRIVATE_NAMES (removed after the test)."""
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        path = os.path.join(d, "client_names.txt")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
        patcher = mock.patch.dict(os.environ, {rules.PRIVATE_NAMES_ENV: path})
        patcher.start()
        self.addCleanup(patcher.stop)
        return path

    def test_a_client_name_is_found_however_it_is_written(self):
        """The private list's names match one word or several joined; the name is never in
        this repository, and a hit names the entry number, not the name."""
        self.private_list("# fake names, for this test only\nQuillonfake\n\nVan Quill  # two words\n")
        names, why = rules.client_list()
        self.assertEqual(len(names), 2, why)
        for s in ("the Quillonfake lab", "QUILLONFAKE", "2026-10-01_Quillonfake_Bat_HT",
                  "QuillonfakeBat", "Quillonfake08172026", "Qüillonfake", "Van Quill",
                  "VanQuill", "van-quill", "VANQUILL"):
            self.assertTrue(rules.line_hits("references/x.md", s), s)
        for s in ("Quillonfakery", "unquillonfake", "Van", "Quill", "Vanquillon"):
            self.assertEqual(rules.client_names(s), [], s)
        hit = rules.line_hits("references/x.md", "ask the Van Quill lab")
        self.assertEqual(hit, ["a client name (private list entry 4)"])
        self.assertNotIn("quill", " ".join(hit).lower())

    def test_without_the_private_list_only_the_client_part_is_skipped(self):
        """Public CI has no list: the client part says it was skipped and where the list goes;
        the submission-number and staff checks still run."""
        with mock.patch.dict(os.environ, {rules.PRIVATE_NAMES_ENV: os.path.join(
                tempfile.gettempdir(), "no-such-dir-c2", "client_names.txt")}):
            names, why = rules.client_list()
            self.assertIsNone(names)
            self.assertIn(rules.SKIPPED, why)
            self.assertIn("no-such-dir-c2", why)
            self.assertIn(rules.PRIVATE_NAMES_ENV, why)
            self.assertEqual(rules.client_names("the Quillonfake lab"), [])
            self.assertTrue(rules.line_hits("scripts/x.py", "# PROT_" + "0756 2026-09-23"))
            self.assertTrue(rules.line_hits("scripts/x.py", "# brett" + "sp 2026-09-23"))

    def test_a_list_that_is_there_but_unusable_is_an_error_not_a_skip(self):
        """Rule 4: an empty list, or an entry with no letters, fails loudly."""
        for text in ("", "# only a comment\n", "Quillonfake\n0022\n"):
            self.private_list(text)
            with self.assertRaises(ValueError, msg=repr(text)):
                rules.client_list()

    def test_the_default_list_is_outside_the_repository(self):
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop(rules.PRIVATE_NAMES_ENV, None)
            path = rules.private_list_path()
        self.assertEqual(path, os.path.expanduser(
            os.path.join("~", ".config", "ucdavis-proteomics", "client_names.txt")))
        self.assertFalse(os.path.abspath(path).startswith(os.path.abspath(SKILL) + os.sep))

    def test_no_digest_or_salt_of_a_client_name_in_the_public_files(self):
        """C2 security B1: salted digests of names were reversed with a word list in seconds, so
        the guard holds no form of a client name: no digest, no salt, no hash of a name. The one
        hash left is baseline_key, of a (path, line) already in the tree, never of a name."""
        for rel in (os.path.join("tests", "identifier_rules.py"),
                    os.path.join("tests", "test_no_new_identifiers.py")):
            with open(os.path.join(SKILL, rel), encoding="utf-8") as fh:
                src = fh.read()
            self.assertIsNone(re.search(r"['\"][0-9a-f]{16,}['\"]", src), rel)
            for word in ("CLIENT_" + "DIGESTS", "name_" + "digest", "_" + "SALT"):
                self.assertNotIn(word, src, rel)
        with open(os.path.join(HERE, "identifier_rules.py"), encoding="utf-8") as fh:
            src = fh.read()
        a = src.index("def baseline_key(")
        key_fn = src[a:src.index("\ndef ", a + 10)]
        self.assertEqual(src.count("sha" + "256("), key_fn.count("sha" + "256("), "a hash outside "
                         "baseline_key")
        self.assertEqual(src.count("hash" + "lib."), key_fn.count("hash" + "lib."))

    def test_the_allowlist_is_exact(self):
        rel, line = rules.ALLOWED[0]
        self.assertEqual(rules.problems(rel, line), [])
        self.assertTrue(rules.problems(rel, line + " # and more"))
        for rel, line in rules.ALLOWED:
            with open(os.path.join(SKILL, rel), encoding="utf-8") as fh:
                self.assertIn(line, [x.strip() for x in fh], f"{rel}: an allowlist entry that is gone")

    def test_the_fragpipe_licensor_is_named_where_a_user_is_sent_to_it(self):
        """C2 docs N2 (Ruling 82): the scrub had replaced the commercial licensor's name and
        contact with "the tools' commercial licensor", so "tell them to ask ..." named no one.
        It is the tools' public licence text, so it is back, in these four places only."""
        lic = "frag" + "matics"
        want = {"SKILL.md": (f"license via {lic}.com", f"send them to info@{lic}.com"),
                os.path.join("references", "fragpipe-diatracer.md"): (f"commercial via {lic}.com",),
                os.path.join("references", "search-engines.md"): (
                    f"({lic}.com, info@{lic}.com)", f"tools via {lic.title()}.",
                    f"re-confirmed with {lic.title()}.", f"info@{lic}.com and their own licensing")}
        for rel, phrases in want.items():
            with open(os.path.join(SKILL, rel), encoding="utf-8") as fh:
                text = " ".join(fh.read().split())
            for p in phrases:
                self.assertIn(p, text, rel)
            self.assertNotIn("the tools' commercial licensor", text, rel)

    def test_the_licensor_is_allowed_only_in_its_licence_text(self):
        """PUBLIC_TEXT takes out the exact phrases, in their own files: every line that names
        the licensor there passes, the phrase elsewhere still hits, and anything else on a
        licence line is still checked."""
        lic = "frag" + "matics"
        seen = 0
        for rel in sorted({r for r, _p, _w in rules.PUBLIC_TEXT}):
            with open(os.path.join(SKILL, rel), encoding="utf-8") as fh:
                for line in fh:
                    if lic in line.lower():
                        seen += 1
                        self.assertEqual(rules.line_hits(rel, line), [], f"{rel}: {line.strip()}")
        self.assertEqual(seen, 7)
        for rel, phrase, why in rules.PUBLIC_TEXT:
            self.assertIn("licence", why)
            self.assertTrue(rules.line_hits(os.path.join("scripts", "x.py"), f"# {phrase}"), phrase)
            self.assertTrue(rules.line_hits(rel, f"{phrase} -- " + "brett" + "sp"), (rel, phrase))
            self.assertTrue(rules.line_hits(rel, f"{phrase} PROT_" + "0756"), (rel, phrase))
            self.assertTrue(rules.line_hits(rel, f"ask {lic.title()} Inc"), (rel, "bare name"))


class NoNewIdentifiers(unittest.TestCase):
    def test_the_baseline_holds_no_identifier_text(self):
        """Keys and counts only: no line of the file is a staff or submission-number hit, its
        text outside the keys holds no client name, and every entry is a 64-hex key."""
        path = os.path.join(SKILL, BASELINE_REL)
        with open(path, encoding="utf-8") as fh:
            raw = fh.read()
        for n, line in enumerate(raw.splitlines(), 1):
            self.assertEqual(rules.problems(BASELINE_REL, line), [], f"line {n}")
        self.assertEqual(rules.client_names(KEY.sub(" ", raw)), [])
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        self.assertEqual(set(data), {"about", "names", "lines"})
        self.assertIn("baseline_key", data["about"])
        for key in data["names"]:
            self.assertRegex(key, KEY.pattern + "$")
        for entry in data["lines"]:
            self.assertEqual(len(entry), 2, entry)
            self.assertRegex(entry[0], KEY.pattern + "$")
            self.assertIsInstance(entry[1], int)
            self.assertGreater(entry[1], 0)
        self.assertEqual(len({k for k, _n in data["lines"]}), len(data["lines"]), "a key twice")

    def test_the_baseline_only_shrinks(self):
        lines, names = load_baseline()
        self.assertLessEqual(len(lines), FROZEN_LINES)
        self.assertLessEqual(sum(lines.values()), FROZEN_OCCURRENCES)
        self.assertLessEqual(len(names), FROZEN_NAMES)

    def test_a_new_identifier_is_still_found(self):
        """A tree of its own: a new staff login, submission number or named file is found; the
        same line excused by its key (and count) passes; a client name is found even when its
        line's key is in the baseline."""
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        os.makedirs(os.path.join(d, "scripts"))
        sub = "# from PROT_" + "0756, see the notes"
        login = "# asked by brett" + "sp"
        client = "# for the Quillonfake lab"
        with open(os.path.join(d, "scripts", "x.py"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(("import os", sub, login, sub, client, "")))
        named = "notes_brett" + "sp.md"
        open(os.path.join(d, named), "w").close()
        private = {"quillonfake": 1}
        none = (collections.Counter(), set())
        found = new_identifiers(d, none, private)
        self.assertEqual(len(found), 5, found)            # 2 x sub, login, client, file name
        self.assertTrue(any("in the file name" in f for f in found), found)
        k = rules.baseline_key
        rel = os.path.join("scripts", "x.py")
        once = (collections.Counter({k(rel, sub): 1, k(rel, login): 1, k(rel, client): 1}),
                {k(named, named)})
        found = new_identifiers(d, once, private)
        self.assertEqual(found, [f"{rel}:4: submission number 'PROT_" + "0756'",
                                 f"{rel}:5: a client name (private list entry 1)"])
        twice = (collections.Counter({k(rel, sub): 2, k(rel, login): 1}), {k(named, named)})
        self.assertEqual(new_identifiers(d, twice, {}), [])
        self.assertEqual(k(rel, "  " + sub + "  "), k(rel, sub))          # the line, stripped
        self.assertEqual(k(rel, sub), k("scripts/x.py", sub))            # "/" in the path
        self.assertNotEqual(k(rel, sub), k(os.path.join("tests", "x.py"), sub))

    def test_the_client_name_check_ran(self):
        """The client-name part needs the Core's private list; without it this test is SKIPPED
        with the reason, so a run that never checked a client name says so."""
        names, why = rules.client_list()
        if names is None:
            self.skipTest(why)
        self.assertGreater(len(names), 0, why)

    def test_no_identifier_beyond_the_baseline(self):
        new = new_identifiers()
        self.assertEqual(new, [], "the repository is public: write 'a Core submission' or a fake "
                                  "PROT_9xxx, `<your HIVE user>` or a generic login, 'a Core member', "
                                  "and 'a client lab' or the study for a client (tests/identifier_rules.py)")


if __name__ == "__main__":
    unittest.main()
