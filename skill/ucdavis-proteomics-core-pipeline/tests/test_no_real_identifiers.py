#!/usr/bin/env python3
"""No real client or staff identifier in the published skill (2.13 the skill slim-down; 2.14.2
the C2 scrub; 2.15 one guard).

The repository is public. The skill used to cite real CoreOmics submissions by number, a client
company and lab by name, and Core staff by their HIVE logins and first names -- in references,
script comments and test fixtures. The rule is tests/identifier_rules.py, the one definition
(this file imports it; 2.15 review LOW 1: the two guards had kept nine copied definitions):
submission numbers in every spelling, the fixed names, and two PRIVATE lists that never enter
the repository -- the client names and the Core's staff list. What may appear instead is in its
docstring: PROT_0807 or a fake PROT_9xxx, `<your HIVE user>` or a generic login (staff1,
tester), "a Core member", "a client lab" or the study.

It reads the whole skill (every text line, every file and folder name) and what the new
repository ships outside it (identifier_rules.shipped_outside_skill: marketplace.json, the docs
pages, the CHANGELOG's skill sections, and $SKILL_SCAN_EXTRA -- the README draft at a merge).

The private lists, and SKILL_REQUIRE_CLIENT_NAMES:
  * the client names: $SKILL_PRIVATE_NAMES, else ~/.config/ucdavis-proteomics/client_names.txt;
  * the staff: $SKILL_PRIVATE_STAFF, else ~/.config/ucdavis-proteomics/core_staff.txt (the Core's
    staff file, scripts/staff.py's format, with the names after each login);
  * without a list its part is SKIPPED, and test_the_client_name_check_ran /
    test_the_staff_name_check_ran say so -- CI and HIVE never have them (Ruling 93: they never
    leave the Core director's Mac; no GitHub secret, no HIVE copy);
  * SKILL_REQUIRE_CLIENT_NAMES=1 turns that skip into a FAILURE, for both lists. The controller
    sets it on the Mac at every merge into the slim or release branch and before the release
    push, so a merge is never checked without them.

A line a script needs as it is (a real admin login the code reads) goes into ALLOWED with the
reason. The patterns are built from parts so this file does not match itself.
"""
import os
import re
import shutil
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(SKILL))        # <repo>/skill/<the skill>
sys.path.insert(0, HERE)
import identifier_rules as rules  # noqa: E402  the one definition
from identifier_rules import (TEXT, PROT, SUBMISSION, FAKE_NUMBER, SYNTHETIC, NAMES,  # noqa: E402,F401
                              ALLOWED, _number_ok, problems)

FIX = ("the repository is public: write 'a Core submission' or a fake PROT_9xxx, `<your HIVE "
       "user>` or a generic login, 'a Core member', and 'a client lab' or the study for a client "
       "(tests/identifier_rules.py)")


def private(kind):
    """(list, message) of the private `kind` list ("client" or "staff")."""
    return rules.client_list() if kind == "client" else rules.staff_list()


def list_or_skip(test, kind):
    """The private list, or: SKIP with the reason -- FAIL instead under
    SKILL_REQUIRE_CLIENT_NAMES=1 (identifier_rules.REQUIRE_ENV)."""
    lst, why = private(kind)
    if lst is None:
        if rules.required():
            test.fail(f"{why} -- and {rules.REQUIRE_ENV}=1 requires it")
        test.skipTest(why)
    return lst


class TheRules(unittest.TestCase):
    def test_the_patterns_catch_the_real_shapes(self):
        n = "08" + "12"                      # a real number, in parts so this file is clean
        real_numbers = ["PROT_" + "07" + "56", "PROT" + "07" + "93", "prot-" + "756",
                        "prot_" + "08" + "03", "PROT-" + n, "coreomics_prot" + "07" + "56",
                        "PROT_" + "10234", "PROT " + "08" + "01", "PROT#" + "00" + "02",
                        "submission " + "756", "submission #" + "07" + "93",
                        "Submission number " + "03" + "80",
                        # 2.15 review LOW 15: the spellings that slipped through
                        "PROT #" + n, "PROT." + n, "PROT\u2013" + n, "PROT \u2014 " + n,
                        "submission: " + "812", "submission ID " + n,
                        "the " + n + " submission", "Core" + "Omics " + n,
                        "submissions " + "07" + "56 and " + "07" + "93"]
        for s in real_numbers:
            self.assertTrue(problems("references/x.md", f"see {s} here"), s)
        two = problems("references/x.md", "submissions " + "07" + "56 and " + "07" + "93")
        self.assertEqual(len(two), 2, two)
        names = ["brett" + "sp", "AD3+" + "gab" + "rig", "C:\\Users\\" + "msal" + "emi",
                 "20260929T231000Z_" + "brett" + "sp" + "_x.md", "x.read." + "msal" + "emi",
                 "info@" + "frag" + "matics.com", "Mich" + "elle's Claude", "Gabri" + "ela",
                 "the " + "Dick" + "son lab", "".join(reversed(("va08172026", "Sil"))),
                 " ".join(reversed(("Schaal", "Adam"))), "MSAL" + "EMI", "0022066" + "cd85f",
                 "the " + "Dupan" + "loup search", "HIVE, " + "ga" + "rg cohort"]
        for s in names:
            self.assertTrue(problems("scripts/x.py", f"# {s} 2026-09-23"), s)

    def test_the_allowed_fakes_pass(self):
        for s in ("PROT_0807", "PROT0807", "prot-807", "PROT_807", "PROT 807", "PROT_08070",
                  "PROT_0000", "prot-0000", "PROT_9001", "PROT_9999", "prot9756", "PROT_91234",
                  "PROT_99999", "submission 807", "submission #9756", "Total_prot_10ug",
                  "WT_prot1", "Uniprot_1234", "PROT_123456", "PROT_####",
                  "PROT #0807", "PROT.9756", "PROT\u20139001", "submission: 9756",
                  "the 0807 submission", "CoreOmics 9756",
                  "bstaff", "mstaff", "gstaff", "staff1", "tester", "bsphinney@ucdavis.edu",
                  "Silva Campitelli de Barros BC", "amschaal/coreomics_fs"):
            self.assertEqual(problems("tests/x.py", f"x {s} y"), [], s)
        podcast = os.path.join("tests", "fixtures", "podcast", "x.md")
        self.assertEqual(problems(podcast, "Submission " + "04" + "21."), [])
        self.assertTrue(problems("references/x.md", "Submission " + "04" + "21."))

    def test_a_bare_number_counts_only_beside_a_submission_word(self):
        """LOW 15: a bare 4-digit number is a submission number only with a leading 0 and at most
        NEAR_WORDS word between it and "submission(s)" / "CoreOmics" (a list's other numbers do
        not count) -- years, dates, paths and counts beside those words are not (any 4-digit
        number within 3 words: 32 lines of the 2.15 tree, none a submission)."""
        for s in ("a 2026 submission", "submissions from 2025-10-01", "1060: Several submissions",
                  "coreomics/projects/2026/08/<id>/share", "umask 0002), so the submission's",
                  "0812 runs", "PROT_9756, run 0812 of the plate"):
            self.assertEqual(problems("references/x.md", s), [], s)
        self.assertEqual(rules.NEAR_WORDS, 1)
        self.assertTrue(problems("references/x.md", "submission for " + "08" + "12"))

    def test_the_allowlist_is_exact(self):
        rel, line = ALLOWED[0]
        self.assertEqual(problems(rel, line), [])
        self.assertTrue(problems(rel, line + " # and more"))
        self.assertTrue(problems(os.path.join("scripts", "other.txt"), line))
        for rel, line in ALLOWED:
            with open(os.path.join(SKILL, rel), encoding="utf-8") as fh:
                self.assertIn(line, [x.strip() for x in fh], f"{rel}: an allowlist entry that is gone")

    def test_one_definition(self):
        """LOW 1: the names this guard uses are identifier_rules' own objects, not copies."""
        g = globals()
        for name in ("TEXT", "PROT", "SUBMISSION", "FAKE_NUMBER", "SYNTHETIC", "NAMES", "ALLOWED",
                     "_number_ok", "problems"):
            self.assertIs(g[name], getattr(rules, name), name)
        with open(__file__, encoding="utf-8") as fh:
            src = fh.read()
        self.assertNotIn("re." + "compile(", src.split("class TheRules")[0])
        self.assertFalse(os.path.exists(os.path.join(HERE, "test_no_new_" + "identifiers.py")))
        self.assertFalse(os.path.exists(os.path.join(HERE, "identifier_" + "baseline.json")))


class PrivateLists(unittest.TestCase):
    def private_file(self, env, text):
        """A private list holding `text`, as $<env> (removed after the test)."""
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        path = os.path.join(d, "list.txt")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
        patcher = mock.patch.dict(os.environ, {env: path})
        patcher.start()
        self.addCleanup(patcher.stop)
        return path

    def test_a_client_name_is_found_however_it_is_written(self):
        """The private list's names match one word or several joined; the name is never in
        this repository, and a hit names the entry number, not the name."""
        self.private_file(rules.PRIVATE_NAMES_ENV,
                          "# fake names, for this test only\nQuillonfake\n\nVan Quill  # two words\n")
        names, why = rules.client_list()
        self.assertEqual(len(names), 2, why)
        for s in ("the Quillonfake lab", "QUILLONFAKE", "2026-10-01_Quillonfake_Bat_HT",
                  "QuillonfakeBat", "Quillonfake08172026", "Qüillonfake", "Van Quill",
                  "VanQuill", "van-quill", "VANQUILL"):
            self.assertTrue(rules.line_hits("references/x.md", s, staff={}), s)
        for s in ("Quillonfakery", "unquillonfake", "Van", "Quill", "Vanquillon"):
            self.assertEqual(rules.client_names(s), [], s)
        hit = rules.line_hits("references/x.md", "ask the Van Quill lab", staff={})
        self.assertEqual(hit, ["a client name (private list entry 4)"])
        self.assertNotIn("quill", " ".join(hit).lower())

    def test_a_staff_login_or_name_is_found_by_entry_number(self):
        """M8: the Core's staff list, private like the client names -- a login as a whole word,
        each name after it as a client name is matched; a hit names the entry, never the person."""
        self.private_file(rules.PRIVATE_STAFF_ENV,
                          "# fake staff, for this test only\nqstafffake   Quilla, Quilla Fakestaff\n"
                          "zz.fake2\n")
        staff, why = rules.staff_list()
        self.assertEqual((len(staff["logins"]), len(staff["names"])), (2, 2), why)
        for s in ("asked by qstafffake", "C:\\Users\\QSTAFFFAKE", "~zz.fake2/x", "Quilla's Claude",
                  "Quilla Fakestaff", "QuillaFakestaff"):
            self.assertTrue(rules.line_hits("scripts/x.py", f"# {s}", names={}), s)
        for s in ("qstafffakes", "xqstafffake", "Quillan", "zz.fake23"):
            self.assertEqual(rules.staff_names("scripts/x.py", s), [], s)
        hit = rules.line_hits("scripts/x.py", "# from Quilla Fakestaff", names={})
        self.assertEqual(hit, ["a Core staff login or name (private staff list entry 2)"])
        self.assertNotIn("quilla", " ".join(hit).lower())

    def test_without_a_private_list_only_its_part_is_skipped(self):
        """Public CI has no list: its part says it was skipped and where the list goes; the
        submission-number and fixed-name checks still run."""
        gone = os.path.join(tempfile.gettempdir(), "no-such-dir-c2", "x.txt")
        with mock.patch.dict(os.environ, {rules.PRIVATE_NAMES_ENV: gone,
                                          rules.PRIVATE_STAFF_ENV: gone}):
            for kind, env, skipped in (("client", rules.PRIVATE_NAMES_ENV, rules.SKIPPED),
                                       ("staff", rules.PRIVATE_STAFF_ENV, rules.STAFF_SKIPPED)):
                lst, why = private(kind)
                self.assertIsNone(lst, kind)
                self.assertIn(skipped, why)
                self.assertIn("no-such-dir-c2", why)
                self.assertIn(env, why)
            self.assertEqual(rules.client_names("the Quillonfake lab"), [])
            self.assertEqual(rules.staff_names("x", "qstafffake"), [])
            self.assertTrue(rules.line_hits("scripts/x.py", "# PROT_" + "07" + "56 2026-09-23"))
            self.assertTrue(rules.line_hits("scripts/x.py", "# brett" + "sp 2026-09-23"))

    def test_require_turns_the_skip_into_a_failure(self):
        """Ruling 93: SKILL_REQUIRE_CLIENT_NAMES=1 -- a missing list fails, never skips."""
        gone = os.path.join(tempfile.gettempdir(), "no-such-dir-c2", "x.txt")
        for kind, env in (("client", rules.PRIVATE_NAMES_ENV), ("staff", rules.PRIVATE_STAFF_ENV)):
            with mock.patch.dict(os.environ, {env: gone, rules.REQUIRE_ENV: "1"}):
                with self.assertRaises(AssertionError, msg=kind) as cm:
                    list_or_skip(self, kind)
                self.assertIn(rules.REQUIRE_ENV, str(cm.exception))
            with mock.patch.dict(os.environ, {env: gone, rules.REQUIRE_ENV: "0"}):
                with self.assertRaises(unittest.SkipTest, msg=kind):
                    list_or_skip(self, kind)

    def test_a_list_that_is_there_but_unusable_is_an_error_not_a_skip(self):
        """Rule 4: an empty list, or an entry with no letters, fails loudly."""
        for text in ("", "# only a comment\n", "Quillonfake\n0022\n"):
            self.private_file(rules.PRIVATE_NAMES_ENV, text)
            with self.assertRaises(ValueError, msg=repr(text)):
                rules.client_list()
        for text in ("", "# only a comment\n", "Doe, Jane\n", "qstafffake  1234\n"):
            self.private_file(rules.PRIVATE_STAFF_ENV, text)
            with self.assertRaises(ValueError, msg=repr(text)):
                rules.staff_list()

    def test_the_default_lists_are_outside_the_repository(self):
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop(rules.PRIVATE_NAMES_ENV, None)
            os.environ.pop(rules.PRIVATE_STAFF_ENV, None)
            paths = (rules.private_list_path(), rules.private_staff_path())
        self.assertEqual(paths, tuple(os.path.expanduser(os.path.join(
            "~", ".config", "ucdavis-proteomics", n)) for n in ("client_names.txt",
                                                                 "core_staff.txt")))
        for path in paths:
            self.assertFalse(os.path.abspath(path).startswith(os.path.abspath(REPO) + os.sep))

    def test_no_digest_or_salt_of_a_name_in_the_public_files(self):
        """C2 security B1: salted digests of names were reversed with a word list in seconds, so
        the guard holds no form of a name: no digest, no salt, no hash."""
        for rel in (os.path.join("tests", "identifier_rules.py"),
                    os.path.join("tests", "test_no_real_identifiers.py")):
            with open(os.path.join(SKILL, rel), encoding="utf-8") as fh:
                src = fh.read()
            self.assertIsNone(re.search(r"['\"][0-9a-f]{16,}['\"]", src), rel)
            for word in ("CLIENT_" + "DIGESTS", "name_" + "digest", "_" + "SALT", "sha" + "256(",
                         "hash" + "lib"):
                self.assertNotIn(word, src, rel)

    def test_the_client_name_check_ran(self):
        """The client-name part needs the Core's private list; without it this test is SKIPPED
        with the reason (FAILS under SKILL_REQUIRE_CLIENT_NAMES=1), so a run that never checked a
        client name says so."""
        names = list_or_skip(self, "client")
        self.assertGreater(len(names), 0)

    def test_the_staff_name_check_ran(self):
        """The same for the Core's staff list (M8)."""
        staff = list_or_skip(self, "staff")
        self.assertGreater(len(staff["logins"]), 0)


class TheLicensor(unittest.TestCase):
    def test_the_fragpipe_licensor_is_named_where_a_user_is_sent_to_it(self):
        """C2 docs N2 (Ruling 82): the scrub had replaced the commercial licensor's name and
        contact with "the tools' commercial licensor", so "tell them to ask ..." named no one.
        It is the tools' public licence text, so it is back, in these places only."""
        lic = "frag" + "matics"
        # 2.15: SKILL.md's two lines (old step 7) live in fragpipe-diatracer.md in the slim tree
        want = {os.path.join("references", "fragpipe-diatracer.md"): (
                    f"commercial via {lic}.com", f"license via {lic}.com",
                    f"send them to info@{lic}.com"),
                os.path.join("references", "search-engines.md"): (
                    f"({lic}.com, info@{lic}.com)", f"tools via {lic.title()}.",
                    f"re-confirmed with {lic.title()}.", f"info@{lic}.com and their own licensing")}
        for rel, phrases in want.items():
            with open(os.path.join(SKILL, rel), encoding="utf-8") as fh:
                text = " ".join(fh.read().split())
            for p in phrases:
                self.assertIn(p, text, rel)
            self.assertNotIn("the tools' commercial licensor", text, rel)
        with open(os.path.join(SKILL, "SKILL.md"), encoding="utf-8") as fh:
            self.assertNotIn("the tools' commercial licensor", fh.read())

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
                        self.assertEqual(rules.line_hits(rel, line, {}, {}), [],
                                         f"{rel}: {line.strip()}")
        self.assertEqual(seen, 7)
        for rel, phrase, why in rules.PUBLIC_TEXT:
            self.assertIn("licence", why)
            self.assertTrue(rules.line_hits(os.path.join("scripts", "x.py"), f"# {phrase}", {}, {}),
                            phrase)
            self.assertTrue(rules.line_hits(rel, f"{phrase} -- " + "brett" + "sp", {}, {}),
                            (rel, phrase))
            self.assertTrue(rules.line_hits(rel, f"{phrase} PROT_" + "07" + "56", {}, {}),
                            (rel, phrase))
            self.assertTrue(rules.line_hits(rel, f"ask {lic.title()} Inc", {}, {}),
                            (rel, "bare name"))


class OutsideTheSkill(unittest.TestCase):
    """LOW 15: what the new repository ships outside skill/ is public too."""

    def test_the_files_scanned(self):
        rels = [r for r, _p, _l in rules.shipped_outside_skill(REPO)]
        for rel in (".claude-plugin/marketplace.json", "CHANGELOG.md",
                    os.path.join("docs", "skill-install.html"),
                    os.path.join("docs", "how-reports-are-made.html")):
            if os.path.isfile(os.path.join(REPO, rel)):     # an installed plugin has no repo root
                self.assertIn(rel, rels)

    def test_a_data_uri_and_the_public_space_url_are_not_text(self):
        img = "data:image/png;base64,iVBORw0K" + "brett" + "spGgoAAAAN+/PROT_" + "07" + "56=="
        self.assertEqual(rules.outside_hits("docs/x.html", f'<img src="{img}"> ok', {}, {}), [])
        self.assertTrue(rules.outside_hits("docs/x.html", "<p>PROT_" + "07" + "56</p>", {}, {}))
        space = "https://huggingface.co/spaces/" + "brett" + "sp" + "/de-limp-proteomics"
        self.assertEqual(rules.outside_hits("README.md", f"[try it]({space})", {}, {}), [])
        self.assertTrue(rules.outside_hits("README.md", "ask " + "brett" + "sp", {}, {}))

    def test_the_extra_paths_are_scanned(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        draft = os.path.join(d, "README_draft.md")
        with open(draft, "w", encoding="utf-8") as fh:
            fh.write("# x\n\nfrom PROT #" + "08" + "12\n")
        with mock.patch.dict(os.environ, {rules.SCAN_EXTRA_ENV: draft}):
            found = [(rel, path) for rel, path, _l in rules.shipped_outside_skill(REPO)]
        self.assertIn((os.path.abspath(draft), draft), found)
        self.assertEqual(scan_outside(REPO, extra=draft, names={}, staff={}),
                         [f"{os.path.abspath(draft)}:3: submission number 'PROT #" + "08" + "12'"])


def scan_skill(skill=SKILL, names=None, staff=None):
    """Every hit under the skill: each text line and each file and folder name."""
    hits = []
    for rel, i, text in rules.walk(skill):
        where = f"{rel}:{i}" if i else f"{rel}: in the file name"
        for what in rules.line_hits(rel, text, names, staff):
            hits.append(f"{where}: {what}")
    return hits


def scan_outside(repo=REPO, extra=None, names=None, staff=None):
    """Every hit in what the new repository ships outside the skill (shipped_outside_skill)."""
    hits = []
    env = {rules.SCAN_EXTRA_ENV: extra} if extra else {}
    with mock.patch.dict(os.environ, env):
        files = rules.shipped_outside_skill(repo)
    for rel, path, lines in files:
        text = rules._read_text(path)
        if text is None:
            continue
        for i, line in enumerate(text.splitlines(), 1):
            if lines is None or i in lines:
                for what in rules.outside_hits(rel, line, names, staff):
                    hits.append(f"{rel}:{i}: {what}")
    return hits


class NoneInTheRepository(unittest.TestCase):
    def test_none_in_the_skill(self):
        names, staff = private("client")[0] or {}, private("staff")[0] or {}
        self.assertEqual(scan_skill(names=names, staff=staff), [], FIX)

    def test_none_in_what_ships_outside_the_skill(self):
        names, staff = private("client")[0] or {}, private("staff")[0] or {}
        self.assertEqual(scan_outside(names=names, staff=staff), [], FIX)


if __name__ == "__main__":
    unittest.main()
