#!/usr/bin/env python3
"""No real client or staff identifier in the published skill (2.13, the skill slim-down).

The repository is public. The skill used to cite real CoreOmics submissions by number, a client
company and lab by name, and Core staff by their HIVE logins and first names -- in references,
script comments and test fixtures. What may appear instead:

* submission numbers: PROT_0807 (the one public example, also as 807 / PROT_807), PROT_08070
  (the malformed test value), PROT_0000 (the placeholder) and the reserved fake ranges
  PROT_9000-PROT_9999 and PROT_90000-PROT_99999, in any spelling (PROT_, PROT-, PROT, prot_ ...);
  "submission <n>" follows the same rule. In docs, say "a submission" or "an earlier submission"
  when a second one is meant;
* people: `<your HIVE user>` in docs, generic logins in tests (bstaff, mstaff, gstaff, staff1,
  tester, ...), "a Core member" in comments.

A line a script needs as it is (a real admin login the code reads) goes into ALLOWED with the
reason. The patterns are built from parts so this file does not match itself.
"""
import os
import re
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)

TEXT = (".py", ".md", ".R", ".sh", ".json", ".txt", ".yaml", ".yml", ".tsv", ".csv", ".html",
        ".cfg", ".toml", ".js", ".css")

# A submission number written with its PROT prefix (any case, separator or none), and one
# introduced by the word "submission" -- the two forms core_submission.py's identify reads.
PROT = re.compile(r"(?<![A-Za-z])" + "pr" + r"ot[ _#-]?([0-9]{3,5})(?![0-9])", re.I)
SUBMISSION = re.compile(r"\b" + "submis" + r"sion\s*(?:number|no\.?)?\s*#?\s*([0-9]{3,5})(?![0-9])", re.I)
FAKE_NUMBER = re.compile(r"0807|807|08070|0000|9[0-9]{3}|9[0-9]{4}")
# The podcast fixtures' study is synthetic and says so ("no such study"); its number stays.
SYNTHETIC = {(os.path.join("tests", "fixtures", "podcast"), "04" + "21")}

# Staff logins and names, the client company, a client lab and dataset, a person outside the
# Core, and the real CoreOmics id of a client submission.
NAMES = re.compile("|".join((
    "brett" + "sp", "msal" + "emi", "sal" + "emi", r"(?<![a-z])gab" + "rig", "gabri" + "ela",
    "mich" + "elle", "frag" + "matics", "dick" + "son", "sil" + r"va(?=[0-9_])", "adam " + "schaal",
    "0022066" + "cd85f")), re.I)

# (path under the skill, the line exactly, why it must stay)
ALLOWED = (
    # core_admins.txt is the trust anchor: skill_version.py/.sh, notes.py, slack_collab.py and
    # make_podcast.py trust only files owned by a HIVE account listed here, so the real admin
    # login is data the code reads (tests pin it by hash, not by name).
    (os.path.join("scripts", "core_admins.txt"), "brett" + "sp"),
)


def _number_ok(rel, digits):
    if FAKE_NUMBER.fullmatch(digits):
        return True
    return any(rel.startswith(folder + os.sep) and digits == n for folder, n in SYNTHETIC)


def problems(rel, line):
    """What is wrong with one line of the file at `rel` (a path under the skill): [] when none."""
    found = []
    for rx in (PROT, SUBMISSION):
        for m in rx.finditer(line):
            if not _number_ok(rel, m.group(1)):
                found.append(f"submission number {m.group(0)!r}")
    if NAMES.search(line) and (rel, line.strip()) not in ALLOWED:
        found.append(f"a real name or login {NAMES.search(line).group(0)!r}")
    return found


class NoRealIdentifiers(unittest.TestCase):
    def test_the_patterns_catch_the_real_shapes(self):
        real_numbers = ["PROT_" + "07" + "56", "PROT" + "0793", "prot-" + "756", "prot_" + "0803",
                        "PROT-" + "0812", "coreomics_prot" + "0756", "PROT_" + "10234",
                        "PROT " + "0801", "PROT#" + "0002", "submission " + "756",
                        "submission #" + "0793", "Submission number " + "0380"]
        for s in real_numbers:
            self.assertTrue(problems("references/x.md", f"see {s} here"), s)
        names = ["brett" + "sp", "AD3+" + "gab" + "rig", "C:\\Users\\" + "msal" + "emi",
                 "20260929T231000Z_" + "brett" + "sp" + "_x.md", "x.read." + "msal" + "emi",
                 "info@" + "frag" + "matics.com", "Mich" + "elle's Claude", "Gabri" + "ela",
                 "the " + "Dick" + "son lab", "Sil" + "va08172026", "Adam " + "Schaal",
                 "MSAL" + "EMI", "0022066" + "cd85f"]
        for s in names:
            self.assertTrue(problems("scripts/x.py", f"# {s} 2026-09-23"), s)

    def test_the_allowed_fakes_pass(self):
        for s in ("PROT_0807", "PROT0807", "prot-807", "PROT_807", "PROT 807", "PROT_08070",
                  "PROT_0000", "prot-0000", "PROT_9001", "PROT_9999", "prot9756", "PROT_91234",
                  "PROT_99999", "submission 807", "submission #9756", "Total_prot_10ug",
                  "WT_prot1", "Uniprot_1234", "PROT_123456", "PROT_####",
                  "bstaff", "mstaff", "gstaff", "staff1", "tester", "bsphinney@ucdavis.edu",
                  "Silva Campitelli de Barros BC", "amschaal/coreomics_fs"):
            self.assertEqual(problems("tests/x.py", f"x {s} y"), [], s)
        podcast = os.path.join("tests", "fixtures", "podcast", "x.md")
        self.assertEqual(problems(podcast, "Submission " + "04" + "21."), [])
        self.assertTrue(problems("references/x.md", "Submission " + "04" + "21."))

    def test_the_allowlist_is_exact(self):
        rel, line = ALLOWED[0]
        self.assertEqual(problems(rel, line), [])
        self.assertTrue(problems(rel, line + " # and more"))
        self.assertTrue(problems(os.path.join("scripts", "other.txt"), line))
        for rel, line in ALLOWED:
            with open(os.path.join(SKILL, rel), encoding="utf-8") as fh:
                self.assertIn(line, [x.strip() for x in fh], f"{rel}: an allowlist entry that is gone")

    def test_none_in_the_skill(self):
        hits = []
        for dp, dns, fns in os.walk(SKILL):
            dns[:] = [d for d in dns if d not in ("__pycache__", ".pytest_cache")]
            for name in dns + fns:
                rel = os.path.relpath(os.path.join(dp, name), SKILL)
                if PROT.search(name) and not all(_number_ok(rel, m.group(1)) for m in PROT.finditer(name)):
                    hits.append(f"{rel}: a submission number in the name")
                if NAMES.search(name):
                    hits.append(f"{rel}: a real name or login in the name")
            for fn in fns:
                path = os.path.join(dp, fn)
                rel = os.path.relpath(path, SKILL)
                with open(path, "rb") as fh:
                    data = fh.read()
                if not fn.endswith(TEXT) and b"\0" in data[:8192]:
                    continue                      # binary; any text file is read, named or not
                for i, line in enumerate(data.decode("utf-8", "replace").splitlines(), 1):
                    for what in problems(rel, line):
                        hits.append(f"{rel}:{i}: {what}")
        self.assertEqual(hits, [], "the repository is public: use PROT_0807 or a fake number "
                                   "(PROT_9000-9999), and `<your HIVE user>` or a generic login")


if __name__ == "__main__":
    unittest.main()
