#!/usr/bin/env python3
"""
The client report's "Audit & caveats" is a plain-language table, one row per check: what was
checked, the result, what it means for the results. The technical AUDIT.md stays the Core's
record (audit_results.client_section_md).

A staff report (PROT_0000, skill 2.11.2): the report's Audit & caveats section was AUDIT.md
as written -- check codes, the normalisation check's B1-B4 strings, Normalisation.Instability
values, "implausibly high ... not real biology". Neither the staff member nor a client could
read it, and a client could take "not real biology" to mean the results are wrong. The plain
words are written where each finding is made (audit_results.py; a pipeline's own checks carry
theirs, run_silac.R), never by the model writing the report.

Stdlib only.
"""
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL, "scripts")
sys.path.insert(0, SCRIPTS)

import audit_results as ar  # noqa: E402

# A label-free DE as a staff member's AUDIT.json had it (names are placeholders).
FINDINGS = [
    {"check": "replication", "status": "PASS", "message": "All groups have ≥3 replicates."},
    {"check": "id_depth", "status": "PASS", "message": "3183 proteins quantified."},
    {"check": "de_signal", "status": "WARN", "message":
        "DE_dpc_A.B.csv: 3076/3183 (97%) proteins significant — implausibly high. This usually "
        "means a batch effect, a normalization problem, or confounded/mislabeled groups, not real "
        "biology.", "plain": ar.PLAIN_DE_SIGNAL_HIGH},
    {"check": "de_signal", "status": "WARN", "message":
        "DE_dpc_A.C.csv: 2937/3183 (92%) proteins significant — implausibly high.",
        "plain": ar.PLAIN_DE_SIGNAL_HIGH},
    {"check": "normalization_check", "status": "WARN", "message":
        "Experiment type whole-proteome (the user): the data check did not agree with its default "
        "(B1 normalisation factors 415.01-fold apart between groups (A vs B), separating "
        "completely; B4 Normalisation.Instability above 0.3 in 3 run(s)); Core staff chose "
        "non-normalised (raw) quantities because A is loaded low by design.",
        "plain": "The usual normalisation did not fit these samples, so Core staff chose which "
                 "quantities to use: Core staff chose non-normalised (raw) quantities because A is "
                 "loaded low by design."},
    {"check": "normalisation_stability", "status": "INFO", "message":
        "Not relevant: the DE read DIA-NN's non-normalised quantities (de_provenance.json "
        "quantities: raw)."},
]
TECHNICAL = ("B1 ", "B4 ", "Normalisation.Instability", "not real biology", "implausibly high",
             "de_provenance.json", "DE_dpc_", "adj.P", "de_signal", "normalization_check",
             "normalisation_stability", "id_depth", "415.01")


def table_rows(md):
    return [ln for ln in md.splitlines() if ln.startswith("|") and not set(ln) <= set("|-: ")]


class PlainTable(unittest.TestCase):
    def test_one_plain_row_per_check(self):
        md = ar.client_section_md(FINDINGS)
        rows = table_rows(md)
        self.assertEqual(rows[0], "| Check | Result | What it means for the results |")
        self.assertEqual(len(rows) - 1, 5, md)                 # de_signal twice -> one row
        for t in TECHNICAL:
            self.assertNotIn(t, md)
        self.assertIn("AUDIT.md", md)                          # where the technical record is
        self.assertIn("⚠️ Caution", md)
        self.assertIn("✅ Passed", md)
        self.assertIn("Core staff chose non-normalised (raw) quantities", md)

    def test_the_worst_result_wins_a_row(self):
        md = ar.client_section_md([
            {"check": "de_signal", "status": "PASS", "message": "x"},
            {"check": "de_signal", "status": "WARN", "message": "y",
             "plain": ar.PLAIN_DE_SIGNAL_NONE}])
        row = table_rows(md)[1]
        self.assertIn("⚠️ Caution", row)
        self.assertIn(ar.PLAIN_DE_SIGNAL_NONE, row)

    def test_a_blocking_problem_reads_as_one(self):
        md = ar.client_section_md([{"check": "replication", "status": "FAIL", "message": "x"}])
        self.assertIn("⛔ Problem", md)

    def test_every_check_and_result_the_audit_makes_has_plain_words(self):
        src = open(os.path.join(SCRIPTS, "audit_results.py"), encoding="utf-8").read()
        made = set(re.findall(r'add\(findings,\s*"([a-z_]+)",\s*"(PASS|WARN|FAIL|INFO)"', src))
        made |= {("contamination", "PASS"), ("contamination", "WARN"),   # status computed
                 ("sage_lfq", "PASS"), ("sage_lfq", "WARN"), ("sage_lfq", "INFO")}
        self.assertGreater(len(made), 30)
        for check, status in sorted(made):
            with self.subTest(check=check, status=status):
                label, meaning = ar.client_row({"check": check, "status": status, "message": "x"})
                self.assertNotIn("_", label, "a plain label, not the check id")
                self.assertTrue(meaning and meaning != ar.GENERIC_MEANING[status],
                                f"{check} {status} has no plain meaning of its own")

    def test_a_pipeline_check_with_no_plain_words_is_never_quoted(self):
        label, meaning = ar.client_row({"check": "newpipe_balance", "status": "WARN",
                                        "message": "raw technical text 0.37 log2"})
        self.assertEqual(label, "Newpipe balance")
        self.assertEqual(meaning, ar.GENERIC_MEANING["WARN"])
        self.assertNotIn("0.37", meaning)

    def test_an_old_audit_md_reads_the_same(self):
        md = ("# Results audit — common proteomics pitfalls\n\n**Overall: ⚠️ WARN** — 0 blocking, "
              "2 warning(s).\n\n- ✅ **replication** — All groups have ≥3 replicates.\n"
              "- ⚠️ **de_signal** — DE_dpc_A.B.csv: 3076/3183 (97%) proteins significant — "
              "implausibly high. This usually means a batch effect, a normalization problem, or "
              "confounded/mislabeled groups, not real biology.\n"
              "- ℹ️ **target_contaminants** — Not assessed: no <fasta>.meta.json.\n")
        found = ar.findings_from_md(md)
        self.assertEqual([(f["check"], f["status"]) for f in found],
                         [("replication", "PASS"), ("de_signal", "WARN"),
                          ("target_contaminants", "INFO")])
        out = ar.client_section_md(found)
        self.assertEqual(len(table_rows(out)) - 1, 3)
        self.assertNotIn("not real biology", out)
        self.assertIn(ar.PLAIN_DE_SIGNAL_GENERIC, out)          # no site words: the generic one


class PipelineChecksSpeakForThemselves(unittest.TestCase):
    def test_every_silac_check_carries_its_plain_words(self):
        src = open(os.path.join(SCRIPTS, "run_silac.R"), encoding="utf-8").read()
        calls = [m.start() for m in re.finditer(r'\bchk\("silac_', src)]
        self.assertGreater(len(calls), 15)
        for at in calls:
            # the call runs to the next chk( or the next section rule
            end = min([x for x in (src.find('chk("', at + 5), src.find("\n# ----", at)) if x > 0])
            with self.subTest(call=src[at:at + 40]):
                self.assertIn("plain =", src[at:end])
        self.assertIn("SILAC_CHECK_LABELS", src)

    def test_the_audit_passes_them_through(self):
        f = []
        ar.audit_pipeline_checks(f, {"checks": [
            {"check": "silac_replicates", "status": "WARN", "message": "KO.WT: 2 replicates",
             "plain": "A comparison has the minimum number of replicates.",
             "label": "SILAC replicates"}]})
        self.assertEqual(ar.client_row(f[0]), ("SILAC replicates",
                                               "A comparison has the minimum number of replicates."))


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


class TheReport(unittest.TestCase):
    def build(self, report=None, audit_json=True):
        self._tmp = tempfile.TemporaryDirectory()
        s = os.path.join(self._tmp.name, "session")
        o = os.path.join(s, "output")
        write(os.path.join(o, "AUDIT.md"), "# Results audit\n\n" + "".join(
            f"- ⚠️ **{f['check']}** — {f['message']}\n" for f in FINDINGS))
        if audit_json:
            write(os.path.join(o, "AUDIT.json"), json.dumps({"overall": "WARN", "findings": FINDINGS}))
        write(os.path.join(o, "AI_Analysis_Report.md"),
              report or "# A report\n\n## Key Findings\n\nThings changed.\n")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"),
                            "--session", s, "--out", os.path.join(o, "Analysis_Report.html"),
                            "--no-pdf"], capture_output=True, text=True, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(os.path.join(o, "Analysis_Report.html"), encoding="utf-8") as fh:
            page = fh.read()
        with open(os.path.join(o, "Analysis_Report.md"), encoding="utf-8") as fh:
            twin = fh.read()
        return page, twin

    def tearDown(self):
        if hasattr(self, "_tmp"):
            self._tmp.cleanup()

    def section(self, twin):
        m = re.search(r"## Audit & caveats\n(.*?)(?=\n## |\Z)", twin, re.S)
        self.assertIsNotNone(m, twin)
        return m.group(1)

    def test_the_client_reads_the_plain_table(self):
        page, twin = self.build()
        sec = self.section(twin)
        self.assertIn("| Check | Result | What it means for the results |", sec)
        for t in TECHNICAL:
            self.assertNotIn(t, sec)
        self.assertNotIn("not real biology", page)
        self.assertIn("What it means for the results", page)

    def test_an_older_session_with_only_audit_md(self):
        _page, twin = self.build(audit_json=False)
        sec = self.section(twin)
        self.assertIn("| Check | Result | What it means for the results |", sec)
        self.assertNotIn("not real biology", sec)

    def test_a_report_of_its_own_still_wins(self):
        _page, twin = self.build(report="# A report\n\n## Audit & caveats\n\n| Check | Result | "
                                        "What it means |\n|---|---|---|\n| Replicates | ok | fine |\n")
        sec = self.section(twin)
        self.assertIn("| Replicates | ok | fine |", sec)
        self.assertNotIn("What it means for the results", sec)
        self.assertEqual(twin.count("## Audit & caveats"), 1)


class Documented(unittest.TestCase):
    def test_the_writer_leaves_the_audit_section_to_the_tool(self):
        skill = " ".join(open(os.path.join(SKILL, "SKILL.md"), encoding="utf-8").read().split())
        self.assertNotIn('plus an "Audit & caveats" section from `AUDIT.md`', skill)
        self.assertIn("plain-language", skill)
        src = open(os.path.join(SCRIPTS, "analysis_prompt.py"), encoding="utf-8").read()
        self.assertIn("Audit & caveats", src)


if __name__ == "__main__":
    unittest.main(verbosity=2)
