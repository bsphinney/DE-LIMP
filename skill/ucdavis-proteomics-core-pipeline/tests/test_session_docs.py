#!/usr/bin/env python3
"""The session's README.html / README.md / AGENTS.md, and "where this lives on HIVE".

Brett, 2026-09-24: collaborators cannot read a .md, an AI agent handed the folder needs a guide
to it, and nothing said where the files and the raw data are on HIVE. A real hive_remote session
(Silva08172026) had a README that named no HIVE path and claimed an input/raw_files.txt that did
not exist -- the raw paths were only in output/search/file_list.txt and search_provenance.json.

What these tests pin:
  * README.html is valid, self-contained (inline CSS, no script, no external asset) and has the
    same sections as README.md -- both rendered from one source
  * AGENTS.md is generated from the records: a DPC and a MaxLFQ provenance each come out in
    their own words, and neither names the other pipeline
  * HIVE locations resolve for an R:-style drive, /Volumes/proteomics, a UNC path and /quobyte,
    through hive_shares.tsv (the table hive_path.sh reads too); unknown -> "not recorded"
  * input/raw_files.txt is written from the search's record when missing, and the README never
    claims it exists when it does not
  * finalize puts README.html + AGENTS.md in the zip and MANIFEST.txt, and rewrites them with the
    registry record's path when the run-log hook creates it
stdlib only, no network.
"""
import csv
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
import zipfile
from contextlib import redirect_stderr, redirect_stdout
from html.parser import HTMLParser
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import session                                      # noqa: E402
import session_docs                                 # noqa: E402
import share_map                                    # noqa: E402
import make_deposit                                 # noqa: E402
import test_deposit_package as tdp                  # noqa: E402
from job_env import job_env                         # noqa: E402

PY = sys.executable
DPC = dict(tdp.DE_PROV, significance_rule="adj.P.Val < adjp (BH); no fold-change filter",
           significant_per_contrast={"Treated-Control": 7}, groups={"Control": 2, "Treated": 1},
           n_samples=3)
# build_maxlfq.R's descriptor, as run_de.R writes it into de_provenance.json
MAXLFQ = {"pipeline_id": "maxlfq", "display_label": "MaxLFQ + limma",
          "rollup_method": "DIA-NN PG.MaxLFQ",
          "de_engine": "limma::lmFit -> contrasts.fit -> eBayes (NA-tolerant per row)",
          "missing_policy": "NAs left in place; limma drops them per row. All-missing-in-one-"
                            "condition proteins are on/off calls.",
          "citation": "Quantification: DIA-NN MaxLFQ (Demichev et al. 2020, Nat Methods 17:41). "
                      "DE: limma (Ritchie et al. 2015, NAR 43:e47).",
          "method": "maxlfq", "adjp": 0.05, "logfc": 1, "logfc_role": "reference_line_only",
          "significance_rule": "adj.P.Val < adjp (BH); no fold-change filter",
          "design": "~ 0 + groups + Batch", "contrasts": ["Treated-Control"],
          "significant_per_contrast": {"Treated-Control": 3}}
CORE_ONLY = ("Proteomics Core storage: only Core members can open it -- ask the Core for a copy "
             "or for access")
ROWS = [{"server": "128.120.208.24", "share": "proteomics",
         "hive": "/nfs/lssc0/flinders/proteomics", "mac": "/Volumes/proteomics",
         "windows": "\\\\128.120.208.24\\proteomics", "access": None},
        {"server": "*", "share": "proteomics-grp", "hive": "/quobyte/proteomics-grp",
         "mac": "/Volumes/proteomics-grp", "windows": None, "access": CORE_ONLY}]
FLINDERS = "/nfs/lssc0/flinders/proteomics"


def write(path, text):
    tdp.write(path, text)


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def sections_md(text):
    return [ln[3:].strip() for ln in text.splitlines() if ln.startswith("## ")]


class H2s(HTMLParser):
    """Collects <h2> texts and checks the page is well-formed and self-contained."""

    def __init__(self):
        super().__init__()
        self.stack, self.bad, self.h2, self.cur, self.external, self.tags = [], 0, [], None, [], set()

    def handle_starttag(self, tag, attrs):
        self.tags.add(tag)
        a = dict(attrs)
        if tag in ("script", "link", "img", "iframe") or (a.get("src") or "").startswith("http"):
            self.external.append((tag, a))
        if tag not in ("meta", "br"):
            self.stack.append(tag)
        if tag == "h2":
            self.cur = ""

    def handle_endtag(self, tag):
        if self.stack and self.stack[-1] == tag:
            self.stack.pop()
        else:
            self.bad += 1
        if tag == "h2":
            self.h2.append(self.cur.strip())
            self.cur = None

    def handle_data(self, data):
        if self.cur is not None:
            self.cur += data


class Locations(unittest.TestCase):
    """share_map: the one table, both directions."""

    def loc(self, path, resolver=None):
        return share_map.locate(path, ROWS, resolver=resolver, here_is_hive=False)

    def test_the_table_file_is_what_both_sides_read(self):
        rows = share_map.load_table()
        self.assertEqual(rows, ROWS)
        with open(os.path.join(SCRIPTS, "hive_path.sh")) as fh:
            sh = fh.read()
        self.assertIn('SHARES="$HERE/hive_shares.tsv"', sh)
        self.assertNotIn("128.120.208.24/proteomics)", sh, "the table must not live in the script")

    def test_mac_mount(self):
        r = self.loc("/Volumes/proteomics/Data/lab/service/x")
        self.assertEqual(r["hive"], f"{FLINDERS}/Data/lab/service/x")
        self.assertEqual(r["windows"], "\\\\128.120.208.24\\proteomics\\Data\\lab\\service\\x")
        self.assertEqual(r["mac"], "/Volumes/proteomics/Data/lab/service/x")

    def test_drive_letter_through_hive_path_sh(self):
        calls = []

        def resolver(p):                  # hive_path.sh --no-verify: R: is \\128.120.208.24\proteomics
            calls.append(p)
            return "128.120.208.24", "proteomics", "Data/lab/service/P1"
        r = self.loc("R:\\Data\\lab\\service\\P1", resolver)
        self.assertEqual(calls, ["R:\\Data\\lab\\service\\P1"])
        self.assertEqual(r["hive"], f"{FLINDERS}/Data/lab/service/P1")
        self.assertEqual(r["windows"], "\\\\128.120.208.24\\proteomics\\Data\\lab\\service\\P1")

    def test_unc_and_quobyte(self):
        self.assertEqual(self.loc("\\\\128.120.208.24\\proteomics\\Data\\x")["hive"],
                         f"{FLINDERS}/Data/x")
        q = self.loc("/quobyte/proteomics-grp/SERVICE/a")
        self.assertEqual(q["hive"], "/quobyte/proteomics-grp/SERVICE/a")
        self.assertEqual(q["mac"], "/Volumes/proteomics-grp/SERVICE/a")
        self.assertIsNone(q["windows"], "the proteomics-grp Windows server is not recorded")

    def test_unknown_is_not_recorded_never_a_guess(self):
        r = self.loc("/Users/someone/data", resolver=lambda p: None)
        self.assertIsNone(r["hive"])
        self.assertTrue(r["how"].startswith("not recorded"))
        # hive_path.sh's fallback would GUESS /nfs/lssc0/flinders/<share>; share_map does not
        r = self.loc("S:\\x", resolver=lambda p: ("10.0.0.9", "labshare", "x"))
        self.assertIsNone(r["hive"])
        self.assertIn("not recorded", r["how"])

    def test_hive_path_sh_reads_the_table(self):
        """A copy of hive_path.sh next to an edited table maps the new share -- no ssh."""
        with tempfile.TemporaryDirectory() as d:
            shutil.copy(os.path.join(SCRIPTS, "hive_path.sh"), d)
            with open(os.path.join(d, "hive_shares.tsv"), "w") as fh:
                fh.write("#server\tshare\thive_path\tmac_mount\twindows_unc\n"
                         "srv9\tlabshare\t/quobyte/lab\t/Volumes/labshare\t\\\\srv9\\labshare\n")
            r = subprocess.run(["bash", os.path.join(d, "hive_path.sh"), "--no-verify",
                                "\\\\srv9\\labshare\\runs\\a"], capture_output=True, text=True,
                               timeout=60, env=dict(os.environ, HIVE_EXEC="false"))
            j = json.loads(r.stdout)
            self.assertEqual(j["candidates"], ["/quobyte/lab/runs/a"])
            self.assertIs(j["known_share"], True)
            self.assertIs(j["verified"], False)
            self.assertIn("--no-verify", j["how"])


class RawList(unittest.TestCase):
    def test_written_from_search_provenance_when_missing(self):
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            os.remove(p["raw_list"])
            files = ["/quobyte/proteomics-grp/SERVICE/lab/raw/a.d",
                     "/quobyte/proteomics-grp/SERVICE/lab/raw/b.d"]
            sp = json.loads(read(p["search_prov"]))
            write(p["search_prov"], json.dumps(dict(sp, files=files)))
            level, note = session_docs.ensure_raw_list(p)
            self.assertEqual(level, "OK")
            self.assertIn("search_provenance.json", note)
            self.assertEqual(session.read_raw_list(p["session_dir"]), files)

    def test_written_from_file_list_and_never_claimed_when_absent(self):
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            os.remove(p["raw_list"])
            level, note = session_docs.ensure_raw_list(p)
            self.assertEqual(level, "SKIPPED")               # no record anywhere
            f = session_docs.gather(p["session_dir"])
            md = session_docs.readme_md(f)
            self.assertNotIn("raw_files.txt (where the raw data are)", md)
            self.assertRegex(md, r"\| Raw data \| not recorded \|")
            write(os.path.join(p["search_out"], "file_list.txt"), "/quobyte/x/r1.d\n/quobyte/x/r2.d\n")
            self.assertEqual(session_docs.ensure_raw_list(p)[0], "OK")
            self.assertEqual(session.read_raw_list(p["session_dir"]),
                             ["/quobyte/x/r1.d", "/quobyte/x/r2.d"])


class LegacyEncodedRawList(unittest.TestCase):
    """Skill 2.7 and older wrote input/raw_files.txt in the computer's own encoding: on Windows
    cp1252, where even the header's dash (0x97) is not UTF-8. A strict UTF-8 read raised in
    gather(), and finalize lost README.html, AGENTS.md and the Methods."""

    def legacy(self, d, extra=""):
        p = tdp.dia_session(d)
        paths = session.read_raw_list(p["session_dir"])
        text = ("# Raw MS files used in this analysis (not copied — too large).\n"
                + "".join(r + "\n" for r in paths) + extra)
        with open(p["raw_list"], "wb") as fh:
            fh.write(text.encode("cp1252"))
        return p, paths

    def test_a_cp1252_path_is_read_and_said_to_be_damaged(self):
        with tempfile.TemporaryDirectory() as d:
            p, paths = self.legacy(d, "C:\\Daten\\Müller\\HeLa_extra.d\n")
            got = session.read_raw_list(p["session_dir"])
            self.assertEqual(got[:3], paths)
            self.assertEqual(got[3], "C:\\Daten\\M\ufffdller\\HeLa_extra.d")
            note = session.raw_list_encoding_note(p["session_dir"])
            self.assertIn("not UTF-8", note)
            self.assertIn("1 path(s)", note)
            self.assertIn("M\ufffdller", note)
            level, line = session_docs.ensure_raw_list(p)
            self.assertEqual(level, "OK")
            self.assertIn("present (4 files); not UTF-8", line)

    def test_a_header_only_problem_leaves_every_path_as_written(self):
        with tempfile.TemporaryDirectory() as d:
            p, paths = self.legacy(d)
            self.assertEqual(session.read_raw_list(p["session_dir"]), paths)
            self.assertIn("only in its comment lines",
                          session.raw_list_encoding_note(p["session_dir"]))
            write(p["raw_list"], "# utf-8 — fine\n" + "".join(x + "\n" for x in paths))
            self.assertIsNone(session.raw_list_encoding_note(p["session_dir"]))

    def test_finalize_keeps_the_readme_agents_and_methods(self):
        with tempfile.TemporaryDirectory() as d:
            p, _ = self.legacy(d, "C:\\Daten\\Müller\\HeLa_extra.d\n")
            with open(p["conditions"], "ab") as fh:                # a legacy sample name too
                fh.write("C:\\Daten\\Müller\\HeLa_extra.d,Treated\n".encode("cp1252"))
            res = Finalize.run_finalize(Finalize(), p)
            self.assertEqual(res["zip_docs"], "added")
            for rel in ("README.html", "AGENTS.md", "output/methods.md"):
                self.assertGreater(os.path.getsize(os.path.join(p["session_dir"], rel)), 0, rel)
            man = read(p["manifest_txt"])
            self.assertRegex(man, r"\[OK\]\s+input/raw_files\.txt .* not UTF-8 .*1 path\(s\)")
            self.assertRegex(man, r"\[OK\]\s+Publication methods \(output/methods\.md\)")
            self.assertNotIn("UnicodeDecodeError", man)


class Agents(unittest.TestCase):
    def build(self, d, prov, qc=False):
        p = tdp.dia_session(d)
        os.remove(os.path.join(p["de_dir"], "DE_Treated_vs_Control.csv"))
        write(os.path.join(p["de_dir"], "de_provenance.json"), json.dumps(prov))
        write(os.path.join(p["de_dir"], "DE_x_Treated.Control.csv"),
              '"Protein.Group","Genes","logFC","AveExpr","t","P.Value","adj.P.Val","B","Mystery"\n')
        write(os.path.join(p["de_dir"], "Expression_Matrix.csv"),
              '"Protein.Group","Genes","s1","s2","s3"\n')
        if qc:
            write(os.path.join(p["de_dir"], "QC_detected_vs_inferred.csv"),
                  "Sample,Group,Detected,Inferred,Total,PctDetected,PctInferred\n"
                  "s1,Control,80,20,100,80,20\ns2,Treated,60,40,100,60,40\n")
        return p

    def test_the_records_the_other_branches_added(self):
        """Detection_Matrix.csv, run_de.R's contaminant record and the CoreOmics submission are
        all named in AGENTS.md, each in its own record's words."""
        cont = {"policy": "removed", "removed": True, "tag": "Cont_", "id_column": "Protein.Ids",
                "n_precursors": 812, "n_protein_groups": 23, "n_sample_groups_sharing": 4,
                "n_sample_groups_all_shared": 1, "share_table": "QC_contaminant_share.csv",
                "removed_table": "contaminants_removed.csv", "database_risk": True,
                "database_note": "CANARY database note"}
        det = {"file": "Detection_Matrix.csv", "zero_means": "inferred",
               "values": "precursors observed for the protein in that run",
               "na_means": "protein not in the precursor matrix"}
        with tempfile.TemporaryDirectory() as d:
            p = self.build(d, dict(DPC, contaminants=cont, detection_matrix=det), qc=True)
            for n in ("Detection_Matrix.csv", "contaminants_removed.csv",
                      "QC_contaminant_share.csv"):
                write(os.path.join(p["de_dir"], n), "x\n")
            text = session_docs.agents_md(session_docs.gather(p["session_dir"]))
            self.assertIn("| Which values were measured vs inferred, per protein and sample | "
                          "`output/tables/Detection_Matrix.csv` |", text)
            self.assertIn("0 = inferred", text)
            self.assertIn("`Detection_Matrix.csv` marks every cell", text)
            self.assertIn("`output/tables/contaminants_removed.csv`", text)
            self.assertIn("`output/tables/QC_contaminant_share.csv`", text)
            self.assertIn("812 precursors", text)                  # make_methods' sentence
            self.assertIn("**Caveat:** CANARY database note", text)
            self.assertNotIn("kept them out of quantification", text)
            self.assertNotIn("input/submission.json", text)        # no record, no line

            import submission_report
            submission_report.attach(p["session_dir"], {"internal_id": "PROT_0001",
                                                        "sample_prep": "lab",
                                                        "prot_or_pep": "peptides"})
            text = session_docs.agents_md(session_docs.gather(p["session_dir"]))
            self.assertIn("`input/submission.json` |", text)
            self.assertIn("the submitting lab prepared the samples and sent peptides", text)
            self.assertIn("- CoreOmics submission: PROT_0001", text)

    def test_dpc_in_its_own_words(self):
        with tempfile.TemporaryDirectory() as d:
            p = self.build(d, DPC, qc=True)
            text = session_docs.agents_md(session_docs.gather(p["session_dir"]))
            self.assertIn("DPC-Quant + limma (limpa)", text)
            self.assertIn(DPC["missing_policy"], text)
            self.assertNotIn("MaxLFQ", text)
            self.assertIn("Inferred is not measured", text)
            self.assertIn("(20–40% of values per sample)", text)
            self.assertIn("reference line on the volcano only", text)
            self.assertIn("Treated-Control = 7", text)
            self.assertIn("`Mystery` — no description recorded", text)
            self.assertIn("`adj.P.Val` — Benjamini-Hochberg adjusted p-value", text)

    def test_maxlfq_in_its_own_words(self):
        with tempfile.TemporaryDirectory() as d:
            p = self.build(d, MAXLFQ)
            text = session_docs.agents_md(session_docs.gather(p["session_dir"]))
            self.assertIn("MaxLFQ + limma", text)
            self.assertIn("DIA-NN PG.MaxLFQ", text)
            self.assertIn("on/off calls", text)
            self.assertIn("covariates: Batch", text)
            for wrong in ("DPC", "limpa", "Inferred is not measured"):
                self.assertNotIn(wrong, text, wrong)

    def test_no_threshold_claimed_that_was_not_recorded(self):
        with tempfile.TemporaryDirectory() as d:
            prov = {k: v for k, v in DPC.items() if k not in ("significance_rule", "logfc_role")}
            p = self.build(d, prov)
            text = session_docs.agents_md(session_docs.gather(p["session_dir"]))
            self.assertIn("whether a fold-change filter was applied is not recorded", text)
            self.assertNotIn("reference line on the volcano only", text)

    def test_pull_down_controls_by_name(self):
        with tempfile.TemporaryDirectory() as d:
            prov = dict(DPC, groups={"Old_IgG": 3, "Old_RyR": 3},
                        contrasts=["Old_RyR-Old_IgG"])
            p = self.build(d, prov)
            text = session_docs.agents_md(session_docs.gather(p["session_dir"]))
            self.assertIn("Pull-down controls, by group name: Old_IgG", text)
            self.assertIn("enrichment over that control", text)
            p2 = self.build(os.path.join(d, "b"), DPC)          # Control/Treated: no IP claim
            self.assertNotIn("enrichment", session_docs.agents_md(
                session_docs.gather(p2["session_dir"])))

    def test_significant_counts_from_an_old_repr_record(self):
        with tempfile.TemporaryDirectory() as d:
            p = self.build(d, dict(DPC, significant_per_contrast="{'Treated-Control': 11}"))
            md = session_docs.readme_md(session_docs.gather(p["session_dir"]))
            self.assertIn("| Treated-Control | 11 |", md)


def run_brief(test, d, prov):
    """analysis_prompt.py's brief for a DE folder holding `prov` as de_provenance.json."""
    with open(os.path.join(d, "de_provenance.json"), "w") as fh:
        json.dump(prov, fh)
    out = os.path.join(d, "brief.md")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "analysis_prompt.py"),
                        "--out", out, "--de-dir", d], capture_output=True, text=True)
    test.assertEqual(r.returncode, 0, r.stderr)
    return read(out)


class DatabaseCaveat(unittest.TestCase):
    """run_de.R's database_note (an older FASTA: legacy or identity rule only) is a caveat in the
    report brief and in AGENTS.md, as the block warnings are."""

    def test_brief_and_agents_carry_the_database_note(self):
        cont = {"policy": "removed", "removed": True, "tag": "Cont_", "n_precursors": 10,
                "n_protein_groups": 3, "database_checked": True, "database_risk": True,
                "database_note": "CANARY identity rule alone: bovine EEF1A1, YWHAZ and TUBA1D"}
        prov = dict(DPC, contaminants=cont)
        with tempfile.TemporaryDirectory() as d:
            p = Agents().build(d, prov)
            agents = session_docs.agents_md(session_docs.gather(p["session_dir"]))
            self.assertIn("**Caveat:** CANARY identity rule alone", agents)
            b = run_brief(self, p["de_dir"], prov)
            self.assertIn("**Contaminant-filter caveat (say this in Data Quality Notes, naming "
                          "the proteins):** CANARY identity rule alone", b)


class BlockRecord(unittest.TestCase):
    """run_de.R --block (feat/skill-de-block): AGENTS.md and the analysis brief state the block in
    make_methods.de_block_sentence()'s words and carry its warnings (review C3); an unblocked
    record ({"applied": false}) adds nothing."""
    BLOCK = {"applied": True, "column": "Mouse", "effect": "random", "n_blocks": 6,
             "consensus_correlation": 0.412, "n_proteins_estimated": 900, "n_proteins": 1000,
             "fit": "limma lmFit(block=, correlation=)",
             "contrast_model": {"Treated-Control": "blocked"},
             "warnings": ["CANARY between-mouse contrasts are anti-conservative"]}

    def brief(self, d, prov):
        return run_brief(self, d, prov)

    def test_agents_and_brief_state_the_block_and_its_warnings(self):
        prov = dict(DPC, block=self.BLOCK, block_column="Mouse")
        with tempfile.TemporaryDirectory() as d:
            p = Agents().build(d, prov)
            text = session_docs.agents_md(session_docs.gather(p["session_dir"]))
            self.assertIn("- Blocking: Samples sharing a Mouse were modelled as correlated", text)
            self.assertIn("**Caveat:** CANARY between-mouse contrasts", text)
            b = self.brief(p["de_dir"], prov)
            self.assertIn("- Blocking: Samples sharing a Mouse were modelled as correlated", b)
            self.assertIn("**Blocking caveat (say this in Data Quality Notes):** CANARY", b)

    def test_an_unblocked_record_adds_nothing(self):
        prov = dict(DPC, block={"applied": False,
                                "note": "no --block: samples modelled as independent"})
        with tempfile.TemporaryDirectory() as d:
            p = Agents().build(d, prov)
            self.assertNotIn("Blocking", session_docs.agents_md(
                session_docs.gather(p["session_dir"])))
            self.assertNotIn("Blocking", self.brief(p["de_dir"], prov))


class Readme(unittest.TestCase):
    def test_html_is_self_contained_and_matches_the_md(self):
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            session_docs.write_docs(p["session_dir"])
            md = read(p["readme"])
            page = read(os.path.join(p["session_dir"], "README.html"))
            h = H2s()
            h.feed(page)
            self.assertEqual((h.stack, h.bad), ([], 0))
            self.assertEqual(h.external, [])
            self.assertIn("style", h.tags)
            # the look is report_style.py's, shared with Analysis_Report.html (one design
            # system), over semantic markup
            import report_style
            self.assertIn(report_style.CSS + report_style.DOC_CSS, page)
            self.assertIn('<main class="doc" id="main">', page)
            for tag in ("header", "nav", "section", "aside", "table"):
                self.assertIn(tag, h.tags, tag)
            self.assertIn('<aside class="callout note">', page)
            self.assertIn("<nav><h2>Start here</h2>", page)
            self.assertIn('<meta charset="utf-8">', page)
            self.assertEqual(h.h2, sections_md(md))
            for s in ("Start here", "Summary", "Where this lives on HIVE", "Reproduce",
                      "Methods", "Deposit the data (PRIDE / MassIVE)"):
                self.assertIn(s, h.h2)
            self.assertIn('href="AGENTS.md"', page)
            self.assertIn('href="output/tables/"', page)
            self.assertNotIn("Open `README.html`", page.replace("<code>", "`").replace(
                "</code>", "`"))
            self.assertIn("Open `README.html`", md)

    def test_links_only_to_files_that_exist(self):
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            write(os.path.join(p["output_dir"], "Analysis_Report.html"), "<html></html>")
            md = session_docs.readme_md(session_docs.gather(p["session_dir"]))
            for target in re.findall(r"\]\(([^)]+)\)", md):
                self.assertTrue(os.path.exists(os.path.join(p["session_dir"], target)), target)
            self.assertIn("(output/Analysis_Report.html)", md)
            write(os.path.join(p["output_dir"], "AI_Analysis_Report.docx"), "old")   # older session
            md = session_docs.readme_md(session_docs.gather(p["session_dir"]))
            self.assertNotIn("AI_Analysis_Report.docx", md)          # kept, never pointed at
            first = [ln for ln in md.splitlines() if ln.startswith("- [")][0]
            self.assertIn("Analysis_Report.html", first)             # the report, linked first

    def test_the_report_pdf_and_text_twin(self):
        """make_analysis_html writes Analysis_Report.pdf and .md beside the HTML: the README
        links both, AGENTS.md names both, and an HTML with no PDF says how to make one."""
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            write(os.path.join(p["output_dir"], "Analysis_Report.html"), "<html></html>")
            md = session_docs.readme_md(session_docs.gather(p["session_dir"]))
            # why it was not made is MANIFEST.txt's to say (finalize), never assumed here
            self.assertIn("`output/Analysis_Report.pdf` — not made: open `Analysis_Report.html`, "
                          "Print, Save as PDF", md)
            self.assertNotIn("no browser", md)
            self.assertNotIn("(output/Analysis_Report.pdf)", md)
            write(p["manifest_txt"], "[INFO]    Report PDF (output/Analysis_Report.pdf) -- x\n")
            md = session_docs.readme_md(session_docs.gather(p["session_dir"]))
            self.assertIn("not made (the \"Report PDF\" line of `MANIFEST.txt` says why)", md)
            write(os.path.join(p["output_dir"], "Analysis_Report.pdf"), "%PDF-1.4")
            write(os.path.join(p["output_dir"], "Analysis_Report.md"), "# twin")
            f = session_docs.gather(p["session_dir"])
            md = session_docs.readme_md(f)
            self.assertIn("[Analysis report (PDF)](output/Analysis_Report.pdf)", md)
            self.assertIn("[Analysis report (plain text)](output/Analysis_Report.md)", md)
            self.assertIn("for NotebookLM or other AI notebooks", md)
            self.assertNotIn("— not made", md)
            agents = session_docs.agents_md(f)
            self.assertIn("`output/Analysis_Report.md`", agents)
            self.assertIn("`output/Analysis_Report.pdf`", agents)

    def test_relative_links_render_and_schemes_do_not(self):
        s = make_deposit._inline("[r](output/a%20b.html) [x](javascript:alert(1)) [w](https://a.b)")
        self.assertIn('<a href="output/a%20b.html">r</a>', s)
        self.assertNotIn('href="javascript', s)
        self.assertIn('<a href="https://a.b">w</a>', s)


class Finalize(unittest.TestCase):
    def run_finalize(self, p, hooks=None):
        args = type("A", (), dict(dir=p["session_dir"], zip=True, no_deposit=True,
                                  no_notify=True, reanalysis_of=""))()
        out = io.StringIO()
        with mock.patch.dict(os.environ, job_env(os.path.dirname(p["session_dir"])),
                             clear=True), redirect_stdout(out), redirect_stderr(io.StringIO()):
            if hooks is None:
                session.do_finalize(args)
            else:
                with mock.patch.object(session, "_finish_hooks", return_value=hooks):
                    session.do_finalize(args)
        return json.loads(out.getvalue())

    def test_docs_in_the_zip_and_the_manifest(self):
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            write(p["output_files_md"], "# Output files — what each one is\n")
            res = self.run_finalize(p)
            self.assertEqual(res["zip_docs"], "added")
            base = os.path.basename(p["session_dir"])
            with zipfile.ZipFile(res["zip"]) as z:
                names = z.namelist()
                for n in ("README.md", "README.html", "AGENTS.md", "MANIFEST.txt"):
                    self.assertEqual(names.count(f"{base}/{n}"), 1, n)
            man = read(p["manifest_txt"])
            for part in ("README.html", "AGENTS.md", "input/raw_files.txt"):
                self.assertRegex(man, rf"\[OK\]\s+{re.escape(part)}")
            files_md = read(p["output_files_md"])
            self.assertIn("`README.html`", files_md)
            self.assertIn("`AGENTS.md`", files_md)
            self.run_finalize(p)                                  # re-finalize: one block
            self.assertEqual(read(p["output_files_md"]).count("written at finalize: BEGIN"),
                             1)
            md = read(p["readme"])
            self.assertIn("[MANIFEST.txt](MANIFEST.txt)", md)

    def test_registry_path_from_the_hook_is_written_in(self):
        reg = "/quobyte/proteomics-grp/skill_runs/sessions/2026-09-24_demo"
        hooks = {"run_log": {"logged": True, "path": reg}, "slack": {"sent": False,
                 "level": "INFO", "detail": "off"},
                 "manifest": [("OK", "Core run log", "logged")]}
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            res = self.run_finalize(p, hooks)
            self.assertIn(reg, read(p["readme"]))
            self.assertIn(reg, read(os.path.join(p["session_dir"], "AGENTS.md")))
            with zipfile.ZipFile(res["zip"]) as z:
                base = os.path.basename(p["session_dir"])
                self.assertIn(reg, z.read(f"{base}/README.html").decode())

    def test_a_raw_list_that_cannot_be_written_keeps_the_readme(self):
        """ensure_raw_list() failing used to set session_docs = None: README.html and AGENTS.md
        were lost, and MANIFEST.txt blamed session_docs.py for it."""
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            with mock.patch.object(session_docs, "ensure_raw_list",
                                   side_effect=OSError("disk full")):
                res = self.run_finalize(p)
            self.assertEqual(res["zip_docs"], "added")
            for n in ("README.html", "AGENTS.md"):
                self.assertGreater(os.path.getsize(os.path.join(p["session_dir"], n)), 0, n)
            man = read(p["manifest_txt"])
            self.assertRegex(man, r"\[SKIPPED\] input/raw_files\.txt .*-- OSError: disk full")
            self.assertNotIn("could not be loaded", man)
            self.assertRegex(man, r"\[OK\]\s+README\.html")

    def test_a_readme_that_cannot_render_is_neither_empty_nor_zipped(self):
        """README.html was opened, then rendered: a failed render left a 0-byte file, and the zip
        took anything that existed."""
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            html = os.path.join(p["session_dir"], "README.html")
            base = os.path.basename(p["session_dir"])
            broke = mock.patch.object(session_docs, "_render_html",
                                      side_effect=RuntimeError("renderer broke"))
            with broke:
                res = self.run_finalize(p)
            self.assertFalse(os.path.exists(html), "a failed render left a file behind")
            self.assertFalse(os.path.exists(html + ".part"))
            with zipfile.ZipFile(res["zip"]) as z:
                names = z.namelist()
            self.assertNotIn(f"{base}/README.html", names)
            for n in ("README.md", "AGENTS.md"):
                self.assertIn(f"{base}/{n}", names)
            self.assertIn("README.html", res["zip_docs"])
            man = read(p["manifest_txt"])
            self.assertRegex(man, r"\[SKIPPED\] README\.html .*RuntimeError: renderer broke")
            self.assertRegex(man, r"\[SKIPPED\] README / AGENTS\.md in the zip")
            # an older README.html stays whole on disk, and is not zipped as this finalize's
            write(html, "<html>older</html>")
            with broke:
                res = self.run_finalize(p)
            self.assertEqual(read(html), "<html>older</html>")
            with zipfile.ZipFile(res["zip"]) as z:
                self.assertNotIn(f"{base}/README.html", z.namelist())

    def test_docs_subcommand_as_a_copy(self):
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            real = "/Volumes/proteomics/Data/lab/service/x/2026-09-24_demo"
            r = subprocess.run([PY, os.path.join(SCRIPTS, "session.py"), "docs", "--dir",
                                p["session_dir"], "--as", real], capture_output=True, text=True,
                               env=job_env(d))
            self.assertEqual(r.returncode, 0, r.stderr)
            md = read(p["readme"])
            self.assertIn(f"`{FLINDERS}/Data/lab/service/x/2026-09-24_demo`", md)
            self.assertFalse(os.path.exists(p["manifest_txt"]), "docs must not write MANIFEST")
            self.assertNotIn("[MANIFEST.txt]", md)


class WhereEverythingIs(unittest.TestCase):
    def test_the_tables_a_reader_needs_are_named(self):
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            values = "precursors observed for the protein in that run; 0 = value inferred"
            de = dict(DPC, detection_matrix={"file": "Detection_Matrix.csv", "values": values},
                      contaminants={"policy": "removed", "tag": "Cont_",
                                    "removed_table": "contaminants_removed.csv"})
            write(os.path.join(p["de_dir"], "de_provenance.json"), json.dumps(de))
            for n in ("Detection_Matrix.csv", "QC_detected_vs_inferred.csv",
                      "contaminants_removed.csv"):
                write(os.path.join(p["de_dir"], n), "Protein.Group\n")
            write(os.path.join(p["output_dir"], "AUDIT.md"), "**Overall: PASS**\n")
            write(os.path.join(p["output_dir"], "AUDIT.json"), "{}")
            write(os.path.join(p["output_dir"], "SAMPLE_QUALITY.md"), "# Sample quality\n")
            md = session_docs.readme_md(session_docs.gather(p["session_dir"]))
            where = md.split("## Where everything is")[1].split("\n## ")[0]
            for rel in ("output/tables/Detection_Matrix.csv",
                        "output/tables/QC_detected_vs_inferred.csv",
                        "output/tables/contaminants_removed.csv", "output/AUDIT.md",
                        "output/SAMPLE_QUALITY.md"):
                self.assertIn(f"| `{rel}` |", where)
            self.assertIn(values, where)                         # the record's own words
            self.assertIn("`output/AUDIT.md` | the pitfall audit: PASS / WARN / FAIL per check "
                          "(and `.json` beside it)", where)
            self.assertNotIn("SAMPLE_QUALITY.md` | sample quality and contamination flags (and",
                             where)

    def test_core_storage_says_ask_the_core_not_the_flinders_server(self):
        """Bio review S6: a raw folder on /quobyte/proteomics-grp read as if a collaborator
        could open it with smb://128.120.208.24 -- which is Flinders."""
        f = {"locations": [
            dict(share_map.locate("/quobyte/proteomics-grp/raw/P1", ROWS), what="Raw data"),
            dict(share_map.locate(FLINDERS + "/Data/x", ROWS), what="Search output")]}
        t = session_docs.locations_table(f)
        raw = next(ln for ln in t.splitlines() if ln.startswith("| Raw data"))
        self.assertIn(CORE_ONLY, raw)                           # the Windows cell
        self.assertIn("`/Volumes/proteomics-grp/raw/P1` (Core members)", raw)
        self.assertIn("`smb://128.120.208.24/proteomics` for paths under `/Volumes/proteomics`",
                      t)
        self.assertIn(f"Paths under `/quobyte/proteomics-grp` (`/Volumes/proteomics-grp`): "
                      f"{CORE_ONLY}.", t)
        flinders_only = session_docs.locations_table({"locations": f["locations"][1:]})
        self.assertNotIn("Core storage", flinders_only)
        self.assertNotIn("Core members", flinders_only)


class RegistryLocate(unittest.TestCase):
    def test_locate_reads_the_index_only(self):
        import record_run
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            runs = os.path.join(d, "skill_runs")
            folder = os.path.join(runs, "sessions", "2026-09-24_demo")
            os.makedirs(folder)
            ident = {"session": os.path.realpath(p["session_dir"])}
            write(os.path.join(folder, "run_record.json"), json.dumps({"identity": ident}))
            os.makedirs(os.path.join(runs, ".index"))
            with mock.patch.dict(os.environ, {"SKILL_RUNS_DIR": runs, "RECORD_RUN": "on"}):
                self.assertIsNone(record_run.locate(p["session_dir"]))
                os.symlink(os.path.relpath(folder, os.path.join(runs, ".index")),
                           os.path.join(runs, ".index",
                                        f"session_{record_run.key16(ident['session'])}"))
                before = sorted(os.listdir(runs))
                self.assertEqual(record_run.locate(p["session_dir"]), folder)
                self.assertEqual(sorted(os.listdir(runs)), before, "locate() never writes")
                r = subprocess.run([PY, os.path.join(SCRIPTS, "record_run.py"), "locate",
                                    "--session", p["session_dir"]], capture_output=True,
                                   text=True, timeout=60)
                self.assertEqual(json.loads(r.stdout), {"located": folder})
            with mock.patch.dict(os.environ, {"SKILL_RUNS_DIR": runs, "RECORD_RUN": "off"}):
                self.assertIsNone(record_run.locate(p["session_dir"]))


if __name__ == "__main__":
    unittest.main()
