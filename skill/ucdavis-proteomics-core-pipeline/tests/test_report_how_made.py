#!/usr/bin/env python3
"""Every report says how it was made: "How this report was made" (make_analysis_html.py).

The Core director, 2026-10-27: the reports never said that an AI (Claude) ran the pipeline and
wrote the report text. Every page build_page() assembles -- the HTML report and its Markdown
twin -- now carries that section, after the pipeline's own sections and before the appendix,
with a link to the public explainer (docs/how-reports-are-made.html). Each sentence comes from a
record of the run (CLAUDE.md rules 1-2): the search engine and version from
search_provenance.json, the DE by its own display_label, reproducibility_log.R when it is there,
the skill version from plugin.json, the podcast from podcast.json and its meaning check, the
expert review from the report's own section. A fact no record holds is left out, never filled
in; a section of the same title in the AI-written report wins.
Stdlib only.
"""
import argparse
import html
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import make_analysis_html as mah  # noqa: E402
import make_podcast as mp  # noqa: E402
import skill_version as sv  # noqa: E402

TITLE = "How this report was made"
REPORT = "# Study\n\n*A standfirst.*\n\n## Key findings\n\nSome text.\n"
PNG = b"\x89PNG\r\n\x1a\n" + b"x" * 40


def write(path, text, mode="w"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, mode, **({} if "b" in mode else {"encoding": "utf-8"})) as fh:
        fh.write(text)


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def podcast(backend="kokoro", check="PASS", meaning="PASS", reviewer="a reviewer agent (Claude)",
            **extra):
    """podcast.json as render writes it (the keys make_podcast.cmd_render records)."""
    tts = {"backend": backend, "sent_to_cloud": backend == "gemini", "fallback": None}
    if backend == "kokoro":
        tts.update(model="kokoro", models_used=["kokoro"], engine={"name": "Kokoro-82M"})
    elif backend == "gemini":
        tts.update(model="gemini-2.5-pro-preview-tts", models_used=["gemini-2.5-pro-preview-tts"])
    elif backend == "say":
        tts.update(model="say", models_used=["say"])
    man = {"show": "Signal to Noise", "title": "Contact sites", "audio": "podcast.m4a",
           "transcript": "transcript.html", "duration_s": 600.0, "ai_generated": True,
           "script": "podcast_script.md", "script_sha256": "abc123",
           "check": {"status": check, "file": "check.txt", "overridden": check != "PASS"},
           "meaning_check": {"status": meaning, "file": "meaning_check.md",
                             "overridden": False},
           "tts": tts}
    man.update(extra)
    return man, (f"# Meaning check: podcast_script.md\nscript_sha256: abc123\n"
                 f"verdict: {meaning}\nreviewer: {reviewer}\n")


def raw_sections(doc, md):
    """The section as written, markers and all: (inside its HTML <section>, its .md body)."""
    h = re.search(r'<section id="sec-how-this-report-was-made">(.*?)</section>', doc, re.S)
    m = re.search(r"(?ms)^## How this report was made\n(.*?)(?=^## |^---$|\Z)", md)
    return (h.group(1) if h else None), (m.group(1) if m else None)


def blocks_text(blocks):
    """how_made()'s blocks as the .md twin carries them."""
    return "\n\n".join(b[1] if b[0] == "md" else b[2] for b in blocks)


def link(out, *extra):
    return subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_podcast.py"), "link", out,
                           *extra], capture_output=True, text=True)


class HowMade(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.s = self._td.name
        self.out = os.path.join(self.s, "output")
        self.tables = os.path.join(self.out, "tables")
        write(os.path.join(self.out, "AI_Analysis_Report.md"), REPORT)
        write(os.path.join(self.tables, "de_provenance.json"),
              json.dumps({"adjp": 0.05, "display_label": "DPC-Quant + limma (limpa)",
                          "pipeline_id": "dpc", "skill_version": sv.skill_version()}))
        write(os.path.join(self.tables, "reproducibility_log.R"), "# the analysis as R\n")
        write(os.path.join(self.out, "search", "search_provenance.json"),
              json.dumps({"engine": "diann", "version": "2.7.0"}))
        # the Methods (9d), which the podcast's finish gate needs (test_finish_gates.py)
        write(os.path.join(self.out, "methods.md"), "# Methods\n")

    def tearDown(self):
        self._td.cleanup()

    def set_podcast(self, man, meaning_md=None):
        pod = os.path.join(self.out, "podcast")
        write(os.path.join(pod, "podcast.json"), json.dumps(man))
        if meaning_md is not None:
            write(os.path.join(pod, "meaning_check.md"), meaning_md)

    def build(self):
        """-> (html page, Markdown twin, the section's text in the HTML, the section's text in
        the .md) -- the text as a reader sees it, the markers left out"""
        page = os.path.join(self.out, "Analysis_Report.html")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"),
                            "--session", self.s, "--out", page, "--no-pdf"],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        doc, md = read(page), read(os.path.join(self.out, "Analysis_Report.md"))
        h, m = raw_sections(doc, md)
        return (doc, md, html.unescape(re.sub(r"<[^>]+>", " ", h)) if h is not None else None,
                re.sub(r"<!--.*?-->\n?", "", m) if m is not None else None)

    def files(self):
        return (read(os.path.join(self.out, "Analysis_Report.html")),
                read(os.path.join(self.out, "Analysis_Report.md")))

    # -------------------------------------------------------------- the section, both renderers
    def test_in_the_html_and_the_markdown_twin_with_the_link(self):
        doc, md, h, m = self.build()
        self.assertIsNotNone(h)
        self.assertIsNotNone(m)
        self.assertEqual(len(re.findall(r"<h2 [^>]*>How this report was made</h2>", doc)), 1)
        self.assertEqual(len(re.findall(r"(?m)^## How this report was made$", md)), 1)
        self.assertIn(f'<a href="{mah.HOW_MADE_URL}">{mah.HOW_MADE_LINK}</a>', doc)
        self.assertIn(f"[{mah.HOW_MADE_LINK}]({mah.HOW_MADE_URL})", m)
        self.assertEqual(mah.HOW_MADE_URL,
                         "https://bsphinney.github.io/DE-LIMP/how-reports-are-made.html")
        self.assertIsNone(re.search(r"</?[A-Za-z][A-Za-z0-9]*\b[^>]*>", m))   # no tags in the .md
        self.assertIn('href="#how-this-report-was-made"', doc)                  # in the contents

    def test_the_realistic_case_says_each_part_from_its_record(self):
        man, mc = podcast()
        self.set_podcast(man, mc)
        _, _, h, m = self.build()
        version = json.load(open(os.path.join(SCRIPTS, "..", ".claude-plugin", "plugin.json"),
                                 encoding="utf-8"))["version"]
        for text in (h, m):
            self.assertIn("DIA-NN 2.7.0 searched the raw instrument files", text)
            self.assertIn("DPC-Quant + limma (limpa) estimated the protein amounts and ran the "
                          "statistics", text)
            self.assertIn("tables/reproducibility_log.R", text)
            self.assertIn("re-runs the statistics as plain R code, with no AI involved", text)
            self.assertIn("Claude, Anthropic's AI model, ran each step, following the UC Davis "
                          f"Proteomics Core's pipeline (skill {version}), and wrote this "
                          "report's analysis text.", " ".join(text.split()))
            self.assertNotIn("This report was built with skill", text)
            self.assertIn("Claude also wrote the podcast script.", text)
            self.assertIn("A program checked the script's numbers and names against the report.",
                          text)
            self.assertIn("A separate Claude reviewer, not the script's writer, read every "
                          "sentence", text)
            self.assertIn("Kokoro-82M on UC Davis's HIVE computing cluster", text)
            self.assertIn("Nothing was sent off campus to make the voices.", text)
            self.assertNotIn("Expert", text)
            self.assertNotRegex(text, r"Claude (?:Opus|Sonnet|Haiku|Fable)|claude-")  # no model

    def test_after_the_pipeline_sections_and_before_the_appendix(self):
        write(os.path.join(self.tables, "de_provenance.json"),
              json.dumps({"adjp": 0.05, "display_label": "SILAC ratios + limma",
                          "report_sections": [{"title": "Channel mixing",
                                               "markdown": "The pipeline's own section."}]}))
        figs = os.path.join(self.out, "figures")
        write(os.path.join(figs, "qc_pvalue_panel.png"), PNG, "wb")
        write(os.path.join(figs, "figures.json"),
              json.dumps([{"file": "qc_pvalue_panel.png", "caption": "p-values."}]))
        doc, md, _, _ = self.build()
        md_h2 = re.findall(r"(?m)^## (.*)$", md)
        html_h2 = [re.sub(r"<[^>]+>", "", t) for t in re.findall(r'<h2 id="[^"]+">(.*?)</h2>', doc)]
        self.assertEqual(md_h2, html_h2)
        i = md_h2.index(TITLE)
        self.assertEqual(md_h2[i - 1], "Channel mixing")
        self.assertEqual(md_h2[i + 1:], [mah.APPENDIX])

    def test_a_page_without_a_report_still_says_it(self):
        os.remove(os.path.join(self.out, "AI_Analysis_Report.md"))
        figs = os.path.join(self.out, "figures")
        write(os.path.join(figs, "pca.png"), PNG, "wb")
        write(os.path.join(figs, "figures.json"), json.dumps([{"file": "pca.png", "caption": "c"}]))
        _, _, h, m = self.build()
        for text in (h, m):
            self.assertIn("Claude, Anthropic's AI model, ran each step", text)
            self.assertNotIn("wrote this report", text)       # no AI-written text on this page
            self.assertIn(mah.HOW_MADE_LINK, text)

    def test_the_section_alone_is_still_nothing_to_render(self):
        empty = os.path.join(self.s, "empty")
        os.makedirs(empty)
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"),
                            "--figures", empty, "--out", os.path.join(empty, "r.html"),
                            "--no-pdf"], capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("nothing to render", r.stderr)

    # -------------------------------------------------------------- no record, no sentence
    def test_no_search_record_no_engine_or_version(self):
        os.remove(os.path.join(self.out, "search", "search_provenance.json"))
        _, _, h, m = self.build()
        for text in (h, m):
            self.assertNotIn("searched", text)
            self.assertNotIn("DIA-NN", text)
            self.assertNotIn("2.7.0", text)
            self.assertIn("DPC-Quant + limma (limpa) estimated the protein amounts", text)

    def test_a_search_record_without_a_version_names_no_version(self):
        write(os.path.join(self.out, "search", "search_provenance.json"),
              json.dumps({"engine": "diann", "version": None}))
        _, _, h, m = self.build()
        for text in (h, m):
            self.assertIn("DIA-NN searched the raw instrument files", text)
            self.assertNotIn("2.7.0", text)

    def test_no_de_record_no_de_sentence(self):
        os.remove(os.path.join(self.tables, "de_provenance.json"))
        _, _, h, m = self.build()
        for text in (h, m):
            self.assertNotIn("estimated the protein amounts", text)
            self.assertNotIn("limma", text)
            self.assertIn("DIA-NN 2.7.0 searched the raw instrument files.", text)

    def test_no_engine_and_no_de_no_software_sentence(self):
        os.remove(os.path.join(self.tables, "de_provenance.json"))
        os.remove(os.path.join(self.out, "search", "search_provenance.json"))
        os.remove(os.path.join(self.tables, "reproducibility_log.R"))
        _, _, h, m = self.build()
        for text in (h, m):
            self.assertNotIn("established scientific software", text)
            self.assertIn("Claude, Anthropic's AI model", text)

    def test_no_reproducibility_log_no_plain_r_sentence(self):
        os.remove(os.path.join(self.tables, "reproducibility_log.R"))
        _, _, h, m = self.build()
        for text in (h, m):
            self.assertNotIn("plain R", text)
            self.assertNotIn("reproducibility_log", text)

    def test_an_unknown_skill_version_names_none(self):
        a = argparse.Namespace(session=self.s, tables=self.tables,
                               out=os.path.join(self.out, "Analysis_Report.html"))
        with mock.patch.object(sv, "skill_version", return_value=sv.UNKNOWN):
            for prov in ({}, {"skill_version": sv.UNKNOWN}):
                text = blocks_text(mah.how_made(a, prov, REPORT, []))
                self.assertIn("following the UC Davis Proteomics Core's pipeline, and wrote", text)
                self.assertNotIn("skill", text)
                self.assertNotIn("unknown", text)

    def test_the_de_records_version_is_named_and_a_rebuild_says_its_own(self):
        """The version that ran the DE (de_provenance.json skill_version) is the pipeline's; a page
        built by another version says so, and a DE record from before the field names none."""
        a = argparse.Namespace(session=self.s, tables=self.tables,
                               out=os.path.join(self.out, "Analysis_Report.html"))
        flat = lambda prov: " ".join(blocks_text(mah.how_made(a, prov, REPORT, [])).split())  # noqa: E731
        with mock.patch.object(sv, "skill_version", return_value="9.9.9"):
            t = flat({"skill_version": "9.9.9"})
            self.assertIn("pipeline (skill 9.9.9), and wrote", t)
            self.assertNotIn("built with", t)
            t = flat({"skill_version": "9.9.8"})
            self.assertIn("pipeline (skill 9.9.8), and wrote this report's analysis text. "
                          "This report was built with skill 9.9.9.", t)
            t = flat({})                                  # a DE from before the record named it
            self.assertIn("pipeline, and wrote this report's analysis text. This report was built "
                          "with skill 9.9.9.", t)
            t = flat({"skill_version": sv.UNKNOWN})
            self.assertIn("pipeline, and wrote this report's analysis text. This report was built "
                          "with skill 9.9.9.", t)

    def test_the_header_search_fact_is_unchanged(self):
        a = argparse.Namespace(session=self.s)
        self.assertEqual(dict(mah.session_facts(a, {}))["Search"], "DIA-NN 2.7.0")
        self.assertEqual(mah.search_engine(self.s), ("DIA-NN", "2.7.0"))
        self.assertEqual(mah.search_engine(None), (None, None))

    # -------------------------------------------------------------- the podcast's part
    def test_no_podcast_no_podcast_sentence(self):
        _, _, h, m = self.build()
        for text in (h, m):
            self.assertNotIn("podcast", text.lower())
            self.assertNotIn("voices", text)

    def test_a_meaning_check_that_did_not_pass_is_not_claimed(self):
        man, mc = podcast(meaning="FAIL")
        man["meaning_check"]["overridden"] = True
        self.set_podcast(man, mc)
        _, _, h, m = self.build()
        for text in (h, m):
            self.assertIn("Claude also wrote the podcast script.", text)
            self.assertIn("A program checked", text)
            self.assertNotIn("reviewer", text)

    def test_an_overridden_check_is_not_claimed(self):
        man, mc = podcast(check="FAIL", meaning="missing")
        self.set_podcast(man, mc)
        _, _, h, m = self.build()
        for text in (h, m):
            self.assertNotIn("A program checked", text)
            self.assertNotIn("reviewer", text)
            self.assertIn("Kokoro-82M", text)

    def test_a_person_reviewer_is_not_called_claude(self):
        man, mc = podcast(reviewer="Jane Doe (Core staff)")
        self.set_podcast(man, mc)
        _, _, h, m = self.build()
        for text in (h, m):
            self.assertIn("A separate reviewer, not the script's writer", text)
            self.assertNotIn("Claude reviewer", text)
            self.assertNotIn("Jane", text)                  # never a name from the record

    def test_another_scripts_meaning_check_names_no_claude_reviewer(self):
        man, mc = podcast()
        self.set_podcast(man, mc.replace("abc123", "def456"))
        _, _, h, m = self.build()
        self.assertIn("A separate reviewer, not the script's writer", m)
        self.assertNotIn("Claude reviewer", m)

    def test_gemini_says_what_google_received(self):
        man, mc = podcast(backend="gemini",
                          verify={"status": "OK", "cloud_consent": "Core director",
                                  "audio_sent": "AAC 16 kHz"})
        self.set_podcast(man, mc)
        _, _, h, m = self.build()
        for text in (h, m):
            self.assertIn("The voices were made by Google's Gemini (gemini-2.5-pro-preview-tts), "
                          "in Google's cloud.", text)
            self.assertIn("Google received the script's text and the finished audio, to check "
                          "that it was spoken as written, never the report.", text)
            self.assertNotIn("off campus", text)
        man["verify"] = {"status": "NOT RUN", "reason": "--no-verify was given"}
        self.set_podcast(man)
        _, _, _, m = self.build()
        self.assertIn("Google received the script's text, never the report.", m)

    def test_kokoro_after_gemini_failed_says_what_gemini_had(self):
        man, mc = podcast()
        man["tts"]["fallback"] = {"from": "gemini", "to": "kokoro", "why": "billing",
                                  "chunks": 12, "gemini_chunks_sent": 3,
                                  "gemini_chunks_not_in_script": 0, "preflight_requests": 0}
        self.set_podcast(man, mc)
        _, _, h, m = self.build()
        for text in (h, m):
            self.assertIn("Kokoro-82M on UC Davis's HIVE computing cluster", text)
            self.assertIn("Google's Gemini, tried first, failed; it had been sent the text of 3 "
                          "of the 12 parts of the script.", text)
            self.assertNotIn("off campus", text)
        man["tts"]["fallback"] = {"from": "gemini", "to": "kokoro", "why": "nokey"}
        self.set_podcast(man)
        _, _, _, m = self.build()
        self.assertIn("Nothing was sent off campus to make the voices.", m)          # Gemini was never asked
        self.assertNotIn("tried first", m)

    def test_macos_say(self):
        man, mc = podcast(backend="say")
        self.set_podcast(man, mc)
        _, _, _, m = self.build()
        self.assertIn("macOS's built-in speech (say), on the Mac that made the podcast.", m)
        self.assertIn("Nothing was sent anywhere to make the voices.", m)

    # -------------------------------------------------------------- the podcast after the report
    def test_link_puts_the_paragraph_in_as_a_later_build_would(self):
        """The usual order: the report first, then the podcast and `link`. Link fills in the
        podcast's paragraph in place (before the expert review's), in the HTML and the .md; a
        second link changes nothing; and a report made after the podcast says exactly the same."""
        write(os.path.join(self.out, "AI_Analysis_Report.md"),
              REPORT + "\n## Expert Review Notes\n\nThe review checked the names.\n")
        doc, md, h, m = self.build()
        self.assertNotIn(mp.START, doc)
        self.assertNotIn("podcast", h.lower())
        self.assertNotIn("podcast", m.lower())
        self.assertIn(mp.MADE_START + mp.MADE_END, doc)             # its place, empty
        self.assertIn(mp.MADE_START + "\n" + mp.MADE_END, md)
        man, mc = podcast()
        self.set_podcast(man, mc)
        r = link(self.out, "--unchecked")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("podcast paragraph of How this report was made replaced", r.stdout)
        doc, md = self.files()
        self.assertEqual(doc.count(mp.START), 1)                    # the Listen card
        linked = raw_sections(doc, md)
        for text in linked:
            self.assertIn("Claude also wrote the podcast script.", text)
            self.assertIn("Kokoro-82M", text)
            self.assertLess(text.index("Claude also wrote"), text.index("Separate Claude reviewers"))
            self.assertEqual(text.count(mp.MADE_START), 1)
        r = link(self.out, "--unchecked")                           # again: nothing changes
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("podcast paragraph of How this report was made already current", r.stdout)
        self.assertEqual(self.files(), (doc, md))
        doc, md, _, _ = self.build()                                # made after the podcast
        self.assertEqual(raw_sections(doc, md), linked)
        self.assertEqual(doc.count('class="pc-card"'), 1)

    def test_add_how_made_fills_the_place_or_leaves_the_text(self):
        """share's hook (build_share) and link's edit: the one block, or no change."""
        man, mc = podcast()
        self.set_podcast(man, mc)
        empty = "<p>a</p>" + mp.MADE_START + mp.MADE_END + "<p>b</p>"
        got = mp.add_how_made(empty, self.out)
        self.assertEqual(got, "<p>a</p>" + mp.how_made_block(self.out, "html") + "<p>b</p>")
        self.assertIn("<p>Claude also wrote the podcast script.", got)
        self.assertEqual(mp.add_how_made(got, self.out), got)
        self.assertEqual(mp.add_how_made("<p>no place</p>", self.out), "<p>no place</p>")
        os.remove(os.path.join(self.out, "podcast", "podcast.json"))    # no podcast: the place
        self.assertEqual(mp.add_how_made(got, self.out), empty)         # is emptied again
        self.assertEqual(mp.how_made_block(self.out, "md"), mp.MADE_START + "\n" + mp.MADE_END)

    def test_link_reprints_the_pdf_with_the_paragraph(self):
        page = os.path.join(self.out, "Analysis_Report.html")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"),
                            "--session", self.s, "--out", page], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        pdf = os.path.join(self.out, "Analysis_Report.pdf")
        if not (os.path.isfile(pdf) and shutil.which("pdftotext")):
            self.skipTest("no headless browser to print the PDF, or no pdftotext to read it")
        man, mc = podcast()
        self.set_podcast(man, mc)
        os.utime(pdf, (1, 1))                       # older than the HTML link is about to write
        self.assertEqual(link(self.out, "--unchecked").returncode, 0)
        text = subprocess.run(["pdftotext", pdf, "-"], capture_output=True, text=True).stdout
        self.assertIn("Claude also wrote the podcast script.", " ".join(text.split()))

    def test_a_check_on_the_report_itself_holds_across_link_and_a_rebuild(self):
        """A podcast checked against the report's own Analysis_Report.md (references/podcast.md
        allows it): link filling in the paragraph, and a report made again after the podcast,
        change only the section's body, which the check's sources leave out (source_text)."""
        sys.path.insert(0, HERE)
        import test_make_podcast as tpod
        write(os.path.join(self.out, "AI_Analysis_Report.md"), tpod.REPORT)
        self.build()
        twin = os.path.join(self.out, "Analysis_Report.md")
        page = os.path.join(self.out, "Analysis_Report.html")
        script = os.path.join(self.out, "podcast", "podcast_script.md")
        write(script, tpod.script_text())
        rc, _, _ = tpod.run("check", script, "--source", twin)
        self.assertEqual(rc, 0, read(os.path.join(self.out, "podcast", "check.txt")))
        tpod.meaning_pass(script)
        s = mp.parse_script(script)
        man, _ = podcast(script_sha256=s.sha256)
        self.set_podcast(man)

        def holds():
            self.assertEqual(mp.read_check(s)[0], "PASS", mp.read_check(s)[2])
            self.assertEqual(mp.read_meaning(s)[0], "PASS", mp.read_meaning(s)[1])
        holds()
        r = link(self.out)                                          # no --unchecked
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        for text in raw_sections(*self.files()):
            self.assertIn("Claude also wrote the podcast script.", text)
        holds()                                                     # across link ...
        self.build()
        holds()                                                     # ... and a rebuild after it
        self.assertEqual(link(self.out).returncode, 0)              # link still vouches for it
        holds()
        # The HTML report as a source too, checked once its Listen card is in (link's card
        # changes the HTML's hash on its own, in origin/main too): a rebuild still holds
        rc, _, _ = tpod.run("check", script, "--source", twin, page)
        self.assertEqual(rc, 0, read(os.path.join(self.out, "podcast", "check.txt")))
        self.build()
        holds()
        self.assertEqual(link(self.out).returncode, 0)
        holds()
        # ... while a change to the report's own text still counts
        write(os.path.join(self.out, "AI_Analysis_Report.md"), tpod.REPORT + "\nOne more line.\n")
        self.build()
        self.assertEqual(mp.read_check(s)[0], "STALE")

    # -------------------------------------------------------------- the expert review
    def test_expert_review_notes_named_only_when_the_report_has_them(self):
        _, _, h, m = self.build()
        self.assertNotIn("reviewers", m)
        write(os.path.join(self.out, "AI_Analysis_Report.md"),
              REPORT + "\n## Expert Review Notes\n\nThe review checked the names.\n")
        _, _, h, m = self.build()
        for text in (h, m):
            self.assertIn("Separate Claude reviewers — a proteomics expert, a biologist and a "
                          "statistician — checked the analysis; their notes are in the Expert "
                          "Review Notes section.", " ".join(text.split()))

    # -------------------------------------------------------------- the report's own wins
    def test_the_reports_own_section_of_that_title_wins(self):
        write(os.path.join(self.out, "AI_Analysis_Report.md"),
              REPORT + "\n## How this report was made\n\nThe report's own words.\n")
        doc, md, h, m = self.build()
        self.assertEqual(len(re.findall(r"<h2 [^>]*>How this report was made</h2>", doc)), 1)
        self.assertEqual(len(re.findall(r"(?m)^## How this report was made$", md)), 1)
        self.assertIn("The report's own words.", m)
        self.assertNotIn("Anthropic", md)
        self.assertNotIn(mah.HOW_MADE_URL, doc)
        man, mc = podcast()                         # link has no place for the paragraph then
        self.set_podcast(man, mc)
        r = link(self.out, "--unchecked")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("no place for podcast paragraph of How this report was made", r.stdout)
        self.assertNotIn("Claude also wrote", "".join(self.files()))

    def test_the_writer_is_told_to_leave_the_section_to_the_tool(self):
        src = read(os.path.join(SCRIPTS, "analysis_prompt.py"))
        self.assertIn('w(f"Do not write a \\"{HOW_MADE_TITLE}\\" section: the HTML report adds it',
                      src)
        import analysis_prompt
        self.assertEqual(analysis_prompt.HOW_MADE_TITLE, TITLE)


if __name__ == "__main__":
    unittest.main(verbosity=2)
