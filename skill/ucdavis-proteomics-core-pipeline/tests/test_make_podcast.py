#!/usr/bin/env python3
"""make_podcast.py -- the optional audio discussion: check, render, link. Offline, stdlib only.

Why it exists: a first prototype had Gemini write the script. It passed a number check and still
said things the report never did ("limpa is our Core's custom extension to limma"; "Kcnb2, the
Kv2.1 beta subunit"), ran 30 minutes and read numbers badly. Now the agent writes the script
from the report and this tool refuses to voice it until:
  * every number in the transcript is in the sources (commas, unicode minus, 5e-15 vs 5.2e-15,
    "10^-15" as an order of magnitude; integers never rounded) and no number is spelled out;
  * every symbol-like token (Jph3, FKBP12.6, IgG ...) is in the sources or disclosed under
    "Claims beyond the report";
  * the AI disclosure is spoken in the first 3 turns, and no forbidden name appears.
Render caches each chunk (a re-run makes no TTS call; editing one line re-makes one chunk),
sends nothing to the cloud without --cloud-ok, and never lets the key reach any output. Link
puts ONE Listen card near the top of the report however often it runs, and the report generator
and session_docs.py keep it.
"""
import array
import base64
import contextlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import types
import unittest
import wave
import zipfile
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import make_podcast as mp  # noqa: E402

FAKE_KEY = "AIzaSyFAKE0123456789abcdefghijklmnopqrs"   # shaped like a key; not one

REPORT = """# Contact-site interactomes in Old and Young mouse brain

## Overview
30 IPs: 4 baits (JPH3, JPH4, Kv2.1, RyR) plus IgG controls, n = 3 per group. 6,112 proteins
were quantified across 12 comparisons.

## Key findings
Ryr2 is the top hit in the RyR pulldown (adj.P = 5.2e-15; log2FC 10.62). The RyR IP contains
FKBP12.6 (Fkbp1b) and calmodulin. Kv2.1 enriched 215 proteins at both ages, 77 of them
Kv2.1-only. A median 29.4% of values were inferred. Jph3 fell by −1.3 in Old; its partners by a
median of −0.9 (n = 155). C1qa rose in Old Kv2.1 pulldowns. The two thinnest Young IgG runs
detected 2,897 and 3,407 proteins.
"""

SEG1 = [("MAYA", "Welcome to Signal to Noise. Quick note first: this is an AI-generated "
                 "discussion, and our voices are synthetic."),
        ("LEO", "The report is the record. We are just two voices reading it closely."),
        ("MAYA", "30 IPs, 4 baits, and 6,112 proteins across 12 comparisons.")]
SEG2 = [("LEO", "Ryr2 tops the RyR pulldown, adjusted p of 5e-15 and a log2FC of 10.6."),
        ("MAYA", "And FKBP12.6 comes along. Call it 10^-15 if you like round numbers."),
        ("LEO", "Kv2.1 has 215 proteins at both ages; 77 of them are Kv2.1-only."),
        ("MAYA", "Jph3 drops −1.3, its 155 partners about −0.9, C1qa goes up, and Kcnb2 is a "
                 "cousin we will not over-read.")]
CLAIMS = "- Kcnb2 is a Kv2-family subunit (general knowledge, not in the report)."


def script_text(segs=(SEG1, SEG2), claims=CLAIMS, header=None, pron=None, extra=""):
    head = header if header is not None else [
        "# Signal to Noise — Contact-site interactomes",
        "", "Show: Signal to Noise", "Title: Contact-site interactomes, Old vs Young",
        "Hosts: Maya (cell biologist), Leo (statistician)", "Voices: Maya=Kore, Leo=Charon"]
    rows = pron if pron is not None else [("Jph3", "J P H three"), ("C1qa", "C one Q A")]
    L = list(head) + ["", "## Pronunciation", "", "| Written | Spoken |", "|---|---|"]
    L += [f"| {w} | {s} |" for w, s in rows]
    L += ["", "## Claims beyond the report", "", claims, "", "## Transcript", "", mp.T_START, ""]
    for i, seg in enumerate(segs):
        if i:
            L += ["---", ""]
        for who, text in seg:
            L += [f"**{who}:** {text}", ""]
    L += [mp.T_END, extra]
    return "\n".join(L) + "\n"


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def run(*argv):
    """make_podcast.main(argv) -> (rc, stdout, stderr)."""
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        rc = mp.main(list(argv))
    return rc, out.getvalue(), err.getvalue()


class Workspace(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.d = self._td.name
        self.out = os.path.join(self.d, "output")
        self.pod = os.path.join(self.out, "podcast")
        self.report = os.path.join(self.out, "AI_Analysis_Report.md")
        self.script = os.path.join(self.pod, "podcast_script.md")
        write(self.report, REPORT)

    def tearDown(self):
        self._td.cleanup()

    def check(self, text=None, *extra):
        write(self.script, text if text is not None else script_text())
        rc, out, err = run("check", self.script, "--source", self.report, *extra)
        return rc, read(os.path.join(self.pod, "check.txt"))


# ------------------------------------------------------------------------------------ parse
class Parse(Workspace):
    def test_the_script_format(self):
        write(self.script, script_text(header=[
            "# Signal to Noise", "Title: T", "Hosts: Maya (cell biologist), Leo (statistician)",
            "Voices: Leo=Puck", "Say voices: Maya=Samantha, Leo=Daniel"]))
        s = mp.parse_script(self.script)
        self.assertEqual(s.problems, [])
        self.assertEqual(s.hosts, [("Maya", "cell biologist"), ("Leo", "statistician")])
        self.assertEqual(s.gemini_voices, {"Maya": "Kore", "Leo": "Puck"})   # default + override
        self.assertEqual(s.say_voices, {"Maya": "Samantha", "Leo": "Daniel"})
        self.assertEqual(s.pronunciation, [("Jph3", "J P H three"), ("C1qa", "C one Q A")])
        self.assertEqual(len(s.claims), 1)
        self.assertEqual([len(x) for x in s.segments], [3, 4])
        lines = read(self.script).splitlines()
        for t in s.turns():                      # line numbers point at the turn itself
            self.assertIn(t.text[:30], lines[t.line - 1])
        self.assertEqual(s.turns()[3].speaker, "Leo")

    def test_claims_none_is_explicit_and_empty_is_a_problem(self):
        write(self.script, script_text(claims="None"))
        s = mp.parse_script(self.script)
        self.assertEqual((s.claims, s.problems), ([], []))
        # no "## Transcript" heading: the claims section runs into the transcript markers
        write(self.script, script_text(claims="None").replace("## Transcript\n", ""))
        s = mp.parse_script(self.script)
        self.assertEqual((s.claims, s.problems, len(s.turns())), ([], [], 7))
        write(self.script, script_text(claims=""))
        self.assertTrue(any("Claims beyond the report" in m for _, m in
                            mp.parse_script(self.script).problems))

    def test_a_wrapped_turn_or_a_stranger_is_a_problem_not_silently_dropped(self):
        text = script_text().replace("**LEO:** The report is the record.",
                                     "**LEO:** The report is\nthe record.")
        text = text.replace("**MAYA:** 30 IPs", "**NADIA:** 30 IPs")
        write(self.script, text)
        msgs = [m for _, m in mp.parse_script(self.script).problems]
        self.assertTrue(any("not a turn" in m for m in msgs), msgs)
        self.assertTrue(any("NADIA" in m and "not one of the hosts" in m for m in msgs), msgs)


class RealFormat(Workspace):
    """tests/fixtures/podcast/format_v1_script.md copies, line for line, the layout of the first
    real Signal to Noise script (PROT_0756, 2026-09-25): bold header bullets, hosts separated by
    '·' with the voice in the parentheses, a Pronunciation table with a digits row, a Claims
    section with prose, nested bullets and no '## Transcript' heading, blank lines between turns,
    10 segments. The study itself is synthetic (the real one is unpublished client data)."""
    FIX = os.path.join(HERE, "fixtures", "podcast")

    def test_parses_exactly(self):
        s = mp.parse_script(os.path.join(self.FIX, "format_v1_script.md"))
        self.assertEqual(s.problems, [])
        self.assertEqual(s.hosts, [("Maya", "cell biologist"), ("Leo", "statistician")])
        self.assertEqual(s.gemini_voices, {"Maya": "Kore", "Leo": "Charon"})
        self.assertEqual(s.show, "Signal to Noise")
        self.assertTrue(s.title.startswith("Who stays with the chaperone?"))
        self.assertTrue(s.author.startswith("a test author"))
        self.assertEqual(len(s.segments), 10)
        self.assertEqual(s.turns()[0].speaker, "Maya")
        self.assertIn(("0421", "zero four two one"), s.pronunciation)
        self.assertEqual(len(s.pronunciation), 12)
        self.assertEqual(len(s.claims), 8)                    # nested bullets are claims too
        self.assertTrue(any(c.startswith("a 1% FDR means") for c in s.claims))
        self.assertFalse(any("TRANSCRIPT" in c or "**MAYA" in c for c in s.claims))

    def test_passes_check_teaches_and_closes(self):
        os.makedirs(self.pod, exist_ok=True)
        shutil.copy(os.path.join(self.FIX, "format_v1_script.md"), self.script)
        rc, out, err = run("check", self.script, "--source",
                           os.path.join(self.FIX, "synthetic_report.md"))
        txt = read(os.path.join(self.pod, "check.txt"))
        self.assertEqual(rc, 0, txt)
        self.assertNotIn("never taught", txt)
        self.assertNotIn("what to do with this", txt)
        self.assertIn("segments: 10", txt)
        rc, out, err = run("render", self.script, "--tts", "say", "--dry-run")
        self.assertEqual(rc, 0, err)
        self.assertIn("Maya: Leo, picture a yeast cell", out)
        self.assertIn("submission zero four two one", out)
        self.assertIn("adjusted p of 3.1 times 10 to the minus 12.", out)
        self.assertIn("D I A N N 2.6.1", out)

    def test_teaching_and_close_warnings(self):
        rc, txt = self.check()                                # the short script teaches little
        self.assertIn("the collaborator is never taught:", txt)
        self.assertIn("what LC-MS/MS does", txt)
        self.assertIn("the last two segments never say what to do with this", txt)


# ------------------------------------------------------------------------------------ check
class Check(Workspace):
    def test_a_clean_script_passes(self):
        rc, txt = self.check()
        self.assertEqual(rc, 0, txt)
        self.assertIn("status: PASS", txt)
        self.assertIn(f"source: {os.path.abspath(self.report)} sha256=", txt)
        # rounding and order of magnitude are allowed, and said so
        self.assertRegex(txt, r"10\.6 matches a source value rounded")
        self.assertRegex(txt, r"5e-15 matches a source value rounded")
        self.assertRegex(txt, r"10\^-15 is spoken as an order of magnitude")
        self.assertIn("Kcnb2 is not in the sources; it is listed under Claims", txt)

    def test_an_invented_number_fails(self):
        rc, txt = self.check(script_text().replace("215 proteins", "4,555 proteins"))
        self.assertEqual(rc, 1)
        self.assertIn("number 4555 is not in the sources", txt)
        self.assertIn("turn 6", txt)

    def test_integers_are_never_rounded(self):
        # the source says 29.4%; "29%" is a different claim (and the brief says to quote it)
        rc, txt = self.check(script_text().replace("215 proteins at both ages",
                                                   "215 proteins at both ages, 29% inferred"))
        self.assertEqual(rc, 1)
        self.assertIn("number 29 is not in the sources", txt)

    def test_an_invented_symbol_fails_unless_disclosed(self):
        bad = script_text().replace("C1qa goes up", "Kcnb9 goes up")
        rc, txt = self.check(bad)
        self.assertEqual(rc, 1)
        self.assertIn("Kcnb9 looks like a gene/protein symbol", txt)
        rc, txt = self.check(script_text(claims=CLAIMS + "\n- Kcnb9 is named from memory.")
                             .replace("C1qa goes up", "Kcnb9 goes up"))
        self.assertEqual(rc, 0, txt)
        rc, txt = self.check(script_text(claims="None"))          # Kcnb2 no longer disclosed
        self.assertEqual(rc, 1)
        self.assertIn("Kcnb2 looks like", txt)

    def test_symbols_match_the_report_case_insensitively(self):
        rc, txt = self.check(script_text().replace("Ryr2 tops", "RYR2 tops"))
        self.assertEqual(rc, 0, txt)

    def test_the_ai_disclosure_must_be_in_the_first_three_turns(self):
        seg1 = [(w, t.replace("this is an AI-generated discussion", "this is a discussion"))
                for w, t in SEG1]
        rc, txt = self.check(script_text(segs=(seg1, SEG2)))
        self.assertEqual(rc, 1)
        self.assertIn("no AI disclosure in the first 3 turns", txt)

    def test_a_forbidden_name_fails(self):
        text = script_text().replace("We are just two voices", "The Dickson lab asked us, two voices")
        rc, txt = self.check(text, "--forbid-name", "Dickson", "Silva")
        self.assertEqual(rc, 1)
        self.assertIn("forbidden name 'Dickson' appears in: turn 2", txt)
        self.assertIn("forbidden names: Dickson, Silva", txt)

    def test_spelled_out_numbers_fail(self):
        rc, txt = self.check(script_text().replace("30 IPs", "Thirty IPs"))
        self.assertEqual(rc, 1)
        self.assertIn("spelled-out number 'Thirty'", txt)

    def test_required_sections(self):
        text = script_text().replace("## Pronunciation", "## Notes")
        rc, txt = self.check(text)
        self.assertEqual(rc, 1)
        self.assertIn("missing section '## Pronunciation'", txt)
        text = script_text().replace("## Claims beyond the report", "## Extras")
        rc, txt = self.check(text)
        self.assertIn("missing section '## Claims beyond the report'", txt)

    def test_length_is_reported_and_warned(self):
        rc, txt = self.check()
        self.assertRegex(txt, r"words: \d+ · turns: 7 · segments: 2 · about [\d.]+ min at 150 wpm")
        self.assertIn("the brief asks for 2,200-3,600", txt)           # a short test script
        self.assertIn("2 segments", txt)

    def test_numbers_hidden_in_an_embedded_image_do_not_count(self):
        html_src = os.path.join(self.out, "Analysis_Report.html")
        write(html_src, "<html><body><p>" + REPORT.replace("\n", " ") + "</p><img src=\"data:image"
                        "/png;base64,iVBORw0KGgo4555AAAA4555==\"><script>var x=4555;</script>"
                        "</body></html>")
        write(self.script, script_text().replace("215 proteins", "4555 proteins"))
        rc, _, _ = run("check", self.script, "--source", html_src)
        self.assertEqual(rc, 1)
        self.assertIn("number 4555 is not in the sources", read(os.path.join(self.pod, "check.txt")))

    def test_number_forms_normalise(self):
        book = mp.NumberBook()
        book.add_text("6,112 proteins; −1.3; 5.2×10⁻¹⁵; 3.407e3; p = 0.004")
        v = lambda s: [book.verdict(n) for n in mp.numbers_in(s)]       # noqa: E731
        self.assertEqual(v("6112"), ["exact"])
        self.assertEqual(v("-1.3 and 1.3 and +1.3"), ["exact"] * 3)
        self.assertEqual(v("5e-15"), ["rounded"])
        self.assertEqual(v("5.2 times 10 to the minus 15"), ["exact"])
        self.assertEqual(v("10 to the minus 15"), ["magnitude"])
        self.assertEqual(v("10^-3"), ["magnitude"])                   # p = 0.004
        self.assertEqual(v("3407"), ["exact"])
        self.assertEqual(v("6113"), [None])
        self.assertEqual(v("7 and 40 and 100"), ["trivial"] * 3)
        self.assertEqual(v("Kv2.1 FKBP12.6 log2"), [])                # symbols, not numbers
        self.assertEqual([(n.kind, str(n.value)) for n in mp.numbers_in("p of 5.2e-15.")],
                         [("sci", "5.2E-15")])                        # not "5.2" and "15"


# ------------------------------------------------------------------------------ pronunciation
class Pronunciation(unittest.TestCase):
    def test_longest_first_whole_tokens(self):
        sp = mp.pronouncer()
        self.assertEqual(sp("Kv2.1 and Kv2 and Kv2.2"),
                         "K V two point one and K V two and K V two point two")
        self.assertEqual(sp("the ER, ERK and hERG"), "the E R, ERK and hERG")
        self.assertEqual(sp("log2FC vs log2"), "log two fold change versus log two")
        self.assertEqual(sp("adj.P of 5e-15"), "adjusted p of 5 times 10 to the minus 15")
        self.assertEqual(sp("1e-9 and 10^-9 and 10⁻⁹"),
                         "10 to the minus 9 and 10 to the minus 9 and 10 to the minus 9")
        self.assertEqual(sp("fell by −1.3 (≥ 2, ~3)"), "fell by minus 1.3 ( at least 2, about 3)")
        self.assertEqual(sp("**really** *big*"), "really big")
        # the end of a sentence is not part of the number (live Gemini test, 2026-09-25)
        self.assertEqual(sp("Adjusted p of 5e-15."), "Adjusted p of 5 times 10 to the minus 15.")

    def test_the_script_table_adds_and_wins(self):
        sp = mp.pronouncer([("Kv2.1", "kay vee two one"), ("Tecr", "teck R")])
        self.assertEqual(sp("Kv2.1, Kv2 and Tecr"), "kay vee two one, K V two and teck R")


class Chunks(unittest.TestCase):
    @staticmethod
    def seg(*words):
        return [mp.Turn("Maya", " ".join(["w"] * n), i, 0) for i, n in enumerate(words)]

    def test_long_segments_split_at_turns_short_ones_do_not(self):
        c = mp.make_chunks([self.seg(*[100] * 12), self.seg(150, 150), self.seg(800)])
        self.assertEqual([x.words for x in c], [400, 400, 400, 300, 800])
        self.assertEqual([x.seg for x in c], [0, 0, 0, 1, 2])


# ------------------------------------------------------------------------------------ render
class FakeTTS(object):
    """A backend that returns a tone and counts its calls."""
    name, cloud = "fake", False
    calls = []

    def __init__(self, s, a):
        self.voices = {h: "tone" for h, _ in s.hosts}
        self.model = "fake-1"
        self.produced, self.models_used, self.warnings = 0, [], []

    def prepare(self, chunks, speak, cache_dir):
        pass

    def payload(self, chunk, speak, model=None):
        return {"backend": "fake", "turns": [[t.speaker, speak(t.text)] for t in chunk.turns]}

    def synth(self, chunk, speak, label):
        FakeTTS.calls.append([speak(t.text) for t in chunk.turns])
        if self.model not in self.models_used:
            self.models_used.append(self.model)
        return array.array("h", [1200, -1200]) * int(mp.RATE * 0.05 * len(chunk.turns))  # 0.1 s/turn


class Render(Workspace):
    def setUp(self):
        super().setUp()
        mp.BACKENDS["fake"] = FakeTTS
        FakeTTS.calls = []
        self.enc = mock.patch.object(mp, "encode_aac", return_value=(None, "test: no encoder"))
        self.enc.start()

    def tearDown(self):
        self.enc.stop()
        mp.BACKENDS.pop("fake", None)
        super().tearDown()

    def test_render_is_gated_on_a_passing_check(self):
        write(self.script, script_text())
        rc, out, err = run("render", self.script, "--tts", "fake", "--dry-run")   # preview: ok
        self.assertEqual(rc, 0, err)
        self.assertIn("Leo: K V two point one has 215 proteins", out)
        rc, out, err = run("render", self.script, "--tts", "fake")
        self.assertEqual(rc, 2)
        self.assertIn("no check.txt", err)
        self.assertEqual(FakeTTS.calls, [])
        self.assertEqual(self.check()[0], 0)
        write(self.script, script_text() + "\n")                      # edited after the check
        rc, out, err = run("render", self.script, "--tts", "fake")
        self.assertEqual(rc, 2)
        self.assertIn("changed after check.txt", err)
        rc, out, err = run("render", self.script, "--tts", "fake", "--unchecked")
        self.assertEqual(rc, 0, err)
        man = json.loads(read(os.path.join(self.pod, "podcast.json")))
        self.assertTrue(man["check"]["overridden"])

    def test_render_outputs_cache_and_one_line_edits(self):
        self.assertEqual(self.check()[0], 0)
        rc, out, err = run("render", self.script, "--tts", "fake")
        self.assertEqual(rc, 0, err)
        self.assertEqual(len(FakeTTS.calls), 2)                       # 2 segments, 2 chunks
        man = json.loads(read(os.path.join(self.pod, "podcast.json")))
        for k in ("show", "title", "hosts", "tts", "script_sha256", "sources", "words",
                  "duration_s", "created", "cloud_tts_consent", "audio", "transcript"):
            self.assertIn(k, man)
        self.assertIs(man["ai_generated"], True)
        self.assertEqual(man["audio"], "podcast.wav")                 # no encoder: WAV kept
        self.assertEqual(man["script_sha256"], mp.parse_script(self.script).sha256)
        self.assertEqual(man["sources"][0]["file"], os.path.abspath(self.report))
        self.assertEqual(man["check"]["status"], "PASS")
        self.assertGreater(man["duration_s"], 3)
        self.assertEqual(len(os.listdir(os.path.join(self.pod, ".cache"))), 2)
        # speech gets the substitutions; the transcript keeps the written form
        spoken = " ".join(x for c in FakeTTS.calls for x in c)
        self.assertIn("K V two point one", spoken)
        self.assertIn("J P H three", spoken)
        self.assertNotIn("Kv2.1", spoken)
        page = read(os.path.join(self.pod, "transcript.html"))
        self.assertIn("Kv2.1", page)
        self.assertNotIn("K V two point one", page.split('id="made"')[0])
        self.assertIn(mp.esc(mp.DISCLOSURE), page)
        self.assertIn("Kcnb2 is a Kv2-family subunit", page)          # the claims ledger
        self.assertIn('src="podcast.wav"', page)

        FakeTTS.calls = []                                            # re-run: all cached
        rc, out, err = run("render", self.script, "--tts", "fake")
        self.assertEqual((rc, FakeTTS.calls), (0, []))
        self.assertEqual(json.loads(out)["chunks_cached"], 2)

        write(self.script, read(self.script).replace("77 of them", "77 of those"))
        self.assertEqual(self.check(read(self.script))[0], 0)
        rc, out, err = run("render", self.script, "--tts", "fake")    # one line -> one chunk
        self.assertEqual(rc, 0, err)
        self.assertEqual(len(FakeTTS.calls), 1)
        self.assertIn("77 of those are K V two point one-only.", " ".join(FakeTTS.calls[0]))

    def test_render_elsewhere_copies_the_script_and_check(self):
        self.assertEqual(self.check()[0], 0)
        dest = os.path.join(self.d, "elsewhere")
        rc, out, err = run("render", self.script, "--tts", "fake", "--out", dest)
        self.assertEqual(rc, 0, err)
        for f in ("podcast_script.md", "check.txt", "podcast.json", "transcript.html",
                  "podcast.wav"):
            self.assertTrue(os.path.isfile(os.path.join(dest, f)), f)

    def test_assembled_audio_has_chimes_gaps_and_normalised_level(self):
        self.assertEqual(self.check()[0], 0)
        run("render", self.script, "--tts", "fake")
        a = mp.read_wav(os.path.join(self.pod, "podcast.wav"))
        speech = sum(len(c) for c in FakeTTS.calls) * 0.1             # seconds of fake speech
        extra = 2 * 1.4 + 0.5 + 0.6 + mp.GAP_SEGMENT                  # chimes + pads + 1 gap
        self.assertAlmostEqual(len(a) / float(mp.RATE), speech + extra, delta=0.05)
        self.assertLess(max(a), 29500)
        self.assertGreater(max(a), 3000)                              # boosted from 1200

    @unittest.skipUnless(shutil.which("afconvert") or shutil.which("ffmpeg"), "no AAC encoder")
    def test_real_encoder_writes_m4a_and_drops_the_wav(self):
        self.enc.stop()
        try:
            self.assertEqual(self.check()[0], 0)
            rc, out, err = run("render", self.script, "--tts", "fake")
            self.assertEqual(rc, 0, err)
            man = json.loads(read(os.path.join(self.pod, "podcast.json")))
            self.assertEqual(man["audio"], "podcast.m4a")
            self.assertTrue(os.path.getsize(os.path.join(self.pod, "podcast.m4a")) > 1000)
            self.assertFalse(os.path.exists(os.path.join(self.pod, "podcast.wav")))
        finally:
            self.enc.start()

    @unittest.skipUnless(shutil.which("say"), "macOS say not installed")
    def test_say_backend_offline(self):
        segs = ([("MAYA", "An AI-generated test, Kv2.1 first."), ("LEO", "Then 215 proteins.")],)
        self.assertEqual(self.check(script_text(segs=segs, claims="None"))[0], 0)
        rc, out, err = run("render", self.script, "--tts", "say")
        self.assertEqual(rc, 0, err)
        man = json.loads(read(os.path.join(self.pod, "podcast.json")))
        self.assertEqual((man["tts"]["backend"], man["tts"]["sent_to_cloud"]), ("say", False))
        self.assertIs(man["cloud_tts_consent"], False)
        self.assertGreater(man["duration_s"], 4)


# ------------------------------------------------------------------------------ gemini (fake)
def fake_genai(behaviour, models=("gemini-2.5-pro-preview-tts", "gemini-2.5-flash-preview-tts",
                                  "gemini-3.8-flash-tts", "gemini-3.8-flash-lite-tts"),
               list_error=None):
    """google.genai stand-ins. `behaviour(model, contents, config)` returns PCM bytes or raises.
    generate_content (2.x) passes the prompt text as `contents`; interactions.create (3.x) passes
    its `input` list, and the fake answers the way the docs describe: base64 audio/wav."""
    google, genai, gt = (types.ModuleType("google"), types.ModuleType("google.genai"),
                         types.ModuleType("google.genai.types"))

    class Cfg(object):
        def __init__(self, **kw):
            self.__dict__.update(kw)
    for n in ("GenerateContentConfig", "SpeechConfig", "MultiSpeakerVoiceConfig",
              "SpeakerVoiceConfig", "VoiceConfig", "PrebuiltVoiceConfig"):
        setattr(gt, n, Cfg)
    log = {"clients": 0, "calls": [], "apis": []}

    class Models(object):
        def list(self):
            if list_error:
                raise list_error
            return [types.SimpleNamespace(name="models/" + m, supported_actions=["generateContent"])
                    for m in models] + [types.SimpleNamespace(name="models/gemini-9-pro-tts",
                                                              supported_actions=["generateContent"])]

        def generate_content(self, model, contents, config):
            log["calls"].append((model, contents, config))
            log["apis"].append("generate_content")
            data = behaviour(model, contents, config)
            part = types.SimpleNamespace(inline_data=types.SimpleNamespace(
                data=data, mime_type="audio/L16;codec=pcm;rate=24000"))
            return types.SimpleNamespace(candidates=[types.SimpleNamespace(
                content=types.SimpleNamespace(parts=[part]), finish_reason="STOP")])

    class Interactions(object):
        def create(self, model, input, response_format, generation_config, **kw):
            log["calls"].append((model, input, {"response_format": response_format,
                                                "generation_config": generation_config}))
            log["apis"].append("interactions")
            pcm = behaviour(model, input, generation_config)
            buf = io.BytesIO()
            with wave.open(buf, "wb") as w:
                w.setnchannels(1)
                w.setsampwidth(2)
                w.setframerate(24000)
                w.writeframes(pcm)
            return types.SimpleNamespace(status="completed", output_audio=types.SimpleNamespace(
                data=base64.b64encode(buf.getvalue()).decode("ascii"), mime_type="audio/wav"))

    class Client(object):
        def __init__(self, api_key):
            log["clients"] += 1
            self.models = Models()
            self.interactions = Interactions()
    genai.Client, genai.types, google.genai = Client, gt, genai
    return {"google": google, "google.genai": genai, "google.genai.types": gt}, log


def pcm_for(contents):
    if isinstance(contents, list):                                # interactions input
        words = sum(len(p["text"].split()) for turn in contents for p in turn["content"])
    else:
        words = len(contents.split("\n\n", 1)[-1].split())
    n = int(mp.RATE * words * 60.0 / mp.WPM / 2)
    return (array.array("h", [900, -900]) * n).tobytes()


class ApiError(Exception):
    """Shaped like the SDK's errors: generate_content's APIError has .code and .details (the
    response JSON); the interactions client's has .status_code and .body."""

    def __init__(self, code, msg, details=None, style="genai"):
        super().__init__(f"{code} {msg}")
        if style == "genai":
            self.code, self.details = code, details
        else:
            self.status_code, self.body = code, details


def rpc_error(quota_id, delay=None, value="10"):
    """A 429 body the way Google sends it: QuotaFailure + RetryInfo in error.details."""
    det = [{"@type": "type.googleapis.com/google.rpc.QuotaFailure",
            "violations": [{"quotaMetric": "generativelanguage.googleapis.com/"
                                           "generate_content_free_tier_requests",
                            "quotaId": quota_id, "quotaValue": value}]}]
    if delay is not None:
        det.append({"@type": "type.googleapis.com/google.rpc.RetryInfo", "retryDelay": delay})
    return {"error": {"code": 429, "status": "RESOURCE_EXHAUSTED",
                      "message": "You exceeded your current quota.", "details": det}}


PER_MINUTE = "GenerateRequestsPerMinutePerProjectPerModel-FreeTier"
PER_DAY = "GenerateRequestsPerDayPerProjectPerModel-FreeTier"


class Gemini(Workspace):
    def setUp(self):
        super().setUp()
        self.env = mock.patch.dict(os.environ, {"GEMINI_API_KEY": FAKE_KEY})
        self.env.start()
        self.nosleep = mock.patch.object(mp, "_sleep", lambda s: None)
        self.nosleep.start()
        self.enc = mock.patch.object(mp, "encode_aac", return_value=(None, "test"))
        self.enc.start()
        self.assertEqual(self.check()[0], 0)

    def tearDown(self):
        for p in (self.env, self.nosleep, self.enc):
            p.stop()
        super().tearDown()

    def render(self, behaviour, *extra, **kw):
        mods, log = fake_genai(behaviour, **kw)
        with mock.patch.dict(sys.modules, mods):
            rc, out, err = run("render", self.script, "--tts", "gemini", *extra)
        return rc, out, err, log

    def assert_no_key_anywhere(self, *texts):
        for t in texts:
            self.assertNotIn(FAKE_KEY, t)
        for root, _, files in os.walk(self.d):
            for f in files:
                with open(os.path.join(root, f), "rb") as fh:
                    self.assertNotIn(FAKE_KEY.encode(), fh.read(), f)

    def test_nothing_is_sent_without_consent(self):
        rc, out, err, log = self.render(lambda m, c, k: pcm_for(c))
        self.assertEqual(rc, 2)
        self.assertEqual((log["clients"], log["calls"]), (0, []))
        self.assertIn(mp.TERMS_URL, err)
        self.assertIn("--cloud-ok", err)

    def test_pro_first_then_fallback_and_the_key_never_leaks(self):
        def behaviour(model, contents, config):
            if "pro" in model:
                raise ApiError(429, "RESOURCE_EXHAUSTED quota limit: 0 for "
                                    f"GenerateRequestsPerDay (key={FAKE_KEY})")
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--cloud-ok", "Brett, 2026-09-25")
        self.assertEqual(rc, 0, err)
        tried = [m for m, _, _ in log["calls"]]
        # read from the API: pro first, then the 3.8 flash model (never a preview or lite one)
        self.assertEqual(tried[:3], ["gemini-9-pro-tts", "gemini-2.5-pro-preview-tts",
                                     "gemini-3.8-flash-tts"])
        self.assertEqual(log["apis"][:3], ["interactions", "generate_content",   # 9 >= 3
                                           "interactions"])
        man = json.loads(read(os.path.join(self.pod, "podcast.json")))
        self.assertEqual(man["tts"]["model"], "gemini-3.8-flash-tts")
        self.assertEqual(man["tts"]["models_used"], ["gemini-3.8-flash-tts"])
        self.assertEqual(man["cloud_tts_consent"], "Brett, 2026-09-25")
        self.assertIn("[key]", err)                                   # scrubbed, not dropped
        self.assert_no_key_anywhere(out, err)
        # what was sent: the transcript's turns (spoken form) as the two hosts, nothing else
        _, contents, cfg = log["calls"][-1]
        parts = contents[0]["content"]
        self.assertEqual(contents[0]["type"], "user_input")
        self.assertIn("Ryr2 tops the RyR pulldown", parts[0]["text"])
        self.assertEqual(parts[0]["annotations"], [{"type": "speech_metadata", "speaker": "Leo",
                                                   "style": "calm, precise and dryly funny"}])
        self.assertTrue(any("K V two point one" in p["text"] for p in parts))
        self.assertNotIn("Contact-site interactomes", json.dumps(contents))   # no report text
        self.assertEqual(cfg["response_format"], {"type": "audio"})
        self.assertEqual(cfg["generation_config"], {"speech_config": {
            "mode": "conversational", "speakers": [{"speaker": "Maya", "voice": "Kore"},
                                                   {"speaker": "Leo", "voice": "Charon"}]}})

    def test_a_failure_names_the_chunk_keeps_the_cache_and_resumes_on_the_same_model(self):
        state = {"n": 0}

        def flaky(model, contents, config):
            state["n"] += 1
            if state["n"] >= 2:
                raise ApiError(429, f"RESOURCE_EXHAUSTED PerDay limit reached key={FAKE_KEY}")
            return pcm_for(contents)
        rc, out, err, log = self.render(flaky, "--cloud-ok", models=("gemini-2.5-pro-preview-tts",
                                                                     "gemini-2.5-flash-preview-tts"))
        self.assertEqual(rc, 1)
        self.assertIn("chunk 2/2", err)
        self.assertIn("out of its daily quota", err)
        self.assert_no_key_anywhere(out, err)
        self.assertEqual(len(os.listdir(os.path.join(self.pod, ".cache"))), 1)
        first = log["calls"][0][0]
        # next day: every model works again; the cached chunk is reused, chunk 2 on the SAME model
        rc, out, err, log = self.render(lambda m, c, k: pcm_for(c), "--cloud-ok",
                                        models=("gemini-2.5-pro-preview-tts",
                                                "gemini-2.5-flash-preview-tts"))
        self.assertEqual(rc, 0, err)
        self.assertEqual([m for m, _, _ in log["calls"]], [first])
        self.assertIn("resuming with", err)

    def test_an_error_carrying_the_key_is_scrubbed_everywhere(self):
        def boom(model, contents, config):
            raise RuntimeError(f"500 INTERNAL upstream said x-goog-api-key={FAKE_KEY} {FAKE_KEY}")
        rc, out, err, log = self.render(boom, "--cloud-ok", "--model", "gemini-2.5-flash-preview-tts")
        self.assertEqual(rc, 1)
        self.assertIn("chunk 1/2", err)
        self.assertEqual(len(log["calls"]), 4)                          # retried, then gave up
        self.assert_no_key_anywhere(out, err)

    def test_model_listing_failure_falls_back_to_known_names(self):
        rc, out, err, log = self.render(lambda m, c, k: pcm_for(c), "--cloud-ok",
                                        list_error=RuntimeError(f"offline {FAKE_KEY}"))
        self.assertEqual(rc, 0, err)
        self.assertEqual(log["calls"][0][0], mp.KNOWN_TTS_MODELS[0])
        self.assert_no_key_anywhere(out, err)

    def clock(self):
        """A fake clock that sleeping advances: waits and pacing become exact numbers."""
        t, slept = [1000.0], []

        def sleep(sec):
            slept.append(round(sec, 3))
            t[0] += sec
        self.nosleep.stop()
        patches = [mock.patch.object(mp, "_sleep", sleep), mock.patch.object(mp, "_now", lambda: t[0])]
        for q in patches:
            q.start()
        self.addCleanup(lambda: [q.stop() for q in patches] and None)
        self.addCleanup(self.nosleep.start)
        return slept

    def test_per_minute_429_waits_retry_delay_plus_5_and_keeps_the_model(self):
        slept = self.clock()
        state = {"n": 0}

        def limited(model, contents, config):
            state["n"] += 1
            if state["n"] in (1, 2):
                raise ApiError(429, f"RESOURCE_EXHAUSTED key={FAKE_KEY}",
                               rpc_error(PER_MINUTE, "7s"))
            return pcm_for(contents)
        rc, out, err, log = self.render(limited, "--cloud-ok", "--model",
                                        "gemini-2.5-flash-preview-tts")
        self.assertEqual(rc, 0, err)
        self.assertEqual(err.count("[render] rate-limited, waiting 12s"), 2)
        self.assertIn(PER_MINUTE, err)
        self.assertEqual(slept.count(12.0), 2)                        # retryDelay 7 s + 5 s
        self.assertEqual({m for m, _, _ in log["calls"]}, {"gemini-2.5-flash-preview-tts"})
        self.assertEqual(len(log["calls"]), 4)                        # 2 waits, then 2 chunks
        self.assert_no_key_anywhere(out, err)

    def test_the_interactions_error_shape_is_read_too(self):
        slept = self.clock()
        state = {"n": 0}

        def limited(model, contents, config):
            state["n"] += 1
            if state["n"] == 1:
                raise ApiError(429, "Too Many Requests", rpc_error(PER_MINUTE, "30s"),
                               style="interactions")
            return pcm_for(contents)
        rc, out, err, log = self.render(limited, "--cloud-ok", "--model", "gemini-3.8-flash-tts")
        self.assertEqual(rc, 0, err)
        self.assertIn("rate-limited, waiting 35s", err)
        self.assertIn(35.0, slept)

    def test_rate_limit_waits_are_capped_per_chunk(self):
        self.clock()

        def always(model, contents, config):
            raise ApiError(429, "RESOURCE_EXHAUSTED", rpc_error(PER_MINUTE, "10s"))
        rc, out, err, log = self.render(always, "--cloud-ok", "--model",
                                        "gemini-2.5-flash-preview-tts")
        self.assertEqual(rc, 1)
        self.assertEqual(len(log["calls"]), mp.RATE_WAITS + 1)
        self.assertIn(f"still rate-limited on gemini-2.5-flash-preview-tts after {mp.RATE_WAITS} "
                      "waits", err)
        self.assertIn("re-run the same command to resume", err)

    def test_a_daily_quota_stops_with_the_resume_message_and_no_wait(self):
        slept = self.clock()
        state = {"n": 0}

        def daily(model, contents, config):
            state["n"] += 1
            if state["n"] == 2:
                raise ApiError(429, "RESOURCE_EXHAUSTED", rpc_error(PER_DAY, "20s"))
            return pcm_for(contents)
        rc, out, err, log = self.render(daily, "--cloud-ok", "--model",
                                        "gemini-2.5-flash-preview-tts")
        self.assertEqual(rc, 1)
        self.assertIn("chunk 2/2", err)
        self.assertIn("out of its daily quota (quota " + PER_DAY + "; retry in 20s;", err)
        self.assertIn("re-run the same command later to resume", err)
        self.assertNotIn("rate-limited", err)
        self.assertEqual([x for x in slept if x != 20.0], [])         # pacing only, no 25 s wait
        self.assertEqual(len(os.listdir(os.path.join(self.pod, ".cache"))), 1)

    def test_the_quota_id_leads_every_error_line(self):
        long = "You exceeded your current quota, please check your plan and billing details. " * 5

        def pro_limited(model, contents, config):
            if "pro" in model:
                raise ApiError(429, long, rpc_error(PER_DAY, "20s"))
            return pcm_for(contents)
        rc, out, err, log = self.render(pro_limited, "--cloud-ok", "--model",
                                        "gemini-2.5-pro-preview-tts", "--model",
                                        "gemini-3.8-flash-tts")
        self.assertEqual(rc, 0, err)
        self.assertIn(f"gemini-2.5-pro-preview-tts: no quota (quota {PER_DAY}; retry in 20s; "
                      "ApiError: 429 You exceeded", err)

        def all_limited(model, contents, config):
            raise ApiError(429, long, rpc_error(PER_DAY))
        shutil.rmtree(os.path.join(self.pod, ".cache"))
        rc, out, err, log = self.render(all_limited, "--cloud-ok", "--model",
                                        "gemini-2.5-pro-preview-tts", "--model",
                                        "gemini-3.8-flash-tts")
        self.assertEqual(rc, 1)
        self.assertIn(f"last: quota {PER_DAY}; ApiError: 429", err)

    def test_a_text_only_429_uses_its_retry_in(self):
        slept = self.clock()
        state = {"n": 0}

        def limited(model, contents, config):
            state["n"] += 1
            if state["n"] == 1:
                raise RuntimeError("429 RESOURCE_EXHAUSTED. Please retry in 36.15s.")
            return pcm_for(contents)
        rc, out, err, log = self.render(limited, "--cloud-ok", "--model",
                                        "gemini-2.5-flash-preview-tts")
        self.assertEqual(rc, 0, err)
        self.assertIn("rate-limited, waiting 41s", err)

    def test_requests_to_one_model_are_paced(self):
        slept = self.clock()
        rc, out, err, log = self.render(lambda m, c, k: pcm_for(c), "--cloud-ok", "--model",
                                        "gemini-2.5-flash-preview-tts")
        self.assertEqual(rc, 0, err)
        self.assertEqual(slept, [20.0])                               # before the 2nd request
        self.assertIn("pacing: 20s before the next request", err)
        slept[:] = []
        shutil.rmtree(os.path.join(self.pod, ".cache"))
        rc, out, err, log = self.render(lambda m, c, k: pcm_for(c), "--cloud-ok", "--model",
                                        "gemini-2.5-flash-preview-tts", "--min-interval", "0")
        self.assertEqual((rc, slept), (0, []))

    def test_3x_models_use_interactions_with_per_turn_speakers_and_styles(self):
        write(self.script, read(self.script).replace(
            "Voices: Maya=Kore, Leo=Charon",
            "Voices: Maya=Kore, Leo=Charon\nStyles: Maya: excited, fast; Leo: deadpan"))
        self.assertEqual(self.check(read(self.script))[0], 0)
        rc, out, err, log = self.render(lambda m, c, k: pcm_for(c), "--cloud-ok", "--model",
                                        "gemini-3.8-flash-tts")
        self.assertEqual(rc, 0, err)
        self.assertEqual(log["apis"], ["interactions", "interactions"])
        parts = [p for _, inp, _ in log["calls"] for p in inp[0]["content"]]
        self.assertEqual(len(parts), 7)                               # one part per turn
        self.assertEqual({p["annotations"][0]["speaker"]: p["annotations"][0]["style"]
                          for p in parts}, {"Maya": "excited, fast", "Leo": "deadpan"})
        self.assertTrue(all(set(p) == {"type", "text", "annotations"} for p in parts))
        man = json.loads(read(os.path.join(self.pod, "podcast.json")))
        self.assertEqual(man["tts"]["model"], "gemini-3.8-flash-tts")
        # the WAV header is read, not spoken: audio length = the PCM the fake made
        speech = sum(len(pcm_for([{"content": [p]}])) // 2 for p in parts) / float(mp.RATE)
        self.assertGreater(man["duration_s"], speech)
        self.assertLess(man["duration_s"], speech + 6)

    def test_a_rejected_request_before_the_first_chunk_hands_over(self):
        def reject_3x(model, contents, config):
            if model.startswith("gemini-3"):
                raise ApiError(400, "INVALID_ARGUMENT Multi-speaker generation requests must "
                                    "specify speaker names for each part in the contents.")
            return pcm_for(contents)
        rc, out, err, log = self.render(reject_3x, "--cloud-ok", "--model", "gemini-3.8-flash-tts",
                                        "--model", "gemini-2.5-flash-preview-tts")
        self.assertEqual(rc, 0, err)
        self.assertIn("gemini-3.8-flash-tts: request rejected", err)
        man = json.loads(read(os.path.join(self.pod, "podcast.json")))
        self.assertEqual(man["tts"]["models_used"], ["gemini-2.5-flash-preview-tts"])

    def test_2x_models_keep_generate_content(self):
        rc, out, err, log = self.render(lambda m, c, k: pcm_for(c), "--cloud-ok", "--model",
                                        "gemini-2.5-flash-preview-tts")
        self.assertEqual(rc, 0, err)
        self.assertEqual(log["apis"], ["generate_content", "generate_content"])
        _, contents, cfg = log["calls"][-1]
        self.assertTrue(contents.startswith("TTS the following conversation between Maya and Leo, "
                                            "two hosts of a lively, natural science podcast. Maya "
                                            "sounds curious and energetic; Leo sounds calm, "
                                            "precise and dryly funny:\n\nLeo: Ryr2 tops"))
        voices = {v.speaker: v.voice_config.prebuilt_voice_config.voice_name
                  for v in cfg.speech_config.multi_speaker_voice_config.speaker_voice_configs}
        self.assertEqual(voices, {"Maya": "Kore", "Leo": "Charon"})
        self.assertEqual(cfg.response_modalities, ["AUDIO"])

    def test_error_info(self):
        e = mp.error_info(ApiError(429, "x", rpc_error(PER_MINUTE, "37s")))
        self.assertEqual((e["kind"], e["delay"], e["quota"]), ("rate", 37.0, PER_MINUTE))
        e = mp.error_info(ApiError(429, "x", json.dumps(rpc_error(PER_DAY, "20s")), "interactions"))
        self.assertEqual((e["kind"], e["quota"]), ("exhausted", PER_DAY))
        self.assertEqual(mp.error_info(ApiError(429, "x", rpc_error(PER_MINUTE, value="0")))["kind"],
                         "exhausted")                                  # a zero quota never refills
        self.assertEqual(mp.error_info(ApiError(429, "x", rpc_error(PER_MINUTE)))["kind"], "rate")
        self.assertEqual(mp.error_info(ApiError(404, "NOT_FOUND"))["kind"], "missing")
        self.assertEqual(mp.error_info(ApiError(400, "INVALID_ARGUMENT"))["kind"], "rejected")
        self.assertEqual(mp.error_info(ApiError(403, "PERMISSION_DENIED"))["kind"], "fatal")
        self.assertEqual(mp.error_info(RuntimeError("503 UNAVAILABLE"))["kind"], "transient")
        self.assertEqual(mp.GeminiTTS.api("gemini-3.8-flash-tts"), "interactions")
        self.assertEqual(mp.GeminiTTS.api("models/gemini-3.1-flash-tts-preview"), "interactions")
        self.assertEqual(mp.GeminiTTS.api("gemini-2.5-pro-preview-tts"), "generate_content")

    def test_scrub(self):
        mp._SECRETS.append("sekrit-value")
        s = mp.scrub(f"a sekrit-value b {FAKE_KEY} ?key=abc&x=1 (key=def)")
        self.assertEqual(s, "a [key] b [key] ?key=[key]&x=1 (key=[key])")


# -------------------------------------------------------------------------------------- link
MANIFEST = {"show": "Signal to Noise", "title": "Contact sites", "audio": "podcast.m4a",
            "transcript": "transcript.html", "duration_s": 1140.2, "ai_generated": True}

REPORT_HTML = ('<!doctype html><html><head><title>R</title></head><body>'
               '<header class="band"><div class="band-in"><h1>R</h1></div></header>'
               '<div class="layout"><details class="toc"></details><main id="main">'
               '<section id="sec-x"><h2 id="x">X</h2><p>body</p></section></main></div>'
               '</body></html>')


class Link(Workspace):
    def setUp(self):
        super().setUp()
        write(os.path.join(self.pod, "podcast.json"), json.dumps(MANIFEST))
        write(os.path.join(self.out, "Analysis_Report.html"), REPORT_HTML)
        write(os.path.join(self.d, "README.html"),
              '<main class="doc"><nav><h2>Start here</h2><ul><li><a href="output/Analysis_Report'
              '.html">Analysis report</a> — figures</li><li><a href="output/methods.md">M</a></li>'
              '</ul></nav></main>')
        write(os.path.join(self.d, "README.md"),
              "# S\n\n## Start here\n\n- [Analysis report](output/Analysis_Report.html) — x\n"
              "- [M](output/methods.md)\n")
        write(os.path.join(self.d, "AGENTS.md"), "# AGENTS\n\n## Where this lives on HIVE\n\nt\n\n"
                                                 "## Do not\n\n- x\n")

    def files(self):
        return {p: read(p) for p in (os.path.join(self.out, "Analysis_Report.html"), self.report,
                                     os.path.join(self.d, "README.html"),
                                     os.path.join(self.d, "README.md"),
                                     os.path.join(self.d, "AGENTS.md"))}

    def test_link_is_idempotent_and_near_the_top(self):
        rc, out, err = run("link", self.out)
        self.assertEqual(rc, 0, err)
        first = self.files()
        for path, text in first.items():
            self.assertEqual((text.count(mp.START), text.count(mp.END)), (1, 1), path)
        for _ in range(2):
            self.assertEqual(run("link", self.out)[0], 0)
        self.assertEqual(self.files(), first)                         # byte-identical re-runs
        page = first[os.path.join(self.out, "Analysis_Report.html")]
        self.assertIn('<main id="main">' + mp.START, page)           # top of the reading column
        self.assertIn('<audio controls preload="none" src="podcast/podcast.m4a">', page)
        self.assertIn('href="podcast/transcript.html"', page)
        self.assertIn("about 19 min", page)
        self.assertIn("@media print{.pc-card audio{display:none}.pc-card .pc-print{display:block}",
                      page)
        self.assertIn("Audio: podcast/podcast.m4a", page)
        md = first[self.report]
        self.assertTrue(md.startswith("# Contact-site interactomes in Old and Young mouse brain\n\n"
                                      + mp.START + "\n> **Listen:**"), md[:200])
        self.assertIn("[podcast/podcast.m4a](podcast/podcast.m4a)", md)
        readme = first[os.path.join(self.d, "README.html")]
        self.assertRegex(readme, r"Analysis report</a> — figures</li>" + re.escape(mp.START)
                         + r'<li><a href="output/podcast/podcast.m4a">')
        rmd = first[os.path.join(self.d, "README.md")]
        self.assertIn("— x\n" + mp.START + "\n- [Audio discussion of these results]"
                      "(output/podcast/podcast.m4a)", rmd)
        agents = first[os.path.join(self.d, "AGENTS.md")]
        self.assertLess(agents.index(mp.START), agents.index("## Do not"))
        self.assertIn("NOT authoritative", agents)
        self.assertIn("`output/podcast/podcast_script.md`", agents)

    def test_md_hook_matches_link(self):
        md = "# Title\n\nStandfirst.\n\n## A\n"
        once = mp.add_listen_md(md, self.out)
        self.assertEqual(mp.add_listen_md(once, self.out), once)
        run("link", self.out)
        self.assertEqual(read(self.report).split("## Overview")[0].split(mp.START)[1],
                         mp.add_listen_md(REPORT, self.out).split("## Overview")[0].split(mp.START)[1])
        self.assertEqual(mp.strip_block(once).replace("\n\n\n", "\n\n"), md)
        self.assertEqual(mp.add_listen_md(md, self.d), md)             # no podcast: unchanged

    def test_an_older_pdf_is_reprinted_or_flagged(self):
        html_path = os.path.join(self.out, "Analysis_Report.html")
        pdf = os.path.join(self.out, "Analysis_Report.pdf")
        write(pdf, "%PDF old")
        os.utime(pdf, (1, 1))                                 # older than the HTML
        with mock.patch.dict(sys.modules, {"html_to_pdf": None}):     # not on this branch
            rc, out, err = run("link", self.out)
        self.assertEqual(rc, 0, err)
        self.assertIn("[INFO] output/Analysis_Report.pdf: older than the HTML, so it has no "
                      "Listen card; to reprint it, open Analysis_Report.html in a browser", out)
        calls = []

        def convert(h, p, **kw):
            calls.append((h, p))
            self.assertIn(mp.START, read(h))                  # printed AFTER the card went in
            write(p, "%PDF new")
            return True, "3 pages, printed by fake"
        fake = types.ModuleType("html_to_pdf")
        fake.convert = convert
        os.utime(pdf, (1, 1))
        with mock.patch.dict(sys.modules, {"html_to_pdf": fake}):
            rc, out, err = run("link", self.out)
            self.assertIn("[OK] output/Analysis_Report.pdf: reprinted with the Listen card", out)
            self.assertEqual(calls, [(html_path, pdf)])
            rc, out, err = run("link", self.out)              # up to date now: not reprinted
        self.assertEqual(len(calls), 1)
        self.assertIn("up to date with the HTML", out)
        os.remove(pdf)
        self.assertIn("no PDF beside the report", run("link", self.out)[1])

    def test_missing_files_are_info_and_no_podcast_is_an_error(self):
        for f in ("README.html", "README.md", "AGENTS.md"):
            os.remove(os.path.join(self.d, f))
        rc, out, err = run("link", self.out)
        self.assertEqual(rc, 0, err)
        self.assertIn("[INFO] README.html: not there; skipped", out)
        self.assertIn("[OK] output/Analysis_Report.html: Listen card added", out)
        os.remove(os.path.join(self.pod, "podcast.json"))
        self.assertEqual(run("link", self.out)[0], 2)


class ReportGenerator(Workspace):
    """The real make_analysis_html.py output: link adds one card, and regenerating the report
    keeps exactly one (its hook), with the Markdown twin's Listen line never printed as text."""

    def make_report(self):
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"),
                            "--session", self.d, "--out",
                            os.path.join(self.out, "Analysis_Report.html")],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        return read(os.path.join(self.out, "Analysis_Report.html"))

    def test_generated_report(self):
        write(self.report, "# Study\n\nA one-line standfirst.\n\n## Overview\n\nText 6,112.\n")
        page = self.make_report()
        self.assertNotIn(mp.START, page)                              # no podcast, no card
        write(os.path.join(self.pod, "podcast.json"), json.dumps(MANIFEST))
        self.assertEqual(run("link", self.out)[0], 0)
        self.assertEqual(run("link", self.out)[0], 0)
        page = read(os.path.join(self.out, "Analysis_Report.html"))
        self.assertEqual(page.count(mp.START), 1)
        self.assertIn('<main id="main">' + mp.START, page)
        self.assertIn(mp.START, read(self.report))                    # the md twin has its line
        page = self.make_report()                                     # regenerate: card kept
        self.assertEqual(page.count(mp.START), 1)
        self.assertEqual(page.count('class="pc-card"'), 1)
        self.assertNotIn("podcast:start --&gt;", page)                # never printed as text
        self.assertNotIn("Listen:**", page)
        self.assertIn('<p class="subtitle">A one-line standfirst.</p>', page)   # still the standfirst


class SessionFiles(unittest.TestCase):
    """session_docs.py lists the podcast itself (so finalize keeps it), link then adds nothing,
    and the session zip carries the podcast but not its TTS cache."""

    def test_readme_agents_and_zip(self):
        import session_docs
        import test_deposit_package as tdp
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            pod = os.path.join(p["output_dir"], "podcast")
            write(os.path.join(p["output_dir"], "Analysis_Report.html"), REPORT_HTML)
            write(os.path.join(pod, "podcast.json"), json.dumps(MANIFEST))
            for f in ("podcast.m4a", "transcript.html", "podcast_script.md", "check.txt"):
                write(os.path.join(pod, f), "x" * 2000)
            for i in range(3):
                write(os.path.join(pod, ".cache", f"{i:064x}.wav"), "w" * 100)
            session_docs.write_docs(p["session_dir"])
            readme = read(p["readme"])
            self.assertEqual(readme.count("(output/podcast/podcast.m4a)"), 1)
            starts = [ln for ln in readme.splitlines() if ln.startswith("- [")]
            self.assertIn("Analysis_Report.html", starts[0])
            self.assertIn("Audio discussion of these results", starts[1])   # right after it
            agents = read(os.path.join(p["session_dir"], "AGENTS.md"))
            self.assertIn("a derivative, not a record", agents)
            self.assertLess(agents.index("Audio discussion (podcast)"), agents.index("## Do not"))
            rc, out, err = run("link", p["output_dir"])
            self.assertEqual(rc, 0, err)
            self.assertIn("already listed (session_docs.py)", out)
            for f in ("README.html", "README.md", "AGENTS.md"):
                self.assertNotIn(mp.START, read(os.path.join(p["session_dir"], f)), f)

            r = tdp.finalize(p["session_dir"], "--zip")
            self.assertEqual(r.returncode, 0, r.stderr)
            res = json.loads(r.stdout)
            names = zipfile.ZipFile(res["zip"]).namelist()
            self.assertTrue(any(n.endswith("output/podcast/podcast.m4a") for n in names))
            self.assertFalse([n for n in names if "/.cache/" in n], names)
            self.assertEqual(res["zip_excluded"]["output/podcast/.cache (TTS scratch, kept on "
                                                  "disk)"], 3)
            readme = read(p["readme"])                                   # finalize rewrote it
            self.assertEqual(readme.count("(output/podcast/podcast.m4a)"), 1)

    def test_output_files_catalog(self):
        import make_report
        self.assertEqual(make_report.describe("podcast.m4a")[0], "Analysis report")
        self.assertIn("Claims beyond the report", make_report.describe("podcast_script.md")[1])
        with tempfile.TemporaryDirectory() as d:
            write(os.path.join(d, "podcast", "podcast.json"), "{}")
            write(os.path.join(d, "podcast", ".cache", "a.wav"), "w")
            self.assertEqual([os.path.basename(f) for f in make_report.collect([d])],
                             ["podcast.json"])


if __name__ == "__main__":
    unittest.main()
