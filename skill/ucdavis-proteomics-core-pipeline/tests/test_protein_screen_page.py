#!/usr/bin/env python3
"""The protein screen page (scripts/protein_screen.tmpl.html), rendered in headless Chrome.

Pages are built with make_protein_screen.py from the hand-made mini screen data
(tests/fixtures/screen_data_mini: groups A, B, C; contrasts B-A and C-A) and, when R with limpa is
installed, from the synthetic analysis that make_screen_data.R runs on (screen_fixture.py). Each page
is loaded with --headless=new --virtual-time-budget=8000 --dump-dom, the browser found the way
html_to_pdf.py finds it. Every test skips when there is no Chrome, Chromium or Edge.
"""
import copy
import html
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import html_to_pdf  # noqa: E402
import make_protein_screen as mps  # noqa: E402

BROWSER = html_to_pdf.find_browser()
SCRIPT = os.path.join(SCRIPTS, "make_protein_screen.py")
FIXTURE = os.path.join(HERE, "fixtures", "screen_data_mini")

# Words, names, files and dates of the one dataset the page was first drawn for. None may reach a page.
MOCKUP_LITERALS = [r"\bFeb\b", r"\bMar\b", "February", "March", "HeLa", "FL020226", "FL100326", "BASP1", "SURF6",
                   "POTEE", "P80723", "DE_dpc_Mar.Feb", "2026-10-26", "mockup", "Mockup"]
MISSING_BLOCKS = ("xic", "check_routes", "fran", "report_text")

SEQ = "MK" + "AEEATEAQEVVEATPEGACTEPR" + "GG" + "LGSEQGGGK" + "AAAA"
UNIPROT = {"acc": "P10001", "id": "AAA_HUMAN", "reviewed": "UniProtKB reviewed (Swiss-Prot)", "name": "Protein A",
           "gene": "GENEA", "organism": "Homo sapiens", "length": len(SEQ), "sequence": SEQ,
           "function": "Does a thing, as UniProt words it.", "location": ["Cytoplasm"], "uniprot_version": 3,
           "uniprot_updated": "2026-09-02", "release": "2026_03", "fetched_at": "2026-10-28"}


def dump_dom(html_path, hash_=""):
    """The page's DOM after its scripts ran, from headless Chrome with a throwaway profile.

    `hash_` (e.g. "#open=profile") opens folded sections, as the page's own #open= option does.
    """
    profile = tempfile.mkdtemp(prefix="screen_chrome_")
    try:
        p = subprocess.run([BROWSER, "--headless=new", "--disable-gpu", "--no-first-run", "--no-default-browser-check",
                            "--user-data-dir=" + profile, "--virtual-time-budget=8000", "--dump-dom",
                            "file://" + os.path.abspath(html_path) + hash_],
                           capture_output=True, text=True, timeout=180)
    finally:
        shutil.rmtree(profile, ignore_errors=True)
    if not p.stdout.strip():
        raise AssertionError("Chrome returned no DOM:\n" + p.stderr[-2000:])
    return p.stdout


def render_state(dom):
    m = re.search(r'<body[^>]*data-render="([^"]*)"', dom)
    return "data-render=" + (m.group(1) if m else "(not set)")


ALL_OPEN = "#open=weights,profile,sankey,health,picture"


def make_page(*args, screen_data=FIXTURE):
    """make_protein_screen.py --offline -> the page's path."""
    out = os.path.join(tempfile.mkdtemp(prefix="screen_page_"), "page.html")
    p = subprocess.run([sys.executable, SCRIPT, "--screen-data", screen_data, "--out", out, "--offline", *args],
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise AssertionError("make_protein_screen.py failed (%d):\n%s" % (p.returncode, p.stderr[-3000:]))
    return out


def build_page(*args, screen_data=FIXTURE, hash_=""):
    """make_protein_screen.py --offline -> (html text, dom)."""
    out = make_page(*args, screen_data=screen_data)
    with open(out, encoding="utf-8") as fh:
        html = fh.read()
    return html, dump_dom(out, hash_)


def page_from_D(D, hash_=""):
    """Render D into the template directly (no CLI) -> dom."""
    with open(mps.TEMPLATE, encoding="utf-8") as fh:
        page = mps.render(fh.read(), D)
    out = os.path.join(tempfile.mkdtemp(prefix="screen_page_"), "page.html")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(page)
    return dump_dom(out, hash_)


def mini_D(protein="P10001", contrast=None, uniprot=None):
    man = mps.load_manifest(FIXTURE)
    idx = mps.parse_index(os.path.join(FIXTURE, "index.tsv"))
    row = next(r for r in idx if r["Protein.Group"] == protein)
    rec = mps.read_record(FIXTURE, row, man)
    return mps.page_data(man, rec, contrast, uniprot, index_rows=idx,
                         uniprot_reason=None if uniprot else "offline: UniProt was not consulted")


def rendered(dom):
    """The DOM's markup without <script> and <style> contents.

    The dumped DOM still holds the page's script, and with it the data D and every string the script
    could write, so checks for rendered output must not look there.
    """
    return re.sub(r"<(script|style)\b.*?</\1>", " ", dom, flags=re.S | re.I)


def visible_text(dom):
    """The rendered markup without tags, whitespace collapsed: roughly what a reader sees."""
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", rendered(dom)))


def story(dom):
    """The visible text of the short answer."""
    m = re.search(r'<div class="qt-story".*?</ol>', rendered(dom), re.S)
    return visible_text(m.group(0)) if m else ""


class PageChecks:
    """The brief's checks, for a page of an analysis with groups A, B and C, showing the comparison B-A."""

    def test_renders_ok(self):
        self.assertIn('data-render="ok"', self.dom, render_state(self.dom))

    def test_no_mockup_literals(self):
        for pat in MOCKUP_LITERALS:
            m = re.search(pat, self.dom)
            self.assertIsNone(m, "%s: ...%s..." % (pat, self.dom[max(0, m.start() - 80):m.end() + 80] if m else ""))

    def test_names_this_analysis(self):
        text = visible_text(self.dom)
        self.assertIn("B vs A", text)
        self.assertIn("DE_dpc_B.A.csv", text)
        # where the page names the table its numbers were read from (the footer also lists sources by name)
        self.assertIn("Comparison B vs A (tables/DE_dpc_B.A.csv)", text)
        self.assertIn("Read from tables/DE_dpc_B.A.csv \u2192 logFC, adj.P.Val", text)

    def test_missing_blocks_say_so(self):
        for block in MISSING_BLOCKS:
            self.assertRegex(rendered(self.dom), r'data-missing="%s">[^<]*<b>[^<]+</b> Not in this version of the screen data: '
                             % block, block)

    def test_three_groups_render(self):
        text = visible_text(self.dom)
        self.assertIn("not in this comparison", text)            # group C's samples
        self.assertIn("C 1 not in this comparison", text)        # the precursor table labels each of C's samples
        self.assertIn("Samples of C are not in this comparison", text)

    def test_interval_is_not_said_to_come_from_the_de_table(self):
        """The DE tables have no CI.L/CI.R (run_de.R's topTable has no confint); the page shows no interval."""
        text = visible_text(self.dom)
        self.assertNotIn("CI.L", text)
        self.assertNotIn("95% interval", text)

    def test_csp_no_network(self):
        self.assertIn("default-src 'none'", self.html)

    def test_shares_are_not_called_a_nudge(self):
        text = visible_text(self.dom)
        self.assertIsNone(re.search(r"nudg", text, re.I), "the page describes shares as a nudge")
        self.assertIn("checked against direct recalculation", text)

    def test_wording_rules(self):
        """Spec section 7: never 'no measurement' for a gap; 'inferred' only for protein values."""
        text = visible_text(self.dom)
        self.assertNotIn("no measurement", text.lower())
        for m in re.finditer(r"inferred", text):
            ctx = text[max(0, m.start() - 60):m.end() + 60].lower()
            self.assertNotRegex(ctx, r"(peptide|precursor)s? (was|were|is|are) inferred", ctx)


@unittest.skipUnless(BROWSER, "no Chrome, Chromium or Edge found")
class MiniPage(PageChecks, unittest.TestCase):
    """GENEA (P10001) from the mini fixture, offline: the default comparison B-A."""

    @classmethod
    def setUpClass(cls):
        path = make_page("--protein", "GENEA")
        with open(path, encoding="utf-8") as fh:
            cls.html = fh.read()
        cls.dom, cls.dom_open = dump_dom(path), dump_dom(path, ALL_OPEN)

    def test_title_names_the_gene(self):
        self.assertIn("<title>GENEA: where its number comes from</title>", self.dom)

    def test_offline_card_says_uniprot_is_missing(self):
        self.assertRegex(rendered(self.dom), r'data-missing="uniprot">Not on this page: offline: UniProt was not consulted')

    def test_cutoffs_come_from_the_analysis(self):
        text = visible_text(self.dom)
        self.assertIn("Q.Value 1%", text)
        self.assertIn("Lib.Q.Value 1%", text)

    def test_colours_by_role(self):
        """Drawn marks, not the script: the protein-value dots of B (numerator), A (denominator), C (other)."""
        r = rendered(self.dom_open)
        for role in ("num", "den", "other"):
            self.assertRegex(r, r'<circle [^>]*r="6" fill="var\(--%s\)"' % role, role)

    def test_weights_said_even_when_they_are(self):
        self.assertIn("about evenly", story(self.dom))           # fixture: 50% each, 2 per group


@unittest.skipUnless(BROWSER, "no Chrome, Chromium or Edge found")
class MiniPageVariants(unittest.TestCase):
    def test_other_comparison(self):
        html, dom = build_page("--protein", "GENEA", "--contrast", "C-A")
        self.assertIn('data-render="ok"', dom, render_state(dom))
        text = visible_text(dom)
        self.assertIn("C vs A", text)
        self.assertIn("DE_dpc_C.A.csv", text)
        self.assertIn("Read from tables/DE_dpc_C.A.csv \u2192 logFC, adj.P.Val", text)
        self.assertIn("B 1 not in this comparison", text)        # now group B's samples
        self.assertNotIn("DEFAULT", text)                         # chosen, and the quantity is recorded

    def test_default_comparison_is_tagged(self):
        # repo rule 2: 2 comparisons and no --contrast -> the first, and the page says nobody chose it
        text = visible_text(build_page("--protein", "GENEA")[1])
        self.assertIn("Comparison B vs A (tables/DE_dpc_B.A.csv) · the first comparison, chosen by default "
                      "(DEFAULT — not user-confirmed)", text)
        self.assertEqual(text.count("DEFAULT"), 1)

    def test_unrecorded_quantity_is_tagged(self):
        # repo rule 2: no `quantities` in the analysis record -> the page falls back to Precursor.Normalised, tagged
        D = mini_D(contrast="B-A")
        D["analysis"]["quantities"] = None
        text = visible_text(page_from_D(D, ALL_OPEN))
        tag = ("Precursor.Normalised (DEFAULT — not user-confirmed: this screen data does not record which "
               "quantity the analysis read)")
        self.assertIn("log2 of DIA-NN’s " + tag, text)
        self.assertIn("limpa uses the log2 of the result (DEFAULT — not user-confirmed", text)
        D["contrib"] = {}                  # no shares: the Sankey opens on its intensity view, which names the quantity
        self.assertIn("the quantity limpa read, " + tag, visible_text(page_from_D(D, ALL_OPEN)))
        self.assertNotIn("DEFAULT", visible_text(page_from_D(mini_D(contrast="B-A"), ALL_OPEN)))

    def test_untested_protein_with_withheld_shares(self):
        html, dom = build_page("--protein", "P10004")
        self.assertIn('data-render="ok"', dom, render_state(dom))
        self.assertIn("not tested in this comparison", visible_text(dom))
        self.assertRegex(rendered(dom), r'data-missing="shares"><b>Shares\.</b> Withheld for this protein: '
                                        r'refit differs from the delivered value by 5\.15e-05')

    def test_cutoffs_follow_the_manifest(self):
        d = os.path.join(tempfile.mkdtemp(), "screen_data")
        shutil.copytree(FIXTURE, d, ignore=shutil.ignore_patterns("make.py", "__pycache__"))
        mp = os.path.join(d, "manifest.json")
        with open(mp) as fh:
            m = json.load(fh)
        m["analysis"]["q_cutoffs"] = [0.05, 0.005]
        with open(mp, "w") as fh:
            json.dump(m, fh)
        html, dom = build_page("--protein", "GENEA", screen_data=d)
        text = visible_text(dom)
        self.assertIn("Q.Value 5%", text)
        self.assertIn("Lib.Q.Value 0.5%", text)

    def test_uniprot_card_is_verbatim(self):
        dom = page_from_D(mini_D(uniprot=copy.deepcopy(UNIPROT)))
        self.assertIn('data-render="ok"', dom, render_state(dom))
        self.assertNotIn('data-missing="uniprot"', rendered(dom))
        text = visible_text(dom)
        for fact in ("Does a thing, as UniProt words it.", "Location: Cytoplasm", "UniProtKB reviewed (Swiss-Prot)",
                     "AAA_HUMAN", "Homo sapiens", "release 2026_03", "fetched 2026-10-28", "annotation updated 2026-09-02"):
            self.assertIn(fact, text, fact)
        self.assertIn('href="https://www.uniprot.org/uniprotkb/P10001/entry"', rendered(dom))

    @staticmethod
    def group_D(pg, facts=None):
        """GENEA's page data, its record filed as protein group `pg`, with UniProt answering for it."""
        man = mps.load_manifest(FIXTURE)
        idx = mps.parse_index(os.path.join(FIXTURE, "index.tsv"))
        rec = mps.read_record(FIXTURE, next(r for r in idx if r["Protein.Group"] == "P10001"), man)
        rec["pg"] = pg
        return mps.page_data(man, rec, None, copy.deepcopy(UNIPROT), index_rows=idx, facts=facts)

    def test_isoform_page_has_no_coverage_and_says_the_canonical_entry_is_shown(self):
        note = "P10001-2 is an isoform, which has no UniProt entry of its own; the canonical entry, P10001, is shown"
        fact = {"source": "UniProt", "status": "ok", "url": "https://rest.uniprot.org/uniprotkb/P10001.json",
                "release": "2026_03", "fetched_at": "2026-10-28", "licence": "CC BY 4.0",
                "fields": [["Accession", "P10001"]], "note": note}
        dom = page_from_D(self.group_D("P10001-2", [fact]))
        self.assertIn('data-render="ok"', dom, render_state(dom))
        r = rendered(dom)
        self.assertNotIn('aria-label="Sequence coverage"', r)
        self.assertRegex(r, r'data-missing="coverage"><b>Sequence coverage and peptide positions\.</b> Not shown: '
                            r"this protein group's accession is the isoform P10001-2")
        self.assertIn(note, visible_text(dom))

    def test_multi_accession_card_says_the_first_accessions_entry_is_shown(self):
        dom = page_from_D(self.group_D("P10001;P10009"))
        self.assertIn('data-render="ok"', dom, render_state(dom))
        self.assertIn("This protein group holds 2 accessions (P10001, P10009); the facts on this card are for the "
                      "first, P10001.", visible_text(dom))
        self.assertIn('aria-label="Sequence coverage"', rendered(dom))
        self.assertNotIn("protein group holds", visible_text(page_from_D(self.group_D("P10001"))))

    def test_xic_block_draws_when_present(self):
        D = mini_D()
        D["missing"].pop("xic", None)
        pid = D["proteins"]["P10001"]["prec_ids"][0]
        trace = {"rt": [12.40, 12.45, 12.50, 12.55, 12.60], "v": [10, 300, 900, 280, 12]}
        D["xic"] = {pid: {s: {"y7": trace, "b3": trace, "ms1": trace} for s in D["samples"]}}
        dom = page_from_D(D)
        self.assertIn('data-render="ok"', dom, render_state(dom))
        self.assertNotIn('data-missing="xic"', rendered(dom))
        self.assertEqual(len(re.findall(r'<div class="qt-xpanel', rendered(dom))), 6)

    def test_the_weights_matter_cut_off_is_stated(self):
        # GENEA: a plain average of A and B gives +0.20; limma's answer is +0.50
        text = visible_text(page_from_D(mini_D()))
        self.assertIn("Here the weights matter (the plain average and limma’s answer differ by more than 0.05 log2). "
                      "A plain average of the 4 values gives +0.20.", text)
        D = mini_D()
        D["proteins"]["P10001"]["logFC"] = 0.22
        text = visible_text(page_from_D(D))
        self.assertIn("Here the weights barely matter (the plain average and limma’s answer are within 0.05 log2): "
                      "a plain average of the 4 values gives +0.200, limma +0.220.", text)

    def test_uneven_weights_are_not_called_even(self):
        D = mini_D()
        D["routes"]["P10001"]["wshare"] = [0.8, 0.2, 0.5, 0.5, 0.5, 0.5]      # no sample inferred
        text = story(page_from_D(D))
        self.assertNotIn("about evenly", text)
        self.assertIn("A 1 carries 80% of the A mean, where an equal share would be 50%", text)

    def test_inferred_values_with_more_weight_are_not_said_to_get_less(self):
        D = mini_D()
        D["proteins"]["P10001"]["nobs"] = [0, 2, 2, 2, 2, 2]                 # A 1 inferred ...
        D["routes"]["P10001"]["wshare"] = [0.6, 0.4, 0.5, 0.5, 0.5, 0.5]      # ... yet weighted more than A 2
        dom = page_from_D(D)
        self.assertNotIn("less weight", story(dom))
        self.assertNotIn("about evenly", story(dom))
        self.assertIn("A 1: 60% of the A mean", story(dom))
        self.assertNotIn("here the inferred value", visible_text(dom))      # step 1 makes no such claim either

    def test_inferred_values_with_less_weight_are_said_to_get_less(self):
        D = mini_D()
        D["proteins"]["P10001"]["nobs"] = [0, 2, 2, 2, 2, 2]
        D["routes"]["P10001"]["wshare"] = [0.1, 0.9, 0.5, 0.5, 0.5, 0.5]
        dom = page_from_D(D)
        self.assertIn("limma gave the inferred values less weight (A 1: 10% of the A mean, against 90% for its found sample)",
                      story(dom))
        self.assertIn("here the inferred value in A 1", visible_text(dom))

    def test_a_broken_page_says_so(self):
        dom = page_from_D(None)
        self.assertRegex(dom, r'data-render="error: [^"]+"')


def share_D(other_b1=0.05, total_b1=0.95):
    """GENEA's page data with shares that add up: 40% + 50% own, 5% other samples, 5% limpa's prior and the gaps
    (c2, with one gap: 40% + 5% + 55%); in B 1 (sample 2) other `other_b1` and total `total_b1`."""
    D = mini_D()
    cb, n = D["contrib"]["P10001"], len(D["samples"])
    cb["own"] = [[0.4] * n, [0.5 if i != 5 else None for i in range(n)]]
    cb["other"], cb["total"] = [0.05] * n, [0.95] * n
    cb["total"][5] = 0.45
    cb["other"][2], cb["total"][2] = other_b1, total_b1
    return D


def share_bars(dom):
    """{sample index: [(width %, title)]}: the slices of each share bar, as drawn."""
    out = {}
    for m in re.finditer(r'<div class="lab qt-hit" data-s="(\d+)".*?<div class="bar">(.*?)<div class="max">', rendered(dom), re.S):
        out[int(m.group(1))] = [(float(w), html.unescape(t)) for w, t in
                                re.findall(r'<div class="seg"[^>]*? style="width:([\d.]+)%;[^"]*" title="([^"]*)"', m.group(2))]
    return out


@unittest.skipUnless(BROWSER, "no Chrome, Chromium or Edge found")
class ShareTotals(unittest.TestCase):
    """A sample's shares add up to 100% only when none is negative; a figure must never misstate them."""

    def check_widths(self, bars, scale):
        """Every slice is as wide as its % label on one scale (100% = 1/scale of the track); no bar overflows."""
        for si, segs in bars.items():
            self.assertLessEqual(sum(w for w, _ in segs), 100.01, si)
            for w, title in segs:
                pct = int(re.search(r"(\d+)%", title).group(1))
                self.assertLessEqual(abs(w - pct / scale), 0.5 / scale + 0.01, (si, w, title))

    def test_a_total_above_100_percent_is_drawn_to_scale_and_said(self):
        dom = page_from_D(share_D(other_b1=0.18, total_b1=1.08))        # limpa's prior and the gaps: -8%
        self.assertIn('data-render="ok"', dom, render_state(dom))
        bars = share_bars(dom)
        self.assertEqual([round(w, 2) for w, _ in bars[2]], [round(100 * x / 1.08, 2) for x in (0.5, 0.4, 0.18)])
        self.check_widths(bars, 1.08)
        self.assertIn('class="mark" style="left:92.59%"', rendered(dom))            # where 100% is
        text = visible_text(dom)
        self.assertIn("Its slices add up to 108%: limpa’s prior and the gaps pull the other way, by 8%.", text)
        self.assertEqual(text.count("pull the other way"), 1)                          # B 1 only
        self.assertNotIn("adding up to 100%", text)
        self.assertIn("its slices add up to that sample’s protein value (100%), except where part of the value "
                      "pulls the other way", text)

    def test_a_negative_other_samples_share_is_not_hidden(self):
        dom = page_from_D(share_D(other_b1=-0.1, total_b1=0.8))         # other -10%, prior and gaps +20%
        bars = share_bars(dom)
        self.assertEqual([round(w, 2) for w, _ in bars[2]], [round(100 * x / 1.1, 2) for x in (0.5, 0.4, 0.2)])
        self.check_widths(bars, 1.1)
        self.assertIn("Its slices add up to 110%: the same peptides found in other samples pull the other way, "
                      "by 10%.", visible_text(dom))
        # B 1 has no light slice (its "other samples" part is negative): the key says when one is drawn
        self.assertNotIn("other samples", " ".join(t for _, t in bars[2]))
        self.assertIn("The light slice, drawn where that part is positive, comes from the same peptides found in "
                      "other samples, and the amber slice, likewise, rests on peptides not found in this sample",
                      visible_text(dom))

    def test_shares_that_add_up_need_no_mark(self):
        dom = page_from_D(share_D())
        self.check_widths(share_bars(dom), 1.0)
        self.assertNotIn('class="mark"', rendered(dom))
        text = visible_text(dom)
        self.assertIn("One bar per sample: its slices add up to that sample’s protein value (100%).", text)
        self.assertNotIn("pull the other way", text)

    def test_the_short_answer_has_no_negative_rest(self):
        D = mini_D("P10002")                                    # one peptide
        cb = D["contrib"]["P10002"]
        self.assertIn("the rest comes from other samples and limpa’s prior", story(page_from_D(D)))
        cb["own"], cb["total"] = [[1.05] * 6], [1.05] * 6
        s = story(page_from_D(D))
        self.assertIn("Its only peptide carries 105% of each sample’s value; other samples and limpa’s prior make up "
                      "the difference, pulling the other way where it is above 100%.", s)
        self.assertNotIn("the rest comes", s)

    def test_the_sankey_says_where_bands_add_up_to_more(self):
        text = visible_text(page_from_D(share_D(other_b1=0.18, total_b1=1.08), "#open=sankey"))
        self.assertIn("Each sample’s bands add up to 100% of its protein value, except in B 1, where they add up to "
                      "108%: limpa’s prior and the gaps pull the other way, by 8%, which a band cannot show.", text)
        text = visible_text(page_from_D(share_D(), "#open=sankey"))
        self.assertIn("Each sample’s bands add up to 100% of its protein value.", text)
        self.assertNotIn("pull the other way", text)

    def test_the_sankey_summarises_more_than_two_samples(self):
        D = share_D(other_b1=0.15, total_b1=1.05)                    # B 1: 105%
        cb = D["contrib"]["P10001"]
        cb["other"][0], cb["total"][0] = 0.18, 1.08                  # A 1: 108%
        cb["other"][3], cb["total"][3] = 0.12, 1.02                  # B 2: 102%
        text = visible_text(page_from_D(D, "#open=sankey"))
        self.assertIn("Each sample’s bands add up to 100% of its protein value, except in 3 samples, where part of "
                      "the value pulls the other way, by 2–8%, which a band cannot show; each sample’s bar in the "
                      "figure above says which part.", text)
        self.assertNotIn("except in A 1", text)


@unittest.skipUnless(BROWSER, "no Chrome, Chromium or Edge found")
class RFixturePage(PageChecks, unittest.TestCase):
    """P00001 (up in B) from make_screen_data.R's own output on the synthetic limpa analysis."""

    @classmethod
    def setUpClass(cls):
        import screen_fixture as fx
        if not fx.r_has(*fx.R_NEEDS):
            raise unittest.SkipTest("needs R with " + "/".join(fx.R_NEEDS))
        tmp = tempfile.mkdtemp(prefix="screen_r_")
        s = fx.make_session(tmp)
        p = subprocess.run(["Rscript", os.path.join(SCRIPTS, "make_screen_data.R"), "--session", s["session"]],
                           capture_output=True, text=True)
        if p.returncode != 0:
            raise AssertionError("make_screen_data.R failed:\n" + p.stderr[-3000:])
        cls.html, cls.dom = build_page("--protein", "P00001",
                                       screen_data=os.path.join(s["output"], "screen_data"))

    def test_report_checksum_was_verified(self):
        self.assertIn("search/report.parquet\u2019s sha256 matches the one recorded when the screen data was made",
                      visible_text(self.dom))


def build_list(*args, screen_data=FIXTURE):
    """make_protein_screen.py --top ... --offline -> (out dir, the list page's DOM)."""
    out = os.path.join(tempfile.mkdtemp(prefix="screen_list_"), "top")
    p = subprocess.run([sys.executable, SCRIPT, "--screen-data", screen_data, "--out-dir", out, "--offline", *args],
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise AssertionError("make_protein_screen.py --top failed (%d):\n%s" % (p.returncode, p.stderr[-3000:]))
    return out, dump_dom(os.path.join(out, "index.html"))


@unittest.skipUnless(BROWSER, "no Chrome, Chromium or Edge found")
class ListPage(unittest.TestCase):
    """The top-N list page, B-A of the mini data: two proteins are called (P10001 up, P10003 down)."""

    @classmethod
    def setUpClass(cls):
        cls.out, cls.dom = build_list("--top", "5")
        cls.text = html.unescape(visible_text(cls.dom))

    def test_renders_ok(self):
        self.assertIn('data-render="ok"', self.dom, render_state(self.dom))

    def test_header_states_comparison_rule_and_k_of_n(self):
        self.assertIn("B vs A", self.text)
        self.assertIn("2 of 5 requested", self.text)
        self.assertIn("Called changed at adjusted p < 0.05; ranked by p-value; no fold-change cut.", self.text)

    def test_says_only_k_when_fewer_are_called(self):
        self.assertIn("Only 2 proteins are called changed in B vs A", self.text)

    def test_no_only_sentence_when_enough_are_called(self):
        _, dom = build_list("--top", "2")
        t = visible_text(dom)
        self.assertNotRegex(t, r"Only \d+ proteins? (is|are) called")
        self.assertIn("2 of 2 requested", t)

    def test_rows_in_p_value_order_with_their_columns(self):
        rows = re.findall(r"<tr>.*?</tr>", re.search(r"<tbody>.*?</tbody>", rendered(self.dom), re.S).group(0), re.S)
        self.assertEqual(len(rows), 2)
        first, second = (visible_text(r).strip() for r in rows)
        self.assertTrue(first.startswith("1 GENEA P10001 AAA_HUMAN +0.50 log2 1.41× 0.03 A 2/2 · B 2/2 no yes"), first)
        self.assertIn("GENEX P10003", second)
        self.assertIn("−0.10 log2 0.93× 0.04", second)
        self.assertNotIn("P10002", self.text)             # sits at adj.P.Val == adjp: not called, never a filler
        self.assertNotIn("P10004", self.text)             # untested

    def test_an_adjusted_p_just_under_the_cutoff_does_not_print_as_the_cutoff(self):
        d = copy.deepcopy(mini_list_D())
        d["proteins"][0]["adjP"] = 0.04973
        out = os.path.join(tempfile.mkdtemp(prefix="screen_list_"), "x.html")
        with open(mps.LIST_TEMPLATE, encoding="utf-8") as fh:
            open(out, "w", encoding="utf-8").write(mps.render(fh.read(), d))
        self.assertIn("0.0497", visible_text(dump_dom(out)))

    def test_default_comparison_is_tagged_on_the_list(self):
        self.assertIn("Comparison B vs A · the first comparison, chosen by default (DEFAULT — not user-confirmed)",
                      self.text)
        _, dom = build_list("--top", "5", "--contrast", "B-A")
        self.assertNotIn("DEFAULT", visible_text(dom))

    def test_offline_list_has_no_uniprot_link(self):
        self.assertNotRegex(rendered(self.dom), r'href="https?://[^"]*uniprot')
        self.assertNotIn("UniProt release", self.text)

    def _with_names(self, edit=None):
        d = copy.deepcopy(mini_list_D())
        for r in d["proteins"]:
            r["fact"], r["fact_why"] = "Name of " + r["pg"], None
            r["fact_src"] = {"source": "UniProt", "release": "2026_03", "fetched_at": "2026-10-28",
                             "url": "https://rest.uniprot.org/uniprotkb/%s.json" % r["pg"], "acc": "P12345"}
        if edit:
            edit(d)
        out = os.path.join(tempfile.mkdtemp(prefix="screen_list_"), "x.html")
        with open(mps.LIST_TEMPLATE, encoding="utf-8") as fh:
            open(out, "w", encoding="utf-8").write(mps.render(fh.read(), d))
        return dump_dom(out)

    def test_names_carry_source_release_fetch_date_and_link(self):
        dom = self._with_names()
        t = html.unescape(visible_text(dom))
        self.assertEqual(t.count("Protein names: UniProt release 2026_03, fetched 2026-10-28"), 1)
        self.assertEqual(re.findall(r'<a href="([^"]+)">Name of', rendered(dom)),
                         ["https://www.uniprot.org/uniprotkb/P12345/entry"] * 2)      # the human page, not the JSON
        self.assertEqual(re.findall(r'<a class="muted"[^>]*href="([^"]+)">record', rendered(dom)),
                         ["https://rest.uniprot.org/uniprotkb/P10001.json", "https://rest.uniprot.org/uniprotkb/P10003.json"])
        self.assertIn("“record” links the exact record the value was fetched from", t)

    def test_differing_releases_are_shown_per_row(self):
        def edit(d):
            d["proteins"][1]["fact_src"]["release"] = "2026_04"
        t = html.unescape(visible_text(self._with_names(edit)))
        self.assertIn("UniProt release 2026_03, fetched 2026-10-28", t)
        self.assertIn("UniProt release 2026_04, fetched 2026-10-28", t)
        self.assertNotIn("Protein names: UniProt", t)

    def test_a_groups_entry_note_is_shown_under_its_name(self):
        def edit(d):
            d["proteins"][0]["fact_note"] = "UniProt entry of P10001, the group's first accession"
        t = html.unescape(visible_text(self._with_names(edit)))
        self.assertIn("Name of P10001 record UniProt entry of P10001, the group's first accession", t)
        self.assertEqual(t.count("the group's first accession"), 1)

    def test_a_bad_accession_is_not_linked(self):
        def edit(d):
            d["proteins"][0]["fact_src"]["acc"] = "x/../evil"
        self.assertNotIn("evil", rendered(self._with_names(edit)))

    def test_a_foreign_url_is_not_linked(self):
        def edit(d):
            d["proteins"][0]["fact_src"]["url"] = "https://evil.example/x"
        dom = self._with_names(edit)
        self.assertNotIn("evil.example", rendered(dom))

    def test_fact_line_says_why_when_nothing_was_looked_up(self):
        self.assertIn("not reached", self.text)

    def test_every_row_links_its_page(self):
        links = re.findall(r'<a href="([^"]+)"', rendered(self.dom))
        self.assertEqual(links, ["GENEA.html", "GENEX.html"])
        for l in links:
            self.assertTrue(os.path.isfile(os.path.join(self.out, l)))

    def test_a_linked_page_renders(self):
        self.assertIn('data-render="ok"', dump_dom(os.path.join(self.out, "GENEA.html")))

    def test_names_no_dataset_of_the_mockup(self):
        for pat in MOCKUP_LITERALS:
            self.assertIsNone(re.search(pat, rendered(self.dom)), pat)

    def test_zero_called_still_builds_and_says_so(self):
        _, dom = build_list("--top", "5", "--contrast", "B-A", "--direction", "up")
        self.assertIn("Only 1 protein is called changed in B vs A", visible_text(dom))
        d = copy.deepcopy(mini_list_D())
        d["proteins"], d["called"] = [], 0
        out = os.path.join(tempfile.mkdtemp(prefix="screen_list_"), "x.html")
        with open(mps.LIST_TEMPLATE, encoding="utf-8") as fh:
            open(out, "w", encoding="utf-8").write(mps.render(fh.read(), d))
        t = visible_text(dump_dom(out))
        self.assertIn("Only 0 proteins are called changed in B vs A", t)
        self.assertNotIn("<table", rendered(dump_dom(out)))

    def test_hostile_values_are_escaped(self):
        d = copy.deepcopy(mini_list_D())
        d["proteins"][0]["names"] = "<img src=x onerror=alert(1)>"
        d["proteins"][0]["fact"] = "</script><b id=pwned>x</b>"
        out = os.path.join(tempfile.mkdtemp(prefix="screen_list_"), "x.html")
        with open(mps.LIST_TEMPLATE, encoding="utf-8") as fh:
            open(out, "w", encoding="utf-8").write(mps.render(fh.read(), d))
        dom = dump_dom(out)
        self.assertIn('data-render="ok"', dom)
        self.assertNotIn("<b id=pwned", rendered(dom))
        self.assertIn("&lt;img src=x onerror=alert(1)&gt;", rendered(dom))


def mini_list_D():
    import json as _json
    out, _ = None, None
    tmp = os.path.join(tempfile.mkdtemp(prefix="screen_list_"), "top")
    subprocess.run([sys.executable, SCRIPT, "--screen-data", FIXTURE, "--out-dir", tmp, "--offline", "--top", "5"],
                   capture_output=True, text=True, check=True)
    with open(os.path.join(tmp, "index.html"), encoding="utf-8") as fh:
        m = re.search(r"const D = (\{.*?\});\nconst CT", fh.read(), re.S)
    return _json.loads(m.group(1))


if __name__ == "__main__":
    unittest.main()
