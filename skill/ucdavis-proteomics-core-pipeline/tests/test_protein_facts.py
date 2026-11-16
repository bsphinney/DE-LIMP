#!/usr/bin/env python3
"""protein_facts.py: protein and gene facts from UniProt, NCBI, HGNC, Ensembl and the Human Protein Atlas.

Nothing here touches the network. The tests serve responses recorded once from the live services
(tests/fixtures/protein_facts/, recorded by record.py there) by patching protein_facts._get.

Brett's rule is the reason this module exists: look the protein up, never rely on memory. So the
central test, test_every_card_value_is_in_a_fetched_record, requires each value a fact shows to
appear verbatim in a recorded response from the source it names (after JSON parsing; numbers are
compared as numbers).
"""
import io
import json
import os
import sys
import unittest
import urllib.error
import urllib.parse
import urllib.request
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import protein_facts as pf  # noqa: E402
import test_protein_screen_page as tpp  # noqa: E402  (module import only: its test classes must not be collected here)

FIX = os.path.join(HERE, "fixtures", "protein_facts")
REFSEQ_ACC = "NP_006308.3"          # BASP1's RefSeq protein (checked live 2026-10-28); P80723 is its UniProt entry
SOURCE_HOST = {"UniProt": "rest.uniprot.org", "NCBI Gene": "eutils.ncbi.nlm.nih.gov",
               "NCBI Protein": "eutils.ncbi.nlm.nih.gov", "HGNC": "rest.genenames.org",
               "Ensembl": "rest.ensembl.org", "Human Protein Atlas": "www.proteinatlas.org"}
SHAPE = {"source", "status", "url", "release", "fetched_at", "licence", "fields", "note"}


def load_index():
    with open(os.path.join(FIX, "index.json")) as fh:
        return json.load(fh)


class Served:
    """Stands in for protein_facts._get: serves the recorded responses, remembers what it served."""

    def __init__(self, dead=()):
        self.index, self.dead = load_index(), set(dead)
        self.urls, self.bodies, self.unrecorded, self.by_url = [], {}, [], {}

    def __call__(self, url, headers=None, timeout=20):
        self.urls.append(url)
        host = urllib.parse.urlparse(url).hostname
        if host in self.dead:
            raise urllib.error.URLError("%s is unreachable (test)" % host)
        rec = self.index["responses"].get(url)
        if rec is None:
            self.unrecorded.append(url)
            raise urllib.error.URLError("no recorded response for " + url)
        with open(os.path.join(FIX, rec["file"]), "rb") as fh:
            body = fh.read()
        self.bodies.setdefault(host, []).append(body)
        self.by_url[url] = body
        hdrs = dict(rec["headers"])
        hdrs.setdefault("x-final-url", url)
        return rec["status"], hdrs, body


def by_source(facts, source):
    got = [f for f in facts if f["source"] == source]
    assert len(got) == 1, (source, [f["source"] for f in facts])
    return got[0]


def _at(doc, path):
    for k in path:
        doc = doc[k]
    return doc


# Where each shown (label, value) must come from: the record of ONE call (the Fact's own url, or the
# GenBank efetch for the NCBI Protein gene), at one fixed key path. A swapped label or a value taken
# from another record of the same host fails. HPA labels are the record's own keys.
UNIPROT_LIVE = {"Accession": ["primaryAccession"], "Entry name": ["uniProtkbId"], "Status": ["entryType"],
                "Protein name": ["proteinDescription", "recommendedName", "fullName", "value"],
                "Gene": ["genes", 0, "geneName", "value"], "Organism": ["organism", "scientificName"],
                "Protein existence": ["proteinExistence"], "Entry version": ["entryAudit", "entryVersion"],
                "Annotation updated": ["entryAudit", "lastAnnotationUpdateDate"]}
UNIPROT_INACTIVE = {"Accession": ["primaryAccession"], "Status": ["entryType"],
                    "Reason": ["inactiveReason", "inactiveReasonType"],
                    "Merged or demerged into": ["inactiveReason", "mergeDemergeTo"],
                    "UniParc": ["extraAttributes", "uniParcId"]}
PATHS = {
    "NCBI Gene": {"Symbol": ["name"], "Description": ["description"], "Also known as": ["otheraliases"],
                  "Organism": ["organism", "scientificname"], "Chromosome": ["chromosome"],
                  "Map location": ["maplocation"], "Summary": ["summary"]},
    "NCBI Protein": {"Title": ["title"], "Organism": ["organism"], "Accession version": ["accessionversion"],
                     "Length": ["slen"], "Updated": ["updatedate"], "Database": ["sourcedb"]},
    "HGNC": {"Symbol": ["symbol"], "Name": ["name"], "HGNC ID": ["hgnc_id"], "Status": ["status"],
             "Locus group": ["locus_group"], "Locus type": ["locus_type"], "Location": ["location"],
             "Alias symbols": ["alias_symbol"], "Alias names": ["alias_name"], "Previous symbols": ["prev_symbol"],
             "Previous names": ["prev_name"], "Gene groups": ["gene_group"]},
    "Ensembl": {"Gene ID": ["id"], "Name": ["display_name"], "Description": ["description"],
                "Biotype": ["biotype"], "Species": ["species"], "Assembly": ["assembly_name"],
                "Chromosome": ["seq_region_name"], "Start": ["start"], "End": ["end"], "Strand": ["strand"],
                "Source": ["source"]},
}
GENBANK = {"Gene": r'/gene="([^"]+)"', "GeneID": r'/db_xref="GeneID:(\d+)"'}


def check_fact(fact, served):
    """Every (label, value) of `fact` against its one place in the recorded response -> list of mismatches."""
    import re
    src, bad = fact["source"], []
    if not fact["fields"]:
        return bad
    body = served.by_url.get(fact["url"])
    if body is None:
        return ["no recorded response was served for %s" % fact["url"]]
    for label, value in fact["fields"]:
        if src == "NCBI Protein" and label in GENBANK:
            gp = [b for u, b in served.by_url.items() if "efetch" in u]
            m = re.search(GENBANK[label], gp[0].decode("utf-8")) if gp else None
            want = m.group(1) if m else "(no such qualifier)"
        else:
            doc = json.loads(body.decode("utf-8"))
            try:
                if src == "UniProt":
                    path = (UNIPROT_INACTIVE if doc.get("entryType") == "Inactive" else UNIPROT_LIVE)[label]
                elif src == "Human Protein Atlas":
                    path = [label]
                else:
                    path = PATHS[src][label]
                if src == "NCBI Gene":
                    doc = doc["result"][re.search(r"[?&]id=(\d+)", fact["url"]).group(1)]
                elif src == "NCBI Protein":
                    doc = doc["result"][doc["result"]["uids"][0]]
                elif src == "HGNC":
                    doc = doc["response"]["docs"][0]
                want = _at(doc, path)
            except (KeyError, IndexError, TypeError) as e:
                bad.append((src, label, value, "no such place in the record: %r" % (e,)))
                continue
        if want != value or type(want) is not type(value):
            bad.append((src, label, value, "the record says %r" % (want,)))
    return bad


class FactsCase(unittest.TestCase):
    def setUp(self):
        self.served = Served()
        p = mock.patch.object(pf, "_get", self.served)
        p.start()
        self.addCleanup(p.stop)

    def facts_with_dead(self, acc, host):
        served = Served(dead=[host])
        with mock.patch.object(pf, "_get", served):
            return pf.facts_for(acc)

    def tearDown(self):
        self.assertEqual(self.served.unrecorded, [], "a URL was requested that has no recorded response")


class TestAllowList(unittest.TestCase):
    def test_host_allow_list(self):
        self.assertEqual(pf.ALLOWED_HOSTS, {"rest.uniprot.org", "eutils.ncbi.nlm.nih.gov", "rest.genenames.org",
                                            "rest.ensembl.org", "www.proteinatlas.org"})
        for url in ("https://example.org/x", "http://rest.uniprot.org/uniprotkb/P80723.json",
                    "https://rest.uniprot.org.evil.example/x", "https://evil.example/rest.uniprot.org",
                    "file:///etc/passwd", "ftp://rest.uniprot.org/x"):
            with self.assertRaises(ValueError, msg=url):
                pf._get(url)


class FakeResp(io.BytesIO):
    status = 200

    def __init__(self, body, url, headers=None):
        super().__init__(body)
        self._url, self.headers = url, headers or {}

    def geturl(self):
        return self._url

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


class FakeOpener:
    def __init__(self, result):
        self.result, self.requests = result, []

    def open(self, req, timeout=None):
        self.requests.append((req, timeout))
        if isinstance(self.result, Exception):
            raise self.result
        return self.result


class TestGet(unittest.TestCase):
    def get(self, url, result, **kw):
        op = FakeOpener(result)
        with mock.patch.object(pf, "_OPENER", op):
            return pf._get(url, **kw), op

    def test_sends_user_agent_and_times_out_after_20_s(self):
        url = "https://rest.genenames.org/info"
        (status, hdrs, body), op = self.get(url, FakeResp(b"{}", url, {"Content-Type": "application/json"}),
                                            headers={"Accept": "application/json"})
        req, timeout = op.requests[0]
        self.assertEqual(req.get_header("User-agent"), "ucdavis-proteomics-skill (+https://github.com/bsphinney/UCD-Proteomics-AI-Skill)")
        self.assertEqual(req.get_header("Accept"), "application/json")
        self.assertEqual(timeout, 20)
        self.assertEqual((status, body, hdrs["content-type"]), (200, b"{}", "application/json"))

    def test_an_http_error_is_a_status_not_an_exception(self):
        url = "https://rest.uniprot.org/uniprotkb/P80723.json"
        err = urllib.error.HTTPError(url, 503, "busy", {"Retry-After": "5"}, io.BytesIO(b"busy"))
        (status, hdrs, body), _ = self.get(url, err)
        self.assertEqual((status, body, hdrs["retry-after"]), (503, b"busy", "5"))

    def test_final_host_is_checked(self):
        with self.assertRaises(urllib.error.URLError):
            self.get("https://rest.uniprot.org/uniprotkb/P80723.json",
                     FakeResp(b"{}", "https://evil.example/P80723.json"))

    def test_redirect_handler_follows_only_allowed_hosts(self):
        h = pf._AllowedRedirect()
        req = urllib.request.Request("https://rest.uniprot.org/uniprotkb/O43596.json")
        with self.assertRaises(urllib.error.URLError):
            h.redirect_request(req, None, 303, "See Other", {}, "https://evil.example/x")
        with self.assertRaises(urllib.error.URLError):
            h.redirect_request(req, None, 303, "See Other", {}, "http://rest.uniprot.org/x")
        self.assertIsNotNone(h.redirect_request(req, None, 303, "See Other", {},
                                                "https://rest.uniprot.org/uniprotkb/P80723?from=O43596"))

    def test_ncbi_calls_are_at_least_0_34_s_apart(self):
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=1"
        pf._LAST_CALL.clear()
        slept = []
        with mock.patch.object(pf.time, "monotonic", return_value=100.0), \
                mock.patch.object(pf.time, "sleep", side_effect=slept.append):
            self.get(url, FakeResp(b"{}", url))
            self.get(url, FakeResp(b"{}", url))
        self.assertEqual(len(slept), 1)
        self.assertGreaterEqual(slept[0], 0.34 - 1e-9)

    def test_other_hosts_are_not_slowed(self):
        url = "https://rest.genenames.org/info"
        pf._LAST_CALL.clear()
        slept = []
        with mock.patch.object(pf.time, "monotonic", return_value=100.0), \
                mock.patch.object(pf.time, "sleep", side_effect=slept.append):
            self.get(url, FakeResp(b"{}", url))
            self.get(url, FakeResp(b"{}", url))
        self.assertEqual(slept, [])


class TestFetchUniprot(FactsCase):
    """fetch_uniprot moved here from make_protein_screen.py with its hardening intact."""

    def test_fields_extracted_verbatim(self):
        got = pf.fetch_uniprot("P80723")
        self.assertEqual((got["acc"], got["id"], got["gene"], got["organism"]),
                         ("P80723", "BASP1_HUMAN", "BASP1", "Homo sapiens"))
        self.assertEqual(got["reviewed"], "UniProtKB reviewed (Swiss-Prot)")
        self.assertEqual(got["length"], len(got["sequence"]))
        self.assertRegex(got["release"], r"^\d{4}_\d{2}$")
        self.assertRegex(got["fetched_at"], r"^\d{4}-\d{2}-\d{2}$")
        self.assertEqual(self.served.urls, ["https://rest.uniprot.org/uniprotkb/P80723.json"])

    def test_inactive_entry_is_no_page_but_says_why(self):
        errs = []
        self.assertIsNone(pf.fetch_uniprot("P04277", errors=errs))
        self.assertIn("DEMERGED", errs[0])

    def test_http_error(self):
        errs = []
        with mock.patch.object(pf, "_get", lambda *a, **k: (503, {}, b"busy")):
            self.assertIsNone(pf.fetch_uniprot("P80723", errors=errs))
        self.assertIn("HTTP 503", errs[0])

    def test_incomplete_read_is_caught(self):
        import http.client
        errs = []

        def boom(*a, **k):
            raise http.client.IncompleteRead(b"x")
        with mock.patch.object(pf, "_get", boom):
            self.assertIsNone(pf.fetch_uniprot("P80723", errors=errs))
        self.assertTrue(errs)

    def test_not_json_is_caught(self):
        errs = []
        with mock.patch.object(pf, "_get", lambda *a, **k: (200, {}, b"<html>")):
            self.assertIsNone(pf.fetch_uniprot("P80723", errors=errs))
        self.assertTrue(errs)

    def test_bad_accession_not_requested(self):
        for acc in ("../etc/passwd", "P10001/../x", "ABC", "P8072", "Cont_P99999", REFSEQ_ACC):
            errs = []
            self.assertIsNone(pf.fetch_uniprot(acc, errors=errs), acc)
            self.assertIn("not a UniProt accession", errs[0])
        self.assertEqual(self.served.urls, [])
        for acc in ("P80723", "A0A023GPI8", "Q9Y223-2"):
            self.assertTrue(pf.UNIPROT_ACCESSION.fullmatch(acc.split("-")[0]), acc)

    def test_isoform_id_asks_for_the_entry(self):
        pf.fetch_uniprot("P80723-2")
        self.assertEqual(self.served.urls, ["https://rest.uniprot.org/uniprotkb/P80723.json"])

    def test_isoform_note_says_the_canonical_entry_is_shown(self):
        self.assertEqual((pf.entry_accession("P80723-2"), pf.entry_accession("P80723")), ("P80723", "P80723"))
        self.assertEqual(pf.entry_accession("ABC123-ALT"), "ABC123-ALT")       # protein_ids.py's rule: -<digits> only
        up = by_source(pf.facts_for("P80723-2"), "UniProt")
        self.assertEqual(up["status"], "ok")
        self.assertIn(["Accession", "P80723"], up["fields"])
        self.assertIn("P80723-2 is an isoform", up["note"])
        self.assertIn("the canonical entry, P80723, is shown", up["note"])
        self.assertIsNone(by_source(pf.facts_for("P80723"), "UniProt")["note"])      # a canonical id: no note


class TestFactsFor(FactsCase):
    def test_every_card_value_is_in_a_fetched_record(self):      # Brett's rule: nothing from memory
        for acc in ("P80723", "P02768", "P15108", "O43596", "P04277", REFSEQ_ACC):
            n = 0
            for fact in pf.facts_for(acc):
                n += len(fact["fields"])
                self.assertEqual(check_fact(fact, self.served), [], (acc, fact["source"]))
            self.assertGreater(n, 3, acc)

    def test_the_check_itself_can_fail(self):
        facts = pf.facts_for("P80723")
        ens, hpa = by_source(facts, "Ensembl"), by_source(facts, "Human Protein Atlas")
        self.assertEqual(check_fact(ens, self.served), [])
        swapped = json.loads(json.dumps(ens))
        for f in swapped["fields"]:
            if f[0] in ("Start", "End"):
                f[1] = ens["fields"][[x[0] for x in ens["fields"]].index("End" if f[0] == "Start" else "Start")][1]
        self.assertEqual({b[1] for b in check_fact(swapped, self.served)}, {"Start", "End"})   # two labels swapped
        typed = json.loads(json.dumps(by_source(facts, "HGNC")))
        typed["fields"][1][1] = "a value typed from memory"
        self.assertEqual(len(check_fact(typed, self.served)), 1)
        alb = json.loads(json.dumps(hpa))
        alb["fields"][0][1] = 40000000000000                    # ALB's number shown under BASP1
        self.assertEqual(len(check_fact(alb, self.served)), 1)
        mixed = json.loads(json.dumps(by_source(facts, "NCBI Gene")))
        mixed["url"] = by_source(facts, "Human Protein Atlas")["url"]      # a value looked up in another record
        self.assertTrue(check_fact(mixed, self.served))

    def test_fact_shape_and_order_for_a_human_protein(self):
        facts = pf.facts_for("P80723")
        self.assertEqual([f["source"] for f in facts],
                         ["UniProt", "NCBI Gene", "HGNC", "Ensembl", "Human Protein Atlas"])
        for f in facts:
            self.assertEqual(set(f), SHAPE, f["source"])
            self.assertEqual(f["status"], "ok", (f["source"], f["note"]))
            self.assertRegex(f["fetched_at"], r"^\d{4}-\d{2}-\d{2}$")
            self.assertTrue(f["licence"], f["source"])
            self.assertTrue(f["url"].startswith("https://"), f["source"])
            self.assertTrue(f["fields"], f["source"])
            for pair in f["fields"]:
                self.assertEqual(len(pair), 2)
                self.assertIsInstance(pair[0], str)
        json.dumps(facts, allow_nan=False)

    def test_releases_are_the_ones_the_sources_state(self):
        facts = pf.facts_for("P80723")
        up = by_source(facts, "UniProt")
        self.assertEqual(up["release"], load_index()["responses"][up["url"]]["headers"]["x-uniprot-release"])
        self.assertRegex(by_source(facts, "Ensembl")["release"], r"^\d+$")
        self.assertRegex(by_source(facts, "HGNC")["release"], r"^\d{4}-\d{2}-\d{2}T")
        self.assertIsNone(by_source(facts, "Human Protein Atlas")["release"])      # HPA states none; never typed in
        self.assertIsNone(by_source(facts, "NCBI Gene")["release"])

    def test_ids_come_from_uniprot_crossrefs_not_a_gene_name_search(self):
        pf.facts_for("P02768")
        joined = "\n".join(self.served.urls)
        self.assertIn("https://www.proteinatlas.org/ENSG00000163631.json", self.served.urls)
        self.assertNotIn("search_download", joined)
        self.assertNotIn("[Gene Name]", joined)
        self.assertNotIn("esearch", joined)
        self.assertIn("id=213", joined)

    def test_hpa_blood_concentrations_verbatim(self):
        hpa = by_source(pf.facts_for("P02768"), "Human Protein Atlas")
        self.assertIn(["Blood concentration - Conc. blood IM [pg/L]", 40000000000000], hpa["fields"])
        self.assertIn(["Blood concentration - Conc. blood MS [pg/L]", None], hpa["fields"])        # null = no HPA value
        self.assertIn(["Secretome location", "Secreted to blood"], hpa["fields"])
        self.assertIn(["Evidence", "Evidence at protein level"], hpa["fields"])
        self.assertIn(["RNA blood cell specificity", "Not detected in immune cells"], hpa["fields"])
        self.assertEqual(len(hpa["fields"]), 5)
        self.assertIsNone(hpa["release"])

    def test_hpa_keeps_the_sources_own_types(self):
        hpa = by_source(pf.facts_for("P80723"), "Human Protein Atlas")
        self.assertIn(["Blood concentration - Conc. blood MS [pg/L]", "49000000"], hpa["fields"])   # a JSON string there
        self.assertIn(["Blood concentration - Conc. blood IM [pg/L]", None], hpa["fields"])
        self.assertIn(["RNA blood cell specificity", "Immune cell enriched"], hpa["fields"])

    def test_ncbi_gene_fields(self):
        g = by_source(pf.facts_for("P80723"), "NCBI Gene")
        d = dict((k, v) for k, v in g["fields"])
        self.assertEqual(d["Description"], "brain abundant membrane attached signal protein 1")
        self.assertEqual(d["Also known as"], "CAP-23, CAP23, NAP-22, NAP22")
        self.assertEqual(d["Map location"], "5p15.1")
        self.assertIn("[provided by RefSeq,", d["Summary"])

    def test_yeast_gets_no_human_only_sources(self):
        facts = pf.facts_for("P15108")
        st = {f["source"]: f["status"] for f in facts}
        self.assertEqual(st["HGNC"], "not_applicable")
        self.assertEqual(st["Human Protein Atlas"], "not_applicable")
        self.assertEqual(st["Ensembl"], "ok")
        self.assertEqual(st["NCBI Gene"], "ok")
        joined = "\n".join(self.served.urls)
        self.assertNotIn("genenames.org", joined)
        self.assertNotIn("proteinatlas.org", joined)
        self.assertIn("YMR186W", by_source(facts, "Ensembl")["fields"][0][1])
        self.assertIn("human", by_source(facts, "HGNC")["note"])

    def test_merged_accession_follows_the_target(self):
        facts = pf.facts_for("O43596")
        up = by_source(facts, "UniProt")
        self.assertEqual(up["status"], "ok")
        self.assertIn("P80723", up["note"])
        self.assertIn("O43596", up["note"])
        self.assertIn(["Accession", "P80723"], up["fields"])
        self.assertEqual(by_source(facts, "NCBI Gene")["status"], "ok")

    def test_demerged_accession_lists_targets_and_picks_none(self):
        facts = pf.facts_for("P04277")
        self.assertEqual([f["source"] for f in facts], ["UniProt"])         # no gene sources: there is no one gene to ask about
        up = facts[0]
        self.assertIn("P02768", up["note"])
        self.assertIn("P02769", up["note"])
        self.assertIn("DEMERGED", up["note"])
        self.assertIn(["Merged or demerged into", ["P02769", "P02768"]], up["fields"])
        self.assertTrue(any(k == "UniParc" and v.startswith("UPI") for k, v in up["fields"]))
        self.assertNotIn("NCBI", "\n".join(self.served.urls))

    def test_not_in_uniprot_uses_ncbi_protein(self):
        facts = pf.facts_for(REFSEQ_ACC)
        self.assertEqual([f["source"] for f in facts], ["UniProt", "NCBI Protein", "NCBI Gene"])
        self.assertEqual(by_source(facts, "UniProt")["status"], "not_applicable")
        self.assertEqual(by_source(facts, "NCBI Protein")["status"], "ok")
        self.assertEqual(by_source(facts, "NCBI Gene")["status"], "ok")
        prot = dict(by_source(facts, "NCBI Protein")["fields"])
        self.assertEqual(prot["Title"], "brain acid soluble protein 1 [Homo sapiens]")
        self.assertEqual(prot["Gene"], "BASP1")
        self.assertEqual(prot["GeneID"], "10409")
        self.assertEqual(dict(by_source(facts, "NCBI Gene")["fields"])["Symbol"], "BASP1")

    def test_ncbi_calls_name_the_tool_and_send_no_address_or_key(self):
        env = {k: v for k, v in os.environ.items() if k not in ("NCBI_EMAIL", "NCBI_API_KEY")}
        with mock.patch.dict(os.environ, env, clear=True):
            pf.facts_for("P80723")
        ncbi = [u for u in self.served.urls if "eutils" in u]
        self.assertTrue(ncbi)
        for u in ncbi:
            self.assertIn("tool=ucdavis_proteomics_skill", u)
            self.assertNotIn("email=", u)
            self.assertNotIn("api_key=", u)

    def test_ncbi_email_and_key_only_from_the_environment(self):
        with mock.patch.dict(os.environ, {"NCBI_EMAIL": "x@example.org", "NCBI_API_KEY": "k123"}):
            url = pf._ncbi_url("esummary.fcgi", db="gene", id="10409", retmode="json")
        q = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)
        self.assertEqual((q["email"], q["api_key"], q["tool"]), (["x@example.org"], ["k123"], ["ucdavis_proteomics_skill"]))
        self.served.unrecorded.clear()

    def test_unreachable_source_says_so(self):
        for host, source in (("rest.genenames.org", "HGNC"), ("rest.ensembl.org", "Ensembl"),
                             ("www.proteinatlas.org", "Human Protein Atlas"),
                             ("eutils.ncbi.nlm.nih.gov", "NCBI Gene")):
            facts = self.facts_with_dead("P80723", host)
            f = by_source(facts, source)
            self.assertEqual(f["status"], "not_reached", source)
            self.assertIn(host, f["note"])
            self.assertEqual(f["fields"], [])
            self.assertEqual(by_source(facts, "UniProt")["status"], "ok")

    def test_uniprot_unreachable_stops_there_and_says_so(self):
        facts = self.facts_with_dead("P80723", "rest.uniprot.org")
        self.assertEqual([(f["source"], f["status"]) for f in facts], [("UniProt", "not_reached")])

    def test_a_server_error_is_not_reached_and_a_404_is_not_found(self):
        real = self.served

        def get(url, headers=None, timeout=20):
            if "genenames" in url:
                return 503, {}, b"busy"
            if "ensembl" in url and "lookup" in url:
                return 404, {}, b'{"error":"ID not found"}'
            return real(url, headers, timeout)
        with mock.patch.object(pf, "_get", get):
            facts = pf.facts_for("P80723")
        self.assertEqual(by_source(facts, "HGNC")["status"], "not_reached")
        self.assertIn("HTTP 503", by_source(facts, "HGNC")["note"])
        self.assertEqual(by_source(facts, "Ensembl")["status"], "not_found")

    def test_an_answer_that_is_not_json_is_not_reached(self):
        real = self.served

        def get(url, headers=None, timeout=20):
            if "proteinatlas" in url:
                return 200, {}, b"<html>maintenance</html>"
            return real(url, headers, timeout)
        with mock.patch.object(pf, "_get", get):
            hpa = by_source(pf.facts_for("P80723"), "Human Protein Atlas")
        self.assertEqual(hpa["status"], "not_reached")

    def test_offline_facts(self):
        facts = pf.offline_facts()
        self.assertEqual([f["source"] for f in facts],
                         ["UniProt", "NCBI Gene", "HGNC", "Ensembl", "Human Protein Atlas"])
        for f in facts:
            self.assertEqual((f["status"], f["note"], f["fields"]), ("not_reached", "--offline", []))
            self.assertEqual(set(f), SHAPE)
        self.assertEqual(self.served.urls, [])


class TestFixRound1(FactsCase):
    def page_text(self, facts):
        """The page for `facts`, built the way the CLI builds it (page_data + render): its text, plus the data."""
        D = tpp.mini_D(uniprot=dict(tpp.UNIPROT))
        D["info"]["P10001"]["facts"] = facts
        with open(mps_template()) as fh:
            return mps_render(fh.read(), D), D

    # 1. the user's NCBI address and key must never leave the machine's request line
    def test_email_and_key_are_sent_but_never_shown_or_dumped(self):
        served, seen = self.served, []

        def get(url, headers=None, timeout=20):
            seen.append(url)
            for k in ("email", "api_key"):                       # the recordings were made without them
                url = "&".join(kv for kv in url.split("&") if not kv.startswith(k + "="))
                url = url.replace("?" + k + "=", "?")
            return served(url.rstrip("?"), headers, timeout)
        marks = {"NCBI_EMAIL": "secret.person@example.invalid", "NCBI_API_KEY": "KEYMARKER0123456789"}
        with mock.patch.dict(os.environ, marks), mock.patch.object(pf, "_get", get):
            facts = pf.facts_for("P80723") + pf.facts_for(REFSEQ_ACC)
        ncbi = [u for u in seen if "eutils" in u]
        self.assertTrue(ncbi and all("api_key=KEYMARKER0123456789" in u for u in ncbi))      # the request carried them
        self.assertTrue(by_source(facts[:5], "NCBI Gene")["url"].startswith("https://eutils"))
        page, D = self.page_text(facts)
        blob = json.dumps(facts) + json.dumps(D) + page
        for m in ("secret.person", "example.invalid", "KEYMARKER", "api_key", "email="):
            self.assertNotIn(m, blob)
        self.assertIn("tool=ucdavis_proteomics_skill", by_source(facts[:5], "NCBI Gene")["url"])

    # 3. valid JSON of the wrong shape is a reported failure, never a crash
    def facts_with_body(self, acc, host_part, body, path_part=""):
        real = self.served

        def get(url, headers=None, timeout=20):
            if host_part in url and path_part in url:
                return 200, {}, body
            return real(url, headers, timeout)
        with mock.patch.object(pf, "_get", get):
            return pf.facts_for(acc)

    def test_wrong_shape_bodies_are_reported_and_the_page_still_builds(self):
        cases = [
            ("P80723", "NCBI Gene", "eutils", b"[]"), ("P80723", "NCBI Gene", "eutils", b'{"result":[]}'),
            ("P80723", "NCBI Gene", "eutils", b'{"result":{"10409":["x"]}}'),
            (REFSEQ_ACC, "NCBI Protein", "esummary.fcgi?db=protein", b"[1]"),
            (REFSEQ_ACC, "NCBI Protein", "esummary.fcgi?db=protein", b'{"result":{"uids":[7],"7":"a string"}}'),
            ("P80723", "HGNC", "fetch/hgnc_id", b"[]"),
            ("P80723", "HGNC", "fetch/hgnc_id", b'{"response":{"docs":["a string"]}}'),
            ("P80723", "HGNC", "fetch/hgnc_id", b'{"response":["x"]}'),
            ("P80723", "Ensembl", "lookup/id", b"[]"), ("P80723", "Ensembl", "lookup/id", b'"a string"'),
            ("P80723", "Human Protein Atlas", "proteinatlas", b"[]"),
            ("P80723", "Human Protein Atlas", "proteinatlas", b'"a string"'),
            ("P80723", "UniProt", "rest.uniprot.org", b'{"primaryAccession":"P80723","sequence":"abc"}'),
            ("P80723", "UniProt", "rest.uniprot.org", b'{"entryType":"Inactive","inactiveReason":"gone"}'),
            ("P80723", "UniProt", "rest.uniprot.org", b"[]"),
        ]
        for acc, source, host_part, body in cases:
            facts = self.facts_with_body(acc, host_part, body)
            f = by_source(facts, source)
            self.assertIn(f["status"], ("not_reached", "not_found"), (source, body))
            self.assertTrue(f["note"], (source, body))
            self.assertEqual(f["fields"], [], (source, body))
            json.dumps(facts)
            page, _ = self.page_text(facts)
            self.assertNotIn("__SCREEN_DATA__", page)

    def test_one_wrong_shape_source_leaves_the_others_alone(self):
        facts = self.facts_with_body("P80723", "genenames", b"[1]")
        self.assertEqual([(f["source"], f["status"]) for f in facts],
                         [("UniProt", "ok"), ("NCBI Gene", "ok"), ("HGNC", "not_reached"), ("Ensembl", "ok"),
                          ("Human Protein Atlas", "ok")])
        self.assertIn("unexpected response shape", by_source(facts, "HGNC")["note"])

    # 4. a release lookup that fails is recorded, not passed off as "the source states no version"
    def test_failed_release_lookup_is_recorded(self):
        real = self.served
        for source, part, mode in (("HGNC", "genenames.org/info", "503"), ("Ensembl", "info/data", "503"),
                                   ("HGNC", "genenames.org/info", "dead"), ("Ensembl", "info/data", "dead"),
                                   ("HGNC", "genenames.org/info", "shape"), ("Ensembl", "info/data", "shape"),
                                   ("Ensembl", "info/data", "empty")):
            def get(url, headers=None, timeout=20, part=part, mode=mode):
                if part in url:
                    if mode == "503":
                        return 503, {}, b"busy"
                    if mode == "dead":
                        raise urllib.error.URLError("unreachable (test)")
                    return 200, {}, (b"[]" if mode == "shape" else b'{"releases":[]}')
                return real(url, headers, timeout)
            with mock.patch.object(pf, "_get", get):
                f = by_source(pf.facts_for("P80723"), source)
            self.assertEqual(f["status"], "ok", (source, mode))
            self.assertTrue(f["fields"], (source, mode))              # the record itself is still shown
            self.assertTrue(f["release"].startswith("not retrieved: "), (source, mode, f["release"]))
            if mode == "503":
                self.assertIn("HTTP 503", f["release"])
            if mode == "dead":
                self.assertIn("could not be reached", f["release"])


def mps_template():
    import make_protein_screen
    return make_protein_screen.TEMPLATE


def mps_render(text, D):
    import make_protein_screen
    return make_protein_screen.render(text, D)


@unittest.skipUnless(tpp.BROWSER, "no Chrome, Chromium or Edge found")
class TestPageCard(unittest.TestCase):
    def test_page_card_shows_sources(self):
        served = Served()
        with mock.patch.object(pf, "_get", served):
            facts = pf.facts_for("P02768") + pf.facts_for("P15108")
        D = tpp.mini_D(uniprot=dict(tpp.UNIPROT))
        D["info"]["P10001"]["facts"] = facts
        dom = tpp.page_from_D(D)
        self.assertIn('data-render="ok"', dom, tpp.render_state(dom))
        text = tpp.visible_text(dom)
        for want in ("Human Protein Atlas", "pg/L", "NCBI Gene", "HGNC", "Ensembl", "40,000,000,000,000",
                     "Secreted to blood", "version not stated by the source", "CC BY 4.0",
                     "Not applicable", "human genes only"):
            self.assertIn(want, text, want)
        self.assertNotIn("nextprot", text.lower())                       # retired: no statement on every page
        self.assertIn('href="https://www.proteinatlas.org/ENSG00000163631.json"', tpp.rendered(dom))

    def test_hostile_values_render_safely(self):
        facts = pf.offline_facts()
        facts[1] = dict(facts[1], status="ok", note=None, release="<b>9</b>",
                        fields=[["Summary", "<img src=x onerror=alert(1)>"]], licence="<i>x</i>")
        D = tpp.mini_D(uniprot=dict(tpp.UNIPROT))
        D["info"]["P10001"]["facts"] = facts
        dom = tpp.page_from_D(D)
        self.assertIn('data-render="ok"', dom, tpp.render_state(dom))
        r = tpp.rendered(dom)
        self.assertNotIn("<img src=x", r)
        self.assertNotIn("<b>9</b>", r)
        text = tpp.visible_text(dom)
        self.assertIn("Not reached: --offline", text)
        self.assertIn("&lt;img src=x onerror=alert(1)&gt;", text)      # shown as text, not run

    def test_offline_page_says_the_sources_were_not_reached(self):
        out = tpp.make_page("--protein", "GENEA")
        text = tpp.visible_text(tpp.dump_dom(out))
        self.assertIn("Not reached: --offline", text)

    def test_a_failed_release_lookup_is_not_shown_as_no_version(self):
        facts = pf.offline_facts()
        facts[3] = dict(facts[3], status="ok", note=None, release="not retrieved: rest.ensembl.org answered HTTP 503",
                        fields=[["Gene ID", "ENSG00000176788"]], licence="x")
        D = tpp.mini_D(uniprot=dict(tpp.UNIPROT))
        D["info"]["P10001"]["facts"] = facts
        text = tpp.visible_text(tpp.page_from_D(D))
        self.assertIn("release not retrieved: rest.ensembl.org answered HTTP 503", text)
        self.assertNotIn("version not stated by the source", text)


if __name__ == "__main__":
    unittest.main()
