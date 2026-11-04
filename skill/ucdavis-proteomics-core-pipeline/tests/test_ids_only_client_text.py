#!/usr/bin/env python3
"""
An identifications-only search is never described to a client as having quantified anything
(2.14.2 C2 review: biologist B1 = errors_ux B1).

A Sage --ids-only search, or a FragPipe workflow that runs no quantification, writes
identifications.tsv and no report.parquet; the DE gate refuses it (sage_lfq_check.json status
not_applicable). The delivery README said "searched to identify and quantify proteins" and
listed "report.parquet and the matrices" for its folder; the interim README and email said "the
database search that identified and quantified the proteins". Now ONE predicate,
sage_lfq_check.search_quantified (identifications_only() of the record the gate reads, else its
copies in search_provenance.json), decides every one of those sentences:
  - the analysis / no-de README's claim and its search-folder row;
  - the interim README's sentence, folder rows and "still to come" line;
  - delivery.json `search_quantified`, which the interim email draft reads.
Hermetic: temp dirs, the CoreOmics fake of test_core_submission, no network.
"""
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import core_submission as cs                   # noqa: E402
import sage_lfq_check as slc                   # noqa: E402
import test_core_submission as tcs             # noqa: E402
# the fixtures only: a module-level TestCase with tests would run its tests again here
import test_deliver_interim as tdi             # noqa: E402

IDS_ONLY_REC = slc.not_applicable()           # Sage --ids-only's record, the gate's file
FRAGPIPE_IDS_REC = slc.not_applicable(engine="fragpipe", why="the workflow runs no quantification",
                                      ran="FragPipe ran no quantification")


def folder(root, name, files):
    d = os.path.join(root, name)
    os.makedirs(d, exist_ok=True)
    for rel, obj in files.items():
        with open(os.path.join(d, rel), "w") as fh:
            fh.write(obj if isinstance(obj, str) else json.dumps(obj))
    return d


class ThePredicate(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.t = self._t.name

    def tearDown(self):
        self._t.cleanup()

    def test_the_gates_record_says_identifications_only(self):
        d = folder(self.t, "search_sage", {slc.RECORD: IDS_ONLY_REC,
                                            "search_provenance.json": {"engine": "sage"}})
        self.assertFalse(slc.search_quantified(d))
        self.assertTrue(slc.identifications_only(slc.load(d)))

    def test_a_copy_with_only_search_provenance_is_read_too(self):
        """A delivered or laptop copy has search_provenance.json and not the gate's file."""
        for prov in ({"engine": "sage", "quantified": False},
                     {"engine": "sage", "sage_lfq_check": IDS_ONLY_REC},
                     {"engine": "fragpipe", "quantification": FRAGPIPE_IDS_REC}):
            with self.subTest(prov=sorted(prov)):
                d = folder(self.t, "s", {"search_provenance.json": prov})
                self.assertFalse(slc.search_quantified(d))

    def test_a_quantifying_search_is_quantified(self):
        for files in ({"search_provenance.json": {"engine": "diann", "version": "2.7.0"}},
                      {"search_provenance.json": {"engine": "sage", "quantified": True},
                       slc.RECORD: {"status": "ok"}},
                      {"search_provenance.json": {"engine": "fragpipe", "quantification": {
                          "status": "quantified", "why": "LFQ-MBR"}}},
                      {}):
            with self.subTest(files=sorted(files)):
                d = folder(tempfile.mkdtemp(dir=self.t), "s", files)
                self.assertTrue(slc.search_quantified(d))


SUMMARY = {"internal_id": "PROT_TEST", "id": "x"}


class TheReadme(unittest.TestCase):
    IDS = {"search_sage/identifications.tsv", "search_sage/sage_identifications.json",
           "search_sage/search_provenance.json", "Analysis_Report.html", "methods.md"}

    def readme(self, mode, delivered, quantified):
        return cs.build_readme(SUMMARY, set(delivered), [], mode, search_dir="search_sage",
                               quantified=quantified)

    def test_no_de_readme_of_an_ids_only_search_never_says_quantified(self):
        text = self.readme("no-de", self.IDS, {"search_sage": False})
        self.assertNotIn("quantif", text.lower())
        self.assertNotIn("report.parquet", text)
        self.assertIn("identifications only: no quantities", text)
        self.assertIn("`search_sage/identifications.tsv`", text)

    def test_a_quantifying_search_keeps_its_words(self):
        """(pin) a quantifying search reads as before."""
        delivered = {"search_sage/report.parquet", "search_sage/search_provenance.json",
                     "Analysis_Report.html"}
        text = self.readme("no-de", delivered, {"search_sage": True})
        self.assertIn("searched to identify and quantify proteins", text)
        self.assertIn("report.parquet and the matrices delivered", text)
        # no record at all (an older caller): as before
        self.assertIn("identify and quantify", self.readme("no-de", delivered, None))

    def test_interim_readme_of_an_ids_only_search(self):
        delivered = {"search_sage/identifications.tsv", "search_sage/search_provenance.json"}
        text = self.readme("interim", delivered, {"search_sage": False})
        self.assertNotIn("quantif", text.lower())
        self.assertNotIn("statistical comparison", text)
        self.assertIn("identifications only", text)
        self.assertIn("There is no report in this folder yet.", text)
        q = self.readme("interim", {"search_dda/report.parquet",
                                    "search_dda/search_provenance.json"}, {"search_dda": True})
        self.assertIn("identified and quantified the proteins", q)
        self.assertIn("its report and matrices", q)
        self.assertIn("There is no report or statistical comparison in this folder yet.", q)

    def test_qc_panels_are_promised_only_with_figures(self):
        """(nit) a no-DE report has no figures (step 8b is skipped), so no QC panels."""
        self.assertNotIn("quality-control panels", self.readme("no-de", self.IDS,
                                                               {"search_sage": False}))
        self.assertIn("quality-control panels",
                      self.readme("no-de", self.IDS | {"figures/pca.png"}, {"search_sage": False}))


class TheE2EShape(unittest.TestCase):
    """The two-run DDA end-to-end session (E2E low list): Sage --ids-only in the default
    output/search/, delivered --mode no-de. The predicate is asked of the session, as deliver
    asks it, and the README never names report.parquet."""

    def test_no_de_readme_of_the_default_search_folder(self):
        with tempfile.TemporaryDirectory() as t:
            out = os.path.join(t, "output")
            folder(out, "search", {"identifications.tsv": "run\tpeptide\n",
                                   "sage_identifications.json": {"quantified": False},
                                   slc.RECORD: IDS_ONLY_REC,
                                   "search_provenance.json": {"engine": "sage",
                                                              "quantified": False}})
            delivered = {"search/identifications.tsv", "search/sage_identifications.json",
                         "search/search_provenance.json", "Analysis_Report.html", "methods.md"}
            q = cs.searches_quantified(out, cs.delivered_search_folders(delivered, "no-de",
                                                                        "search"))
            self.assertEqual(q, {"search": False})
            text = cs.build_readme(SUMMARY, delivered, [], "no-de", search_dir="search",
                                   quantified=q)
            self.assertNotIn("report.parquet", text)
            self.assertNotIn("quantif", text.lower())
            self.assertIn("`search/identifications.tsv`", text)


class InterimOfAnIdsOnlySearch(tdi.InterimBase):
    """deliver --mode interim, end to end, of a search for identifications only."""

    def setUp(self):
        super().setUp()
        for rel in ("report.parquet", "report.pg_matrix.tsv"):
            os.remove(os.path.join(self.dda, rel))
        tcs.write(self.dda, "identifications.tsv", "run\tpeptide\nr1\tPEPTIDE\n")
        tcs.write(self.dda, "sage_identifications.json", json.dumps({"quantified": False}))
        tcs.write(self.dda, slc.RECORD, json.dumps(IDS_ONLY_REC))

    def test_readme_and_record_say_identifications_only(self):
        rc, out, p = self.interim("--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        readme = tcs.read(os.path.join(self.delivery, "README.md"))
        self.assertNotIn("quantif", readme.lower())
        self.assertIn("identifications only", readme)
        self.assertEqual(self.delivery_json()["search_quantified"], {"search_dda": False})


class InterimOfAQuantifyingSearch(tdi.InterimBase):
    def test_the_record_says_quantified(self):
        rc, out, p = self.interim("--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertEqual(self.delivery_json()["search_quantified"], {"search_dda": True})
        self.assertIn("identified and quantified the proteins",
                      tcs.read(os.path.join(self.delivery, "README.md")))


class TheEmailDraft(unittest.TestCase):
    """email-draft from delivery.json alone (test_deliver_interim.InterimOutward's fixture)."""
    setUp = tdi.InterimOutward.setUp
    tearDown = tdi.InterimOutward.tearDown
    draft = tdi.InterimOutward.draft

    def test_an_ids_only_interim_draft_never_says_quantified(self):
        text = self.draft(search_quantified={"search_dda": False})
        self.assertNotIn("identified and quantified", text)
        self.assertIn("identified the peptides and proteins in each sample (identifications "
                      "only: no quantities)", text)

    def test_a_quantifying_or_older_record_keeps_its_words(self):
        for sq in ({"search_dda": True}, None):
            with self.subTest(sq=sq):
                self.assertIn("identified and quantified the proteins in each sample",
                              self.draft(search_quantified=sq))


if __name__ == "__main__":
    unittest.main()
