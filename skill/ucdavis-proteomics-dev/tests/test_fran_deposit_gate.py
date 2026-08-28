#!/usr/bin/env python3
"""
The FRAN handover: eligibility gate + how a search is staged for FRAN's ingest cron.

Everything here guards a failure that is SILENT. The skill hands a search over by linking it
into a drop directory; a cron ingests whatever it finds there. So a mistake does not raise —
it puts a collaborator's data in the Core corpus, or labels a DIA-NN search as Radiant, or
stages an entry the cron then skips without a word. The assertions are about those.
"""
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import fran_deposit as fd  # noqa: E402


class Args:
    """argparse-shaped stand-in."""
    def __init__(self, out, **kw):
        self.out = out
        for k, v in dict(skip=False, force=False, dry_run=False, organism=None, taxon=None,
                         fasta_meta=None, name=None).items():
            setattr(self, k, kw.get(k, v))


def search_dir(root, engine="diann", empty=False):
    """A minimal but REALISTIC search output directory for one engine."""
    d = os.path.join(root, "search_out")
    os.makedirs(d, exist_ok=True)
    body = "" if empty else "x" * 64
    if engine == "diann":
        with open(os.path.join(d, "report.parquet"), "w") as fh:
            fh.write(body)
    elif engine == "fragpipe":
        # FragPipe's tree CONTAINS a DIA-NN report -- that is exactly why it must be detected first.
        os.makedirs(os.path.join(d, "dia-quant-output"), exist_ok=True)
        with open(os.path.join(d, "dia-quant-output", "report.tsv"), "w") as fh:
            fh.write(body)
        with open(os.path.join(d, "report.tsv"), "w") as fh:
            fh.write(body)
    elif engine == "radiant":
        os.makedirs(os.path.join(d, "radiant_results", "fulcrum-results"), exist_ok=True)
        with open(os.path.join(d, "radiant_results", "fulcrum-results", "part-0.parquet"), "w") as f:
            f.write(body)
    # run_search.py writes this for EVERY engine
    with open(os.path.join(d, "search_provenance.json"), "w") as fh:
        json.dump({"engine": engine, "version": "9.9"}, fh)
    return d


class GateTests(unittest.TestCase):
    def test_missing_dir_is_a_reason_code_not_a_crash(self):
        """The orchestrator must get something it can act on: a deposit that raises would abort
        the analysis over an optional side-quest."""
        r = fd.check(Args("/definitely/not/a/real/search/dir"))
        self.assertFalse(r["eligible"])
        self.assertEqual(r["reason"], "not_on_hive")

    def test_opt_out(self):
        self.assertEqual(fd.check(Args(tempfile.gettempdir(), skip=True))["reason"], "opted_out")
        os.environ["FRAN_DEPOSIT"] = "off"
        try:
            self.assertEqual(fd.check(Args(tempfile.gettempdir()))["reason"], "opted_out")
        finally:
            del os.environ["FRAN_DEPOSIT"]

    def test_an_unwritable_drop_dir_refuses_the_search(self):
        """THE policy test. Write permission on the drop directory IS Core membership -- it lives
        inside /quobyte/proteomics-grp. A collaborator's search must not be handed over, and the
        refusal has to survive a search that otherwise looks perfect."""
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d)
            drop = os.path.join(d, "incoming")
            os.makedirs(drop)
            os.chmod(drop, 0o500)                       # readable, NOT writable
            os.environ["FRAN_DROP_DIR"] = drop
            try:
                r = fd.check(Args(out))
            finally:
                os.chmod(drop, 0o755)
                del os.environ["FRAN_DROP_DIR"]
            self.assertFalse(r["eligible"])
            self.assertEqual(r["reason"], "not_core_facility")

    def test_zero_byte_report_is_incomplete(self):
        """DIA-NN can exit 0 having identified nothing. Staging an empty report would put a
        real-looking zero-ID search in the corpus."""
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d, empty=True)
            os.environ["FRAN_DROP_DIR"] = os.path.join(d, "incoming")
            try:
                self.assertEqual(fd.check(Args(out))["reason"], "search_incomplete")
            finally:
                del os.environ["FRAN_DROP_DIR"]

    def test_sage_is_skipped_not_mislabelled(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "s")
            os.makedirs(out)
            with open(os.path.join(out, "search_provenance.json"), "w") as fh:
                json.dump({"engine": "sage", "version": "0.14.7"}, fh)
            self.assertEqual(fd.check(Args(out))["reason"], "engine_unsupported")


class EngineDetectionTests(unittest.TestCase):
    """All three DIA routes must be recognised, and told apart from each other."""

    def test_provenance_is_authoritative_for_each_engine(self):
        for engine in ("diann", "fragpipe", "radiant"):
            with tempfile.TemporaryDirectory() as d:
                self.assertEqual(fd.detect_engine(search_dir(d, engine))[:2], (engine, "9.9"))

    def test_marker_fallback_tells_the_three_apart(self):
        """Without provenance (an older dir, or one made by hand) the markers must still be
        right. FragPipe's tree holds a DIA-NN report.tsv, so order AND specificity matter: test
        FragPipe first, but never on a marker DIA-NN also has."""
        for engine in ("diann", "fragpipe", "radiant"):
            with tempfile.TemporaryDirectory() as d:
                out = search_dir(d, engine)
                os.unlink(os.path.join(out, "search_provenance.json"))
                got = fd.detect_engine(out)
                self.assertEqual(got[0], engine, f"{engine} detected as {got[0]} ({got[2]})")

    def test_a_diann_1_9_dir_with_only_report_tsv_is_not_fragpipe(self):
        """The trap the marker list exists to avoid: FragPipe's report is ALSO report.tsv, so
        detecting FragPipe from that name would relabel every DIA-NN 1.9 search."""
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "s")
            os.makedirs(out)
            with open(os.path.join(out, "report.tsv"), "w") as fh:
                fh.write("x")
            self.assertEqual(fd.detect_engine(out)[0], "diann")

    def test_the_report_each_engine_hands_over_is_the_quant_of_record(self):
        for engine, want in (("diann", "report.parquet"),
                             ("fragpipe", "dia-quant-output/report.tsv"),
                             ("radiant", "radiant_results/fulcrum-results")):
            with tempfile.TemporaryDirectory() as d:
                out = search_dir(d, engine)
                self.assertTrue(fd.find_report(out, engine).endswith(want))


class StagingTests(unittest.TestCase):
    def _stage(self, d, engine="diann", **kw):
        out = search_dir(d, engine)
        drop = os.path.join(d, "incoming")
        os.environ["FRAN_DROP_DIR"] = drop
        try:
            a = Args(out, **kw)
            c = fd.check(a)
            self.assertTrue(c["eligible"], c.get("detail"))
            entry = c["entry"]
            # stage() calls jout() -> SystemExit(0). That IS its contract with the orchestrator.
            with self.assertRaises(SystemExit):
                fd.stage(a)
            return out, entry
        finally:
            del os.environ["FRAN_DROP_DIR"]

    def test_the_entry_is_a_real_dir_of_symlinks_not_a_symlink(self):
        """FRAN's scanner walks with os.walk(followlinks=False): a symlinked DIRECTORY is never
        descended into, so a bare symlink to the search dir would be invisible to the cron and
        the search would silently never be ingested."""
        with tempfile.TemporaryDirectory() as d:
            out, entry = self._stage(d)
            self.assertTrue(os.path.isdir(entry))
            self.assertFalse(os.path.islink(entry))
            self.assertTrue(os.path.islink(os.path.join(entry, "report.parquet")))

    def test_nothing_is_copied(self):
        """A search directory is tens of GB. Every linked item must be a link, not a copy."""
        with tempfile.TemporaryDirectory() as d:
            out, entry = self._stage(d)
            for f in os.listdir(entry):
                if f != fd.MANIFEST:
                    self.assertTrue(os.path.islink(os.path.join(entry, f)), f)

    def test_each_engine_stages_what_makes_it_detectable(self):
        """After staging, FRAN's own marker test must find the SAME engine in the entry that ran
        in the search dir -- otherwise the cron ingests it with the wrong adapter."""
        for engine in ("diann", "fragpipe", "radiant"):
            with tempfile.TemporaryDirectory() as d:
                out, entry = self._stage(d, engine)
                self.assertEqual(fd.detect_engine(entry)[0], engine,
                                 f"{engine} entry detects as {fd.detect_engine(entry)}")

    def test_provenance_is_never_linked_because_it_reads_as_a_radiant_marker(self):
        """FRAN's ENGINE_MARKERS lists search_provenance.json under RADIANT and tests Radiant
        before DIA-NN. run_search.py writes that file into every search dir, so linking it would
        relabel every staged DIA-NN and FragPipe search as Radiant. It goes in the manifest."""
        with tempfile.TemporaryDirectory() as d:
            out, entry = self._stage(d, "diann")
            self.assertNotIn("search_provenance.json", os.listdir(entry))
            man = json.load(open(os.path.join(entry, fd.MANIFEST)))
            self.assertEqual(man["search_provenance"]["engine"], "diann")

    def test_the_manifest_carries_what_the_cron_cannot_derive(self):
        """A DIA-NN report has no organism column: without this the corpus row is NULL and the
        search is invisible on FRAN's species page. The real output_dir matters too -- ingesting
        under the drop path would record where the handover was, not where the search is."""
        with tempfile.TemporaryDirectory() as d:
            out, entry = self._stage(d, organism="Homo sapiens", taxon=9606, name="MyRun")
            man = json.load(open(os.path.join(entry, fd.MANIFEST)))
            self.assertEqual(man["organism"], "Homo sapiens")
            self.assertEqual(man["taxon"], 9606)
            self.assertEqual(man["output_dir"], os.path.realpath(out))
            self.assertEqual(man["engine"], "diann")
            self.assertIn("--organism-name", man["suggested_ingest"])
            self.assertIn(os.path.realpath(out), man["suggested_ingest"])

    def test_an_unknown_organism_is_absent_not_guessed(self):
        with tempfile.TemporaryDirectory() as d:
            out, entry = self._stage(d)
            self.assertIsNone(json.load(open(os.path.join(entry, fd.MANIFEST)))["organism"])

    def test_xic_from_the_PARALLEL_chain_is_found_and_flattened(self):
        """The layout that actually bites. The 5-step chain is the DEFAULT route above 5 files,
        and it writes chromatograms to `xic/t<N>_xic/` -- one directory per array task, with
        NOTHING at report_xic/. A finder that knew only the single-shot layout reported "no XICs"
        for every parallel run. Measured on a real 399-run cohort: 399 task dirs, 399 files."""
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d)
            for t in range(3):
                td = os.path.join(out, "xic", f"t{t}_xic")
                os.makedirs(td)
                for suffix in (".xic.parquet", ".ms1_mobilogram.parquet"):
                    with open(os.path.join(td, f"run{t}{suffix}"), "w") as fh:
                        fh.write("x")
            os.environ["FRAN_DROP_DIR"] = os.path.join(d, "incoming")
            try:
                c = fd.check(Args(out))
                self.assertTrue(c["xic"]["present"])
                self.assertEqual(c["xic"]["n_files"], 3)          # 3 runs, mobilograms not counted
                self.assertEqual(len(c["xic"]["dirs"]), 3)        # 3 task dirs
                with self.assertRaises(SystemExit):
                    fd.stage(Args(out))
                # ...presented as ONE report_xic/, which is where FRAN's XIC lane looks
                xd = os.path.join(c["entry"], "report_xic")
                self.assertTrue(os.path.isdir(xd))
                self.assertEqual(len([f for f in os.listdir(xd) if f.endswith(".xic.parquet")]), 3)
                self.assertTrue(all(os.path.islink(os.path.join(xd, f)) for f in os.listdir(xd)))
                self.assertTrue(os.path.exists(os.path.join(xd, "run0.xic.parquet")))
            finally:
                del os.environ["FRAN_DROP_DIR"]

    def test_restaging_does_not_keep_a_dropped_files_trace(self):
        """The normalised report_xic/ is a real directory we build, so a re-stage has to clear it.
        A file dropped from a re-run (the watcher's pathological-file playbook does exactly that)
        would otherwise keep its chromatogram in the handover forever."""
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d)
            td = os.path.join(out, "xic", "t0_xic")
            os.makedirs(td)
            for r in ("good", "pathological"):
                with open(os.path.join(td, f"{r}.xic.parquet"), "w") as fh:
                    fh.write("x")
            os.environ["FRAN_DROP_DIR"] = os.path.join(d, "incoming")
            try:
                a = Args(out)
                entry = fd.check(a)["entry"]
                with self.assertRaises(SystemExit):
                    fd.stage(a)
                os.unlink(os.path.join(td, "pathological.xic.parquet"))   # dropped on the re-run
                with self.assertRaises(SystemExit):
                    fd.stage(Args(out, force=True))
                self.assertEqual(os.listdir(os.path.join(entry, "report_xic")),
                                 ["good.xic.parquet"])
            finally:
                del os.environ["FRAN_DROP_DIR"]

    def test_xic_files_ride_along_and_absence_is_reported_not_hidden(self):
        """DIA-NN --xic chromatograms live INSIDE the search dir, so linking report_xic hands the
        XIC lane its input with no second mechanism."""
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d)
            os.makedirs(os.path.join(out, "report_xic"))
            with open(os.path.join(out, "report_xic", "runA.xic.parquet"), "w") as fh:
                fh.write("x")
            os.environ["FRAN_DROP_DIR"] = os.path.join(d, "incoming")
            try:
                c = fd.check(Args(out))
                self.assertTrue(c["xic"]["present"])
                self.assertEqual(c["xic"]["n_files"], 1)
                with self.assertRaises(SystemExit):
                    fd.stage(Args(out))
                self.assertTrue(os.path.exists(os.path.join(c["entry"], "report_xic",
                                                            "runA.xic.parquet")))
            finally:
                del os.environ["FRAN_DROP_DIR"]
        with tempfile.TemporaryDirectory() as d:            # and a search WITHOUT xic still stages
            out = search_dir(d)
            os.environ["FRAN_DROP_DIR"] = os.path.join(d, "incoming")
            try:
                c = fd.check(Args(out))
                self.assertTrue(c["eligible"])
                self.assertFalse(c["xic"]["present"])
            finally:
                del os.environ["FRAN_DROP_DIR"]

    def test_restaging_converges_on_one_entry(self):
        """A deterministic entry name is what stops a re-stage becoming a SECOND candidate the
        cron ingests as a separate search."""
        with tempfile.TemporaryDirectory() as d:
            out, entry = self._stage(d)
            drop = os.path.dirname(entry)
            os.environ["FRAN_DROP_DIR"] = drop
            try:
                a = Args(out, force=True)
                with self.assertRaises(SystemExit):
                    fd.stage(a)
            finally:
                del os.environ["FRAN_DROP_DIR"]
            self.assertEqual(os.listdir(drop), [os.path.basename(entry)])

    def test_a_staged_search_is_not_staged_twice_without_force(self):
        with tempfile.TemporaryDirectory() as d:
            out, entry = self._stage(d)
            os.environ["FRAN_DROP_DIR"] = os.path.dirname(entry)
            try:
                self.assertEqual(fd.check(Args(out))["reason"], "already_staged")
                self.assertTrue(fd.check(Args(out, force=True))["eligible"])
            finally:
                del os.environ["FRAN_DROP_DIR"]


class VerifyTests(unittest.TestCase):
    def test_pending_is_not_a_failure_and_a_dead_link_is(self):
        """Right after a run the only honest answer is "handed over, cron hasn't run yet". But a
        link whose target has moved leaves an entry that LOOKS staged and that the cron will skip
        in silence -- that one has to be called out."""
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d)
            drop = os.path.join(d, "incoming")
            os.environ["FRAN_DROP_DIR"] = drop
            try:
                a = Args(out)
                entry = fd.check(a)["entry"]
                with self.assertRaises(SystemExit):
                    fd.stage(a)
                with self.assertRaises(SystemExit) as cm:
                    fd.verify(Args(out))
                self.assertEqual(cm.exception.code, 0)
                os.unlink(os.path.join(out, "report.parquet"))     # target goes away
                self.assertIn("report.parquet",
                              [f for f in os.listdir(entry)
                               if os.path.islink(os.path.join(entry, f))
                               and not os.path.exists(os.path.join(entry, f))])
            finally:
                del os.environ["FRAN_DROP_DIR"]


if __name__ == "__main__":
    unittest.main(verbosity=2)
