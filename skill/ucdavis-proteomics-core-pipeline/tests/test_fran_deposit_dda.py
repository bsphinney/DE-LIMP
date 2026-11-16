#!/usr/bin/env python3
"""
FRAN refuses a DDA search, with an honest reason (2.15, Ruling 57.2, review M8).

FRAN's corpus is DIA-only. Before 2.15 a FragPipe DDA search reached `check` as an engine FRAN
supports and was turned away as `search_incomplete` (it has no dia-quant-output/report.tsv): a
FAILED search, which it was not -- and a DIA-NN --dda search (report.parquet) was staged, so its
rows would have entered the corpus under a DIA label. Now both are `dda_search`, read from the
search's own files (fran_deposit.search_acquisition): FragPipe's fp-manifest data types, the
DIA-NN cfg's --dda (estimate_params.DIANN_DDA_FLAG, the one definition), DIA-NN's log, else the
bundle. A DIA-NN DIA search is staged exactly as before (pin), and so is a FragPipe DIA one.
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

import fran_deposit as fd  # noqa: E402
from test_fran_deposit_gate import Args, search_dir  # noqa: E402

# the QC name rule reads the temp folder names; keep "_" out of them (test_fran_deposit_gate)
tempfile._RandomNameSequence.characters = "abcdefghijklmnopqrstuvwxyz0123456789"


def _write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    return path


class FranRefusesDda(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(__import__("shutil").rmtree, self.d, True)
        self.saved = {k: os.environ.get(k) for k in ("FRAN_DROP_DIR", "FRAN_DEPOSIT",
                                                     "FRAN_CORPUS_QUERY", "FRAN_HEALTH")}
        os.environ["FRAN_DROP_DIR"] = os.path.join(self.d, "incoming")
        os.environ["FRAN_CORPUS_QUERY"] = "off"
        os.environ["FRAN_HEALTH"] = "off"
        os.environ.pop("FRAN_DEPOSIT", None)

        def restore():
            for k, v in self.saved.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v
        self.addCleanup(restore)

    def provenance(self, out, **rec):
        import fasta_fixture
        db = fasta_fixture.write_real_db(os.path.dirname(out))
        with open(os.path.join(out, "search_provenance.json"), "w") as fh:
            json.dump(dict({"version": "9.9", "fasta": db}, **rec), fh)

    def fragpipe_dda(self):
        """A finished FragPipe DDA search as run_search.py leaves it: report.parquet adapted from
        combined_protein.tsv, its fp-manifest (DDA rows), FragPipe's log ending ALL JOBS DONE."""
        out = os.path.join(self.d, "search_out")
        _write(os.path.join(out, "report.parquet"), "x" * 64)
        _write(os.path.join(out, "combined_protein.tsv"), "Protein\n")
        _write(os.path.join(out, "fragpipe.fp-manifest"),
               f"{out}/raw_stage/a.d\ta\t\tDDA\n{out}/raw_stage/b.d\tb\t\tDDA\n")
        _write(os.path.join(out, "log_2026-11-04_10-00-00.txt"),
               "=====ALL JOBS DONE IN 6.5 MINUTES=====\n")
        self.provenance(out, engine="fragpipe")
        return out

    def test_a_fragpipe_dda_search_is_dda_search_not_incomplete(self):
        out = self.fragpipe_dda()
        r = fd.check(Args(out))
        self.assertFalse(r["eligible"])
        self.assertEqual(r["reason"], "dda_search", r.get("detail"))
        self.assertIn("DIA corpus", r["detail"])
        self.assertIn("Nothing is wrong with the search", r["detail"])
        self.assertIn("fragpipe.fp-manifest", r["detail"])
        self.assertEqual(r["acquisition"]["value"], "DDA")

    def test_stage_never_hands_it_over(self):
        out = self.fragpipe_dda()
        with self.assertRaises(SystemExit):
            fd.stage(Args(out))
        self.assertFalse(os.path.exists(os.path.join(self.d, "incoming", fd.entry_name(out))))

    def test_a_diann_dda_search_is_refused_from_its_cfg(self):
        out = search_dir(self.d, "diann")
        cfg = _write(os.path.join(self.d, "wf", "diann.cfg"),
                     "--dda\n--cut K*,R*\n--missed-cleavages 1\n")
        self.provenance(out, engine="diann", params_file=cfg)
        r = fd.check(Args(out))
        self.assertEqual(r["reason"], "dda_search", r.get("detail"))
        self.assertIn("DIA-NN with --dda", r["detail"])
        self.assertIn("diann.cfg has --dda", r["detail"])

    def test_a_diann_dda_search_is_refused_from_its_log(self):
        """no cfg left to read: DIA-NN's own log echoes its command line"""
        out = search_dir(self.d, "diann")
        _write(os.path.join(out, "report.log.txt"),
               "DIA-NN 2.7.0\ndiann --f a.raw --dda --fasta x.fasta --out report.parquet\n")
        self.provenance(out, engine="diann")
        self.assertEqual(fd.check(Args(out))["reason"], "dda_search")

    def test_a_diann_dia_search_is_still_staged(self):
        """(pin) a DIA-NN DIA search: its cfg has no --dda, and it is handed over as before"""
        out = search_dir(self.d, "diann")
        cfg = _write(os.path.join(self.d, "wf", "diann.cfg"), "--cut K*,R*\n--missed-cleavages 1\n")
        self.provenance(out, engine="diann", params_file=cfg)
        r = fd.check(Args(out))
        self.assertTrue(r["eligible"], r.get("detail"))
        self.assertEqual(r["acquisition"]["value"], "DIA")
        with self.assertRaises(SystemExit):
            fd.stage(Args(out))
        self.assertTrue(os.path.isdir(os.path.join(self.d, "incoming", fd.entry_name(out))))

    def test_a_fragpipe_dia_search_is_still_eligible(self):
        """(pin) FragPipe's diaTracer reuse lists the pseudo-spectra as DDA rows beside the
        DIA-Quant .d: any DIA row makes the search DIA"""
        out = search_dir(self.d, "fragpipe")
        _write(os.path.join(out, "fragpipe.fp-manifest"),
               f"{out}/diatracer_stage/a_diatracer.mzML\ta\t\tDDA\n"
               f"{out}/diatracer_stage/a.d\ta\t\tDIA-Quant\n")
        self.assertTrue(fd.check(Args(out))["eligible"])

    def test_the_bundle_decides_when_nothing_else_does(self):
        out = search_dir(self.d, "diann")
        bundle = _write(os.path.join(self.d, "wf", "workflow.manifest.json"),
                        json.dumps({"acquisition": "DDA"}))
        self.provenance(out, engine="diann", bundle=bundle)
        self.assertEqual(fd.search_acquisition(out, {"bundle": bundle})[0], "DDA")
        self.assertEqual(fd.check(Args(out))["reason"], "dda_search")

    def test_nothing_said_is_not_a_refusal(self):
        out = search_dir(self.d, "diann")
        self.assertEqual(fd.search_acquisition(out, {})[0], None)
        self.assertTrue(fd.check(Args(out))["eligible"])


if __name__ == "__main__":
    unittest.main()
